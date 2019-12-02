program Depth3D;
//Compute each voxels distance from edge of cluster
//   fpc -CX -Xs -XX -O3 Depth3D

{$mode Delphi} //{$mode objfpc}
{$H+}
uses 
	{$ifdef unix}
	cthreads, cmem, // the c memory manager is on some systems much faster for multi-threading
	{$endif}
	SimdUtils, VectorMath,
	distance_field,  dateutils, StrUtils, sysutils, Classes, nifti_types, 
	nifti_loadsave, nifti_foreign, resize, math;

const
    kEXIT_SUCCESS = 0;
    kEXIT_FAIL = 1;
    kEXIT_PARTIALSUCCESS = 2; //processed some but not all input files
    kText_No = 0; //do not create text report
    kText_Yes = 1; //create text report and image
    kText_Screen = 2; //show text report on screen but do not create files
    kSave_Depth = 0;
    kSave_Intensity = 1;
    kSave_Both = 2;
    kSave_None = 3;
    kVers = 'v1.0.20191126';
    
    
function isIsotropic(var hdr: TNIFTIhdr): boolean;
var
	mn, mx: single;
begin
	mn := min(abs(hdr.pixdim[1]), min(abs(hdr.pixdim[2]), abs(hdr.pixdim[3])));
	mx := max(abs(hdr.pixdim[1]), max(abs(hdr.pixdim[2]), abs(hdr.pixdim[3])));
	if (mn = 0) or ((mx/mn) > 1.02) then
		exit(false);
	exit(true);
end;

function distanceFieldAll(fnm, outName: string; isGz, is3D: boolean; isImg, isTxt: integer; threshold: single = 0.5; maxthreads: integer = 0; superSample: integer = 1; clusterType: integer = 26; smallestClusterVox: integer = 1; maskName: string = ''): boolean;
var
	hdr: TNIFTIhdr;
	img, intensityImg: TUInt8s;
	isInputNIfTI: boolean;
	txtNam, ext: string;
	superSampleXYZ: TVec3;
begin
	result := false;
	intensityImg := nil;
	superSampleXYZ := Vec3(superSample, superSample, superSample);
	if (threshold = 0) and ((isImg = kSave_Intensity) or (isImg = kSave_Both)) then begin
		writeln('Saving intensity maps is not designed for atlases (where threshold = 0).');
		exit;
	end;
	if not loadVolumes(fnm, hdr, img, isInputNIfTI) then exit;
	if (hdr.dim[1] < 2) or (hdr.dim[2] < 2) then begin
		writeln('File dimensions too small');
		exit;
	end;
	if outName = '' then begin
		if isImg = kSave_Intensity then
			outName := extractfilepath(fnm)+'thresh_'+extractfilename(changefileext(fnm,''))
		else
			outName := extractfilepath(fnm)+'depth_'+extractfilename(changefileext(fnm,''))
	end else
		outName := changefileext(outName,'');
	//handle double extensions: img.nii.gz and img.BRIK.gz
	ext := upcase(extractfileext(outName));
	if (ext = '.NII') or (ext = '.BRIK') then
		outName := changefileext(outName,'');
	txtNam := '';
	if (isTxt = kText_Yes) then
		txtNam := outName+'.1D';
	if (isTxt = kText_Screen) then
		txtNam := '-';
	hdr.intent_code := kNIFTI_INTENT_NONE; //just in case input is labelled map
	//if clusterType > 0 then begin
	//	if not makeClusters(hdr, img, threshold, clusterType, smallestClusterMM3) then exit;	
	//end;
	
	if (superSample > 1) or (not isIsotropic(hdr)) then begin
		if (threshold = 0) or (hdr.dim[4] > 1) then begin 
			if superSample <= 1 then
				writeln('Error: atlases and 4D images must be isotropic.')
			else
				writeln('Error: super sampling not supported for atlases or 4D images.');
			exit;
		end;
		if (isTxt <> kText_No) then
			writeln('Be aware the report is based on up-sampled voxels');	
		if (smallestClusterVox > 1) then
			writeln('Be aware that minimum cluster extent based on up-sampled voxels');
		changeDataType(hdr, img, kDT_FLOAT32);
		result := ShrinkOrEnlarge(hdr, img, superSampleXYZ, maxthreads);
		if not result then begin
			writeln('Supersampling error');
			exit;
		end;
		writeln(format('Depth estimated on %g*%g*%g supersampling (%d*%d*%d voxels).',[superSampleXYZ.X, superSampleXYZ.Y, superSampleXYZ.Z, hdr.dim[1], hdr.dim[2], hdr.dim[3] ]));
		result := distanceFieldVolume(hdr, img, intensityImg, txtNam, threshold, maxthreads, clusterType, smallestClusterVox, maskName);
		superSampleXYZ := 1.0 / superSampleXYZ;
		ShrinkOrEnlarge(hdr, img, superSampleXYZ, maxthreads);
	end else if (threshold = 0) then begin
		if (smallestClusterVox > 1) then
			writeln('Warning: minimum cluster size not used for atlases');
		result := distanceFieldAtlas(hdr, img, txtNam, maxthreads, maskName)
	end else
		result := distanceFieldVolume(hdr, img, intensityImg, txtNam, threshold, maxthreads, clusterType, smallestClusterVox, maskName);
	if not result then exit;
	if (isImg = kSave_None) then exit;
	hdr.descrip := 'Depth3D '+kVers;
	if (isImg = kSave_Depth) or (isImg = kSave_Both) then 
		result := saveNii(outName, hdr, img, isGz, is3D);
	if not result then  exit;
	if (intensityImg = nil) and ((isImg = kSave_Intensity) or (isImg = kSave_Both)) then begin
		result := false;
		writeln('Unable to create intensity image');
		exit;
	end;
	if (isImg = kSave_Intensity) then 
		result := saveNii(outName, hdr, intensityImg, isGz, is3D);
	if (isImg = kSave_Both) then
		result := saveNii(outName+'_intensity', hdr, intensityImg, isGz, is3D);	
end;

procedure showhelp;
var
    exeName, outDir, inDir: string;
begin
    exeName := extractfilename(ParamStr(0));
    {$IFDEF WINDOWS}
    exeName := ChangeFileExt(exeName, ''); //i2nii.exe -> i2nii
    {$ENDIF}
    writeln('Chris Rorden''s '+exeName+' '+kVers);
    writeln(format('usage: %s [options] <in_file(s)>', [exeName]));
	writeln('Reads volume and computes distance fields');
	writeln('OPTIONS');
    writeln(' -3 : save 4D data as 3D files (y/n, default n)');
    writeln(' -c : connectivity neighbors (6=faces, 18=edges, 26=corners, default 26)');
    writeln(' -h : show help');
    writeln(' -i : inverted mask image (ignore voxels with value of non-zero in mask)');
    writeln(' -k : mask image (ignore voxels with value of zero in mask)');
    writeln(' -o : output name (omit to save as input name with "depth_" prefix)');
    writeln(' -r : report table (y/n/s: yes, no, screen, default n)');
    writeln(' -t : threshold, less extreme values treated as outside (default 0.5)');
    writeln('       set to 0 for separate field for each region of an atlas');
    writeln(' -m : minimum cluster extent in voxels (default 1)');
    writeln(' -p : parallel threads (0=optimal, 1=one, 5=five, default 0)');
    writeln(' -s : save images (d,i,b,n: depth, intensity, both, none, default t) ');
    writeln(' -u : upsample for continuous images (1=x1, 2=x2, 5=x5, default 1)');
    writeln(' -z : gz compress images (y/n, default n)');
    writeln(' Examples :');
    {$IFDEF WINDOWS}
    OutDir := 'c:\out dir\';
    InDir := '.\test\';
    {$ELSE}
    OutDir := '~/out dir/';
    InDir := './test/';
    {$ENDIF} 
    writeln(format('  %s -t 0 %sAICHA.nii.gz', [exeName, InDir]));
    writeln(format('  %s -t 0 %sinia19-NeuroMaps.nii.gz', [exeName, InDir]));
    writeln(format('  %s -t 0.5 %s4Dgraywhite.nii.gz', [exeName, InDir]));
    writeln(format('  %s -t 0.5 %sisotropic.nii.gz', [exeName, InDir]));
    writeln(format('  %s -t 0.25 -u 3 %savg152T1_gray.nii.gz', [exeName, InDir]));
    writeln(format('  %s -t 3 -m 5 -r s -s n %smotor.nii.gz', [exeName, InDir]));
    writeln(format('  %s -k %sisotropic_mask.nii.gz %sisotropic.nii.gz', [exeName, InDir, InDir]));
    writeln(format('  %s -o "%sOutputImg" -t 0.5 "%sname with spaces"', [exeName, OutDir, InDir]));
end;

function doRun: integer;
var
	i, nOK, nAttempt: integer;
	clusterType: integer = 26;//AFNI:1=faces,2=edges,3=corners or 6=faces,18=edges,26=corners
    is3D: boolean = false;
    isGz: boolean = false;
    isTxt: integer = kText_No;
    isImg: integer = kSave_Depth;
    isShowHelp: boolean = false;
    smallestClusterVox: integer = 1;
    startTime: TDateTime;
    maxthreads: integer = 0;
    threshold: single = 0.5; 
    superSample: integer = 1;
    outName: string = '';
    maskName: string = '';
    s, v: string;
    c: char;
begin
	startTime := Now;
	result := kEXIT_SUCCESS;
	nOK := 0;
	nAttempt := 0;
	i := 1;
	while i <= ParamCount do begin
        s := ParamStr(i);
        i := i + 1;
        if length(s) < 1 then continue; //possible?
        if s[1] <> '-' then begin
            nAttempt := nAttempt + 1;
            if distanceFieldAll(s, outName, isGz, is3D, isImg, isTxt, threshold, maxthreads, superSample, clusterType, smallestClusterVox, maskName) then
                nOK := nOK + 1;
            continue;
        end;
        if length(s) < 2 then continue; //e.g. 'i2nii -'
        c := upcase(s[2]);
		if c = 'H' then showhelp;
		if i = ParamCount then continue;
        v := ParamStr(i);
        i := i + 1;
        if length(v) < 1 then continue; //e.g. 'i2nii -o ""'
        if c =  '3' then
            is3D := upcase(v[1]) = 'Y'; 
        if c =  'C' then
            clusterType := strtointdef(v, 26); 
        if c = 'I' then
        	maskName := '*'+v;
        if c = 'K' then
        	maskName := v;
        if c =  'M' then
            smallestClusterVox := strtointdef(v, 1);   
        if c =  'O' then
            outName := v;
        if c =  'P' then
            maxthreads := strtointdef(v, 0); 
        if c =  'R' then begin
            if upcase(v[1]) = 'N' then
            	isTxt := kText_No; //yes, save text table	
            if upcase(v[1]) = 'Y' then
            	isTxt := kText_Yes; //yes, save text table
            if upcase(v[1]) = 'S' then
            	isTxt := kText_Screen; //yes, save text table           
        end;
        if c =  'S' then begin
            if upcase(v[1]) = 'T' then
            	isImg := kSave_Depth;	
            if upcase(v[1]) = 'I' then
            	isImg := kSave_Intensity;
            if upcase(v[1]) = 'B' then
            	isImg := kSave_Both;   
            if upcase(v[1]) = 'N' then
            	isImg := kSave_None;          
        end;
        if c =  'T' then
            threshold := strtofloatdef(v, 0.5);    
        if c =  'U' then
            superSample := strtointdef(v, 1);    
        if c =  'Z' then
            isGz := upcase(v[1]) = 'Y';			
	end;
    if (ParamCount = 0) or (isShowHelp) then
        ShowHelp;
    if nOK > 0 then 
        writeln(format('Filter required %.3f seconds.', [MilliSecondsBetween(Now,startTime)/1000.0]));
    if (nOK = nAttempt) then
        ExitCode := kEXIT_SUCCESS
    else if (nOK = 0) then
        ExitCode := kEXIT_FAIL
    else
        ExitCode := kEXIT_PARTIALSUCCESS;
end;

begin
    ExitCode := doRun;
end.