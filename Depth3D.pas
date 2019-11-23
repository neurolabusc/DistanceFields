program Depth3D;
//Compute each voxels distance from edge of cluster
//   fpc -CX -Xs -XX -O3 Depth3D

{$mode Delphi} //{$mode objfpc}
{$H+}
uses 
	{$ifdef unix}
	cthreads, cmem, // the c memory manager is on some systems much faster for multi-threading
	{$endif}
	distance_field, SimdUtils, dateutils, StrUtils, sysutils, Classes, nifti_types, 
	nifti_loadsave, nifti_foreign, resize;

const
    kEXIT_SUCCESS = 0;
    kEXIT_FAIL = 1;
    kEXIT_PARTIALSUCCESS = 2; //processed some but not all input files

function distanceFieldAll(fnm, outName: string; isGz, is3D, isImg, isTxt: boolean; threshold: single = 0.5; maxthreads: integer = 0; superSample: integer = 1): boolean;
var
	hdr: TNIFTIhdr;
	img: TUInt8s;
	isInputNIfTI: boolean;
	txtNam, ext: string;
begin
	result := false;
	if not loadVolumes(fnm, hdr, img, isInputNIfTI) then exit;
	if (hdr.dim[1] < 2) or (hdr.dim[2] < 2) then begin
		writeln('File dimensions too small');
		exit;
	end;
	if outName = '' then 
		outName := extractfilepath(fnm)+'x'+extractfilename(changefileext(fnm,''))
	else
		outName := changefileext(outName,'');
	//handle double extensions: img.nii.gz and img.BRIK.gz
	ext := upcase(extractfileext(outName));
	if (ext = '.NII') or (ext = '.BRIK') then
		outName := changefileext(outName,'');
	txtNam := '';
	if isTxt then
		txtNam := outName+'.1D';
	hdr.intent_code := kNIFTI_INTENT_NONE; //just in case input is labelled map
	if (superSample > 1) then begin
		if (threshold = 0) or (hdr.dim[4] > 1) then begin 
			writeln('Error: super sampling not supported for atlases or 4D images.');
			exit;
		end;
		changeDataType(hdr, img, kDT_FLOAT32);
		result := ShrinkOrEnlarge(hdr, img, superSample, maxthreads);
		if not result then begin
			writeln('Supersampling error');
			exit;
		end;
		writeln(format('Depth estimated on x%d supersampling (%d*%d*%d).',[superSample, hdr.dim[1], hdr.dim[2], hdr.dim[3] ]));
		result := distanceFieldVolume(hdr, img, txtNam, threshold, maxthreads);
		ShrinkOrEnlarge(hdr, img, 1.0/superSample, maxthreads);
	end else if (threshold = 0) then
		result := distanceFieldAtlas(hdr, img, txtNam, maxthreads)
	else
		result := distanceFieldVolume(hdr, img, txtNam, threshold, maxthreads);
	if not result then
		exit;
	ext := upcase(ExtractFileExt(fnm));
	if (ext = '.GZ') then
		fnm := changefileext(fnm,''); //e.g. file.raw.gz -> file.nii not file.raw.nii
	if saveNii(outName, hdr, img, isGz, is3D) then
		result := true;
end;

procedure showhelp;
var
    exeName, outDir, inDir: string;
begin
    exeName := extractfilename(ParamStr(0));
    {$IFDEF WINDOWS}
    exeName := ChangeFileExt(exeName, ''); //i2nii.exe -> i2nii
    {$ENDIF}
    writeln('Chris Rorden''s '+exeName+' v1.0.20191123');
    writeln(' see https://prideout.net/blog/distance_fields/');
	writeln(format('usage: %s [options] <in_file(s)>', [exeName]));
	writeln('Reads volume and computes distance fields');
	writeln('OPTIONS');
    writeln(' -3 : save 4D data as 3D files (y/n, default n)');
    writeln(' -h : show help');
    writeln(' -o : output name (omit to save as input name with "x" prefix)');
    writeln(' -r : generate text report (y/n/o[only, no image], default n)');
    writeln(' -t : threshold, less extreme values treated as outside (default 0.5)');
    writeln('       set to 0 for separate field for each region of an atlas');
    writeln(' -p : parallel threads (0=optimal, 1=one, 5=five, default 0)');
    writeln(' -s : supersample for continuous images (1=x1, 2=x2, 5=x5, default 1)');
    writeln(' -z : gz compress images (y/n, default n)');
    writeln(' Examples :');
    {$IFDEF WINDOWS}
    OutDir := 'c:\out dir\';
    InDir := 'c:\in\';
    {$ELSE}
    OutDir := '~/out dir/';
    InDir := './test/';
    {$ENDIF} 
    writeln(format('  %s -t 0 %sAICHA.nii.gz', [exeName, InDir]));
    writeln(format('  %s -t 0 %sinia19-NeuroMaps.nii.gz', [exeName, InDir]));
    writeln(format('  %s -t 0.5 %s4Dgraywhite.nii.gz', [exeName, InDir]));
    writeln(format('  %s -t 0.5 %sisotropic.nii.gz', [exeName, InDir]));
    writeln(format('  %s -t 0.25 -s 3 %savg152T1_gray.nii.gz', [exeName, InDir]));
    writeln(format('  %s -o "%sOutputImg" -t 0.5 "%sname with spaces"', [exeName, OutDir, InDir]));     	
end;

function doRun: integer;
var
	i, nOK, nAttempt: integer;
    is3D: boolean = false;
    isGz: boolean = false;
    isImg: boolean = true;
    isTxt: boolean = false;
    isShowHelp: boolean = false;
    startTime: TDateTime;
    maxthreads: integer = 0;
    threshold: single = 0.5; 
    superSample: integer = 1;
    outName: string = '';
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
            if distanceFieldAll(s, outName, isGz, is3D, isImg, isTxt, threshold, maxthreads, superSample) then
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
        if c =  'O' then
            outName := v;
        if c =  'P' then
            maxthreads := strtointdef(v, 0); 
        if c =  'R' then begin
            isTxt := upcase(v[1]) = 'Y'; //yes, save text table
            if (upcase(v[1]) = 'O') then begin //text only, no image data
           		isTxt := true;
           		isImg := false; 
            end;
        end;
        if c =  'S' then
            superSample := strtointdef(v, 1);    
        if c =  'T' then
            threshold := strtofloatdef(v, 0.5);    
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