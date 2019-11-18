unit nifti_loadsave;
//written by William ('Bill') G. Miller 2004, released under BSD 2-Clause License
interface
{$mode Delphi}
{$DEFINE GZIP}
{$DEFINE BZIP2}
uses
  {$IFDEF BZIP2}bzip2stream,{$ENDIF}
  {$IFDEF GZIP}zstream, gziputils, {$ENDIF}
  {$IFNDEF UNIX}system, {$ENDIF}
  VectorMath, SimdUtils, nifti_types, 
  DateUtils, Classes, SysUtils, StrUtils, Math, 
  nifti_foreign;
  
const
    kDT_input = 0; //save data output as same datatype as input

function saveNii(fnm: string; var lHdr: TNIFTIhdr; var rawData: TUInt8s; isGz, is3D: boolean): boolean; overload;
function saveNii(fnm: string; var oHdr: TNIFTIhdr; var orawVolBytes: TUInt8s; isGz: boolean): boolean; overload;
function loadVolumes(var Filename: string; out lHdr: TNIFTIhdr; out rawData: TUInt8s; out isInputNIfTI: boolean): boolean;
procedure printf(s: string);
procedure changeDataType(var lHdr: TNIFTIhdr; var rawData: TUInt8s; outDT: integer; isVerbose: boolean = false);

implementation

procedure printf(s: string);
begin
	{$IFDEF UNIX}
     writeln(s);
	{$ELSE}
	if IsConsole then //uses System
		writeln(s);	
	{$ENDIF}    
end;

procedure changeDataType(var lHdr: TNIFTIhdr; var rawData: TUInt8s; outDT: integer; isVerbose: boolean = false);
var
   i8s: TInt8s;
   i16s: TInt16s;
   i32s: TInt32s;
   u8s: TUInt8s;
   u16s: TUInt16s;
   u32s: TUInt32s;
   f32s, x32s: TFloat32s;
   i, n: int64;
   inDT: word;
   omx, omn, mn, mx, intercept, slope, islope: double;
begin
	if (outDT = kDT_input) then exit;
	if ((outDT <> kDT_FLOAT32) and (outDT <> kDT_INT8) and (outDT <> kDT_INT16) and (outDT <> kDT_UINT8) and (outDT <> kDT_UINT16) and (outDT <> kDT_INT32) ) then begin
		printf('Unsupported output datatype');
		exit;
	end;
	inDT := lHdr.datatype; 
	//see if output format is input
	if (inDT = kDT_INT8) and (outDT = kDT_INT8) then exit;
	if (inDT = kDT_INT16) and (outDT = kDT_INT16) then exit;
	if (inDT = kDT_UINT8) and (outDT = kDT_UINT8) then exit;
	if (inDT = kDT_UINT16) and (outDT = kDT_UINT16) then exit;
	if (inDT = kDT_INT32) and (outDT = kDT_INT32) then exit;
	if (inDT = kDT_FLOAT32) and (outDT = kDT_FLOAT32) then exit;
	if ((lHdr.bitpix mod 8) <> 0) or (lHdr.bitpix < 8) then exit;
	n := length(rawData) div (lHdr.bitpix div 8);
	if n < 1 then exit;
	setlength(x32s, n);
	if (inDT = kDT_INT8) then begin
		i8s := TInt8s(rawData);
		for i := 0 to (n - 1) do
			x32s[i] := i8s[i];
	end else if (inDT = kDT_INT16) then begin
		i16s := TInt16s(rawData);
		for i := 0 to (n - 1) do
			x32s[i] := i16s[i];
	end else if (inDT = kDT_INT32) then begin	
		i32s := TInt32s(rawData);
		for i := 0 to (n - 1) do
			x32s[i] := i32s[i];
	end else if (inDT = kDT_UINT8) then begin
		u8s := TUInt8s(rawData);
		for i := 0 to (n - 1) do
			x32s[i] := u8s[i];
	end else if (inDT = kDT_UINT16) then begin
		u16s := TUInt16s(rawData);
		for i := 0 to (n - 1) do
			x32s[i] := u16s[i];
	end else if (inDT = kDT_UINT32) then begin
		u32s := TUInt32s(rawData);
		for i := 0 to (n - 1) do
			x32s[i] := u32s[i];
	end else if (inDT = kDT_FLOAT32) then begin
		f32s := TFloat32s(rawData);
		for i := 0 to (n - 1) do
			x32s[i] := f32s[i];		
	end else begin
		printf('Unable to swap datatype '+inttostr(lHdr.bitpix ));
		setlength(x32s, 0);	
		exit;
	end;
	if outDT = kDT_INT8 then lHdr.datatype := kDT_INT8;
	if outDT = kDT_INT16 then lHdr.datatype := kDT_INT16;
	if outDT = kDT_UINT8 then lHdr.datatype := kDT_UINT8;
	if outDT = kDT_UINT16 then lHdr.datatype := kDT_UINT16;
	if outDT = kDT_INT32 then lHdr.datatype := kDT_INT32;
	if outDT = kDT_FLOAT32 then lHdr.datatype := kDT_FLOAT32;
	if (outDT = kDT_INT8) or (outDT = kDT_UINT8) then
		lHdr.bitpix := 8
	else if (outDT = kDT_INT16) or (outDT = kDT_UINT16) then
		lHdr.bitpix := 16
	else //if kDT_FLOAT32, kDT_INT32) then
		lHdr.bitpix := 32;
	setlength(rawData, n * (lHdr.bitpix div 8));
	if lHdr.scl_slope = 0 then lHdr.scl_slope := 1;
	if (outDT = kDT_FLOAT32) then begin //output float - no need to scale range
		if (isVerbose) then printf('Converted to FLOAT32');
		f32s := TFloat32s(rawData);
		for i := 0 to (n - 1) do
				f32s[i] := (x32s[i] * lHdr.scl_slope)+ lHdr.scl_inter;
	   lHdr.scl_slope := 1.0;
	   lHdr.scl_inter := 0;
	   setlength(x32s, 0);
	   exit;
	end;
	mn := x32s[0];
	mx := mn;
	for i := 0 to (n - 1) do begin
		if x32s[i] < mn then mn := x32s[i];
		if x32s[i] > mx then mx := x32s[i];
	end;
	omn := 0; //unsigned
	omx := 255; //avoid compiler waring
	if (outDT = kDT_INT8) then omn := -128;
	if (outDT = kDT_INT16) then omn := -32768;
	if (outDT = kDT_INT32) then omn := -2147483648;
	if (outDT = kDT_INT8) then omx := 127;
	if (outDT = kDT_INT16) then omx := 32767;	
	if (outDT = kDT_UINT8) then omx := 255;
	if (outDT = kDT_UINT16) then omx := 65535;
	if (outDT = kDT_INT32) then omx := 2147483647;
	if mn > 0 then mn := 0;
	if (inDT <> kDT_FLOAT32) and (omn <= mn) and (omx >= mx) then begin //integer with no need to downsample
		//printf(format('Lossless Conversion: input range %g..%g output range %g..%g', [mn,mx, omn, omx]));
		intercept := 0;
		slope := 1;
	end else begin
		for i := 0 to (n - 1) do
			x32s[i] := (x32s[i] * lHdr.scl_slope)+ lHdr.scl_inter;
		mn := (mn * lHdr.scl_slope)+ lHdr.scl_inter;
		mx := (mx * lHdr.scl_slope)+ lHdr.scl_inter;
		intercept := 0;
		if (omn < 0) and (mn < 0) then //output is signed integer
			slope := max(abs(mx/omx),abs(mn/omn))
		else begin//output is UNSIGNED integer
			if mn < 0 then begin//some negative values: offset values
				intercept := mn;
				slope := (mx-mn)/omx;
			end else
				slope := mx/omx;
		end;
		if slope = 0 then slope := 1;	
		printf(format('scaled input range %g..%g integer output range %g..%g slope %g intercept %g', [mn,mx, omn, omx, slope, intercept]));
		lHdr.scl_slope := slope;
		lHdr.scl_inter := intercept;	
	end;
	islope := 1/slope;
	if (outDT = kDT_INT8) then begin
		i8s := TInt8s(rawData);
		for i := 0 to (n - 1) do
			i8s[i] := round((x32s[i]-intercept)*islope);
	end else if (outDT = kDT_INT16) then begin
		i16s := TInt16s(rawData);
		for i := 0 to (n - 1) do
			i16s[i] := round((x32s[i]-intercept)*islope);
	end else if (outDT = kDT_UINT8) then begin
		u8s := TUInt8s(rawData);
		for i := 0 to (n - 1) do
			u8s[i] := round((x32s[i]-intercept)*islope);
	end else if (outDT = kDT_UINT16) then begin
		u16s := TUInt16s(rawData);
		for i := 0 to (n - 1) do
			u16s[i] := round((x32s[i]-intercept)*islope);
	end else if (outDT = kDT_INT32) then begin
		i32s := TInt32s(rawData);
		for i := 0 to (n - 1) do
			i32s[i] := round((x32s[i]-intercept)*islope);
	end;
	setlength(x32s, 0);	
end;

procedure SwapImg(var rawData: TUInt8s; bitpix: integer);
var
   i16s: TInt16s;
   i32s: TInt32s;
   f64s: TFloat64s;
   i, n: int64;
begin
     if bitpix < 15 then exit;
     if bitpix = 16 then begin
        n := length(rawData) div 2;
        i16s := TInt16s(rawData);
        for i := 0 to (n-1) do
            i16s[i] := swap2(i16s[i]);
     end;
     if bitpix = 32 then begin
        n := length(rawData) div 4;
        i32s := TInt32s(rawData);
        for i := 0 to (n-1) do
            swap4(i32s[i]);
     end;
     if bitpix = 64 then begin
        n := length(rawData) div 8;
        f64s := TFloat64s(rawData);
        for i := 0 to (n-1) do
            Xswap8r(f64s[i]);
     end;
end;

function HdrVolumes(hdr: TNIfTIhdr): integer;
var
  i: integer;
begin
     result := 1;
     for i := 4 to 7 do
         if hdr.dim[i] > 1 then
            result := result * hdr.dim[i];
end;

function HdrVolBytes4D(hdr: TNIfTIhdr): integer;
begin
	result := hdr.Dim[1]*hdr.Dim[2]*hdr.Dim[3] * (hdr.bitpix div 8) *  HdrVolumes(hdr);
end;

{$IFDEF BZIP2}
function LoadImgBZ(FileName : AnsiString; var  rawData: TUInt8s; var lHdr: TNIFTIHdr; gzBytes: integer): boolean;
//foreign: both image and header compressed
label
	123;
var
  Decompressed: TDecompressBzip2Stream;
  InFile: TFileStream;
  i, volBytes, offset: int64;
begin
  result := false;
  volBytes := HdrVolBytes4D(lHdr);
  InFile:=TFileStream.Create(FileName, fmOpenRead);
  offset := round(lHdr.vox_offset);
  if (offset > 0) and (gzBytes = K_bz2Bytes_onlyImageCompressed) then
    Infile.Seek(offset, soFromBeginning);
  Decompressed:=TDecompressBzip2Stream.Create(InFile);
  if (offset > 0) and (gzBytes = K_bz2Bytes_headerAndImageCompressed) then begin
      SetLength (rawData, offset);
      i:=Decompressed.Read(Pointer(@rawData[0])^,volBytes);
      if i <> offset then begin
	      printf(format('BZip2 error: unable to skip header for %s', [FileName]));
	      goto 123;
      end;
  end;
  SetLength (rawData, volBytes);
  i:=Decompressed.Read(Pointer(@rawData[0])^,volBytes);
  result := (i = volBytes);
  if not result then
  	printf(format('BZip2 error: read %d but expected %d bytes for %s', [i, volBytes, FileName]));
  123:
  Decompressed.Free;
  InFile.Free;
end; 
{$ENDIF}//BZIP2

{$IFDEF GZIP}
function LoadHdrRawImgGZ(FileName : AnsiString; var  rawData: TUInt8s; var lHdr: TNIFTIHdr; gzBytes: integer): boolean;
var
   fStream: TFileStream;
   inStream, outStream: TMemoryStream;
   volBytes: int64;
   skip: int64 = 0;
label
     123;
begin
 result := false;
 if not fileexists(Filename) then exit;
 if (lHdr.bitpix <> 8) and (lHdr.bitpix <> 16) and (lHdr.bitpix <> 24) and (lHdr.bitpix <> 32) and (lHdr.bitpix <> 64) then begin
   printf('Unable to load '+Filename+' - this software can only read 8,16,24,32,64-bit image data.');
   exit;
 end;
 volBytes := HdrVolBytes4D(lHdr);
 outStream := TMemoryStream.Create();
 fStream := TFileStream.Create(Filename, fmOpenRead);
 fStream.seek(0, soFromBeginning);
 inStream := TMemoryStream.Create();
 if (gzBytes = K_gzBytes_onlyImageCompressed) then begin
 	fStream.seek(round(lHdr.vox_offset), soFromBeginning);
 	inStream.CopyFrom(fStream, fStream.Size - round(lHdr.vox_offset));
 end else begin
	inStream.CopyFrom(fStream, fStream.Size);
	skip := round(lHdr.vox_offset) 
 end;
 result := unzipStream(inStream, outStream);
 fStream.Free;
 inStream.Free;
 if (not result) and ((volBytes+skip) >=outStream.size) then begin
   printf('unzipStream error but sufficient bytes extracted (perhaps GZ without length in footer)');
   result := true;
 end;
 if not result then goto 123;
 if outStream.Size < (volBytes+skip) then begin
    result := false;
    printf(format('GZ error expected %d found %d bytes: %s',[volBytes+skip,outStream.Size, Filename]));
    goto 123;
 end;
 SetLength (rawData, volBytes);
 outStream.Position := skip; 
 outStream.ReadBuffer (rawData[0], volBytes);
 123:
   outStream.Free;
end;
{$ENDIF}//GZIP

procedure planar3D2RGB8(var  rawData: TUInt8s; var lHdr: TNIFTIHdr);
var
   img: TUInt8s;
   xy, xys, i, j, s, xyz, nBytesS, SamplesPerPixel: int64;
begin
  if lHdr.datatype <> kDT_RGBplanar3D then exit;
  lHdr.datatype := kDT_RGB;
  SamplesPerPixel := 3;
  xy := lHdr.dim[1] * lHdr.dim[2];
  xyz := xy * lHdr.dim[3];
  xys := xy *  SamplesPerPixel;
  nBytesS := xys * lHdr.dim[3] ;
  setlength(img, nBytesS);
  img := Copy(rawData, Low(rawData), Length(rawData));
  j := 0;
  for i := 0 to (xyz-1) do begin
      for s := 0 to (SamplesPerPixel-1) do begin
          rawData[j] := img[i+(s * xyz)] ;
          j := j + 1;
      end;
  end;
  img := nil;
end;

procedure DimPermute2341(var  rawData: TUInt8s; var lHdr: TNIFTIHdr);
//NIfTI demands first three dimensions are spatial, NRRD often makes first dimension non-spatial (e.g. DWI direction)
// This function converts NRRD TXYZ to NIfTI compatible XYZT
var
   i, x,y,z,t,xInc,yInc,zInc,tInc, nbytes: int64;
   i8, o8: TUint8s;
   i16, o16: TUInt16s;
   i32, o32: TUInt32s;
begin
     if HdrVolumes(lHdr) < 2 then exit;
     if (lHdr.bitpix mod 8) <> 0 then exit;
     if (lHdr.bitpix <> 8) and (lHdr.bitpix <> 16) and (lHdr.bitpix <> 32) then exit;
     nbytes := lHdr.Dim[1] * lHdr.Dim[2] * lHdr.Dim[3] * HdrVolumes(lHdr) * (lHdr.bitpix div 8);
     if nbytes < 4 then exit;
     setlength(i8, nbytes);
     i8 := copy(rawData, low(rawData), high(rawData));
     i16 := TUInt16s(i8);
     i32 := TUInt32s(i8);
     o8 := TUInt8s(rawData);
     o16 := TUInt16s(rawData);
     o32 := TUInt32s(rawData);
     t :=  lHdr.Dim[1];
     x :=  lHdr.Dim[2];
     y :=  lHdr.Dim[3];
     z :=  HdrVolumes(lHdr);
     lHdr.Dim[1] := x;
     lHdr.Dim[2] := y;
     lHdr.Dim[3] := z;
     lHdr.Dim[4] := t;
     lHdr.Dim[5] := 1;
     lHdr.Dim[6] := 1;
     lHdr.Dim[7] := 1;
     tInc := 1;
     xInc := t;
     yInc := t * x;
     zInc := t * x * y;
     i := 0;
     if (lHdr.bitpix = 8) then
        for t := 0 to (lHdr.Dim[4]-1) do
            for z := 0 to (lHdr.Dim[3]-1) do
                for y := 0 to (lHdr.Dim[2]-1) do
                    for x := 0 to (lHdr.Dim[1]-1) do begin
                        o8[i] := i8[(x*xInc)+(y*yInc)+(z*zInc)+(t*tInc)];
                        i := i + 1;
                    end;
     if (lHdr.bitpix = 16) then
        for t := 0 to (lHdr.Dim[4]-1) do
            for z := 0 to (lHdr.Dim[3]-1) do
                for y := 0 to (lHdr.Dim[2]-1) do
                    for x := 0 to (lHdr.Dim[1]-1) do begin
                        o16[i] := i16[(x*xInc)+(y*yInc)+(z*zInc)+(t*tInc)];
                        i := i + 1;
                    end;
     if (lHdr.bitpix = 32) then
        for t := 0 to (lHdr.Dim[4]-1) do
            for z := 0 to (lHdr.Dim[3]-1) do
                for y := 0 to (lHdr.Dim[2]-1) do
                    for x := 0 to (lHdr.Dim[1]-1) do begin
                        o32[i] := i32[(x*xInc)+(y*yInc)+(z*zInc)+(t*tInc)];
                        i := i + 1;
                    end;
     i8 := nil;
end;

function saveNii(fnm: string; var oHdr: TNIFTIhdr; var orawVolBytes: TUInt8s; isGz: boolean): boolean; overload;
var
  mStream : TMemoryStream;
  zStream: TGZFileStream;
  oPad32: Uint32; //nifti header is 348 bytes padded with 4
  lExt,NiftiOutName: string;
begin
 result := true;
 fnm :=  fnm + '.nii';
 if isGz then
 	fnm := fnm + '.gz';
 if fileexists(fnm) then
 	printf('Overwriting "'+fnm+'"')
 else
 	printf('Converted "'+fnm+'"');
 mStream := TMemoryStream.Create;
 oHdr.vox_offset :=  sizeof(oHdr) + 4;
 mStream.Write(oHdr,sizeof(oHdr));
 oPad32 := 4;
 mStream.Write(oPad32, 4);
 mStream.Write(oRawVolBytes[0], length(orawVolBytes));
 oRawVolBytes := nil;
 mStream.Position := 0;
 FileMode := fmOpenWrite;
 NiftiOutName := fnm;
 lExt := uppercase(extractfileext(NiftiOutName));
 if (lExt = '.GZ') or (lExt = '.VOI') then begin  //save gz compressed
    //if (lExt = '.GZ') then
    //   NiftiOutName := ChangeFileExt(NiftiOutName,'.nii.gz'); //img.gz -> img.nii.gz
    zStream := TGZFileStream.Create(NiftiOutName, gzopenwrite);
    zStream.CopyFrom(mStream, mStream.Size);
    zStream.Free;
 end else begin
     if (lExt <> '.NII') then
        NiftiOutName := NiftiOutName + '.nii';
     mStream.SaveToFile(NiftiOutName); //save uncompressed
 end;
 mStream.Free;
 FileMode := fmOpenRead;
end;

function saveNii(fnm: string; var lHdr: TNIFTIhdr; var rawData: TUInt8s; isGz, is3D: boolean): boolean; overload;
//handles saving 4D data as 3D volumes...
var
	v, nvol, pad: integer;
	bytesPerVol: int64;
	rawDataVol: TUInt8s;
begin
	 nvol := lHdr.dim[4];
	 if (nvol < 2) or (not is3D) then begin
	 	result := saveNii(fnm, lHdr, rawData, isGz);
	 	exit;
	 end; 
	 //only 4D->3D datasets follow...	
	 pad := trunc(log10(nvol))+1;
	 lHdr.dim[4] := 1;
	 lHdr.dim[0] := 3; //3D
	 bytesPerVol := length(rawData) div nvol;
	 if (bytesPerVol * nVol) <> length(rawData) then exit(false); //should never happen
	 for v := 0 to (nvol-1) do begin
		setlength(rawDataVol, bytesPerVol);
		rawDataVol := copy(rawData, v * bytesPerVol, bytesPerVol);
		result := saveNii(fnm+strutils.Dec2Numb(v+1,pad,10), lHdr, rawDataVol, isGz);	
	end; //for v: each volume
	rawData := nil;
end;

function readNiftiHdr (var lFilename: string; var lHdr: TNIFTIhdr; var gzBytes: int64; var swapEndian, isDimPermute2341: boolean): boolean; overload;
var
    lExt, lExt2GZ: string;
    Stream : TFileStream;
    zStream: TGZFileStream;
begin
    result := false;
    isDimPermute2341 := false;
    swapEndian := false;
    gzBytes := 0;
    lExt := upcase(ExtractFileExt(lFilename));
    lExt2GZ := '';
    if (lExt = '.GZ') then begin
        lExt2GZ := changefileext(lFilename,'');
        lExt2GZ := upcase(ExtractFileExt(lExt2GZ));
    end;
    if (lExt2GZ <> '.NII') and (lExt <> '.NII') then exit;
    lHdr.HdrSz := 0;
    if (lExt = '.GZ') then begin
        zStream := TGZFileStream.Create (lFilename, gzopenread);
        try
            zStream.ReadBuffer (lHdr, SizeOf (TNIFTIHdr));
        finally
            zStream.Free;
        end;
        gzBytes := K_gzBytes_headerAndImageCompressed;
    end else begin
        Stream := TFileStream.Create (lFilename, fmOpenRead or fmShareDenyWrite);
        try
            Stream.ReadBuffer (lHdr, SizeOf (TNIFTIHdr));
        finally
            Stream.Free;
        end;
    end;
    if lHdr.HdrSz <> SizeOf (TNIFTIHdr) then begin
        NIFTIhdr_SwapBytes(lHdr);
        swapEndian := true;
    end;
    if (lHdr.scl_slope = 0.0) then
    	lHdr.scl_slope := 1.0;
    result := lHdr.HdrSz = SizeOf (TNIFTIHdr) 
end;

function loadVolumes(var FileName: string; out lHdr: TNIFTIhdr; out rawData: TUInt8s; out isInputNIfTI: boolean): boolean;
var
  ifnm: string;
  Stream : TFileStream;
  gzBytes, volBytes, FSz: int64;
  ok, swapEndian, isDimPermute2341: boolean;
begin
     result := false;
     ifnm := FileName;
     if not fileexists(FileName) then begin
        FileName := GetCurrentDir+PathDelim+ FileName;
        if not fileexists(FileName) then begin
           printf('Unable to find "'+ifnm+'"');
           exit;
        end;
     end;
     //expand filenames, optional but makes output easier to parse 'Converted "./test/nrrd.raw.nii"'
     ifnm := ExpandFileName(ifnm);
     //set output directory
     ok := readNiftiHdr(FileName, lHdr, gzBytes, swapEndian, isDimPermute2341);
     isInputNIfTI := ok;
     if not ok then 
        ok := readForeignHeader(FileName, lHdr, gzBytes, swapEndian, isDimPermute2341);
     if not ok then begin
        printf('Unable to interpret header "'+FileName+'"');
        exit;
     end;
     //printf('DEBUGMODE>>>read but do not convert images');  exit;
     if not fileexists(FileName) then begin
        printf('Unable to find image data "'+FileName+'" described by header "'+ifnm+'"');
        exit;
     end;
     {$IFDEF BZIP2}
     if ((gzBytes = K_bz2Bytes_headerAndImageCompressed) or (gzBytes = K_bz2Bytes_onlyImageCompressed)) then
     	result := LoadImgBZ(FileName,  rawData, lHdr, gzBytes)
     {$ELSE}
     if (gzBytes = K_bz2Bytes_headerAndImageCompressed) then begin
        printf('Not compiled to read BZip2 files: '+FileName);
        exit;
     end     
     {$ENDIF}      
     {$IFDEF GZIP}
     else if (gzBytes = K_gzBytes_headerAndImageCompressed) or (gzBytes = K_gzBytes_onlyImageCompressed) then
       result := LoadHdrRawImgGZ(FileName,  rawData, lHdr, gzBytes)
     {$ELSE}
     else if (gzBytes = K_gzBytes_onlyImageCompressed) or (gzBytes < 0) then begin
        printf('Not compiled to read GZip files: '+FileName);
        exit;
     end
     {$ENDIF}
     else  begin //not compressed: read raw
       volBytes := HdrVolBytes4D(lHdr);
       FSz := FSize(FileName);
       if (FSz < (round(lHdr.vox_offset)+volBytes)) then begin
        printf(format('Unable to load '+Filename+' - file size (%d) smaller than expected (%d) ', [FSz, round(lHdr.vox_offset)+volBytes]));
        exit(false);
       end;
       Stream := TFileStream.Create (FileName, fmOpenRead or fmShareDenyWrite);
       try
        Stream.Seek(round(lHdr.vox_offset),soFromBeginning);
        SetLength (rawData, volBytes);
        Stream.ReadBuffer (rawData[0], volBytes);
       finally
        Stream.Free;
       end;
       result := true;
     end; //read raw
     if not result then begin
        printf('Unable to read image data for "'+FileName+'"');
       exit;
     end;
     if swapEndian then
     	SwapImg(rawData, lHdr.bitpix);
     planar3D2RGB8(rawData, lHdr);
     if isDimPermute2341 then
        DimPermute2341(rawData, lHdr);
     result := true;
end;

end.
