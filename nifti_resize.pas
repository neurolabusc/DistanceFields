unit nifti_resize;

{$mode delphi}{$H+}
{$inline on}
{$DEFINE MYTHREADS} //multi-threaded
{$IFDEF MYTHREADS}{$ModeSwitch nestedprocvars}{$ENDIF}

interface

uses
{$IFDEF MYTHREADS}mtprocs,mtpcpu,{$ENDIF}
  Classes, SysUtils, nifti_types, SimdUtils, VectorMath;

function ShrinkOrEnlarge(var lHdr: TNIFTIhdr; var lBuffer: TUInt8s; lScale: single; MaxThreads: PtrInt = 0; filterIndex: integer = 0): boolean; overload;
function ShrinkOrEnlarge(var lHdr: TNIFTIhdr; var lBuffer: TUInt8s; var lScale: TVec3; MaxThreads: PtrInt = 0; filterIndex: integer = 0): boolean; overload;//
function ShrinkMax(var hdr: TNIFTIhdr; var img: TUInt8s; scale: integer; maxthreads: integer = 0): boolean;

implementation

uses math;

{$IFDEF MYTHREADS}
function ShrinkMax(var hdr: TNIFTIhdr; var img: TUInt8s; scale: integer; maxthreads: integer = 0): boolean;
//integer reduction, for each output voxel reports peak of all contributing input voxels
var
	tPerThread, tHi, inX, outX, outXY, inY, inXY, outVox: int64;
	inImg, outImg: TFloat32s;	
procedure SubProc(ThreadIndex: PtrInt; Data: Pointer; Item: TMultiThreadProcItem);
var
	tStart, tEnd, i, inSlice, outSlice, inRow, outPos, inPos: int64;	
begin
	tStart := (ThreadIndex * tPerThread);
	if (tStart > tHi) then exit; //more threads than slices in Z direction
	tEnd := tStart + tPerThread - 1; //e.g. if zStart=4 and zPerThread=1 then zEnd=4 
	tEnd := min(tEnd, tHi); //final thread when slices in Z not evenly divisible by number of threads
	for inSlice := tStart to tEnd do begin
		outSlice := (inSlice div scale) * outXY;
		for inRow := 0 to (inY-1) do begin
			inPos := (inSlice * inXY) + (inRow * inX);
			outPos := outSlice + (inRow div scale) * outX;
			for i := 1 to (inX) do begin
				outImg[outPos] := max(inImg[inPos], outImg[outPos]);
				inPos += 1;
				if (i mod scale) = 0 then
					outPos += 1;
			end;
		end;
	end;
end; //nested SubProc()	
var
	outImg8: TUInt8s;
	i: int64;
begin
	result := false;
	if (scale < 2) then exit;
	if ((hdr.dim[1] mod scale) <> 0) then exit;
	if ((hdr.dim[2] mod scale) <> 0) then exit;
	if ((hdr.dim[3] mod scale) <> 0) then exit;
	if (MaxThreads < 1) then MaxThreads :=  GetSystemThreadCount();
	tHi := (hdr.dim[3]*max(hdr.dim[4], 1)) -1;
	tPerThread := ceil((tHi+1)/MaxThreads);
	inX := hdr.dim[1];
	inY := hdr.dim[2];
	inXY := inX * inY;
	//inZT := 
	//inVox := inX * inYZT;
	hdr.dim[1] := hdr.dim[1] div scale;
	hdr.dim[2] := hdr.dim[2] div scale;
	hdr.dim[3] := hdr.dim[3] div scale;
	hdr.pixdim[1] := hdr.pixdim[1] * scale;
	hdr.pixdim[2] := hdr.pixdim[2] * scale;
	hdr.pixdim[3] := hdr.pixdim[3] * scale;
	outX := hdr.dim[1];
	outXY := outX * hdr.dim[2];
	outVox := hdr.dim[1] * hdr.dim[2] * hdr.dim[3]*max(hdr.dim[4], 1);
	inImg := TFloat32s(img);
	setlength(outImg8, outVox * sizeof(single));
	outImg := TFloat32s(outImg8);
	for i := 0 to (outVox-1) do
		outImg[i] := -infinity;
	ProcThreadPool.DoParallelNested(SubProc,0,MaxThreads-1, nil, MaxThreads);
	img := nil; //free input image
	img := outImg8; //return output image
	for i := 0 to 2 do begin
		hdr.srow_x[i] := hdr.srow_x[i] * scale;
		hdr.srow_y[i] := hdr.srow_y[i] * scale;
		hdr.srow_z[i] := hdr.srow_z[i] * scale;
		
	end;
	hdr.qform_code := kNIFTI_XFORM_UNKNOWN;
	result := true;
end;
{$ELSE}
function ShrinkMax(var hdr: TNIFTIhdr; var img: TUInt8s; scale: integer; maxthreads: integer = 0): boolean;
//integer reduction, for each output voxel reports peak of all contributing input voxels
var
	inX, outX, outXY, inY, inXY, outVox: int64;
	inImg, outImg: TFloat32s;	
	i, inSlice, outSlice, inRow, outPos, inPos: int64;	
	outImg8: TUInt8s;
begin
	result := false;
	if (scale < 2) then exit;
	if ((hdr.dim[1] mod scale) <> 0) then exit;
	if ((hdr.dim[2] mod scale) <> 0) then exit;
	if ((hdr.dim[3] mod scale) <> 0) then exit;
	inX := hdr.dim[1];
	inY := hdr.dim[2];
	inXY := inX * inY;
	//inZT := 
	//inVox := inX * inYZT;
	hdr.dim[1] := hdr.dim[1] div scale;
	hdr.dim[2] := hdr.dim[2] div scale;
	hdr.dim[3] := hdr.dim[3] div scale;
	hdr.pixdim[1] := hdr.pixdim[1] * scale;
	hdr.pixdim[2] := hdr.pixdim[2] * scale;
	hdr.pixdim[3] := hdr.pixdim[3] * scale;
	outX := hdr.dim[1];
	outXY := outX * hdr.dim[2];
	outVox := hdr.dim[1] * hdr.dim[2] * hdr.dim[3]*max(hdr.dim[4], 1);
	inImg := TFloat32s(img);
	setlength(outImg8, outVox * sizeof(single));
	outImg := TFloat32s(outImg8);
	for i := 0 to (outVox-1) do
		outImg[i] := -infinity;
	for inSlice := 0 to ((hdr.dim[3]*max(hdr.dim[4], 1)) -1) do begin
		outSlice := (inSlice div scale) * outXY;
		for inRow := 0 to (inY-1) do begin
			inPos := (inSlice * inXY) + (inRow * inX);
			outPos := outSlice + (inRow div scale) * outX;
			for i := 1 to (inX) do begin
				outImg[outPos] := max(inImg[inPos], outImg[outPos]);
				inPos += 1;
				if (i mod scale) = 0 then
					outPos += 1;
			end;
		end;
	end;
	img := nil; //free input image
	img := outImg8; //return output image
	for i := 0 to 2 do begin
		hdr.srow_x[i] := hdr.srow_x[i] * scale;
		hdr.srow_y[i] := hdr.srow_y[i] * scale;
		hdr.srow_z[i] := hdr.srow_z[i] * scale;
		
	end;
	hdr.qform_code := kNIFTI_XFORM_UNKNOWN;
	result := true;
end;
{$ENDIF}

// Extends image shrink code by Anders Melander, anders@melander.dk
// Here's some additional copyrights for you:
//
// The algorithms and methods used in this library are based on the article
// "General Filtered Image Rescaling" by Dale Schumacher which appeared in the
// book Graphics Gems III, published by Academic Press, Inc.
// From filter.c:
// The authors and the publisher hold no copyright restrictions
// on any of these files; this source code is public domain, and
// is freely available to the entire computer graphics community
// for study, use, and modification.  We do request that the
// comment at the top of each file, identifying the original
// author and its original publication in the book Graphics
// Gems, be retained in all programs that use these files.

// Hermite filter

function HermiteFilter(Value: Single): Single;
begin
  // f(t) = 2|t|^3 - 3|t|^2 + 1, -1 <= t <= 1
  if (Value < 0.0) then
    Value := -Value;
  if (Value < 1.0) then
    Result := (2.0 * Value - 3.0) * Sqr(Value) + 1.0
  else
    Result := 0.0;
end;

// Triangle filter
// a.k.a. "Linear" or "Bilinear" filter

function TriangleFilter(Value: Single): Single;
begin
  if (Value < 0.0) then
    Value := -Value;
  if (Value < 1.0) then
    Result := 1.0 - Value
  else
    Result := 0.0;
end;

// Bell filter

function BellFilter(Value: Single): Single;
begin
  if (Value < 0.0) then
    Value := -Value;
  if (Value < 0.5) then
    Result := 0.75 - Sqr(Value)
  else if (Value < 1.5) then
  begin
    Value := Value - 1.5;
    Result := 0.5 * Sqr(Value);
  end
  else
    Result := 0.0;
end;

// B-spline filter

function SplineFilter(Value: Single): Single;
var
  tt: single;
begin
  if (Value < 0.0) then
    Value := -Value;
  if (Value < 1.0) then
  begin
    tt := Sqr(Value);
    Result := 0.5 * tt * Value - tt + 2.0 / 3.0;
  end
  else if (Value < 2.0) then
  begin
    Value := 2.0 - Value;
    Result := 1.0 / 6.0 * Sqr(Value) * Value;
  end
  else
    Result := 0.0;
end;

// Lanczos3 filter

function Lanczos3Filter(Value: Single): Single;
function SinC(Value: Single): Single;
  begin
    if (Value <> 0.0) then
    begin
      Value := Value * Pi;
      Result := sin(Value) / Value
    end
    else
      Result := 1.0;
  end;
begin
  if (Value < 0.0) then
    Value := -Value;
  if (Value < 3.0) then
    Result := SinC(Value) * SinC(Value / 3.0)
  else
    Result := 0.0;
end;


function MitchellFilter(Value: Single): Single;
const
  B = (1.0 / 3.0);
  C = (1.0 / 3.0);
var
  tt: single;
begin
  if (Value < 0.0) then
    Value := -Value;
  tt := Sqr(Value);
  if (Value < 1.0) then
  begin
    Value := (((12.0 - 9.0 * B - 6.0 * C) * (Value * tt))
      + ((-18.0 + 12.0 * B + 6.0 * C) * tt)
      + (6.0 - 2 * B));
    Result := Value / 6.0;
  end
  else if (Value < 2.0) then
  begin
    Value := (((-1.0 * B - 6.0 * C) * (Value * tt))
      + ((6.0 * B + 30.0 * C) * tt)
      + ((-12.0 * B - 48.0 * C) * Value)
      + (8.0 * B + 24 * C));
    Result := Value / 6.0;
  end
  else
    Result := 0.0;
end;

type
  // Contributor for a pixel
  TFilterProc = function(Value: Single): Single;
  TContributor = record
    pixel: integer; // Source pixel
    weight: single; // Pixel weight
  end;
  TContributorList = array[0..0] of TContributor;
  PContributorList = ^TContributorList;
  // List of source pixels contributing to a destination pixel
  TCList = record
    n: integer;
    p: PContributorList;
  end;
  TCListList = array[0..0] of TCList;
  PCListList = ^TCListList;

procedure SetContrib(out contrib: PCListList; SrcPix, DstPix, Delta: integer; xscale, fwidth: single; filter: TFilterProc);
var
  i,j,k: int64;
  width, fscale: single;
  sum, center, weight: single; // Filter calculation variables
  left, right: integer; // Filter calculation variables
begin
  if (DstPix < 1) or (xscale <= 0) then exit;
  if (xscale < 1)and (fwidth > 0.6) then //not for nearest neighbor!
  	fscale := 1.0 / xscale
  else
  	fscale := 1.0;
  width := fwidth * fscale;
  GetMem(contrib, DstPix * sizeof(TCList));
  for i := 0 to DstPix - 1 do begin
      contrib^[i].n := 0;
      GetMem(contrib^[i].p, trunc(width * 2.0 + 1) * sizeof(TContributor));
      center := i / xscale;
      left := floor(center - width);
      left := max(left,0);
      right := ceil(center + width);
      right := min(right, SrcPix - 1);
      sum := 0.0;
      for j := left to right do begin
        weight := filter((center - j) / fscale) / fscale;
        if (weight = 0.0) then
          continue;
        sum := sum + weight;
        k := contrib^[i].n;
        contrib^[i].n := contrib^[i].n + 1;
        contrib^[i].p^[k].pixel := j * Delta;
        contrib^[i].p^[k].weight := weight;
      end;
      for k := 0 to contrib^[i].n - 1 do
          contrib^[i].p^[k].weight := contrib^[i].p^[k].weight/sum;
      (*showmessage(format('n=%d l=%d r=%d c=%g sum=%g',[contrib^[i].n, left, right, center, sum]));
      for k := 0 to contrib^[i].n - 1 do
          showmessage(format('%d %g',[contrib^[i].p^[k].pixel, contrib^[i].p^[k].weight])); *)
    end;
end;

procedure Zoom(var lHdr: TNIFTIhdr; xScale, yScale, zScale: single);
//if we have a 256x256x256 pixel image with scale of 0.5, output is 128x128x128
//if we have a 1x1x1mm pixel image with a scale of 2.0, output is 2x2x2mm
var
   i, inDim: int64;
   scale: array[1..3] of single;
begin
     for i := 1 to 3 do begin
         if i = 1 then
                 scale[i] := xScale
         else if i = 2 then
                 scale[i] := yScale
         else
             scale[i] := zScale;
         if (round(lHdr.dim[i] * scale[i]) < 1) then begin
              scale[i] := 1/ lHdr.dim[i]  //e.g. for reducing 2D images, Z dimension does not change
         end;
         inDim := lHdr.dim[i];
         lHdr.dim[i] := round(lHdr.dim[i] * scale[i]);
         scale[i] := lHdr.dim[i] / inDim; //e.g. rounding error
         lHdr.pixdim[i] := lHdr.pixdim[i] / scale[i];
         //fx(lHdr.srow_x[i] ,lHdr.srow_y[i] ,lHdr.srow_z[i] );
     end;
     for i :=0 to 2 do begin
         lHdr.srow_x[i] := lHdr.srow_x[i]/ scale[i+1];
         lHdr.srow_y[i] := lHdr.srow_y[i]/ scale[i+1];
         lHdr.srow_z[i] := lHdr.srow_z[i]/ scale[i+1];
     end;
     lHdr.qform_code := kNIFTI_XFORM_UNKNOWN;
end;

{$DEFINE LESS_RAM}
{$IFDEF LESS_RAM}
procedure Resize32(var lHdr: TNIFTIhdr; var lImg8: TUInt8s; xScale, yScale, zScale, fwidth: single; filter: TFilterProc; MaxThreads: PtrInt = 0);
//rescales images with any dimension larger than lMaxDim to have a maximum dimension of maxdim...
label
  666;
var
	lXo,lYo,lZo,lXi,lYi,lZi, linePerThread, lineMax: int64;
  	contrib: PCListList;
  	finalImg, img8x, img8y, img8z: TUInt8s;
  	inImg, tempImgX, tempImgY, tempImgZ, out32: TFloat32s;
{$IFDEF MYTHREADS}
procedure SubProcX(ThreadIndex: PtrInt; Data: Pointer; Item: TMultiThreadProcItem);
{$ELSE}
procedure SubProcX(ThreadIndex: PtrInt); inline;
{$ENDIF}
var
	lineStart, lineStartOut, lineEnd, line, x, j : int64;
	sum: double;
begin
	lineStart := (ThreadIndex * linePerThread);
	if (lineStart > lineMax) then exit; //more threads than slices in Z direction
	lineEnd := lineStart + linePerThread - 1; //e.g. if zStart=4 and zPerThread=1 then zEnd=4 
	lineEnd := min(lineEnd, lineMax); //final thread when slices in Z not evenly divisible by number of threads
	for line := lineStart to lineEnd do begin
		lineStart :=  (line*lXi);
		lineStartOut := (line*lXo);
		for x := 0 to (lXo - 1) do begin
			sum := 0.0;
			for j := 0 to contrib^[x].n - 1 do begin
			  sum := sum + (contrib^[x].p^[j].weight * inImg[lineStart +contrib^[x].p^[j].pixel]);
			end;
			tempImgX[lineStartOut+x] := sum;
		end; //for X
	end; //for Z
end; //SubProcX()
{$IFDEF MYTHREADS}
procedure SubProcY(ThreadIndex: PtrInt; Data: Pointer; Item: TMultiThreadProcItem);
{$ELSE}
procedure SubProcY(ThreadIndex: PtrInt); inline;
{$ENDIF}
var
	lineStart, lineEnd, x, y, z, j, i : int64;
	sum: double;
begin
	lineStart := (ThreadIndex * linePerThread);
	if (lineStart > lineMax) then exit; //more threads than slices in Z direction
	lineEnd := lineStart + linePerThread - 1; //e.g. if zStart=4 and zPerThread=1 then zEnd=4 
	lineEnd := min(lineEnd, lineMax); //final thread when slices in Z not evenly divisible by number of threads
	//i := lineStart * lXo * lYo;
	i := lineStart * lYo * lXo;
	for z := lineStart to lineEnd do begin
	  for y := 0 to (lYo - 1) do begin
		  for x := 0 to (lXo-1) do begin
		  	lineStart :=  x+((lXo*lYi) * z);
		  	sum := 0.0;
		  	for j := 0 to contrib^[y].n - 1 do
		  		sum := sum + (contrib^[y].p^[j].weight * tempImgX[lineStart +contrib^[y].p^[j].pixel] );
			tempImgY[i] := sum;
			i := i + 1;
		  end; //for X
		end; //for Y
	end; //for Z
end; //SubProcY()
  
{$IFDEF MYTHREADS}
procedure SubProcZ(ThreadIndex: PtrInt; Data: Pointer; Item: TMultiThreadProcItem);
{$ELSE}
procedure SubProcZ(ThreadIndex: PtrInt); inline;
{$ENDIF}
var
	lineStart, lineEnd, x, y, z, j, i : int64;
	sum: double;
begin
	lineStart := (ThreadIndex * linePerThread);
	if (lineStart > lineMax) then exit; //more threads than slices in Z direction
	lineEnd := lineStart + linePerThread - 1; //e.g. if zStart=4 and zPerThread=1 then zEnd=4 
	lineEnd := min(lineEnd, lineMax); //final thread when slices in Z not evenly divisible by number of threads
	i := lineStart * lXo * lYo;
	for z := lineStart to lineEnd do begin
      for y := 0 to (lYo - 1) do begin
          for x := 0 to (lXo-1) do begin
            lineStart :=  x+(lXo * y);
            sum := 0.0;
            for j := 0 to contrib^[z].n - 1 do begin
              sum := sum + (contrib^[z].p^[j].weight * tempImgY[lineStart +contrib^[z].p^[j].pixel] );
            end;
            tempImgZ[i] := sum;
            i := i + 1;
        end; //for X
    end; //for Y
  end; //for Z
end; //SubProcZ()  
var 
  mx, mn: double;
  i: int64;
begin
  {$IFDEF MYTHREADS}
  if (MaxThreads < 1) then MaxThreads := GetSystemThreadCount();
  {$ELSE}
  MaxThreads := 1;
  {$ENDIF}
  lXi := lHdr.dim[1]; //input X
  lYi := lHdr.dim[2]; //input Y
  lZi := lHdr.dim[3]; //input Z
  lZi := max(lZi, 1); //e.g. 2D image might claim dim[3]=0
  lXo := lXi; lYo := lYi; lZo := lZi; //output initially same as input
  //inBytes := lHdr.dim[1]*lHdr.dim[2]*lHdr.dim[3]*bytesPerVox;
  inImg := TFloat32s(lImg8);
  //find min/max values
  mn := inImg[0];
  mx := mn;
  for i := 0 to (lHdr.dim[1]*lHdr.dim[2]*lHdr.dim[3])-1 do begin
      if inImg[i] < mn then mn := inImg[i];
      if inImg[i] > mx then mx := inImg[i];
  end;
  Zoom(lHdr,xScale, yScale, zScale);
  //shrink in 1st dimension : do X as these are contiguous = faster, compute slower dimensions at reduced resolution
  lXo := lHdr.dim[1]; //input X
  setlength(img8x,lXo*lYi*lZi*4);
  tempImgX := TFloat32s(img8x);
  SetContrib(contrib, lXi, lXo, 1, xScale, fwidth, filter);
  //parallel
  lineMax := (lYi * lZi)-1;
  linePerThread := ceil((lineMax+1)/MaxThreads);
  {$IFDEF MYTHREADS}
  ProcThreadPool.DoParallelNested(SubProcX,0,MaxThreads-1, nil, MaxThreads);
  {$ELSE}
  SubProcX(0);
  {$ENDIF}
  for i := 0 to lXo - 1 do
     FreeMem(contrib^[i].p);
  FreeMem(contrib);
  setlength( lImg8, 0);
  finalImg := img8x;
  out32 := tempImgX;
  if ((lYi = lHdr.dim[2]) and (lZi = lHdr.dim[3])) then goto 666; //e.g. 1D image
  //shrink in 2nd dimension
  lYo := lHdr.dim[2]; //reduce Y output
  setlength(img8y,lXo*lYo*lZi*4);
  tempImgY := TFloat32s(img8y);
  //SetLength( tempImgY,lXo*lYo*lZi); //8
  SetContrib(contrib, lYi, lYo, lXo, yScale, fwidth, filter);
  //parallel
  lineMax := lZi-1;
  linePerThread := ceil((lineMax+1)/MaxThreads);
  {$IFDEF MYTHREADS}
  ProcThreadPool.DoParallelNested(SubProcY,0,MaxThreads-1, nil, MaxThreads);
  {$ELSE}
  SubProcY(0);
  {$ENDIF}
  for i := 0 to lYo - 1 do
     FreeMem(contrib^[i].p);
  FreeMem(contrib);
  SetLength( img8x,0);
  finalImg := img8y;
  out32 := tempImgY;
  if (lZi = lHdr.dim[3]) then goto 666; //e.g. 2D image
  //shrink the 3rd dimension
  lZo := lHdr.dim[3]; //reduce Z output
  setlength(img8z,lXo*lYo*lZo*4);
  tempImgZ := TFloat32s(img8z);
  //SetLength( tempImgZ,lXo*lYo*lZo); //8
  SetContrib(contrib, lZi, lZo, (lXo*lYo), zScale, fwidth, filter);
  //parallel
  lineMax := lZo-1;
  linePerThread := ceil((lineMax+1)/MaxThreads);
  {$IFDEF MYTHREADS}
  ProcThreadPool.DoParallelNested(SubProcZ,0,MaxThreads-1, nil, MaxThreads);
  {$ELSE}
  SubProcZ(0);
  {$ENDIF}
  for i := 0 to lZo - 1 do
     FreeMem(contrib^[i].p);
  FreeMem(contrib);
  Setlength( img8y,0);
  finalImg := img8z;
  out32 := tempImgZ;
666:
  lHdr.dim[1] := lXo;
  lHdr.dim[2] := lYo;
  lHdr.dim[3] := lZo;
  for i := 0 to ((lXo*lYo*lZo)-1) do begin
      //check image range - some interpolation can cause ringing
      // e.g. if input range 0..1000 do not create negative values!
      if out32[i] > mx then out32[i] := mx;
      if out32[i] < mn then out32[i] := mn;
  end;
  lImg8 := finalImg;
end; //Resize32()
{$ELSE} //not LESS_RAM
procedure FindRange(var hdr: TNIFTIhdr; var img: TFloat32s; out mn, mx: single);
var
	i, n: integer;
begin
	n := hdr.dim[1] * hdr.dim[2] * hdr.dim[3];
	if n < 1 then exit;
	mn := img[0];
	mx := mn;
	for i := 0 to (n-1) do begin
		if (img[i] < mn) then mn := img[i];
		if (img[i] > mx) then mx := img[i];		
	end;	
end;

procedure CropRange(var hdr: TNIFTIhdr; var img: TFloat32s; mn, mx: single);
var
	i, n: integer;
begin
	n := hdr.dim[1] * hdr.dim[2] * hdr.dim[3];
	if n < 1 then exit;
	for i := 0 to (n-1) do begin
		if (img[i] < mn) then img[i] := mn;
		if (img[i] > mx) then img[i] := mx;		
	end;	
end;

procedure ResizeX(var lHdr: TNIFTIhdr; var lImg8: TUInt8s; xScale, fwidth: single; filter: TFilterProc; MaxThreads: PtrInt = 0);
//expand/shrink image in 1st dimension
var
	i, lXo, lXi, lYi, lZi, linePerThread, lineMax: int64;
  	contrib: PCListList;
  	inImg, outImg: TFloat32s; 
{$IFDEF MYTHREADS}
procedure SubProcX(ThreadIndex: PtrInt; Data: Pointer; Item: TMultiThreadProcItem);
{$ELSE}
procedure SubProcX(ThreadIndex: PtrInt); inline;
{$ENDIF}
var
	lineStart, lineStartOut, lineEnd, line, x, j : int64;
	sum: double;
begin
	lineStart := (ThreadIndex * linePerThread);
	if (lineStart > lineMax) then exit; //more threads than slices in Z direction
	lineEnd := lineStart + linePerThread - 1; //e.g. if zStart=4 and zPerThread=1 then zEnd=4 
	lineEnd := min(lineEnd, lineMax); //final thread when slices in Z not evenly divisible by number of threads
	for line := lineStart to lineEnd do begin
		lineStart :=  (line*lXi);
		lineStartOut := (line*lXo);
		for x := 0 to (lXo - 1) do begin
			sum := 0.0;
			for j := 0 to contrib^[x].n - 1 do begin
			  sum := sum + (contrib^[x].p^[j].weight * inImg[lineStart +contrib^[x].p^[j].pixel]);
			end;
			outImg[lineStartOut+x] := sum;
		end; //for X
	end; //for Z
end; //SubProcX()
begin
	if (xScale = 1) then exit(); //nothing to do
	{$IFDEF MYTHREADS}
	if (MaxThreads < 1) then MaxThreads := GetSystemThreadCount();
	{$ELSE}
	MaxThreads := 1;
	{$ENDIF}
	lXi := lHdr.dim[1]; //input X
	lYi := lHdr.dim[2]; //input X
	lZi := lHdr.dim[3]; //input X
	Zoom(lHdr,xScale, 1, 1);
	//shrink in 1st dimension : do X as these are contiguous = faster, compute slower dimensions at reduced resolution
	lXo := lHdr.dim[1]; //output X	
	SetContrib(contrib, lXi, lXo, 1, xScale, fwidth, filter);
	inImg := copy(TFloat32s(lImg8), 0, lXi * lYi * lZi);
	setlength(lImg8,lXo*lYi*lZi*sizeof(single));
	outImg := TFloat32s(lImg8);
	lineMax := (lYi * lZi)-1;
	linePerThread := ceil((lineMax+1)/MaxThreads);
	{$IFDEF MYTHREADS}
	ProcThreadPool.DoParallelNested(SubProcX,0,MaxThreads-1, nil, MaxThreads);
	{$ELSE}
	SubProcX(0);
	{$ENDIF}
	for i := 0 to lXo - 1 do
    	FreeMem(contrib^[i].p);
    FreeMem(contrib);
    //writeln(format('x %d %d %d -> %d',[lXi,lYi,lZi, lXo]));
end; //ResizeX()

procedure ResizeY(var lHdr: TNIFTIhdr; var lImg8: TUInt8s; yScale, fwidth: single; filter: TFilterProc; MaxThreads: PtrInt = 0);
//expand/shrink image in 1st dimension
var
	i, lYo, lXi, lYi, lZi, linePerThread, lineMax: int64;
  	contrib: PCListList;
  	inImg, outImg: TFloat32s; 
{$IFDEF MYTHREADS}
procedure SubProcY(ThreadIndex: PtrInt; Data: Pointer; Item: TMultiThreadProcItem);
{$ELSE}
procedure SubProcY(ThreadIndex: PtrInt); inline;
{$ENDIF}
var
	lineStart, lineEnd, x, y, z, j, i : int64;
	sum: double;
begin
	lineStart := (ThreadIndex * linePerThread);
	if (lineStart > lineMax) then exit; //more threads than slices in Z direction
	lineEnd := lineStart + linePerThread - 1; //e.g. if zStart=4 and zPerThread=1 then zEnd=4 
	lineEnd := min(lineEnd, lineMax); //final thread when slices in Z not evenly divisible by number of threads
	//i := lineStart * lXo * lYo;
	i := lineStart * lYo * lXi;
	for z := lineStart to lineEnd do begin
	  for y := 0 to (lYo - 1) do begin
		  for x := 0 to (lXi-1) do begin
		  	lineStart :=  x+((lXi*lYi) * z);
		  	sum := 0.0;
		  	for j := 0 to contrib^[y].n - 1 do
		  		sum := sum + (contrib^[y].p^[j].weight * inImg[lineStart +contrib^[y].p^[j].pixel] );
			outImg[i] := sum;
			//outImg[i] := random();
			i := i + 1;
		  end; //for X
		end; //for Y
	end; //for Z
end; //SubProcY()
begin
	if (yScale = 1) then exit(); //nothing to do
	{$IFDEF MYTHREADS}
	if (MaxThreads < 1) then MaxThreads := GetSystemThreadCount();
	{$ELSE}
	MaxThreads := 1;
	{$ENDIF}
	lXi := lHdr.dim[1]; //input X
	lYi := lHdr.dim[2]; //input X
	lZi := lHdr.dim[3]; //input X
	Zoom(lHdr,1, yScale, 1);
	lYo := lHdr.dim[2]; //output Y	
	SetContrib(contrib, lYi, lYo, lXi, yScale, fwidth, filter);
  	inImg := copy(TFloat32s(lImg8), 0, lXi * lYi * lZi);
	setlength(lImg8,lXi*lYo*lZi*sizeof(single));
	outImg := TFloat32s(lImg8);
	lineMax := lZi-1;
	linePerThread := ceil((lineMax+1)/MaxThreads);
	{$IFDEF MYTHREADS}
	ProcThreadPool.DoParallelNested(SubProcY,0,MaxThreads-1, nil, MaxThreads);
	{$ELSE}
	SubProcY(0);
	{$ENDIF}
	for i := 0 to lYo - 1 do
    	FreeMem(contrib^[i].p);
    FreeMem(contrib);
    //writeln(format('y d %d %d %d -> %d',[lXi,lYi,lZi, lYo]));
end; //ResizeY()

procedure ResizeZ(var lHdr: TNIFTIhdr; var lImg8: TUInt8s; zScale, fwidth: single; filter: TFilterProc; MaxThreads: PtrInt = 0);
//expand/shrink image in 1st dimension
var
	i, lZo, lXi, lYi, lZi, linePerThread, lineMax: int64;
  	contrib: PCListList;
  	inImg, outImg: TFloat32s; 
{$IFDEF MYTHREADS}
procedure SubProcZ(ThreadIndex: PtrInt; Data: Pointer; Item: TMultiThreadProcItem);
{$ELSE}
procedure SubProcZ(ThreadIndex: PtrInt); inline;
{$ENDIF}
var
	lineStart, lineEnd, x, y, z, j, i : int64;
	sum: double;
begin
	lineStart := (ThreadIndex * linePerThread);
	if (lineStart > lineMax) then exit; //more threads than slices in Z direction
	lineEnd := lineStart + linePerThread - 1; //e.g. if zStart=4 and zPerThread=1 then zEnd=4 
	lineEnd := min(lineEnd, lineMax); //final thread when slices in Z not evenly divisible by number of threads
	i := lineStart * lXi * lYi;
	for z := lineStart to lineEnd do begin
      for y := 0 to (lYi - 1) do begin
          for x := 0 to (lXi-1) do begin
            lineStart :=  x+(lXi * y);
            sum := 0.0;
            for j := 0 to contrib^[z].n - 1 do begin
              sum := sum + (contrib^[z].p^[j].weight * inImg[lineStart +contrib^[z].p^[j].pixel] );
            end;
            outImg[i] := sum;
            i := i + 1;
        end; //for X
    end; //for Y
  end; //for Z
end; //SubProcZ()  
begin
	if (zScale = 1) then exit(); //nothing to do
	{$IFDEF MYTHREADS}
	if (MaxThreads < 1) then MaxThreads := GetSystemThreadCount();
	{$ELSE}
	MaxThreads := 1;
	{$ENDIF}
	lXi := lHdr.dim[1]; //input X
	lYi := lHdr.dim[2]; //input X
	lZi := lHdr.dim[3]; //input X
	Zoom(lHdr,1, 1, zScale);
	lZo := lHdr.dim[3]; //output Z	
	inImg := copy(TFloat32s(lImg8), 0, lXi * lYi * lZi);
	setlength(lImg8,lXi*lYi*lZo*sizeof(single));
	outImg := TFloat32s(lImg8);
	SetContrib(contrib, lZi, lZo, (lXi*lYi), zScale, fwidth, filter);
	lineMax := lZo-1;
	linePerThread := ceil((lineMax+1)/MaxThreads);
	{$IFDEF MYTHREADS}
	ProcThreadPool.DoParallelNested(SubProcZ,0,MaxThreads-1, nil, MaxThreads);
	{$ELSE}
	SubProcZ(0);
	{$ENDIF}
	for i := 0 to lZo - 1 do
    	FreeMem(contrib^[i].p);
    FreeMem(contrib);
    //writeln(format('z d %d %d %d -> %d',[lXi,lYi,lZi, lZo]));
end; //ResizeZ()

procedure Resize32(var lHdr: TNIFTIhdr; var lImg8: TUInt8s; xScale, yScale, zScale, fwidth: single; filter: TFilterProc; MaxThreads: PtrInt = 0);
var
	mn, mx: single;
begin
	FindRange(lHdr,  TFloat32s(lImg8), mn, mx);
	if (zScale > 1.0) or (zScale > xScale) then begin
		//z is slowest dimension, do first when up-sampling
		//writeln('rescale zyx');
		ResizeZ(lHdr, lImg8, zScale, fwidth, filter, MaxThreads);
		ResizeY(lHdr, lImg8, yScale, fwidth, filter, MaxThreads);
		ResizeX(lHdr, lImg8, xScale, fwidth, filter, MaxThreads);
	end else begin
		//z is slowest dimension, do last when down-sampling
		//writeln('rescale xyz');
		ResizeX(lHdr, lImg8, xScale, fwidth, filter, MaxThreads);
		ResizeY(lHdr, lImg8, yScale, fwidth, filter, MaxThreads);
		ResizeZ(lHdr, lImg8, zScale, fwidth, filter, MaxThreads);
	end;
	CropRange(lHdr,  TFloat32s(lImg8), mn, mx);
end;

{$ENDIF}

function SetFilter(filterIndex: integer; out filter: TFilterProc): single;
begin
	//Lanczos nice for downsampling
	//Mitchell nice for upsampling
	if filterIndex = 1 then begin
	  filter := @TriangleFilter; result := 1;
	end else if filterIndex = 2 then begin
	  filter := @HermiteFilter; result := 1;
	end else if filterIndex = 3 then begin
	 filter := @BellFilter; result := 1.5;
	end else if filterIndex = 4 then begin
	 filter := @SplineFilter; result := 2;
	end else if filterIndex = 5 then begin
	 filter := @Lanczos3Filter; result := 3;
	end else begin
	 filter := @MitchellFilter; result := 2;
	end;
end;

function ShrinkOrEnlarge(var lHdr: TNIFTIhdr; var lBuffer: TUInt8s; lScale: single; MaxThreads: PtrInt = 0; filterIndex: integer = 0): boolean; overload;
var
   fwidth: single;
   filter : TFilterProc;
begin
  result := false;
  if lHdr.datatype <> kDT_FLOAT then exit;
  if (lScale = 1.0) then exit; //no resize
  if (lHdr.dim[4] > 1) then exit; //not for 4D
  fwidth := SetFilter(filterIndex, filter);
  Resize32(lHdr, lBuffer, lScale, lScale, lScale, fwidth, @filter, MaxThreads);
  result := true;
end;

function ShrinkOrEnlarge(var lHdr: TNIFTIhdr; var lBuffer: TUInt8s; var lScale: TVec3; MaxThreads: PtrInt = 0; filterIndex: integer = 0): boolean; overload;//
var
   fwidth, mx, mn: single;
   filter : TFilterProc;
begin
  result := false;
  if lHdr.datatype <> kDT_FLOAT then exit;
  if (lScale.x = lScale.y) and (lScale.y = lScale.z) then begin
	mn := min(abs(lHdr.pixdim[1]), min(abs(lHdr.pixdim[2]), abs(lHdr.pixdim[3])));
	mx := max(abs(lHdr.pixdim[1]), max(abs(lHdr.pixdim[2]), abs(lHdr.pixdim[3])));
	if (mn <> 0) and ((mx/mn) > 1.02) then begin
		lScale.x := lScale.x* (lHdr.pixdim[1]/mn);
		lScale.y := lScale.y* (lHdr.pixdim[2]/mn);
		lScale.z := lScale.z* (lHdr.pixdim[3]/mn);
		{$ifdef unix} //windows GUI can not write to console
		//writeln(format('Anisotropic image resliced %gx%gx%g', [lScale.x, lScale.y, lScale.z]));
		{$endif}
	end;
  end;
  if (lScale.x = 1.0) and (lScale.y = 1.0) and (lScale.z = 1.0) then exit; //no resize
  if (lHdr.dim[4] > 1) then exit; //not for 4D
  fwidth := SetFilter(filterIndex, filter);
  Resize32(lHdr, lBuffer, lScale.x, lScale.y, lScale.z, fwidth, @filter, MaxThreads);
  result := true;
end;

end.

