unit resize;

{$mode delphi}{$H+}
{$inline on}
{$DEFINE MYTHREADS} //multi-threaded
{$IFDEF MYTHREADS}{$ModeSwitch nestedprocvars}{$ENDIF}

interface

uses
{$IFDEF MYTHREADS}mtprocs,mtpcpu,{$ENDIF}
  Classes, SysUtils, nifti_types, SimdUtils;

function ShrinkOrEnlarge(var lHdr: TNIFTIhdr; var lBuffer: TUInt8s; lScale: single; MaxThreads: PtrInt = 0): boolean; overload;


implementation

uses math;

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
     //showmessage(format('%g -> %g %g %g; %g %g %g; %g %g %g',[iScale, lHdr.srow_x[0],lHdr.srow_y[0],lHdr.srow_z[0], lHdr.srow_x[1],lHdr.srow_y[1],lHdr.srow_z[1], lHdr.srow_x[2],lHdr.srow_y[2],lHdr.srow_z[2]]));
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
end;

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
  //lineStart, 
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

function ShrinkOrEnlarge(var lHdr: TNIFTIhdr; var lBuffer: TUInt8s; lScale: single; MaxThreads: PtrInt = 0): boolean; overload;
var
   fwidth: single;
   filter : TFilterProc;
begin
  result := false;
  if lHdr.datatype <> kDT_FLOAT then exit;
  if (lScale = 1.0) then exit; //no resize
  if (lHdr.dim[4] > 1) then exit; //not for 4D
  filter := @MitchellFilter;
  fwidth := 2;
  Resize32(lHdr, TUInt8s(lBuffer), lScale, lScale, lScale, fwidth, @filter, MaxThreads);
  result := true;
end;

end.

