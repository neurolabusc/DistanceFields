unit nifti_smooth;
{$mode Delphi} 
{$H+}
{$DEFINE MYTHREADS} //multi-threaded
{$IFDEF MYTHREADS}{$ModeSwitch nestedprocvars}{$ENDIF}
interface

uses 
   {$IFDEF MYTHREADS}mtprocs,mtpcpu,{$ENDIF}
	SysUtils,
	SimdUtils, VectorMath, Classes, nifti_types, math;

function nii_smooth_gauss (var hdr: TNIFTIhdr; img: TFloat32s; FWHMmm: single; MaxThreads: integer = 0; isMasked: boolean = false):boolean;

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

type
	TFloatRA = array [0..0] of TScalar;
	TFloatRAp = ^TFloatRA;
	TKernel = record
		Len: integer; //pixels in row
		k: TFloat32s; //Gaussian kernel 
		kWeight: TFloat32s; //Weight for Kernel at each pixel
		kStart,kEnd: TInt32s; //start and end kernel at each pixel, e.g. -3..+3
	end;
	
procedure MakeKernel (FWHMmm, mm: single; len: integer; out k: TKernel; isFaster: boolean = true); inline;
var
  i, j, cutoffvox: integer;
  sigma, expd, wt: double;
begin
	 k.Len := len;
	 sigma  := (FWHMmm/mm)/sqrt(8*ln(2));  //      % FWHM -> sigma
	 if (isFaster) then //mimic faster but lower precision AFNI_BLUR_FIRFAC = 2.5
	 	//https://github.com/afni/afni/blob/699775eba3c58c816d13947b81cf3a800cec606f/src/edt_blur.c
	 	cutoffvox := ceil(  2.5 * sigma)
	 else
     	cutoffvox  := round(6*sigma);       //    % highest / lowest voxel to go out to
     //printf(format('cutoffvox %d', [cutoffvox]));
     setlength(k.k,cutoffvox+1);
     expd := 2*sigma*sigma;
     for i := 0 to cutoffvox do begin
            k.k[i] := exp(-1*(i*i)/expd) ;
     end;
     //calculate start, end for each voxel in 
     setlength(k.kStart,len); //first pixel in kernel, e.g. -3: cutoffvox except edges
     setlength(k.kEnd,len); //last voxel in kernel, e.g. +3: cutoffvox except edges
     setlength(k.kWeight,len); //sum of all voxels in kernel, ~1 except near edges
     for i := 0 to len-1 do begin
     	k.kStart[i] := max(-cutoffvox, -i);//do not read below 0
     	k.kEnd[i] := min(cutoffvox, len-i-1);//do not read above (len-1)
     	if (i > 0) and (k.kStart[i] = (k.kStart[i-1])) and (k.kEnd[i] = (k.kEnd[i-1])) then begin
     		k.kWeight[i] := k.kWeight[i-1];
     		continue;	
     	end;
     	wt := 0;
     	for j := k.kStart[i] to k.kEnd[i] do
     		wt += k.k[abs(j)];
     	k.kWeight[i] := 1 / wt;
     end;
end; //MakeKernel

procedure blur(var row: TFloatRAp; var tmp: TFloat32s; var k: TKernel); inline;
var
	i, j: integer;
	sum: single;
begin
	for i := 0 to (k.Len - 1) do
		tmp[i] := row[i];
	for i := 0 to (k.Len - 1) do begin
		sum := 0;
		for j := k.kStart[i] to k.kEnd[i] do
			sum += tmp[i+j] * k.k[abs(j)];
		row[i] := sum * k.kWeight[i];
	end;	
end;

FUNCTION specialsingle (var s:single): boolean; inline;
//returns true if s is Infinity, NAN or Indeterminate
//4byte IEEE: msb[31] = signbit, bits[23-30] exponent, bits[0..22] mantissa
//exponent of all 1s =   Infinity, NAN or Indeterminate
CONST kSpecialExponent = 255 shl 23;
VAR Overlay: LongInt ABSOLUTE s;
BEGIN
  IF ((Overlay AND kSpecialExponent) = kSpecialExponent) THEN
     RESULT := true
  ELSE
      RESULT := false;
END;

procedure blurMasked(var row: TFloatRAp; var tmp: TFloat32s; var k: TKernel); inline;
var
	i, j: integer;
	sum, wt: single;
begin
	for i := 0 to (k.Len - 1) do
		tmp[i] := row[i];
	for i := 0 to (k.Len - 1) do begin
		if specialsingle(row[i]) then continue;
		sum := 0;
		wt := 0;
		for j := k.kStart[i] to k.kEnd[i] do begin
			if (specialsingle(tmp[i+j])) then continue;
			sum += tmp[i+j] * k.k[abs(j)];
			wt += k.k[abs(j)];
		end;
		if (wt > 0.0) then
			wt := 1.0 / wt;
		row[i] := sum * wt;
	end;	
end;

procedure RemoveNaNs(var hdr: TNIFTIhdr; img: TFloat32s);
var
	i, n: int64;
begin
	n := hdr.dim[1] * hdr.dim[2] * max(hdr.dim[3],1) * max(hdr.dim[4],1);
	for i := 0 to (n-1) do
		if specialsingle(img[i]) then
			img[i] := 0;
end;

{$IFDEF MYTHREADS}
function nii_smooth_gauss (var hdr: TNIFTIhdr; img: TFloat32s; FWHMmm: single; MaxThreads: integer = 0; isMasked: boolean = false):boolean;
//fslmaths ax.nii -s 8 ssax.nii
var
	k: TKernel;
	rows, cols, tHi, tPerThread: int64;
procedure SubProcX(ThreadIndex: PtrInt; Data: Pointer; Item: TMultiThreadProcItem);
var
	i, tStart, tEnd, x: int64;
	f: TFloatRAp;
	tmp: TFloat32s;
begin
	tStart := (ThreadIndex * tPerThread);
	if (tStart > tHi) then exit; //more threads than slices in Z direction
	tEnd := tStart + tPerThread - 1; //e.g. if zStart=4 and zPerThread=1 then zEnd=4 
	tEnd := min(tEnd, tHi); //final thread when slices in Z not evenly divisible by number of threads
	//printf(format('x%d %d..%d', [ThreadIndex, tStart, tEnd]));
	x := hdr.dim[1];
	setlength(tmp, k.Len); //input data before blur
	for i := tStart to tEnd do begin
		f := @img[i * x];
		if isMasked then
			blurMasked(f, tmp, k)	
		else
			blur(f, tmp, k);	
	end; //blur all rows
end; //nested SubProcX()
procedure SubProcY(ThreadIndex: PtrInt; Data: Pointer; Item: TMultiThreadProcItem);
var
	r, s, so, x, y, tStart, tEnd: int64;
	f: TFloatRAp;
	tmp, img2Dy: TFloat32s;
begin
	tStart := (ThreadIndex * tPerThread);
	if (tStart > tHi) then exit; //more threads than slices in Z direction
	tEnd := tStart + tPerThread - 1; //e.g. if zStart=4 and zPerThread=1 then zEnd=4 
	tEnd := min(tEnd, tHi); //final thread when slices in Z not evenly divisible by number of threads
	setlength(img2Dy, rows*cols);
	setlength(tmp, k.Len); //input data before blur
	//printf(format('y%d %d..%d', [ThreadIndex, tStart, tEnd]));
	for s := tStart to tEnd do begin	
		//transpose
		so := s * (rows * cols); //slice offset
		for x := 0 to (cols-1) do begin
			for y := 0 to (rows-1) do begin
				img2Dy[x+(y*cols)] := img[so];
				so := so + 1;
			end;
		end;
		for r := 0 to (rows-1) do begin
			f := @img2Dy[r*cols];
			if isMasked then
				blurMasked(f, tmp, k)	
			else
				blur(f, tmp, k);
		end;
		//transpose
		so := s * (rows * cols); //slice offset
		for x := 0 to (cols-1) do begin
			for y := 0 to (rows-1) do begin
				img[so] := img2Dy[x+(y*cols)];
				so := so + 1;
			end;
		end;
	end;
end; //nested SubProcY()
procedure SubProcZ(ThreadIndex: PtrInt; Data: Pointer; Item: TMultiThreadProcItem);
var
	si, so, r, s, x, y, n, sxy, sxyz, tStart, tEnd, yStart, yEnd: int64;
	f: TFloatRAp;
	tmp, img2Dz: TFloat32s;
begin
	tStart := (ThreadIndex * tPerThread);
	if (tStart > tHi) then exit; //more threads than slices in Z direction
	tEnd := tStart + tPerThread - 1; //e.g. if zStart=4 and zPerThread=1 then zEnd=4 
	tEnd := min(tEnd, tHi); //final thread when slices in Z not evenly divisible by number of threads
	setlength(img2Dz, rows*cols);
	setlength(tmp, k.Len); //input data before blur
	sxy := hdr.dim[1] * hdr.dim[2];	
	sxyz := sxy * hdr.dim[3];
	yStart := 0;
	yEnd := hdr.dim[2]-1;
	if (hdr.dim[4] < MaxThreads) then begin
		yStart := tStart;
		yEnd := tEnd;
		tStart := 0;
		tEnd := hdr.dim[4]-1;
	end;
	//printf(format('z%d %d..%d', [ThreadIndex, tStart, tEnd]));
	for n := tStart to tEnd do begin
		for s := yStart to (yEnd) do begin	
			//transpose
			si := (s * rows) + (n * sxyz);
			so := 0;
			for x := 0 to (rows-1) do begin
				for y := 0 to (cols-1) do begin
					img2Dz[so] := img[x + si + (y*sxy)];
					so += 1;
				end;
			end;
			for r := 0 to (rows-1) do begin
				f := @img2Dz[r*cols];
				if isMasked then
					blurMasked(f, tmp, k)	
				else
					blur(f, tmp, k);
			end;
			//transpose
			so := 0; //slice offset along Y axis
			for x := 0 to (rows-1) do begin
				for y := 0 to (cols-1) do begin
					img[x + si + (y*sxy)] := img2Dz[so];
					so += 1;
				end;
			end;
		end; //for s: slice
	end; //for n: time volume
end; //nested SubProcZ()
begin
	result := false;
	if (MaxThreads < 1) then MaxThreads :=  GetSystemThreadCount();
	hdr.dim[4] := max(hdr.dim[4],1);
	if (hdr.dim[1] < 2) or (hdr.dim[2] < 2) then 
		exit;
	//smooth in X: everything aligned so no need to transpose
	MakeKernel(FWHMmm, hdr.pixdim[1], hdr.dim[1], k);
	tHi := (hdr.dim[2] * hdr.dim[3] * hdr.dim[4])-1;
	tPerThread := ceil((tHi+1)/MaxThreads);
	ProcThreadPool.DoParallelNested(SubProcX,0,MaxThreads-1, nil, MaxThreads);
	MakeKernel(FWHMmm, hdr.pixdim[2], hdr.dim[2], k);
	//setlength(tmp, k.Len); //input data before blur
	cols := hdr.dim[2];
	rows := hdr.dim[1];
	tHi := hdr.dim[3] * hdr.dim[4]-1;
	tPerThread := ceil((tHi+1)/MaxThreads);
	ProcThreadPool.DoParallelNested(SubProcY,0,MaxThreads-1, nil, MaxThreads);
	//smooth in Z: transpose each sagittal slice
	MakeKernel(FWHMmm, hdr.pixdim[3], hdr.dim[3], k);
	cols := hdr.dim[3];
	rows := hdr.dim[1];
	if (hdr.dim[4] < MaxThreads) then
		tHi := hdr.dim[2]-1	//for 3D data parallel process different columns
	else
		tHi := hdr.dim[4]-1; //for 4D data parallel processes different volumes
	tPerThread := ceil((tHi+1)/MaxThreads);
	ProcThreadPool.DoParallelNested(SubProcZ,0,MaxThreads-1, nil, MaxThreads);
	result := true;
	if (isMasked) then
		RemoveNaNs(hdr, img);
end;
{$ELSE}
function nii_smooth_gauss (var hdr: TNIFTIhdr; img: TFloat32s; FWHMmm: single; MaxThreads: integer = 0; isMasked: boolean = false):boolean;
//fslmaths ax.nii -s 8 ssax.nii
var
	k: TKernel;
var
	rows, cols, si, so, r, s, i, x, y, n, sxy, sxyz: int64;
	f: TFloatRAp;
	tmp, img2D: TFloat32s;
begin
	result := false;
	hdr.dim[4] := max(hdr.dim[4],1);
	if (hdr.dim[1] < 2) or (hdr.dim[2] < 2) then 
		exit;
	//smooth in X: everything aligned so no need to transpose
	MakeKernel(FWHMmm, hdr.pixdim[1], hdr.dim[1], k);
	setlength(tmp, k.Len); //input data before blur
	x := hdr.dim[1];
	n := hdr.dim[2] * hdr.dim[3] * hdr.dim[4];
	for i := 0 to (n-1) do begin
		f := @img[i * x];
		if isMasked then
			blurMasked(f, tmp, k)	
		else
			blur(f, tmp, k);	
	end; //blur all rows
	//smooth in Y: transpose each axial slice
	MakeKernel(FWHMmm, hdr.pixdim[2], hdr.dim[2], k);
	setlength(tmp, k.Len); //input data before blur
	cols := hdr.dim[2];
	rows := hdr.dim[1];
	setlength(img2D, rows*cols);
	n := hdr.dim[3] * hdr.dim[4];
	for s := 0 to (n-1) do begin	
		//transpose
		so := s * (rows * cols); //slice offset
		for x := 0 to (cols-1) do begin
			for y := 0 to (rows-1) do begin
				img2D[x+(y*cols)] := img[so];
				so := so + 1;
			end;
		end;
		for r := 0 to (rows-1) do begin
			f := @img2D[r*cols];
			if isMasked then
				blurMasked(f, tmp, k)	
			else
				blur(f, tmp, k);
		end;
		//transpose
		so := s * (rows * cols); //slice offset
		for x := 0 to (cols-1) do begin
			for y := 0 to (rows-1) do begin
				img[so] := img2D[x+(y*cols)];
				so := so + 1;
			end;
		end;
	end;
	//smooth in Z: transpose each sagittal slice
	MakeKernel(FWHMmm, hdr.pixdim[3], hdr.dim[3], k);
	setlength(tmp, k.Len); //input data before blur
	cols := hdr.dim[3];
	rows := hdr.dim[1];
	setlength(img2D, rows*cols);
	sxy := hdr.dim[1] * hdr.dim[2];	
	sxyz := sxy * hdr.dim[3];	
	for n := 0 to (hdr.dim[4]-1) do begin
		for s := 0 to (hdr.dim[2]-1) do begin	
			//transpose
			si := (s * rows) + (n * sxyz);
			so := 0;
			for x := 0 to (rows-1) do begin
				for y := 0 to (cols-1) do begin
					img2D[so] := img[x + si + (y*sxy)];
					so += 1;
				end;
			end;
			for r := 0 to (rows-1) do begin
				f := @img2D[r*cols];
				if isMasked then
					blurMasked(f, tmp, k)	
				else
					blur(f, tmp, k);
			end;
			//transpose
			so := 0; //slice offset along Y axis
			for x := 0 to (rows-1) do begin
				for y := 0 to (cols-1) do begin
					img[x + si + (y*sxy)] := img2D[so];
					so += 1;
				end;
			end;
		end; //for s: slice
	end; //for t: time volume
	result := true;
	if (isMasked) then
		RemoveNaNs(hdr, img);
end;
{$ENDIF}

end.

