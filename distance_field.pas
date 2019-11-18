unit distance_field;
{$DEFINE MYTHREADS}
{$mode Delphi} 
{$IFDEF MYTHREADS}{$ModeSwitch nestedprocvars}{$ENDIF}
{$H+}
interface
uses  
	{$IFDEF MYTHREADS}mtprocs,{$ENDIF}
	{$IFNDEF UNIX}system, {$ENDIF}
	dateutils, Math, VectorMath, SimdUtils, StrUtils, 
	SysUtils, Classes, nifti_types, nifti_loadsave, nifti_foreign;
  

function distanceFieldAtlas(var hdr: TNIFTIhdr; var img: TUInt8s; txtnam: string = ''; MaxThreads: PtrInt = 0): boolean;
function distanceFieldVolume(var hdr: TNIFTIhdr; var img: TUInt8s; txtnam: string = ''; threshold: single = 0.5): boolean;

implementation

type
	FloatRA = array [0..0] of TScalar;
	FloatRAp = ^FloatRA;

procedure printf(s: string);
begin
	{$IFDEF UNIX}
     writeln(s);
	{$ELSE}
	if IsConsole then //uses System
		writeln(s);	
	{$ENDIF}    
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

procedure edt(var f: FloatRAp; var d,z: TFloat32s; var v: TInt32s; n: integer; mm: single);
var
	p, k, q: integer;
	s, dx: TScalar;
function vx(): TScalar; inline;
begin
	if specialsingle(f[q]) or specialsingle(f[p])   then
		result := infinity
	else
		result := ((f[q] + q*q) - (f[p] + p*p)) / (2.0*q - 2.0*p);
end;
begin
    (*# Find the lower envelope of a sequence of parabolas.
    #   f...source data (returns the Y of the parabola vertex at X)
    #   d...destination data (final distance values are written here)
    #   z...temporary used to store X coords of parabola intersections
    #   v...temporary used to store X coords of parabola vertices
    #   i...resulting X coords of parabola vertices
    #   n...number of pixels in "f" to process
    # Always add the first pixel to the enveloping set since it is
    # obviously lower than all parabolas processed so far.*)
    k := 0;
    v[0] := 0;
    z[0] := -infinity;
    z[1] := infinity;
    for q := 1 to n-1 do begin
        (*# If the new parabola is lower than the right-most parabola in
        # the envelope, remove it from the envelope. To make this
        # determination, find the X coordinate of the intersection (s)
        # between the parabolas with vertices at (q,f[q]) and (p,f[p]).
        *)
        p := v[k];
        s := vx();
        while (s <= z[k]) do begin
            k := k - 1;
            p := v[k];
            s := vx();
        end;
        //# Add the new parabola to the envelope.
        k := k + 1;
        v[k] := q;
        z[k] := s;
        z[k + 1] := infinity;
    end;
    (*# Go back through the parabolas in the envelope and evaluate them
    # in order to populate the distance values at each X coordinate.*)
    k := 0;
    for q := 0 to n-1 do begin
        while z[k + 1] < q do
            k := k + 1;
        dx := (q - v[k]);
        d[q] := dx * dx + f[v[k]];
        //i[q] = v[k];
    end;
    for q := 0 to n-1 do
    	f[q] := d[q];
end;

procedure distanceFieldLR(var hdr: TNIFTIhdr; var img: TFloat32s; pixDimX: single = 1.0);
//filter data in the X dimension (Left/Right)
var
	r, rows, cols: integer;
	f: FloatRAp;
	d,z : TFloat32s;
	v: TInt32s;
begin
	cols := hdr.dim[1];
	rows := hdr.dim[2];
	setlength(z, cols+1);
	setlength(d, cols);
	setlength(v, cols);
	if hdr.dim[3] > 1 then rows := rows * hdr.dim[3];
	if hdr.dim[4] > 1 then rows := rows * hdr.dim[4];	
	for r := 0 to (rows-1) do begin		
		f := @img[r*cols];
		edt(f, d, z,v, cols, pixDimX);	
	end;
	//printf(format('LR%g',[pixDimX]));
end;

procedure distanceFieldAP(var hdr: TNIFTIhdr; var img: TFloat32s; pixDimY: single = 1.0);
//filter data in the Y dimension (Anterior/Posterior)
var
	x,y, k, s, r, rows, cols, slices: integer;
	f: FloatRAp;
	d,z, img2D : TFloat32s;
	v: TInt32s;
begin
	cols := hdr.dim[2];
	rows := hdr.dim[1];
	setlength(z, cols+1);
	setlength(d, cols);
	setlength(v, cols);
	setlength(img2D, rows*cols);
	slices := 1;
	if hdr.dim[3] > 1 then slices := slices * hdr.dim[3];
	if hdr.dim[4] > 1 then slices := slices * hdr.dim[4];
	for s := 0 to (slices-1) do begin	
		//transpose
		k := s * (rows * cols); //slice offset
		for x := 0 to (cols-1) do begin
			for y := 0 to (rows-1) do begin
				img2D[x+(y*cols)] := img[k];
				k := k + 1;
			end;
		end;
		for r := 0 to (rows-1) do begin		
			f := @img2D[r*cols];
			edt(f, d, z,v, cols, pixDimY);	
		end;
		//transpose
		k := s * (rows * cols); //slice offset
		for x := 0 to (cols-1) do begin
			for y := 0 to (rows-1) do begin
				img[k] := img2D[x+(y*cols)];
				k := k + 1;
			end;
		end;
	end;
	//printf(format('AP%g',[pixDimY]));
end;

procedure distanceFieldHF(var hdr: TNIFTIhdr; var img: TFloat32s; pixDimZ: single = 1.0);
//filter data in the Z dimension (Head/Foot)
//by far the most computationally expensive pass
// unlike LR and AP, we must process 3rd (Z) and 4th (volume number) dimension separately
var
	vol, sx, sxy, sxyz, x,y,k, s, r, rows, cols, slices, vols: integer;
	f: FloatRAp;
	d,z, img2D : TFloat32s;
	v: TInt32s;
begin
	if (hdr.dim[3] < 2) then exit;
	//we could transpose [3,2,1] or [3,1,2] - latter improves cache?
	vols := min(1,  hdr.dim[4]);
	cols := hdr.dim[3];
	rows := hdr.dim[1];
	sxy := hdr.dim[1] * hdr.dim[2];
	sxyz := sxy * hdr.dim[3];
	slices := hdr.dim[2];;
	setlength(z, cols+1);
	setlength(d, cols);
	setlength(v, cols);
	setlength(img2D, rows*cols);
	for vol := 0 to (vols-1) do begin
		for s := 0 to (slices-1) do begin	
			//transpose
			sx := (vol * sxyz) +(s * rows);
			k := 0; //slice offset along Y axis
			for x := 0 to (rows-1) do begin
				for y := 0 to (cols-1) do begin
					img2D[k] := img[x + sx + (y*sxy)];
					k := k + 1;
				end;
			end;
			for r := 0 to (rows-1) do begin		
				f := @img2D[r*cols];
				edt(f, d, z,v, cols, pixDimZ);	
			end;
			//transpose
			k := 0; //slice offset along Y axis
			for x := 0 to (rows-1) do begin
				for y := 0 to (cols-1) do begin
					img[x + sx + (y*sxy)] := img2D[k];
					k := k + 1;
				end;
			end;
		end; //slice
	end; //vols
	//printf(format('HF%g',[pixDimZ]));
end; //distanceFieldHF()

procedure fixHdr(var hdr: TNIFTIhdr);
begin
	if hdr.pixdim[1] = 0 then hdr.pixdim[1] := 1;
	if hdr.pixdim[2] = 0 then hdr.pixdim[2] := 1;
	if hdr.pixdim[3] = 0 then hdr.pixdim[3] := 1;	
	if hdr.dim[1] < 1 then hdr.dim[1] := 1;
	if hdr.dim[2] < 1 then hdr.dim[2] := 1;
	if hdr.dim[3] < 1 then hdr.dim[3] := 1;
	if hdr.dim[4] < 1 then hdr.dim[4] := 1;
end;

Type
  TCluster = record
		CenterXYZ: TVec3; //Center of cluster as defined by furthest from edge
        mxThick: single; //Peak thickness of cluster
        SzMM3: single; //cluster volume in mm^3
  end;
  TClusters = array of TCluster;
  
procedure saveClusterAsTxt(fnm: string; var c: TClusters);
function f2s(f: single): string;
begin
	result := PadLeft(FloatToStrf(f,ffFixed, 8,2), 8)+' ';
end;
function i2s(f: single): string;
begin
	result := PadLeft(Inttostr(round(f)), 8)+' ';
end;
var
	Txt: TextFile;
	i: integer;  
begin
	if length(c) < 1 then exit;
	AssignFile(Txt, fnm);
	Rewrite(Txt);
	WriteLn(Txt, '# Thick3D interactive cluster table');
	//https://www.slicer.org/wiki/Coordinate_systems
	//https://afni.nimh.nih.gov/afni/community/board/read.php?1,144396,144396
	// NIfTI coordinates go from L>R, P>A, I>S
	WriteLn(Txt, '#Coordinate order = RAS'); 
	WriteLn(Txt, '#VolMM3    Th x     Th y     Th z    Peak     Label');
	WriteLn(Txt, '#------- -------- -------- -------- -------- --------');
	for i := 0 to (length(c)-1) do begin
		if (c[i].SzMM3 <= 0) then continue;
		WriteLn(Txt, i2s(c[i].SzMM3)+f2s(c[i].CenterXYZ.x)+f2s(c[i].CenterXYZ.y)+f2s(c[i].CenterXYZ.z)+f2s(c[i].mxThick) +i2s(i+1));
	end;	
	CloseFile(Txt);
end;

function clusterCenter(var hdr: TNIFTIhdr; var img: TFloat32s; vol: integer = 0): TCluster;
// hdr: NIfTI header specifies image dimensions
// img: Thickness map
// mxThick: reports maximum thickness for entire cluster
// vol: select volume to inspect (for 4D volumes, indexed from 0).
// result: location of thickest point in cluster, for ties this will be near center of mass
var
	x,y,z, i, n: integer;
	s, dx, mnDx, mxThick: single;
	sum, cog, cogDx, center: TVec3;
begin
	result.CenterXYZ :=  Vec3 (0,0,0);
	//first pass: find peak and center of gravity
	sum := Vec3 (0,0,0);
	mxThick := img[0];	
	i := 0;
	n := 0;
	for z := 0 to (max(1, hdr.dim[3])-1) do begin
		for y := 0 to (hdr.dim[2] -1) do begin
			for x := 0 to (hdr.dim[1]-1) do begin
				s := img[i];
				i += 1;
				if (s = 0) then continue;
				sum += Vec3(x,y,z);
				n += 1;
				if (s > mxThick) then mxThick := s;
			end; //x
		end; //y
	end; //z
	result.mxThick := mxThick;
	result.SzMM3 := n * hdr.pixdim[1] * hdr.pixdim[2] * hdr.pixdim[3]; 
	if n < 1 then
		exit;
	cog := sum / n;
	//second pass: find voxel at thickest location, if a tie, choose voxel nearest center of gravity
	i := 0;
	mnDx := infinity;
	center := cog;
	for z := 0 to (max(1, hdr.dim[3])-1) do begin
		for y := 0 to (hdr.dim[2] -1) do begin
			for x := 0 to (hdr.dim[1]-1) do begin
				s := img[i];
				i += 1;
				if (s < mxThick) then continue; //only inspect peaks
				cogDx := cog -  Vec3(x,y,z);
				dx := cogDx.Length;
				if dx > mnDx then continue; //for ties, choose peak nearest CoG
				mnDx := Dx; //
				center := Vec3(x,y,z);
			end; //x
		end; //y
	end; //z
	//convert XYZ from voxels to mm
	cog := center;
	result.CenterXYZ.x := cog.x*hdr.srow_x[0]+cog.y*hdr.srow_x[1]+cog.z*hdr.srow_x[2]+hdr.srow_x[3];
	result.CenterXYZ.y := cog.x*hdr.srow_y[0]+cog.y*hdr.srow_y[1]+cog.z*hdr.srow_y[2]+hdr.srow_y[3];
	result.CenterXYZ.z := cog.x*hdr.srow_z[0]+cog.y*hdr.srow_z[1]+cog.z*hdr.srow_z[2]+hdr.srow_z[3];
end;

function isIsotropic(var hdr: TNIFTIhdr; out voxMM: single): boolean;
var
	mn, mx: single;
begin
	mn := min(abs(hdr.pixdim[1]), min(abs(hdr.pixdim[2]), abs(hdr.pixdim[3])));
	mx := max(abs(hdr.pixdim[1]), max(abs(hdr.pixdim[2]), abs(hdr.pixdim[3])));
	voxMM := mn+(0.5*(mx-mn));
	if (mn = 0) or (specialsingle(mn)) then exit(false);
	if (mx/mn) > 1.02 then exit(false);
	exit(true);
end;

function distanceFieldVolume(var hdr: TNIFTIhdr; var img: TUInt8s;  txtnam: string = '';  threshold: single = 0.5): boolean;
//process a 3D scalar volume
var
	i, vx, vxObj, vol: integer;
	voxMM: single;
	img32: TFloat32s;
	c: TClusters;
begin
	result := false;
	if not isIsotropic(hdr, voxMM) then
		exit;
	fixHdr(hdr);
	changeDataType(hdr, img, kDT_FLOAT32);
	img32 := TFloat32s(img);
	vx := hdr.dim[1] * hdr.dim[2];
	if hdr.dim[3] > 1 then vx *= hdr.dim[3];
	if hdr.dim[4] > 1 then vx *= hdr.dim[4];
	if (threshold > 0) then begin	
		for i := 0 to (vx-1) do
			if img32[i] >= threshold then
				img32[i] := infinity
			else
				img32[i] := 0;
	end else begin
		for i := 0 to (vx-1) do
			if img32[i] <= threshold then
				img32[i] := infinity
			else
				img32[i] := 0;	
	end;
	vxObj := 0; //count how many voxels are inside object
	for i := 0 to (vx-1) do
		if img32[i] <> 0 then
			vxObj += 1;
	if (vxObj = 0) or (vxObj = vx) then begin
		printf(format('Error: the threshold %g does not discriminate between air and object.', [threshold]));
		exit;	
	end;	
	distanceFieldLR(hdr,img32);
	distanceFieldAP(hdr,img32);
	distanceFieldHF(hdr,img32);	
	for i := 0 to (vx-1) do
		img32[i] := sqrt(img32[i]) * voxMM;
	result := true;	
	if  txtnam = '' then exit;
	vol := max(1, hdr.dim[4]); 
	setlength(c, vol);
	for i := 0 to (vol-1) do
		c[i] := clusterCenter(hdr, img32, i);
	saveClusterAsTxt(txtnam, c);
end;

// 4 core machine
//  2.8s +THREADCAPABLE+MYTHREADS runs multi-threaded
// 13.3s +THREADCAPABLE-MYTHREADS runs single-threaded
// 12.5s -THREADCAPABLE-MYTHREADS optimized for single threaded 

function distanceFieldAtlas(var hdr: TNIFTIhdr; var img: TUInt8s;  txtnam: string = '';  MaxThreads: PtrInt = 0): boolean;
//process a 3D indexed integer volume
var
	i, roi, vx, mn, mx: integer;
	out32: TFloat32s;
	imgIdx: TInt32s;
	voxMM: single;
	c: TClusters;
{$IFDEF MYTHREADS}
//https://wiki.freepascal.org/Parallel_procedures
procedure SubProc(Index: PtrInt; Data: Pointer; Item: TMultiThreadProcItem);
{$ELSE}	
procedure SubProc(Index: PtrInt);
{$ENDIF}
  var
  	img32: TFloat32s;
  	j, vxObj: integer;	
  begin
		setlength(img32, vx);
		for j := 0 to (vx-1) do
			if imgIdx[j] = Index then
				img32[j] := infinity
			else
				img32[j] := 0;
		vxObj := 0; //count how many voxels are inside object
		for j := 0 to (vx-1) do
			if img32[j] <> 0 then
					vxObj += 1;
		if (vxObj = 0) or (vxObj = vx) then 
			exit; //sparse atlas, e.g. atlas has region 1 and 3 but not 2
		distanceFieldLR(hdr,img32, hdr.pixdim[1]);
		distanceFieldAP(hdr,img32, hdr.pixdim[2]);
		distanceFieldHF(hdr,img32, hdr.pixdim[3]);
		//only write to voxels INSIDE region, each thread writes to unique set of voxels	
		for j := 0 to (vx -1) do
			if (img32[j] > 0) then //this voxel is inside region "Index"
				out32[j] := img32[j];
		if  txtnam = '' then exit;
		c[Index-1] := clusterCenter(hdr, img32); 
		c[Index-1].mxThick := sqrt(c[Index-1].mxThick) * voxMM;	
  end; //nested SubProc()	
begin
	result := false;
	if not isIsotropic(hdr, voxMM) then
		exit;
	if (hdr.datatype <> kDT_INT8) and (hdr.datatype <> kDT_UINT8) and (hdr.datatype <> kDT_INT16) and (hdr.datatype <> kDT_UINT16) and (hdr.datatype <> kDT_INT32) then begin
		printf('distance fields for atlases (threshold=0) expects integer datatypes.');
		exit;
	end;
	fixHdr(hdr);
	changeDataType(hdr, img, kDT_INT32);  //convert header to float
	vx := hdr.dim[1] * hdr.dim[2];
	if hdr.dim[3] > 1 then vx *= hdr.dim[3];
	if hdr.dim[4] > 1 then vx *= hdr.dim[4];
	if vx < 1 then exit;
	imgIdx := copy(TInt32s(img), 0, vx);
	mn := imgIdx[0];
	mx := mn;
	for i := 0 to (vx-1) do begin
		roi := imgIdx[i];
		if (roi < mn) then mn := roi;
		if (roi > mx) then mx := roi;
	end;
	if (mx < 0) or (mn < 0) then begin
		printf('Invalid atlas: an atlas must have positive labels and must not have negative labels');
		exit;
	end;
	if mn = 0 then mn := 1; //do not compute distance field for "air"
	changeDataType(hdr, img, kDT_FLOAT32);  //convert header to float
	out32 := TFloat32s(img);
	for i := 0 to (vx-1) do
		out32[i] := 0; //assume all voxels outside atlas
	setlength(c, mx);
	for i := 0 to (mx-1) do
		c[i].SzMM3 := 0; //skip for sparse atlases
	
	{$IFDEF MYTHREADS}
	printf(format('Threaded filter of atlas with %d..%d regions',[mn,mx]));
	ProcThreadPool.DoParallelNested(subproc,mn,mx, nil, MaxThreads);
	//ProcThreadPool.DoParallelNested(subproc,mn,mx,nil, 4);
	//DoParallelNested(TMTNestedProcedure, StartIndex, EndIndex, Pointer = nil; MaxThreads: PtrInt = 0);
	{$ELSE}	
	printf(format('Filter of atlas with %d..%d regions',[mn,mx]));
	for roi := mn to mx do	
		subproc(roi);
	{$ENDIF}
	for i := 0 to (vx-1) do
		out32[i] := sqrt(out32[i]) * voxMM;		
	result := true;
	if  txtnam = '' then exit;
	saveClusterAsTxt(txtnam, c);	
end;

end.
 