unit clusters;

{$mode objfpc}{$H+}
interface

uses SimdUtils;


function Clusterize(var lImg: TFloat32s; Xi,Yi,Zi: integer; out clusterNumber: integer; out img32: TInt32s; NeighborMethod: integer; smallestClusterVox: integer = 0): boolean;

implementation

function Clusterize(var lImg: TFloat32s; Xi,Yi,Zi: integer; out clusterNumber: integer; out img32: TInt32s; NeighborMethod: integer; smallestClusterVox: integer = 0): boolean;
label
     123;
var
  i, j, XY, XYZ, qlo, qhi: integer;
  qimg: TInt32s;
procedure checkPixel(vxl: integer); inline;
begin
     if img32[vxl] <> -1 then exit; //already found or not a target
     qhi := qhi + 1;
     img32[vxl] := clusterNumber; //found
     qimg[qhi] := vxl; //location
end;//nested checkPixel()
procedure retirePixel6(); inline;
var
  vxl: integer;
begin
     vxl := qimg[qlo];
     checkPixel(vxl-1);
     checkPixel(vxl+1);
     checkPixel(vxl-Xi);
     checkPixel(vxl+Xi);
     checkPixel(vxl-XY);
     checkPixel(vxl+XY);
     qlo := qlo + 1;
end;//nested retirePixel()
procedure retirePixel18(); inline;
var
  vxl: integer;
begin
     vxl := qimg[qlo];
     //edges in plane
     checkPixel(vxl-Xi-1);
     checkPixel(vxl-Xi+1);
     checkPixel(vxl+Xi-1);
     checkPixel(vxl+Xi+1);
     //edges below
     checkPixel(vxl-1-XY);
     checkPixel(vxl+1-XY);
     checkPixel(vxl-Xi-XY);
     checkPixel(vxl+Xi-XY);
     //edges above
     checkPixel(vxl-1+XY);
     checkPixel(vxl+1+XY);
     checkPixel(vxl-Xi+XY);
     checkPixel(vxl+Xi+XY);
     retirePixel6();
end;//nested retirePixel()
procedure retirePixel26(); inline;
var
  vxl: integer;
begin
     vxl := qimg[qlo];
     //corners below
     checkPixel(vxl-Xi-XY-1);
     checkPixel(vxl-Xi-XY+1);
     checkPixel(vxl+Xi-XY-1);
     checkPixel(vxl+Xi-XY+1);
     //corners above
     checkPixel(vxl-Xi+XY-1);
     checkPixel(vxl-Xi+XY+1);
     checkPixel(vxl+Xi+XY-1);
     checkPixel(vxl+Xi+XY+1);
     retirePixel18();
end;
begin //main RemoveSmallClusters()
  clusterNumber := 0;
  result := false;
  if (Zi < 1) then exit;
  XY := Xi * Yi;
  XYZ := XY * Zi;
  setlength(img32, XYZ);
  setlength(qimg, XYZ);
  //set target voxels
  for i := 0 to (XYZ-1) do begin
      img32[i] := 0;
      if lImg[i] > 0 then
         img32[i] := -1;
  end;
  //clear bottom and top slices
  for i := 0 to (XY-1) do
    img32[i] := 0;
  for i := (XYZ-1-XY) to (XYZ-1) do
    img32[i] := 0;
  //now seed each voxel
  for i := (XY) to (XYZ-1-XY) do begin
      if (img32[i] < 0) then begin //voxels not yet part of any region
         clusterNumber := clusterNumber + 1;
         if (clusterNumber < 1) then goto 123; //more than 2^32 clusters!
         qlo := 0;
         qhi := -1;
         checkPixel(i);
         while qlo <= qhi do begin
             case NeighborMethod of
               3,26:
                 retirePixel26;
               2,18:
                 retirePixel18;
               else //method 1 aka 6 neighbor method
                 retirePixel6;
             end;
           end;
         if (qhi +1) < smallestClusterVox then begin
           for j := 0 to qhi do
              img32[qimg[j]] := 0;//qhi + 1;
           clusterNumber := clusterNumber - 1;
         end;
      end;
  end;
  result := true;
123:
  qimg := nil;
  if not result then clusterNumber := 0;
  //img32 := nil;
end; 

end.

