## About

These example images validate different features of the software.

 - AICHA demonstrates how each depth should be estimated independently for each region in an atlas. The "-t 0" option will process each of the 192 regions independently.
  * https://www.ncbi.nlm.nih.gov/pubmed/26213217
  * http://www.gin.cnrs.fr/en/tools/aicha/
 - 4Dgraywhite is a 4D image with two 3D volumes embedded in a single file. The first volume is a gray matter map and the second is a white matter map. Tissue depth should be estimated independently for each volume. Note that this example is of low resolution and binary, which make a compact, fast test dataset. Since this software works on binary images, and most tissue probability maps are continuous images, in many real world situations you would upsample a probability map prior to binarizing. 
 - Region4 shows a C-shaped object, where center of mass is not an ideal choice for selecting the object.
 - Anisotropic is a sphere saved where the distance between voxels is 2x4x8mm in the row, column and slice directions. The rapid method used here does not handle this corner case. Software implementing this fast algorithm need to detect this so the user can reslice the data to be isotropic or use a different algorithm. Note that other methods like iterative erosion also tend to assume isotropic resolution.
 - Isotropic is a sphere saved where the distance between voxels is 2x2x2mm in the row, column and slice directions. This is a reference image for the anisotropic test. 
 

