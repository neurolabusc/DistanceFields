## About

These example images validate different features of the software.

 - AICHA demonstrates how each depth should be estimated independently for each region in an atlas. The "-t 0" option will process each of the 192 regions independently.
  * https://www.ncbi.nlm.nih.gov/pubmed/26213217
  * http://www.gin.cnrs.fr/en/tools/aicha/
1605
 - inia19-NeuroMaps is an atlas with 1605 independent regions. This number of regions does not fit into a traditional 8-bit datatype. Further, this number of regions would be very slow to process with traditional measures.
  * https://www.ncbi.nlm.nih.gov/pubmed/23230398
 - 4Dgraywhite is a 4D image with two 3D volumes embedded in a single file. The first volume is a gray matter map and the second is a white matter map. Tissue depth should be estimated independently for each volume. Note that this example is of low resolution and binary, which make a compact, fast test dataset. Since this software works on binary images, and most tissue probability maps are continuous images, in many real world situations you would upsample a probability map prior to binarizing. 
 - Region4 shows a C-shaped object, where center of mass is not an ideal choice for selecting the object.
 - Anisotropic is a sphere saved where the distance between voxels is 2x4x8mm in the row, column and slice directions. This algorithm (similar to erosion methods) are not designed for anisotropic images. Therefore, Depth3D should up-sample this image (91*55*23 voxels) to be isotropic (91*110*92 voxels) prior to estimating the depth, and then down-sample the resulting image back to the original resolution.
 - Isotropic is a sphere saved where the distance between voxels is 2x2x2mm in the row, column and slice directions. This is a reference image for the anisotropic test. 
 - avg152T1_gray is a gray matter probability map. Where 4Dgraywhite has a binary map of gray matter, here each voxel stores a probability between zero and one. For these continuous images, it is often better to super-sample images during depth estimates. However, our algorithm requires us to categorically threshold each voxel as being either inside or outside the brain. In these situations, super-sampling can be helpful. For example, `Depth3D -t 0.25 -s 2 ./test/avg152T1_gray.nii.gz` upscales the image by x2 in each dimension prior to estimating depth (and downsamples the results). This tends to look less jagged than estimating on the base resolution `Depth3D -t 0.25 ./test/avg152T1_gray.nii.gz`.
 

