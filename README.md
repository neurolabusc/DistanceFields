## About

Consider an image where some pixels/voxels show cortex and the rest indicate air. This program estimates the distance of each cortical voxel from air. Typical methods to do this use iterative erosion methods that are computationally intensive. This project uses the algorithm described by [Felzenszwalb and Huttenlocher](http://cs.brown.edu/people/pfelzens/dt/) to solve this problem as a separable filtering problem. An excellent introduction to this problem and application is provided by [Philip Rideout](https://prideout.net/blog/distance_fields/). Further, his [Snowy](https://github.com/prideout/snowy) project provides an elegant python implementation.

While Snowy is ideal for most situations, this project is specifically designed for high performance with medical images, reading various popular formats and saving distance fields in the popular NIfTI format, which you can view with many tools including [MRIcroGL](https://www.nitrc.org/projects/mricrogl).

![alt tag](https://github.com/neurolabusc/DistanceFields/blob/master/atlas.jpg)

## Installation

You can get Thick3D using two methods:

 - (Recommended) Download latest compiled release from [Github release web page](https://github.com/rordenlab/DistanceFields/releases).
 - (Recommended) You can also download from the command line for Linux, MacOS and Windows:
   * `curl -fLO https://github.com/rordenlab/DistanceFields/releases/latest/download/Thick3D_lnx.zip`
   * `curl -fLO https://github.com/rordenlab/DistanceFields/releases/latest/download/Thick3D_mac.zip`
   * `curl -fLO https://github.com/rordenlab/DistanceFields/releases/latest/download/Thick3D_win.zip`
 - (Developers) Download the source code from [GitHub](https://github.com/rordenlab/DistanceFields).

## Usage

```
Chris Rorden's Thick3D v1.0.20191118
 see https://prideout.net/blog/distance_fields/
usage: Thick3D [options] <in_file(s)>
Reads volume and computes distance fields
OPTIONS
 -3 : save 4D data as 3D files (y/n, default n)
 -h : show help
 -o : output name (omit to save as input name with "x" prefix)
 -r : generate text report (y/n/o[only, no image], default n)
 -t : threshold, less extreme values treated as outside (default 0.5)
       set to 0 for separate field for each region of an atlas
 -p : parallel threads (0=optimal, 1=one, 5=five, default 0)
 -z : gz compress images (y/n, default n)
 Examples :
  Thick3D -t 0 ./test/AICHA.nii.gz
  Thick3D -t 0.5 ./test/4Dgraywhite.nii.gz
  Thick3D -t 0.5 ./test/anisotropic.nii.gz
  Thick3D -o "~/out dir/OutputImg" -t 0.5 "./test/name with spaces"
```

## Results

By default, this software generates a NIfTI format thickness map. The `-r` option allows you to generate text reports. You can generate these reports along with the NIfTIM image (`-r y`) or only generate the text reports (`-r o`). For standard input images, one line of text is generated per volume. However, for an atlas (`-t 0`) one line is generated for each region in the atlas.
 
```
# Thick3D interactive cluster table
#Coordinate order = RAS
#VolMM3    Th x     Th y     Th z    Peak     Label
#------- -------- -------- -------- -------- --------
    1312   -14.00    66.00    12.00     2.24        1 
   17480     8.00    48.00    44.00     2.83        2 
    1184    20.00    18.00    60.00     1.73        3 
    7096    18.00    62.00   -10.00     3.00        4 
...
```

![alt tag](https://github.com/neurolabusc/DistanceFields/blob/master/graymatter.jpg)

## Performance

This is a very fast algorithm. However, atlas images are necessarily complex, as thickness is computed for each region independently. This tool uses two tricks that accelerate these cases. First, thickness is only computed for a box constrained by the size of the region being computed, rather than for the entire volume. This is useful for atlases with many regions, where each region typically is only a small portion of the overall 3D volume. In general, this trick makes the algorithm about 5 times faster. Further, this algorithm can leverage parallel processing. The included sample [AICHA](http://www.gin.cnrs.fr/en/tools/aicha/) atlas has 384 regions. It requires 1.65 seconds to process in single threaded mode and 0.45 seconds for threaded conversion on a 4-core (8-thread) laptop.

## Limitations

 - This software is provided as is. It uses the BSD license.
 - This algorithm assumes that the spatial resolution is identical in each dimension. Anisotropic data must either be resliced to an isotropic grid or a different method must be used.
 
## Compiling

Most people will want to download a [pre-compiled executable](https://github.com/rordenlab/DistanceFields/releases). Alternatively, [Snowy](https://github.com/prideout/snowy) provides a scriptable Python implementation. However, you can compile your own copy from source code:

 - Download and install [FreePascal for your operating system](https://www.freepascal.org/download.html). For Debian-based unix this may be as easy as `sudo apt-get install fp-compiler`. For other operating systems, you may simply want to install FreePascal from the latest [Lazarus distribution](https://sourceforge.net/projects/lazarus/files/).
 - From the terminal, go inside the directory with the source files and run the following commands to build and validate your compilation:

```
fpc Thick3D
./Thick3D ./test/AICHA.nii.gz
```

## Supported Image Formats

Thick3D supports the same formats as [i2nii](https://github.com/rordenlab/i2nii):

 - [AFNI Brik](https://afni.nimh.nih.gov/pub/dist/doc/program_help/README.attributes.html)(.head).
 - [Analyze](http://imaging.mrc-cbu.cam.ac.uk/imaging/FormatAnalyze)(.hdr).
 - [Bio-Rad PIC](https://docs.openmicroscopy.org/bio-formats/5.8.2/formats/bio-rad-pic.html)(.pic).
 - [Blender Voxel data](http://pythology.blogspot.com/2014/08/you-can-do-cool-stuff-with-manual.html)(.bvox).
 - [BrainVoyager VMR](https://support.brainvoyager.com/brainvoyager/automation-development/84-file-formats/343-developer-guide-2-6-the-format-of-vmr-files)(.vmr, .v16).
 - [DeltaVision](https://docs.openmicroscopy.org/bio-formats/5.8.2/formats/deltavision.html)(.dv).
 - [ECAT](http://nipy.org/nibabel/reference/nibabel.ecat.html)(.v).
 - [FreeSurfer MGH/MGZ Volume](https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/MghFormat)(.mgh/.mgz).
 - [Guys Image Processing Lab](http://rview.colin-studholme.net/rview/rv9manual/fileform.html#GIPL)(.gipl).
 - [ICS Image Cytometry Standard](https://onlinelibrary.wiley.com/doi/epdf/10.1002/cyto.990110502)(.ics).
 - [Interfile](https://www.ncbi.nlm.nih.gov/pubmed/2616095)(.varies, limited support).
 - [ITK MHA/MHD](https://itk.org/Wiki/MetaIO/Documentation)(.mha/.mhd).
 - [MRTrix Volume](https://mrtrix.readthedocs.io/en/latest/getting_started/image_data.html)(.mif/.mih; not all variants supported).
 - [NIfTI](https://brainder.org/2012/09/23/the-nifti-file-format/)(.hdr/.nii/.nii.gz/.voi). Useful for NIfTI-to-NIfTI modifications (e.g. `-3`, `-d`, `-r`, `-z`).
 - [NRRD](http://teem.sourceforge.net/nrrd/format.html)(.nhdr/.nrrd).
 - [POV-Ray Density_File](https://www.povray.org/documentation/view/3.6.1/374/)(.df3).
 - [Spectroscopic Imaging, Visualization and Computing (SIVIC)](https://radiology.ucsf.edu/research/labs/nelson#accordion-software)(.idf).
 - [Stimulate Sdt](https://www.cmrr.umn.edu/stimulate/stimUsersGuide/node57.html)(.spr/.sdt)
 - [Vaa3D](https://github.com/Vaa3D)(.v3draw).
 - [VTK Legacy Voxel Format](https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf)(.vtk).

If your image format is not supported, you may want to see if it is supported by the [Bio-Formats module](https://docs.openmicroscopy.org/bio-formats/5.9.2/supported-formats.html) of [ImageJ/Fiji](https://fiji.sc). If so, you can open the image with the module and save it as NIfTI or NRRD.

## Links

 - [Philip Rideout's distance fields](https://prideout.net/blog/distance_fields/) page was the inspiration for this project.
 - [Philip Rideout's snowy](https://prideout.net/snowy/) provides distance fields (and more) for 2D images with a nice Python wrapper.

