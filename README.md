# Fragment Computation
This is a package extract from the [Lemsvxl](http://vision.lems.brown.edu/content/available-software-and-databases#LEMSVXL) for image fragment.
It use the vxl C++ library.

The input includes contour file(.cem/.cemv), the image and the training data.
## User Guide

### 1. Prerequisites

#### GCC-4.x
You probably need gcc-4.x to build this package.  

Use `$ gcc --version` to check your gcc version.  
Download gcc-4.8
```commandline
sudo apt-get install gcc-4.8
```
Use gcc-4.8 in your environment.
```commandline
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.8 100
```
Then the same as g++-4.x.

#### VXL library

Download the correct version of VXL.

Build the VXL with CMake  
e.g. under VXL root dir:
```commandline
$ mkdir build
$ cd build
$ cmake .. -G "Unix Makefiles"
$ make
```

Set VXL_DIR to your build directory:

```commandline
$ export VXL_DIR="~/YourDirToVXLRoot/build"
```

### 2. Build the fragment_computation

Build the package with CMake

e.g. under fragment_computation root dir:
```commandline
$ mkdir build
$ build
$ cmake .. -G "Unix Makefiles"
$ make 
```

### 3. Using the package

- Original image: It can be `.jpg` or `.png`.
- Contour file: Use [MSEL_contour_extraction](https://github.com/yg13/MSEL_contour_extraction_cxx) to get `.cem` (version 2) from image.  
Then use converter to get `.cem`/`.cemv` (version 1) from `.cem` (version 2).