# Optic nerve single slice dMRI registration

## prerequisite

* python 2 (tested: version 2.7)
  * nibabel >= 2.0 (pip install nibabel)
  * numpy (pip install numpy)
  * scipy (pip install scipy)

## install

#### prerequisite installation

* python
  * test if installed:
    * `python --version`
  * linux, ubuntu or debian:
    * sudo apt-get install python python-pip
  * linux, centos or fedora:
    * sudo yum install python python-pip
  * macos
    * install [homebrew](https://brew.sh/)
    * brew install python
  * Or
    * install [anaconda](https://www.anaconda.com/download/)

* nibabel, numpy, scipy: either of the followings
  * pip install nibabel numpy scipy
  * sudo pip install nibabel numpy scipy
  * conda install nibabel numpy scipy

#### optic_single_slice_dMRI_reg installation

Download `optic_single_slice_dMRI_reg.py`

Add executable permission to `optic_single_slice_dMRI_reg.py`

Place `optic_single_slice_dMRI_reg.py` or its symbolic link in a PATH directory. (e.g. /usr/local/bin in linux or macos)

## usage

Convert dMRI to a NIFTI file and make a mask file of dMRI having

1) value of 1 on the optic nerve (all possible optic nerve locations through all frames)

2) value of 0 on others especially rectus muscle.

Run the following command:

```optic_single_slice_dMRI_reg.py image_filename mask_filename```

