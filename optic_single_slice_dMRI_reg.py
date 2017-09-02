#!/usr/bin/env python
#
# coding: utf-8

# Optic nerve single slice dMRI registratino
#
# Joo-won Kim
# Icahn School of Medicine at Mount Sinai
#
# https://github.com/junqianxulab/optic_nerve_single_slice_dMRI_registration

import nibabel as nib
import numpy as np
import os
import sys
import scipy.ndimage
import scipy.interpolate

def resize_img(dat, mask, d=5):
    nonzero = mask.nonzero()
    x0, xn = nonzero[0].min()-d, nonzero[0].max()+d
    y0, yn = nonzero[1].min()-d, nonzero[1].max()+d - 1
    if x0 < 0:
        x0 = 0
    if xn >= dat.shape[0]:
        xn = dat.shape[0] - 1
    if y0 < 0:
        y0 = 0
    if yn >= dat.shape[1]:
        yn = dat.shape[1] - 1

    return dat[x0:xn+1, y0:yn+1, :, :].copy(), mask[x0:xn+1, y0:yn+1].copy(), (x0, y0)


if len(sys.argv) < 2:
    sys.stderr.write('Usage %s nifti_filename [mask_filename]\n' % sys.argv[0])
    sys.exit(-1)

# set filenames
fn = sys.argv[1]
if fn[-7:] == '.nii.gz':
    bn = fn[:-7]
elif fn[-4:] == '.nii':
    bn = fn[:-4]
else:
    bn = fn

if len(sys.argv) > 2:
    fn_mask = sys.argv[2]
else:
    fn_mask = bn + '_mask.nii.gz'
    if not os.path.isfile(fn_mask):
        fn_mask = bn + '_mask.nii'
        if not os.path.isfile(fn_mask):
            sys.stderr.write('mask file not found\n')
            sys.exit(-1)

img = nib.load(fn)
dat_ori = img.get_data()
mask_ori = nib.load(fn_mask).get_data()

if len(mask_ori.shape) > 2:
    mask_ori = mask_ori[:,:,0]

dat, mask, (x0, y0) = resize_img(dat_ori, mask_ori)

sigma = 1
size = (2,2)
size_out = (5,5)

search_domain = [ (xx, yy) for xx in np.arange((size[0]-1), dat.shape[0]-(size[0]-1)+0.01, 0.2)
                           for yy in np.arange((size[1]-1), dat.shape[1]-(size[1]-1)+0.01, 0.2) ]

dat_mod = dat.astype(np.float)
dat_max = np.zeros(dat.shape, dtype=dat.dtype)
dat_reduced = np.zeros(size_out +(1, dat.shape[-1]), dtype=dat.dtype)

for f in range(dat_mod.shape[-1]):
    dat_mod[:,:,:,f][mask == 0] = 0
    dat_mod[:,:,:,f] = scipy.ndimage.gaussian_filter(dat_mod[:,:,:,f], sigma=sigma)
    dat_mod[:,:,:,f] /= dat_mod[:,:,:,f].max()
    dat_mod[:,:,:,f] *= 100

    f_interp = scipy.interpolate.interp2d(range(dat_mod.shape[1]), range(dat_mod.shape[0]),
                                          dat_mod[:,:,0,f], kind='cubic', fill_value=0.0)

    # use optimization?
    max_sum = 0.0
    max_ind = -1
    for i, (xx, yy) in enumerate(search_domain):
        sum_i = f_interp(
                         np.arange(yy-(size[1]-1), yy+(size[1]-1)+0.1, 1.0),
                         np.arange(xx-(size[0]-1), xx+(size[0]-1)+0.1, 1.0)
        )
        if sum_i.sum() > max_sum:
            max_sum = sum_i.sum()
            max_ind = i
    imax = search_domain[max_ind]

    f_interp = scipy.interpolate.interp2d(range(dat.shape[1]), range(dat.shape[0]),
                                          dat[:,:,0,f], kind='cubic')
    xx, yy = imax
    dat_reduced[:, :, 0, f] = f_interp(
                                       np.arange(yy-(size_out[1]-1)/2.0, yy+(size_out[1]-1)/2.0+0.1, 1.0),
                                       np.arange(xx-(size_out[0]-1)/2.0, xx+(size_out[0]-1)/2.0+0.1, 1.0)
    )
    for dx, dy in [ (t1, t2) for t1 in (-1, 0, 1) for t2 in (-1, 0, 1)]:
        dat_max[int(round(imax[0]))+dx, int(round(imax[1]))+dy, 0, f] = 1

dat_mod_lg = np.zeros(dat_ori.shape, dtype=dat_mod.dtype)
dat_mod_lg[x0:x0+dat_mod.shape[0],y0:y0+dat_mod.shape[1]] = dat_mod
img_out = nib.Nifti1Image(dat_mod_lg, img.affine, img.header)
nib.save(img_out, bn + '_gaussian.nii.gz')

dat_max_lg = np.zeros(dat_ori.shape, dtype=dat_mod.dtype)
dat_max_lg[x0:x0+dat_max.shape[0],y0:y0+dat_max.shape[1]] = dat_max
img_out = nib.Nifti1Image(dat_max_lg, img.affine, img.header)
nib.save(img_out, bn + '_rough_seg.nii.gz')

img_out = nib.Nifti1Image(dat_reduced, img.affine, img.header)
nib.save(img_out, bn + '_reduced.nii.gz')

print 'fslview %s %s -l Red -t 0.3' % (fn, bn + '_rough_seg.nii.gz')
print 'run dMRI fitting on %s' % (bn + '_reduced.nii.gz')

