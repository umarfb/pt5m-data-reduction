import astroalign as aa 
from astropy.io import fits
from scipy.ndimage import zoom
import numpy as np
import matplotlib.pyplot as plt
import os

# script to align image

filename = 'r0500344.fits'
filepath = '/local/php18ufb/backed_up_on_astro3/mini_project01/data20183009/astrometry_images/'

# load reference image and image to be aligned from pt5m
src_img = fits.getdata(filepath + 'r0500344_am.fits')
targ_img = fits.getdata(filepath + 'r0500345_am.fits')

# get fits header for target image
hdr = fits.getheader(filepath + 'r0500345_am.fits')

# increase source image size
#src_img = zoom(src_img, 1.5, order=2)
#fits.writeto('zoom349.fits', src_img)

# known star-to-star correspondence on image
targ_str = np.array([(484.0, 267.0), (878.0, 566.0), (257.0, 443.0), (239.0, 460.0), (899.0, 474.0), (781.0, 374.0)])
src_str = np.array([(334.0, 290.0), (727.0, 589.0), (105.0, 466.0), (88.0, 483.0), (748.0, 498.0), (631.0, 398.0)])

# estimate transform
tform = aa.estimate_transform('affine', targ_str, src_str)
#tform = aa.find_transform(src_img, targ_img)

# apply transformation to image to align
aligned_img = aa.apply_transform(tform, src_img, targ_img)

# update header
hdr['HISTORY'] = 'Aligned'

#aligned_img = aa.register(src_img, targ_img)
os.remove(filename)
fits.writeto(filename, aligned_img, hdr)