import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import mad_std
from astropy.visualization import ZScaleInterval, ImageNormalize
from photutils import datasets
from photutils import DAOStarFinder
from photutils import aperture_photometry, CircularAperture


# load image data
img = fits.getdata('R_final.fits')

# subtract background, calculated using the median of the image
img -= np.median(img, axis=0)

# creates ZScale interval as used by DS9
interval = ZScaleInterval()
interval.get_limits(img)

print(interval.get_limits(img))

# set the detection threshold at the 3-sigma noise level, estimated using 
# the median absolute deviation (mad_std) of the image
# use DAOStarFinder to detect the stars in the image

bkg_sigma = mad_std(img)
daofind = DAOStarFinder(fwhm=5., threshold=4*bkg_sigma)
sources = daofind(img)

for col in sources.colnames:
    sources[col].info.format = '%.8g'  # for consistent table output

print(sources)

# using the list of source locations, compute the sum of pixel values in circular apertures with a radius of 4 pixels

positions = (sources['xcentroid'], sources['ycentroid'])
apertures = CircularAperture(positions, r=4.)  
phot_table = aperture_photometry(img, apertures)  

for col in phot_table.colnames:
    phot_table[col].info.format = '%.8g'  # for consistent table output

print(phot_table)

# create normalized image object (ZScale)
norm = ImageNormalize(img, interval=ZScaleInterval())

plt.figure(figsize=(9,12))
plt.imshow(img, cmap='gray_r', origin='lower', norm=norm)
apertures.plot(color='blue', lw=1.5, alpha=0.5)
plt.show()