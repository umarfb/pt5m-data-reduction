import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy import wcs
from astropy import units as u
from astropy.io import fits
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.visualization import ZScaleInterval, ImageNormalize
from ucac4 import extract_star_array as get_local_ucac4

path = '/local/php18ufb/backed_up_on_astro3/mini_project01/data20183009/'

filter_name = input('Specify filter name: ')
if filter_name.upper() in ('B', 'V', 'R', 'I'):
    filter_name = filter_name.upper()
else:
    raise ValueError('filter does not exist for pt5m')

filename = path + 'aperture_photometry/' + filter_name + '_final_sources.txt'
fitsname = path + filter_name + '_final.fits'

# read in list of source data
sources = pd.read_csv(filename)
source_id = sources['id'].values

# get pixel coordinates from pt5m image
pix_coords = sources[['xcenter', 'ycenter']].values

# Load FITS hdulist
hdulist = fits.open(fitsname)

# parse WCS keywords in primary HDU
w = wcs.WCS(hdulist[0].header)

# determine exposure time
t_exp = hdulist[0].header['EXPTIME']

# Convert pixel coordinates to world coordinates
world = w.wcs_pix2world(pix_coords, 1)
wd_coords = SkyCoord(world, frame='icrs', unit='deg')

# get counts from pt5m image and calculate instrumental magnitude
img_counts = sources['residual_aperture_sum'].values
inst_mag = - 2.5 * np.log10(img_counts/t_exp)

if filter_name.upper() in ('B', 'V'):
    filter_name = 'mag_' + filter_name.upper()
elif filter_name.lower() in ('r', 'i'):
    filter_name = 'mag_' + filter_name.lower()
else:
    raise ValueError('cannot calibrate filter {} against UCAC4'.format(filter_name))

# specify centre of field, and dimensions
ra = hdulist[0].header['CRVAL1']
dec = hdulist[0].header['CRVAL2']
lim_coords = w.wcs_pix2world([[0,0], [1092,736]], 1)
lim_coords = SkyCoord(lim_coords, frame='icrs', unit='deg')
width = (lim_coords.ra.max() - lim_coords.ra.min())
height = (lim_coords.dec.max() - lim_coords.dec.min())

# find corresponding catalog stars and remove catalog stars outside pt5m fov
ucac4 = get_local_ucac4(ra, dec, width.deg, height.deg)
good_mags = (ucac4[filter_name] < 18)
ucac4 = ucac4[good_mags]

ucac4_rm_idx = []
for i in range(len(ucac4)):
    star = ucac4[i]
    pix_pos = w.wcs_world2pix([[star['ra'], star['dec']]], 1)
    x_pix = pix_pos[0,0]
    y_pix = pix_pos[0,1]
    
    if x_pix > 0 and x_pix < 1092:
        pass
    else:
        ucac4_rm_idx.append(i)

ucac4.remove_rows(ucac4_rm_idx)
#print(len(ucac4))

# convert catalog stars RA and Dec to pixel coordinates
ucac4_coords = SkyCoord(ucac4['ra'], ucac4['dec'], unit='deg')
ucac4_pos = []

for pos in ucac4_coords:
    ucac4_pos.append([pos.ra.degree, pos.dec.degree])

ucac4_pix = w.wcs_world2pix(ucac4_pos, 1)

'''

img_data = fits.getdata(fitsname)

# creates ZScale interval as used by DS9
interval = ZScaleInterval()
interval.get_limits(img_data)

# show matched catalog stars
norm = ImageNormalize(img_data, interval=ZScaleInterval())
plt.figure(figsize=(9,12))          # set figure size
plt.imshow(img_data, cmap='Greys', origin='lower', norm=norm)
plt.scatter(ucac4_pix[:,0], ucac4_pix[:,1])
plt.show()

idx, sep, _ = wd_coords.match_to_catalog_sky(ucac4_coords)

# calculate magnitude difference
delta_m = ucac4[filter_name][idx] - inst_mag
#print(inst_mag)
#print(np.median(delta_m), np.std(delta_m)/np.sqrt(len(delta_m)))

# define instrument zero point
zp = np.median(delta_m)

m_obj = inst_mag + zp

for j in range(len(source_id)):
    print(source_id[i], '\t', int(pix_coords[j][0]), int(pix_coords[j][1]), '\t', np.round(m_obj[j], 2))

plt.scatter(ucac4[filter_name][idx], delta_m)
plt.title('R band')
plt.xlabel('$m_{cat}$')
plt.ylabel('$\Delta m$')
plt.savefig('cal_plot.pdf', format='pdf', bbox_inches='tight')
plt.show()
'''