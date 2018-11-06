import numpy as np
import matplotlib.pyplot as plt
from astropy import wcs
from astropy import units as u
from astropy.io import fits
from astropy.io import ascii
from astropy.stats import SigmaClip
from astropy.coordinates import SkyCoord
from astropy.visualization import ZScaleInterval, ImageNormalize
from ucac4 import extract_star_array as get_local_ucac4

# method to find the nearest matching coordinate from a catalog
# - returns nearest matching coordinate
def find_nearest_coord(arr, val):
    arr_RA = np.asarray(arr.ra.deg)
    arr_Dec = np.asarray(arr.dec.deg)
    val_RA = np.asarray(val.ra.deg)
    val_Dec = np.asarray(val.dec.deg)

    RA_sep = np.abs(arr_RA - val_RA)
    Dec_sep = np.abs(arr_Dec - val_Dec)
    match_sep = np.sqrt((RA_sep ** 2) + (Dec_sep ** 2))

    idx = match_sep.argmin()
    return match_sep[idx]

path = '/local/php18ufb/backed_up_on_astro3/mini_project01/data20183009/'

filter_letter = input('Specify filter name: ')
filter_name = filter_letter

if filter_name.upper() in ('B', 'V', 'R', 'I'):
    filter_name = filter_name.upper()
else:
    raise ValueError('filter does not exist for pt5m')

filename = path + 'aperture_photometry/' + filter_name + '_final_sources.txt'
fitsname = path + filter_name + '_final.fits'

# read in list of image sources into astropy table
sources = ascii.read(filename)

# get pixel coordinates from pt5m image
x_pix = np.array(sources['xcenter'].data)
y_pix = np.array(sources['ycenter'].data)
pix_coords = np.vstack((x_pix, y_pix)).T

hdulist = fits.open(fitsname)   # Load FITS hdulist
w = wcs.WCS(hdulist[0].header)  # parse WCS keywords in primary HDU

# Convert pt5m pixel coordinates to world coordinates
world = w.wcs_pix2world(pix_coords, 1)
pt5m_coords = SkyCoord(world, frame='icrs', unit='deg')

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

# define 1 pix seperation in WCS
u_pix_dist = ((width/736) + (height/1092)) / 2

# find corresponding catalog stars and remove catalog stars outside pt5m fov
ucac4 = get_local_ucac4(ra, dec, width.deg, height.deg)
good_mags = (ucac4[filter_name] < 18)
ucac4 = ucac4[good_mags]
ucac4_coords = SkyCoord(ucac4['ra'], ucac4['dec'], unit='deg')  # gets RA and Dec of stars from catalog

# convert catalog stars RA and Dec to pixel coordinates
ucac4_pos = []
for pos in ucac4_coords:
    ucac4_pos.append([pos.ra.degree, pos.dec.degree])

ucac4_pix = w.wcs_world2pix(ucac4_pos, 1)

# auto matching function
mag_inst = - 2.5 * np.log10(sources['residual_aperture_sum'])   # get instrumental magnitudes

radius = 4 * u_pix_dist
idx, sep, _ = pt5m_coords.match_to_catalog_sky(ucac4_coords)
mask = sep < radius

delta_m = ucac4[filter_name][idx] - mag_inst
delta_m = delta_m[mask]
ntotal = len(delta_m)
ucac4_mag = ucac4[filter_name][idx]
ucac4_mag = ucac4_mag[mask]

# sigma clip data
sc = SigmaClip(sigma=3, iters=4)
filtered_delta_m = sc(delta_m)
mask = ~filtered_delta_m.mask

nrej = ntotal - np.sum(mask)

zp = np.ma.median(filtered_delta_m)     # sigma-clipped zero point
zp_err = filtered_delta_m.std()/np.sqrt(np.sum(mask))   # zero point error
print('zero point in ' + filter_name + ' band is: ' + str(zp))

counts = sources['residual_aperture_sum']
sources['mag'] = - 2.5 * np.log10(counts) + zp
print(sources['mag'])

# write magnitudes to file
file = open(path + '/aperture_photometry/' + filter_name + '.txt', 'w')
sources.write(file, format='csv')
file.close()

plt.scatter(ucac4_mag, delta_m, marker='.')
plt.plot([ucac4_mag.min() / 1.05, ucac4_mag.max() * 1.05], [zp, zp], color='r', linestyle='--', lw=1)
plt.ylim(delta_m.min() / 1.025, delta_m.max() * 1.025)
plt.xlim(ucac4_mag.min() / 1.025, ucac4_mag.max() * 1.025)
plt.xlabel('$m_{cat}$')
plt.ylabel('$\Delta m$')
plt.title(filter_letter.upper() + ' band')
plt.show()

'''
img_data = fits.getdata(fitsname)

# creates ZScale interval as used by DS9
interval = ZScaleInterval()
interval.get_limits(img_data)

# show matched catalog stars
norm = ImageNormalize(img_data, interval=ZScaleInterval())
plt.figure(figsize=(9,12))          # set figure size
#plt.xlim(1092)
#plt.ylim(736)
plt.imshow(img_data, cmap='Greys', origin='lower', norm=norm)
plt.scatter(ucac4_pix[:,0], ucac4_pix[:,1], color='g', lw=2, marker='o')
plt.scatter(sources['xcenter'], sources['ycenter'], color='r', lw=2, marker='x')
plt.show()
#plt.savefig('cal_skyplot.pdf', format='pdf', dpi=300, bbox_inches='tight')
'''