from astropy.io import fits
import numpy as np
from scipy import stats
import os

# script to create pixel and mask and correct aligned images

filepath = '/local/php18ufb/backed_up_on_astro3/mini_project01/data20183009/'
#filename = 'B_aligned.fits'
#filename = 'r0500345.fits'

#align_img = fits.getdata(filepath + filename)

# method to return list of all FITS files
def list_files(extension):                  
    f_list = []
    i=0
    for i in range(len(os.listdir(filepath))):
        if os.listdir(filepath)[i].endswith(extension):
            #print(os.listdir()[i])
            f_list.append(os.listdir(filepath)[i])
    return f_list

# method to return list of specific frames by type
def frame_type(f_type, file_list):                     
    f_list = []
    for f in file_list:
        header = fits.getheader(f)['TYPE']
        if header == f_type:
            f_list.append(f)
    return f_list

# method to create list by filter
def filterSort(file_list, bandfilt):
    filt_list = []
    for f in file_list:
        header = fits.getheader(f)
        if header['FILTER'] == bandfilt:
            filt_list.append(f)
    return filt_list

# method to get dimensions of image
def getshape(imgdata):
    y = imgdata.shape[0]
    x = imgdata.shape[1]
    return x, y

# method to create an image mask
def createmask(imgdata, blankval):
    img_x, img_y = getshape(imgdata)

    # initialise new array with same dimensions as image, and set all values to 1
    mask = np.array([[1 for x in range(img_x)] for y in range(img_y)], dtype=np.float)

    # finds location of blank pixels in science image and sets corresponding pixel in mask to 0
    for i in range(img_y):
        for j in range(img_x):
            if imgdata[i][j] == blankval:
                mask[i][j] = 0
                
    return mask


# get list of all fits files
file_list = list_files('fits')
# get list of all science frames from list of fits files
sciframes = frame_type('SCIENCE', file_list)
# get list of all B-band science frames
Bsci_frm = filterSort(sciframes, 'B')

# list of mask frames
masklist = []

for f in Bsci_frm:
    f_hdr = fits.getheader(f)
    f_img = fits.getdata(f)

    # determines most recent change to fits file via HISTORY header
    hist = f_hdr['HISTORY'][len(f_hdr['HISTORY'])-1]

    # gets dimensions of image data
    img_x, img_y = getshape(f_img)

    if hist == 'Aligned':
        # create mask for aligned frame
        maskframe = createmask(f_img, 12.0)
        print(maskframe.shape, 1)
        masklist.append(maskframe)
    else:
        # create mask for non-aligned frame (uniform mask)
        maskframe = np.array([[1 for x in range(img_x)] for y in range(img_y)], dtype=np.float)
        print(maskframe.shape, 0)
        masklist.append(maskframe)

# add all masks to produce a 'master mask' frame
master_mask = np.sum(masklist, 0)
# divide maskframe by number of masks
master_mask = master_mask / (len(masklist))

fits.writeto('B_mask.fits', master_mask)