from astropy.io import fits
import numpy as np 
import os

# script to rebin frames of different bin sizes

# file path for combined frames
filepath = '/local/php18ufb/backed_up_on_astro3/documents/mini_project/data20183009/combinedframes/'

# method to return list of all FITS files
def list_files(extension):                  
    f_list = []
    i=0
    for i in range(len(os.listdir())):
        if os.listdir()[i].endswith(extension):
            #print(os.listdir()[i])
            f_list.append(os.listdir()[i])
    return f_list

# method to return list of specific frames by type
def frame_type(f_type, file_list):                     
    f_list = []
    for f in file_list:
        header = fits.getheader(f)
        if header['TYPE'] == f_type:
            f_list.append(f)
    return f_list

# method to create list by filter
def filterSort(file_list, bandfilt):
    filt_list = []
    i=0
    for f in file_list:
        header = fits.getheader(f)
        if header['FILTER'] == bandfilt:
            filt_list.append(f)
    return filt_list

# method to rebin pixels in FITS image
def rebin(arr, scale):
    oldDim = [arr.shape[0], arr.shape[1]]                   # gets dimensions of image to be rebinned
    newDim = [arr.shape[0]//scale, arr.shape[1]//scale]     # sets dimensions of new rebinned image
    
    new_arr = arr.reshape([newDim[0], oldDim[0]//newDim[0], newDim[1], oldDim[1]//newDim[1]]).sum(3).sum(1)     # rebins old image by summing underlying pixels
    return new_arr

file_list = list_files('.fits')                     # returns list of all FITS files in directory
sky_list = frame_type('SKY', file_list)           # returns list of all flat-field frames

print(len(sky_list))

sky_listB = filterSort(sky_list, 'B')           # defines list of flat-field frames by filter
sky_listV = filterSort(sky_list, 'V')
sky_listR = filterSort(sky_list, 'R')
sky_listI = filterSort(sky_list, 'I')
sky_listH = filterSort(sky_list, 'H')

# rebins 1x1 flat-fields in all filters to 2x2
for frame in sky_list:                             # change for each filter
    f_name = frame.rstrip('.fits')                  # gets filename without extension
    hdr = fits.getheader(frame)
    #print(hdr['XBINNING'])
    newFrame = rebin(fits.getdata(frame), 2)

    # updates FITS headers
    hdr['NAXIS1'] = hdr['NAXIS1']//2
    hdr['NAXIS2'] = hdr['NAXIS2']//2
    hdr['XBINNING'] = 2
    hdr['YBINNING'] = 2
    hdr['XPIXSZ'] = hdr['XPIXSZ'] * 2
    hdr['YPIXSZ'] = hdr['YPIXSZ'] * 2
    hdr['CRPIX1'] = hdr['CRPIX1']//2
    hdr['CRPIX2'] = hdr['CRPIX2']//2
    hdr['HISTORY'] = 'Rebinned to 2x2'

    print('Rebinning ' + f_name)
    os.remove(frame)
    fits.writeto(frame, newFrame, hdr)      # writes rebinned data to FITS file