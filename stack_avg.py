from astropy.io import fits
import numpy as np 
import os

# script to combine frames into a stack and produce a median combined frame

# filepath for FITS files
homepath = '/local/php18ufb/backed_up_on_astro3/documents/mini_project/data20183009/'
# file path for combined frames
filepath = '/local/php18ufb/backed_up_on_astro3/documents/mini_project/data20183009/combinedframes/'

# method to return list of all FITS files
def list_files(extension):                  
    f_list = []
    i=0
    for i in range(len(os.listdir(homepath))):
        if os.listdir(homepath)[i].endswith(extension):
            #print(os.listdir()[i])
            f_list.append(os.listdir(homepath)[i])
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
    for f in file_list:
        header = fits.getheader(f)
        if header['FILTER'] == bandfilt:
            filt_list.append(f)
    return filt_list

# method to create frame stack, specifying the pixel binning
def data_stack(file_list, pixelbin):
    data_stack = []
    hdr_stack = []
    for data_f in file_list:
        hdr = fits.getheader(data_f)
        binpix = fits.getheader(data_f)['XBINNING']     # identifies pixel binning in the frame
        if binpix == pixelbin:
            data_stack.append(fits.getdata(data_f))
            hdr_stack.append(hdr)
    return data_stack, hdr_stack

# method to scale image data by median and combine into stack
def medianStack(file_list):
    data_stack = []
    hdr_stack = []
    for f in file_list:
        data = fits.getdata(f)
        data = data / np.median(data)
        data_stack.append(data)
        hdr_stack.append(fits.getheader(f))
    return data_stack, hdr_stack

file_list = list_files('.fits')                # defines list of all FITS files
bias_list = frame_type('BIAS', file_list)       # defines list of all bias frames
dark_list = frame_type('DARK', file_list)       # defines list of all dark frames
sky_list = frame_type('SKY', file_list)         # defines list of all flat-field frames

darkStack, dark_hdr = data_stack(dark_list, 2)  # creates stack of 2x2 binned dark frames

sky_listB = filterSort(sky_list, 'B')           # defines list of flat-field frames by filter
sky_listV = filterSort(sky_list, 'V')
sky_listR = filterSort(sky_list, 'R')
sky_listI = filterSort(sky_list, 'I')
sky_listH = filterSort(sky_list, 'H')

# create median combined flat-fields and normalise
flatB_stack, flatB_hdr = medianStack(sky_listB)
flatV_stack, flatV_hdr = medianStack(sky_listV)
flatR_stack, flatR_hdr = medianStack(sky_listR)
flatI_stack, flatI_hdr = medianStack(sky_listI)

# creates combined median and normalised flat-fields

medianFlatB = np.median(flatB_stack, axis=0)
mB = np.mean(medianFlatB)
medianFlatB = medianFlatB / mB
fits.writeto(filepath + 'combined_flatB.fits', medianFlatB, flatB_hdr[0])

medianFlatV = np.median(flatV_stack, axis=0)
mV = np.mean(medianFlatV)
medianFlatV = medianFlatV / mV
fits.writeto(filepath + 'combined_flatV.fits', medianFlatV, flatV_hdr[0])

medianFlatR = np.median(flatR_stack, axis=0)
mR = np.mean(medianFlatR)
medianFlatR = medianFlatR / mR
fits.writeto(filepath + 'combined_flatR.fits', medianFlatR, flatR_hdr[0])

medianFlatI = np.median(flatI_stack, axis=0)
mI = np.mean(medianFlatI)
medianFlatI = medianFlatI / mI
fits.writeto(filepath + 'combined_flatI.fits', medianFlatI, flatI_hdr[0])

# create median combined dark stack
#print('Creating master dark frame')
#medianDark = np.median(darkStack, axis=0)
#fits.writeto(filepath + 'combined_dark.fits', medianDark, dark_hdr[0])