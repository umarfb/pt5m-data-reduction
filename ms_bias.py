from astropy.io import fits
import numpy as np 
import os

# script to run bias subtraction

# file path for combined frames
filepath = '/local/php18ufb/backed_up_on_astro3/documents/mini_project/data20183009/combinedframes/'

# list of filters
#filters = ['B', 'V', 'R', 'I']

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
        if header[10] == f_type:
            #print(f, header[10])
            f_list.append(f)
    return f_list

# method to create frame stack
def data_stack(file_list, pixelbin):
    data_stack = []
    hdr_stack = []
    for data_f in file_list:
        hdr = fits.getheader(data_f)
        binpix = hdr['XBINNING']     # identifies pixel binning in the frame
        if binpix == pixelbin:
            data_stack.append(fits.getdata(data_f))
            hdr_stack.append(hdr)
    return data_stack, hdr_stack

# method to median combine data stack
def medianComb(file_stack, pixelbin):
    medianData = np.median(file_stack, axis=0)
    return medianData

# method to create list by filter
def filterSort(file_list, bandfilt):
    filt_list = []
    for f in file_list:
        header = fits.getheader(f)
        if header['FILTER'] == bandfilt:
            filt_list.append(f)
    return filt_list

# method to scale image data by median and combine into stack
def medianStack(file_list, pixelbin):
    data_stack = []
    for f in file_list:
        binpix = fits.getheader(f)['XBINNING']     # identifies pixel binning in the frame
        if binpix == pixelbin:
            data = fits.getdata(f)
            data = data / np.median(data)
            data_stack.append(data)
    return data_stack


file_list = list_files('.fits')                # defines list of all FITS files
bias_list = frame_type('BIAS', file_list)       # defines list of all bias frames
dark_list = frame_type('DARK', file_list)       # defines list of all dark frames
sky_list = frame_type('SKY', file_list)         # defines list of all flat-field frames

sky_listB = filterSort(sky_list, 'B')           # defines list of flat-field frames by filter
sky_listV = filterSort(sky_list, 'V')
sky_listR = filterSort(sky_list, 'R')
sky_listI = filterSort(sky_list, 'I')

n=0
# create median combined bias frame
for n in range(3):
    bias_stack, bias_hdr = data_stack(bias_list, n+1)
    #print(bias_hdr[0])
    medianBias = medianComb(bias_stack, n+1)
    fits.writeto(filepath + 'combined_bias_' + str(n+1) + '.fits', medianBias, bias_hdr[0])

n=0
# create median combined dark frame
#for n in range(3):
#    dark_stack = data_stack(dark_list, n+1)
#    medianDark = medianComb(dark_stack, n+1)
#    fits.writeto(filepath + 'combined_dark_' + str(n+1) + '.fits', medianDark)

n=0
# create median combined flat-fields and normalise
flatB_stack = medianStack(sky_listB, 1)
flatV_stack = medianStack(sky_listV, 1)
flatR_stack = medianStack(sky_listR, 1)
flatI_stack = medianStack(sky_listI, 1)

"""

medianFlatB = np.median(flatB_stack, axis=0)
m = np.mean(medianFlatB)
medianFlatB = medianFlatB / m
fits.writeto(filepath + 'combined_flatB_1.fits', medianFlatB)

medianFlatV = np.median(flatV_stack, axis=0)
m = np.mean(medianFlatV)
medianFlatV = medianFlatV / m
fits.writeto(filepath + 'combined_flatV_1.fits', medianFlatV)

medianFlatR = np.median(flatR_stack, axis=0)
m = np.mean(medianFlatR)
medianFlatR = medianFlatR / m
fits.writeto(filepath + 'combined_flatR_1.fits', medianFlatR)

medianFlatI = np.median(flatI_stack, axis=0)
m = np.mean(medianFlatI)
medianFlatI = medianFlatI / m
fits.writeto(filepath + 'combined_flatI_1.fits', medianFlatI)

"""

