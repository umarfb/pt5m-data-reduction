from astropy.io import fits
import numpy as np 
import os

# script to run bias subtraction

# file path for combined frames
filepath = '/local/php18ufb/backed_up_on_astro3/documents/mini_project/data20183009/combinedframes/'
homepath = '/local/php18ufb/backed_up_on_astro3/documents/mini_project/data20183009/'

# method to return list of all FITS files
def list_files(extension):                  
    f_list = []
    i=0
    for i in range(len(os.listdir(homepath))):
        if os.listdir(homepath)[i].endswith(extension):
            #print(os.listdir()[i])
            f_list.append(os.listdir(homepath)[i])
    return f_list

# method to return list of all FITS files in specific directory
def list_files_dir(directory, extension):                  
    f_list = []
    i=0
    for i in range(len(os.listdir(directory))):
        if os.listdir(directory)[i].endswith(extension):
            #print(os.listdir()[i])
            f_list.append(directory + os.listdir(directory)[i])
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

# method to subtract master bias frame
def bias_sub(biasFrame, frame):
    biasData = fits.getdata(biasFrame)
    fData = fits.getdata(frame)
    dataOut = fData - biasData                          # subtracts bias from frame
    return dataOut

masterBias = filepath + 'combined_bias_2.fits'               # loads master bias frame

file_list = list_files('.fits')                # defines list of all FITS files
print(len(file_list))
dark_list = frame_type('DARK', file_list)       # defines list of all dark frames
sky_list = frame_type('SKY', file_list)         # defines list of all flat-field frames
sci_list = frame_type('SCIENCE', file_list)     # gets list of all science frames 
    
sky_listB = filterSort(sky_list, 'B')           # defines list of flat-field frames by filter
sky_listV = filterSort(sky_list, 'V')
sky_listR = filterSort(sky_list, 'R')
sky_listI = filterSort(sky_list, 'I')

# create list of dark frames with 2x2 binning
dark_list2 = []
for item in dark_list:
    binpix = fits.getheader(item)['XBINNING']
    if binpix == 2:
        dark_list2.append(item)

print(len(dark_list2))

# subtract bias from dark frames
for dark in dark_list2:
    dark_bSub = bias_sub(masterBias, dark)
    f_hdr = fits.getheader(dark)
    f_hdr['HISTORY'] = 'Bias subtracted'            # adds HISTORY header to FITS file
    print('Bias subtracting ' + dark)
    os.remove(dark)
    fits.writeto(dark, dark_bSub, f_hdr)

# subtract bias from flat-fields
for sky in sky_list:
    sky_bSub = bias_sub(masterBias, sky)
    f_hdr = fits.getheader(sky)
    f_hdr['HISTORY'] = 'Bias subtracted'            # adds HISTORY header to FITS file
    print('Bias subtracting ' + sky)
    os.remove(sky)
    fits.writeto(sky, sky_bSub, f_hdr)

# subtract bias from science frames
for sci in sci_list:
    sci_bSub = bias_sub(masterBias, sci)
    f_hdr = fits.getheader(sci)
    f_hdr['HISTORY'] = 'Bias subtracted'            # adds HISTORY header to FITS file
    print('Bias subtracting ' + sci)
    os.remove(sci)
    fits.writeto(sci, sci_bSub, f_hdr)