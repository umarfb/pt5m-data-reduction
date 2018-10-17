from astropy.io import fits
import numpy as np 
import os

# script to run dark subtraction

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

# method to subtract master dark frame
def dark_sub(darkFrame, frame, scale):
    darkData = fits.getdata(darkFrame)
    fData = fits.getdata(frame)
    dataOut = fData - (darkData//scale)                          # subtracts bias from frame
    return dataOut

# method to get scaling factor of exposure time
def getScale(darkFrame, frame):
    t_dark = fits.getheader(darkFrame)['EXPTIME']
    t_frame = fits.getheader(frame)['EXPTIME']
    scale = t_dark/t_frame
    # print(t_dark, t_frame)
    return scale

masterDark = filepath + 'combined_dark.fits'               # loads master dark frame

file_list = list_files('.fits')                # defines list of all FITS files
sky_list = frame_type('SKY', file_list)         # defines list of all flat-field frames
sci_list = frame_type('SCIENCE', file_list)     # gets list of all science frames

# subtract dark from flat-fields
for sky in sky_list:
    scale = getScale(masterDark, sky)
    print('Dark subtracting ' + sky)
    sky_dSub = dark_sub(masterDark, sky, scale)
    sky_hdr = fits.getheader(sky)
    sky_hdr['HISTORY'] = 'Dark subtracted'

    os.remove(sky)
    fits.writeto(sky, sky_dSub, sky_hdr)
    #print(sky_dSub)

# subtract dark from science frames
for sci in sci_list:
    scale = getScale(masterDark, sci)
    print('Dark subtracting ' + sci)
    sci_dSub = dark_sub(masterDark, sci, scale)
    sci_hdr = fits.getheader(sci)
    sci_hdr['HISTORY'] = ' Dark subtracted'

    os.remove(sci)
    fits.writeto(sci, sci_dSub, sci_hdr)