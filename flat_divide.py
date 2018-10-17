from astropy.io import fits
import numpy as np 
import os

# script to divide science frames by flats

# filepath for FITS files
homepath = '/local/php18ufb/backed_up_on_astro3/mini_project01/data20183009/'
# file path for combined frames
filepath = '/local/php18ufb/backed_up_on_astro3/mini_project01/data20183009/combinedframes/'

# list of filters
filters = ['B', 'V', 'R', 'I']

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
        header = fits.getheader(f)
        try:
            if header['TYPE'] == f_type:
                f_list.append(f)
        except KeyError:
            print('File has no type')
    return f_list

# method to create list by filter
def filterSort(file_list, bandfilt):
    filt_list = []
    for f in file_list:
        header = fits.getheader(f)
        if header['FILTER'] == bandfilt:
            filt_list.append(f)
    return filt_list

# assigns master flat-fields for each frame

flatB = filepath + 'combined_flatB.fits'
flatV = filepath + 'combined_flatV.fits'
flatR = filepath + 'combined_flatR.fits'
flatI = filepath + 'combined_flatI.fits'

flat_list = []
flat_list.append(flatB)
flat_list.append(flatV)
flat_list.append(flatR)
flat_list.append(flatI)

file_list = list_files('.fits')                # defines list of all FITS files


sci_list = frame_type('SCIENCE', file_list)     # gets list of all science frames

# sorts science frames by filter
sciB = filterSort(sci_list, 'B')
sciV = filterSort(sci_list, 'V')
sciR = filterSort(sci_list, 'R')
sciI = filterSort(sci_list, 'I')

# gets B filter mask frame
Bmask = fits.getdata(homepath + 'B_mask.fits')

# divide individual science frames by master flat-fields

sci_divB = []       # defines divided science frames for B filter
sci_divV = []       # defines divided science frames for V filter
sci_divR = []       # defines divided science frames for R filter
sci_divI = []       # defines divided science frames for I filter

#sciB.remove('r0500344.fits')         # removes badly aligned science frame in B filter
# divide B filter science frames by B flat fields
for fB in sciB:
    data_in = fits.getdata(fB)
    flat_data = fits.getdata(flatB)
    print('Dividing by flats for ' + fB)
    data_out = data_in / flat_data

    f_hdr = fits.getheader(fB)
    f_hdr['HISTORY'] = 'Divided by flat-field'

    sci_divB.append(data_out)

    #os.remove(fB)
    #fits.writeto(fB, data_out, f_hdr)
    #print(data_out)

combinedSciB = np.sum(sci_divB, axis=0)
#fits.writeto('B_final.fits', combinedSciB, fits.getheader(sciB[0]))

# creates mask divived science image
combinedSciB = combinedSciB / Bmask
fits.writeto('B_final_masked.fits', combinedSciB, fits.getheader(sciB[0]))

# divide V filter science frames by V flat fields
for fV in sciV:
    data_in = fits.getdata(fV)
    flat_data = fits.getdata(flatV)
    print('Dividing by flats for ' + fV)
    data_out = data_in / flat_data

    f_hdr = fits.getheader(fV)
    f_hdr['HISTORY'] = 'Divided by flat-field'

    sci_divV.append(data_out)

    #os.remove(fB)
    #fits.writeto(fB, data_out, f_hdr)
    #print(data_out)

combinedSciV = np.sum(sci_divV, axis=0)
#fits.writeto('V_final.fits', combinedSciV, fits.getheader(sciV[0]))

# divide R filter science frames by R flat fields
for fR in sciR:
    data_in = fits.getdata(fR)
    flat_data = fits.getdata(flatR)
    print('Dividing by flats for ' + fR)
    data_out = data_in / flat_data

    f_hdr = fits.getheader(fR)
    f_hdr['HISTORY'] = 'Divided by flat-field'

    sci_divR.append(data_out)

    #os.remove(fB)
    #fits.writeto(fB, data_out, f_hdr)
    #print(data_out)

combinedSciR = np.sum(sci_divR, axis=0)
#fits.writeto('R_final.fits', combinedSciR, fits.getheader(sciR[0]))

# divide I filter science frames by I flat fields
for fI in sciI:
    data_in = fits.getdata(fI)
    flat_data = fits.getdata(flatI)
    print('Dividing by flats for ' + fI)
    data_out = data_in / flat_data

    f_hdr = fits.getheader(fI)
    f_hdr['HISTORY'] = 'Divided by flat-field'

    sci_divI.append(data_out)

    #os.remove(fB)
    #fits.writeto(fB, data_out, f_hdr)
    #print(data_out)

combinedSciI = np.sum(sci_divI, axis=0)
#fits.writeto('I_final.fits', combinedSciI, fits.getheader(sciI[0]))