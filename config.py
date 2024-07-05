import astropy.units as u

# Config file keeps the filenames across the different python 
# files and notebooks
#
# Integrated intensity maps
#
file_TdV_13co = 'data/NGC1333_13CO_1-0_TdV.fits'
file_TdV_c18o = 'data/NGC1333_SE_C18O_TdV.fits'
file_TdV_h13cop_pad = 'data/NGC1333_H13COp_L17-merged_fix_pad_TdV.fits'
file_TdV_h13cop = 'data/NGC1333_H13COp_L17-merged_fix_TdV.fits'
file_TdV_hnc_pad = 'data/NGC1333_HNC_L23-merged_fix_pad_TdV.fits'
file_TdV_hnc = 'data/NGC1333_HNC_L23-merged_fix_TdV.fits'
#
# Pickle files with each power-spectrum
#
file_pickle_13co = 'data/powerspec_NGC1333_13CO_1-0_TdV.pickle'
file_pickle_13co_apod = 'data/powerspec_NGC1333_13CO_1-0_TdV_apod.pickle'
file_pickle_c18o = 'data/powerspec_NGC1333_SE_C18O_TdV.pickle'
file_pickle_c18o_apod = 'data/powerspec_NGC1333_SE_C18O-TdV_apod.pickle'
file_pickle_h13cop_pad = 'data/powerspec_NGC1333_H13COp_L17-merged_fix_pad_TdV.pickle'
file_pickle_h13cop = 'data/powerspec_NGC1333_H13COp_L17-merged_fix_TdV.pickle'
file_pickle_hnc_pad = 'data/powerspec_NGC1333_HNC_L23-merged_fix_pad_TdV.pickle'
file_pickle_hnc = 'data/powerspec_NGC1333_HNC_L23-merged_fix_TdV.pickle'

distance = 300. * u.pc

file_EMCEE_All = 'data/EMCEE_samples_All.h5'
