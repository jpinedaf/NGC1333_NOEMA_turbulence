from spectral_cube import SpectralCube
import numpy as np
from astropy.io import fits
import astropy.units as u
import matplotlib.pyplot as plt 
plt.ion()



file_list_in = ['NGC1333_H13COp_L17-merged_fix.fits',
             	'NGC1333_HNC_L23-merged_fix.fits',
             	'NGC1333_HCN_L21-merged_fix.fits',
             	'ngc1333_c18o_3-2.fits',
             	'NGC1333_13CO_1-0.fits']

file_list_TdV = ['NGC1333_H13COp_L17-TdV.fits',
             	 'NGC1333_HNC_L23-TdV.fits',
             	 'NGC1333_HCN_L21-TdV.fits',
             	 'NGC1333_SE_C18O-TdV.fits',
             	 'NGC1333_13CO_1-0-TdV.fits']

# fits.setval('NGC1333_SE_C18O-TdV.fits', 'BMAJ', value=14.0/3600., before='SPECSYS')
# fits.setval('NGC1333_SE_C18O-TdV.fits', 'BMIN', value=14.0/3600., after='BMAJ')
# fits.setval('NGC1333_SE_C18O-TdV.fits', 'BPA', value=0.0, after='BMIN')
# fits.setval('ngc1333_c18o_3-2.fits', 'BMAJ', value=14.0/3600., before='SPECSYS')
# fits.setval('ngc1333_c18o_3-2.fits', 'BMIN', value=14.0/3600., after='BMAJ')
# fits.setval('ngc1333_c18o_3-2.fits', 'BPA', value=0.0, after='BMIN')
# fits.setval('ngc1333_c18o_3-2.fits', 'BUNIT', value='K', comment='Ta*', after='BZERO')


thre_list = [4e-3*u.Jy/u.beam,
             4e-3*u.Jy/u.beam, 
             4e-3*u.Jy/u.beam, 
             0.2*u.K, 
             0.2*u.K]


v_list_min = [5.3, 3.5, -5., 5.5, -1.]
v_list_max = [10., 11., 19., 10., 11.5]

for file_in, file_out, v_min, v_max, thre_i in zip(file_list_in, file_list_TdV, v_list_min, v_list_max, thre_list):
    cube = SpectralCube.read(file_in)
    sub_cube = (cube.spectral_slab(v_min * u.km / u.s, v_max * u.km / u.s)).with_spectral_unit(u.km / u.s)
    TdV = sub_cube.moment(order=0)
    TdV.write(file_out, overwrite=True)
    sub_cube_mask = sub_cube.with_mask(sub_cube > thre_i)
    Mom1 = sub_cube_mask.moment(order=1)
    Mom1.write(file_out.replace('TdV', 'Mom1'), overwrite=True)
    Tp = sub_cube_mask.max(axis=0)
    Tp.write(file_out.replace('TdV', 'Tp'), overwrite=True)


# Now remove bad velocity points
# and pad them with the mean of the velocity

# C18O
Mom1, hd = fits.getdata('NGC1333_SE_C18O-Mom1.fits', header=True)
TdV = fits.getdata('NGC1333_SE_C18O-TdV.fits', header=False)
bad = (TdV < 0.75)
Mom1[bad] = np.nan
fits.writeto('NGC1333_SE_C18O-Mom1_QA.fits', Mom1, hd, overwrite=True)
Mom1[bad] = np.nanmean(Mom1)
fits.writeto('NGC1333_SE_C18O-Mom1_QA_pad.fits', Mom1, hd, overwrite=True)

# H13CO+
Mom1, hd = fits.getdata('NGC1333_H13COp_L17-Mom1.fits', header=True)
TdV = fits.getdata('NGC1333_H13COp_L17-TdV.fits', header=False)
bad = (TdV < 100.0e-3)  | np.isnan(TdV)
Mom1[bad] = np.nan
fits.writeto('NGC1333_H13COp_L17-Mom1_QA.fits', Mom1, hd, overwrite=True)
Mom1[bad] = np.nanmean(Mom1)
fits.writeto('NGC1333_H13COp_L17-Mom1_QA_pad.fits', Mom1, hd, overwrite=True)

# HNC
Mom1, hd = fits.getdata('NGC1333_HNC_L23-Mom1.fits', header=True)
TdV = fits.getdata('NGC1333_HNC_L23-TdV.fits', header=False)
bad = (TdV < 300.0e-3) | np.isnan(TdV)
Mom1[bad] = np.nan
fits.writeto('NGC1333_HNC_L23-Mom1_QA.fits', Mom1, hd, overwrite=True)
Mom1[bad] = np.nanmean(Mom1)
fits.writeto('NGC1333_HNC_L23-Mom1_QA_pad.fits', Mom1, hd, overwrite=True)

# HCN
Mom1, hd = fits.getdata('NGC1333_HCN_L21-Mom1.fits', header=True)
TdV = fits.getdata('NGC1333_HCN_L21-TdV.fits', header=False)
bad = (TdV < 1000.0e-3) | np.isnan(TdV)
Mom1[bad] = np.nan
fits.writeto('NGC1333_HCN_L21-Mom1_QA.fits', Mom1, hd, overwrite=True)
Mom1[bad] = np.nanmean(Mom1)
fits.writeto('NGC1333_HCN_L21-Mom1_QA_pad.fits', Mom1, hd, overwrite=True)




# Now pad the TdV map
TdV, hd = fits.getdata('NGC1333_H13COp_L17-TdV.fits', header=True)
bad = np.isnan(TdV)
TdV[bad] = np.random.normal(scale=0.020, size=np.sum(bad))
fits.writeto('NGC1333_H13COp_L17-TdV_pad.fits', TdV, hd, overwrite=True)


# HNC
TdV, hd = fits.getdata('NGC1333_HNC_L23-TdV.fits', header=True)
bad = np.isnan(TdV)
TdV[bad] = np.random.normal(scale=0.030, size=np.sum(bad))
fits.writeto('NGC1333_HNC_L23-TdV_pad.fits', TdV, hd, overwrite=True)

# HCN
TdV, hd = fits.getdata('NGC1333_HCN_L21-TdV.fits', header=True)
bad = np.isnan(TdV)
TdV[bad] = np.random.normal(scale=0.090, size=np.sum(bad))
fits.writeto('NGC1333_HCN_L21-TdV_pad.fits', TdV, hd, overwrite=True)
