from spectral_cube import SpectralCube
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt 
plt.ion()

file_list_in = ['NGC1333_H13COp_L17-merged.fits',
                'NGC1333_HNC_L23-merged.fits',
                'NGC1333_HCN_L21-merged.fits']

file_list_out = ['NGC1333_H13COp_L17-merged_fix.fits',
                'NGC1333_HNC_L23-merged_fix.fits',
                'NGC1333_HCN_L21-merged_fix.fits']

file_list_pad = ['NGC1333_H13COp_L17-merged_fix_pad.fits',
                'NGC1333_HNC_L23-merged_fix_pad.fits',
                'NGC1333_HCN_L21-merged_fix_pad.fits']

for file_in, file_out in zip(file_list_in, file_list_out):
    cube = SpectralCube.read(file_in)
    sub_cube = cube.minimal_subcube()
    sub_cube.write(file_out, overwrite=True)

line_list = ['H13CO+', 'HNC', 'HCN']
i_min = [80, 100, 140]
i_max = [130, 170, 330]

for file_in, file_pad, line_i, ii, jj in zip(file_list_out, file_list_pad, line_list, i_min, i_max):
    cube, hd = fits.getdata(file_in, header=True)
    spec = np.nanmedian(cube, axis=[1,2])
    rms = np.round(np.nanstd(np.append(cube[:ii,:,:], cube[jj:,:,:])), 
                   decimals=3)
    plt.plot(spec)
    plt.text(40, 0.8*np.max(spec), line_i, horizontalalignment='right')
    plt.text(50, 0.8*np.max(spec), rms, horizontalalignment='left')
    print('rms for {0} = {1:.3f} Jy/beam'.format(line_i, rms))
    bad = np.isnan(cube)
    cube[bad] = np.random.normal(scale=rms, size=np.sum(bad))
    fits.writeto(file_pad, cube, hd, overwrite=True)


ii = 56
jj = 110
line_i = 'C18O'
file_in = 'ngc1333_c18o_3-2.fits'

cube, hd = fits.getdata(file_in, header=True)
spec = np.nanmedian(cube, axis=[1,2])
rms = np.round(np.nanstd(np.append(cube[:ii,:,:], cube[jj:,:,:])), 
               decimals=3)
plt.plot(spec)
plt.text(40, 0.8*np.max(spec), line_i, horizontalalignment='right')
plt.text(50, 0.8*np.max(spec), rms, horizontalalignment='left')
print('rms for {0} = {1:.3f} K'.format(line_i, rms))


###
# 13CO file
###
from spectral_cube import SpectralCube
import astropy.units as u
file_13co_raw = 'PerA_13coFCRAO_F_xyv.fits'
# cube, hd = fits.getdata(file_13co_raw, header=True)
# hd['RESTFRQ'] = hd['RESTFRQ']*1e9
# hd['CTYPE3'] = 'VRAD'
# hd['RESTFRQ'] = hd['LINEFREQ']
# hd['OBSERVER'] = 'COMPLETE'
# hd.remove('LINEFREQ')
# hd.add('SPECSYS', 'LSRK', 'Standard of rest for spectral axis')
# fits.writeto(file_13co_raw, cube, hd, overwrite=True)
fits.setval(file_13co_raw, 'SPECSYS', value='LSRK', after='EPOCH')
fits.setval(file_13co_raw, 'BMAJ', value=46.0/3600., before='SPECSYS')
fits.setval(file_13co_raw, 'BMIN', value=46.0/3600., after='BMAJ')
fits.setval(file_13co_raw, 'BPA', value=0.0, after='BMIN')
fits.setval(file_13co_raw, 'BEAM', 
    value='Beam: BMAJ=46.0 arcsec BMIN=46.0 arcsec BPA=0.0 deg', after='BPA')

cube = SpectralCube.read(file_13co_raw)
cube.allow_huge_operations = True

# create a target header to reproject to by making the pixel size 2 times larger
# target_header = cube.wcs.celestial[100:303, 668:843].to_header()
target_header = cube.header #wcs.celestial.to_header()
target_header['CTYPE1'] = 'RA---TAN'
target_header['CTYPE2'] = 'DEC--TAN'
target_header['NAXIS1'] = 181
target_header['NAXIS2'] = 191
target_header['CRVAL1'] = 52.2350623
target_header['CRVAL2'] = 31.3057279
target_header['CRPIX1'] = 90
target_header['CRPIX2'] = 95

ii = 100
jj = 330
line_i = '13CO'
file_in = 'NGC1333_13CO_1-0.fits'

sub_cube = cube.reproject(target_header)

sub_cube_kms = sub_cube.with_spectral_unit(u.km/u.s, 
                velocity_convention='radio', 
                rest_value=110.20135430*u.GHz) # Splatalogue

sub_cube_kms.write(file_in, overwrite=True)

cube, hd = fits.getdata(file_in, header=True)

spec = np.nanmedian(cube, axis=[1,2])
rms = np.round(np.nanstd(np.append(cube[:ii,:,:], cube[jj:,:,:])), 
               decimals=3)
plt.plot(spec)
plt.text(40, 0.8*np.max(spec), line_i, horizontalalignment='right')
plt.text(50, 0.8*np.max(spec), rms, horizontalalignment='left')
print('rms for {0} = {1:.3f} K'.format(line_i, rms))
