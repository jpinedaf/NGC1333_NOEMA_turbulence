from spectral_cube import SpectralCube
import astropy.units as u
import numpy as np
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.io import fits
import matplotlib.pyplot as plt 
import radio_beam



plt.rcParams.update({"text.usetex": True,
                     "font.family": "serif",
                     'xtick.direction': 'in',
                     'ytick.direction': 'in',
                     'xtick.major.size': 6,
                     'ytick.major.size': 6,
                     'xtick.minor.size': 3,
                     'ytick.minor.size': 3})
plt.ion()


def make_TdV(file_in, line='HNC'):
    """
    Function to determine integrated intensity map and to save it.
    """
    if line == 'HNC':
        v_min = 3.5
        v_max = 11.0
    elif line == 'H13CO+':
        v_min = 5.3
        v_max = 10.0
    elif line == 'C18O':
        v_min = 5.5
        v_max = 10.0
    elif line == '13CO':
        v_min = -1.0
        v_max = 11.5
    else:
        return
    file_out = file_in.replace('.fits', '_TdV.fits')
    cube = SpectralCube.read(file_in)
    sub_cube = (cube.spectral_slab(v_min * u.km / u.s, v_max * u.km / u.s)).with_spectral_unit(u.km / u.s)
    TdV = sub_cube.moment(order=0)
    TdV.write(file_out, overwrite=True)

# file names needed to 
file_list_in = ['NGC1333_H13COp_L17-merged.fits',
                'NGC1333_HNC_L23-merged.fits']
file_list_out = ['NGC1333_H13COp_L17-merged_fix.fits',
                 'NGC1333_HNC_L23-merged_fix.fits']
file_list_pad = ['NGC1333_H13COp_L17-merged_fix_pad.fits',
                 'NGC1333_HNC_L23-merged_fix_pad.fits']

for file_in, file_out in zip(file_list_in, file_list_out):
    # Load interferometric cubes and reduce to smalles possible size
    # while also circularizing the beam to 5arcsec
    cube = SpectralCube.read(file_in)
    cube.allow_huge_operations = True
    # hd_target = cube.header
    # target_beam = radio_beam.Beam(major=hd_target['BMAJ']*u.deg, 
    #                               minor=hd_target['BMAJ']*u.deg, pa=0*u.deg)
    target_beam = radio_beam.Beam(major=5*u.arcsec, 
                                  minor=5*u.arcsec, pa=0*u.deg)
    sub_cube = (cube.minimal_subcube()).convolve_to(target_beam)
    sub_cube.write(file_out, overwrite=True)

line_list = ['H13CO+', 'HNC']#
i_min = [80, 100]
i_max = [130, 170]
SIGMA_TO_FWHM = np.sqrt(8*np.log(2))


    
for file_in, file_pad, line_i, ii, jj in zip(file_list_out, file_list_pad, line_list, i_min, i_max):
    cube, hd = fits.getdata(file_in, header=True)
    make_TdV(file_in, line=line_i)
    cube_size = cube.shape
    spec = np.nanmedian(cube, axis=[1,2])
    rms = np.round(np.nanstd(np.append(cube[:ii,:,:], cube[jj:,:,:])), 
                   decimals=3)
    plt.plot(spec)
    plt.text(40, 0.8*np.max(spec), line_i, horizontalalignment='right')
    plt.text(50, 0.8*np.max(spec), rms, horizontalalignment='left')
    print('rms for {0} = {1:.3f} Jy/beam'.format(line_i, rms))
    # obtain pixels outside the coverage
    bad = np.isnan(cube)
    # define beam kernel as 
    kernel = Gaussian2DKernel(x_stddev=np.abs(hd['BMAJ']/ 
                                              (hd['CDELT1'] * SIGMA_TO_FWHM)))
    rms_cube = np.random.normal(scale=rms, size=cube_size)
    # Stella'a trick
    for i in range(cube_size[0]):
        slice = rms_cube[i, :, :]
        slice_sm = convolve(slice, kernel)
        rms_cube[i, :, :] = slice_sm
    renorm = rms / np.std(rms_cube)
    rms_cube = rms_cube * renorm
    cube[bad] = rms_cube[bad] #np.random.normal(scale=rms, size=np.sum(bad))
    fits.writeto(file_pad, cube, hd, overwrite=True)
    make_TdV(file_pad, line=line_i)


###
# C18O file
###
ii = 56
jj = 110
line_i = 'C18O'
Delta_V = 10 - 5.5
# file_in = 'ngc1333_c18o_3-2.fits'
file_in = 'NGC1333_SE_C18O.fits'
make_TdV(file_in, line=line_i)

cube, hd = fits.getdata(file_in, header=True)
spec = np.nanmedian(cube, axis=[1,2])
rms = np.round(np.nanstd(np.append(cube[:ii,:,:], cube[jj:,:,:])), 
               decimals=3)
plt.plot(spec)
plt.text(40, 0.8*np.max(spec), line_i, horizontalalignment='right')
plt.text(50, 0.8*np.max(spec), rms, horizontalalignment='left')
print('rms for {0} = {1:.3f} K'.format(line_i, rms))

# print(np.median(rms_map))
# Median value of 0.185 K
rms_map = np.round(np.nanstd(np.vstack([cube[:ii,:,:], cube[jj:,:,:]]), 
                axis=0), decimals=2)
rms_TdV = rms_map * np.sqrt(np.abs(hd['CDELT3']) * Delta_V)

###
# 13CO file
###
# Load cube and fix header
file_13co_raw = 'PerA_13coFCRAO_F_xyv.fits'
fits.setval(file_13co_raw, 'SPECSYS', value='LSRK', after='EPOCH')
fits.setval(file_13co_raw, 'BMAJ', value=46.0/3600., before='SPECSYS')
fits.setval(file_13co_raw, 'BMIN', value=46.0/3600., after='BMAJ')
fits.setval(file_13co_raw, 'BPA', value=0.0, after='BMIN')
fits.setval(file_13co_raw, 'BEAM', 
    value='Beam: BMAJ=46.0 arcsec BMIN=46.0 arcsec BPA=0.0 deg', after='BPA')

cube = SpectralCube.read(file_13co_raw)
cube.allow_huge_operations = True

# create a target header to reproject
# this is to focus on the NGC 1333 region
target_header = cube.header
target_header['CTYPE1'] = 'RA---TAN'
target_header['CTYPE2'] = 'DEC--TAN'
target_header['NAXIS1'] = 181
target_header['NAXIS2'] = 191
target_header['CRVAL1'] = 52.2350623
target_header['CRVAL2'] = 31.3057279
target_header['CRPIX1'] = 90
target_header['CRPIX2'] = 95

# channel numbers for velocity channel range used to calculate rms
ii = 100
jj = 330
line_i = '13CO'
# file name for output file
file_out_C13O = 'NGC1333_13CO_1-0.fits'
sub_cube = cube.reproject(target_header)
# convert to km/s and save as FITS file
sub_cube_kms = sub_cube.with_spectral_unit(u.km/u.s, 
                velocity_convention='radio', 
                rest_value=110.20135430*u.GHz) # Splatalogue
sub_cube_kms.write(file_out_C13O, overwrite=True)
make_TdV(file_out_C13O, line=line_i)
# Reload data cube
cube, hd = fits.getdata(file_out_C13O, header=True)

spec = np.nanmedian(cube, axis=[1,2])
rms = np.round(np.nanstd(np.append(cube[:ii,:,:], cube[jj:,:,:])), 
               decimals=3)
plt.plot(spec)
plt.text(40, 0.8*np.max(spec), line_i, horizontalalignment='right')
plt.text(50, 0.8*np.max(spec), rms, horizontalalignment='left')
print('rms for {0} = {1:.3f} K'.format(line_i, rms))

#
# Prepare for stacked spectra
# 
file_list_spec = ['NGC1333_13CO_1-0.fits',
                'NGC1333_SE_C18O.fits',
                'NGC1333_H13COp_L17-merged_fix.fits',
                'NGC1333_HNC_L23-merged_fix.fits']
label_list = [r'$^{13}$CO (1--0)',
            r'C$^{18}$O (3--2)',
            r'H$^{13}$CO$^+$ (1--0)',
            r'HNC (1--0)']

plt.close('all')
fig, ax = plt.subplots(figsize=(5, 4))

for file_i, label_i in zip(file_list_spec, label_list):
    cube_i = SpectralCube.read(file_i).to(u.K)
    spec_i = cube_i.median(axis=(1, 2))
    vel_i = cube_i.spectral_axis
    ax.plot(vel_i.to(u.km/u.s).value, spec_i.value, label=label_i, drawstyle='steps-mid')
ax.set_xlim(-5, 20)
plt.legend(labelcolor='linecolor', frameon=False)
ax.set_xlabel(r'Velocity, $V_{LSR}$ (km s$^{-1}$)')
ax.set_ylabel(r'Brightness (K)')

text_pspec = 'Median Spectra'# in the maps'
ax.text(0.05, 0.85, text_pspec, horizontalalignment='left', transform=ax.transAxes, size=12)

ax.yaxis.set_ticks(np.arange(0, 2, 0.5))

fig.savefig('../figs/NGC1333_Median_Spectra.pdf', bbox_inches='tight')
