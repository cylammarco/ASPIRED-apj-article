from astropy.io import fits
from matplotlib import pyplot as plt
import numpy as np
from spectres import spectres


#   ___   __  __    ___     ___  
#  / __| |  \/  |  / _ \   / __| 
# | (_ | | |\/| | | (_) |  \__ \ 
#  \___| |_|__|_|  \___/   |___/ 
#_|"""""|_|"""""|_|"""""|_|"""""|
#"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'
#
fits_r400_1 = fits.open('gemini_gmos-gw170817/gw170817_r400_0_science_0.fits')[1]
fits_r400_2 = fits.open('gemini_gmos-gw170817/gw170817_r400_1_science_0.fits')[1]
fits_r400_3 = fits.open('gemini_gmos-gw170817/gw170817_r400_2_science_0.fits')[1]
fits_r400_4 = fits.open('gemini_gmos-gw170817/gw170817_r400_3_science_0.fits')[1]
fits_b600_1 = fits.open('gemini_gmos-gw170817/gw170817_b600_0_science_0.fits')[1]
fits_b600_2 = fits.open('gemini_gmos-gw170817/gw170817_b600_1_science_0.fits')[1]
fits_b600_3 = fits.open('gemini_gmos-gw170817/gw170817_b600_2_science_0.fits')[1]
fits_b600_4 = fits.open('gemini_gmos-gw170817/gw170817_b600_3_science_0.fits')[1]

wave_r400_1_start = fits_r400_1.header['CRVAL1']
wave_r400_1_bin = fits_r400_1.header['CDELT1']
n_wave_r400_1 = fits_r400_1.header['NAXIS1']
wave_r400_2_start = fits_r400_2.header['CRVAL1']
wave_r400_2_bin = fits_r400_2.header['CDELT1']
n_wave_r400_2 = fits_r400_2.header['NAXIS1']
wave_r400_3_start = fits_r400_3.header['CRVAL1']
wave_r400_3_bin = fits_r400_3.header['CDELT1']
n_wave_r400_3 = fits_r400_3.header['NAXIS1']
wave_r400_4_start = fits_r400_4.header['CRVAL1']
wave_r400_4_bin = fits_r400_4.header['CDELT1']
n_wave_r400_4 = fits_r400_4.header['NAXIS1']

wave_b600_1_start = fits_b600_1.header['CRVAL1']
wave_b600_1_bin = fits_b600_1.header['CDELT1']
n_wave_b600_1 = fits_b600_1.header['NAXIS1']
wave_b600_2_start = fits_b600_2.header['CRVAL1']
wave_b600_2_bin = fits_b600_2.header['CDELT1']
n_wave_b600_2 = fits_b600_2.header['NAXIS1']
wave_b600_3_start = fits_b600_3.header['CRVAL1']
wave_b600_3_bin = fits_b600_3.header['CDELT1']
n_wave_b600_3 = fits_b600_3.header['NAXIS1']
wave_b600_4_start = fits_b600_4.header['CRVAL1']
wave_b600_4_bin = fits_b600_4.header['CDELT1']
n_wave_b600_4 = fits_b600_4.header['NAXIS1']

wave_r400_1 = np.linspace(wave_r400_1_start, wave_r400_1_start + wave_r400_1_bin * n_wave_r400_1, n_wave_r400_1)
wave_r400_2 = np.linspace(wave_r400_2_start, wave_r400_2_start + wave_r400_2_bin * n_wave_r400_2, n_wave_r400_2)
wave_r400_3 = np.linspace(wave_r400_3_start, wave_r400_3_start + wave_r400_3_bin * n_wave_r400_3, n_wave_r400_3)
wave_r400_4 = np.linspace(wave_r400_4_start, wave_r400_4_start + wave_r400_4_bin * n_wave_r400_4, n_wave_r400_4)

wave_b600_1 = np.linspace(wave_b600_1_start, wave_b600_1_start + wave_b600_1_bin * n_wave_b600_1, n_wave_b600_1)
wave_b600_2 = np.linspace(wave_b600_2_start, wave_b600_2_start + wave_b600_2_bin * n_wave_b600_2, n_wave_b600_2)
wave_b600_3 = np.linspace(wave_b600_3_start, wave_b600_3_start + wave_b600_3_bin * n_wave_b600_3, n_wave_b600_3)
wave_b600_4 = np.linspace(wave_b600_4_start, wave_b600_4_start + wave_b600_4_bin * n_wave_b600_4, n_wave_b600_4)


wave_r400 = np.arange(4800, 9500, 0.5)
wave_b600 = np.arange(4000, 9000, 0.5)

flux_r400_1 = spectres(wave_r400, wave_r400_1, fits_r400_1.data)
flux_r400_2 = spectres(wave_r400, wave_r400_2, fits_r400_2.data)
flux_r400_3 = spectres(wave_r400, wave_r400_3, fits_r400_3.data)
flux_r400_4 = spectres(wave_r400, wave_r400_4, fits_r400_4.data)

flux_b600_1 = spectres(wave_b600, wave_b600_1, fits_b600_1.data)
flux_b600_2 = spectres(wave_b600, wave_b600_2, fits_b600_2.data)
flux_b600_3 = spectres(wave_b600, wave_b600_3, fits_b600_3.data)
flux_b600_4 = spectres(wave_b600, wave_b600_4, fits_b600_4.data)



flux_r400 = np.mean((flux_r400_1, flux_r400_2, flux_r400_3, flux_r400_4), axis=0) / 5.
flux_b600 = np.mean((flux_b600_1, flux_b600_2, flux_b600_3, flux_b600_4), axis=0) / 10.

# trim the last few hundred A from the blue and the first few hundred A from
# the red in the combined spectrum
r400_limit = 5000
# This value has to match the
b600_limit = 6000

r400_mask = (wave_r400 > r400_limit) & (wave_r400 < b600_limit)
b600_mask = (wave_b600 >= r400_limit) & (wave_b600 <= b600_limit)

# resample the r400 to match b600 resolution
flux_r400_resampled = spectres(
    wave_b600[b600_mask],
    wave_r400[r400_mask],
    flux_r400[r400_mask],
)

flux_averaged = (flux_r400_resampled + flux_b600[b600_mask])/2.


flux_combined = np.concatenate(
    (
        flux_b600[wave_b600 < r400_limit],
        flux_averaged,
        flux_r400[wave_r400 > b600_limit],
    )
)
wave_combined = np.concatenate(
    (wave_b600[wave_b600 <= b600_limit], wave_r400[wave_r400 > b600_limit])
)


public_data_1 = np.genfromtxt('gemini_gmos-gw170817/AT2017gfo_2017-08-20_01-08-00_Gemini-S_GMOS-S_None.dat')
public_data_2 = np.genfromtxt('gemini_gmos-gw170817/AT2017gfo_2017-08-21_00-15-00_Gemini-S_GMOS-S_None.dat')


plt.figure(1)
plt.clf()
plt.plot(wave_combined, flux_combined + 5e-17)
plt.plot(public_data_2[:,0], public_data_2[:,1])
plt.ion()
plt.show()

