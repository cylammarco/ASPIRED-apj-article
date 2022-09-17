from astropy.io import fits
from matplotlib import pyplot as plt
import numpy as np
from spectres import spectres

plt.ion()

#   ___   __  __    ___     ___
#  / __| |  \/  |  / _ \   / __|
# | (_ | | |\/| | | (_) |  \__ \
#  \___| |_|__|_|  \___/   |___/
# _|"""""|_|"""""|_|"""""|_|"""""|
# "`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'
#
fits_r400_1 = fits.open("gemini_gmos-gw170817/gw170817_r400_0_science_0.fits")[
    1
]
fits_r400_2 = fits.open("gemini_gmos-gw170817/gw170817_r400_1_science_0.fits")[
    1
]
fits_r400_3 = fits.open("gemini_gmos-gw170817/gw170817_r400_2_science_0.fits")[
    1
]
fits_r400_4 = fits.open("gemini_gmos-gw170817/gw170817_r400_3_science_0.fits")[
    1
]
fits_b600_1 = fits.open("gemini_gmos-gw170817/gw170817_b600_0_science_0.fits")[
    1
]
fits_b600_2 = fits.open("gemini_gmos-gw170817/gw170817_b600_1_science_0.fits")[
    1
]
fits_b600_3 = fits.open("gemini_gmos-gw170817/gw170817_b600_2_science_0.fits")[
    1
]
fits_b600_4 = fits.open("gemini_gmos-gw170817/gw170817_b600_3_science_0.fits")[
    1
]

wave_r400_1_start = fits_r400_1.header["CRVAL1"]
wave_r400_1_bin = fits_r400_1.header["CDELT1"]
n_wave_r400_1 = fits_r400_1.header["NAXIS1"]
wave_r400_2_start = fits_r400_2.header["CRVAL1"]
wave_r400_2_bin = fits_r400_2.header["CDELT1"]
n_wave_r400_2 = fits_r400_2.header["NAXIS1"]
wave_r400_3_start = fits_r400_3.header["CRVAL1"]
wave_r400_3_bin = fits_r400_3.header["CDELT1"]
n_wave_r400_3 = fits_r400_3.header["NAXIS1"]
wave_r400_4_start = fits_r400_4.header["CRVAL1"]
wave_r400_4_bin = fits_r400_4.header["CDELT1"]
n_wave_r400_4 = fits_r400_4.header["NAXIS1"]

wave_b600_1_start = fits_b600_1.header["CRVAL1"]
wave_b600_1_bin = fits_b600_1.header["CDELT1"]
n_wave_b600_1 = fits_b600_1.header["NAXIS1"]
wave_b600_2_start = fits_b600_2.header["CRVAL1"]
wave_b600_2_bin = fits_b600_2.header["CDELT1"]
n_wave_b600_2 = fits_b600_2.header["NAXIS1"]
wave_b600_3_start = fits_b600_3.header["CRVAL1"]
wave_b600_3_bin = fits_b600_3.header["CDELT1"]
n_wave_b600_3 = fits_b600_3.header["NAXIS1"]
wave_b600_4_start = fits_b600_4.header["CRVAL1"]
wave_b600_4_bin = fits_b600_4.header["CDELT1"]
n_wave_b600_4 = fits_b600_4.header["NAXIS1"]

wave_r400_1 = np.linspace(
    wave_r400_1_start,
    wave_r400_1_start + wave_r400_1_bin * n_wave_r400_1,
    n_wave_r400_1,
)
wave_r400_2 = np.linspace(
    wave_r400_2_start,
    wave_r400_2_start + wave_r400_2_bin * n_wave_r400_2,
    n_wave_r400_2,
)
wave_r400_3 = np.linspace(
    wave_r400_3_start,
    wave_r400_3_start + wave_r400_3_bin * n_wave_r400_3,
    n_wave_r400_3,
)
wave_r400_4 = np.linspace(
    wave_r400_4_start,
    wave_r400_4_start + wave_r400_4_bin * n_wave_r400_4,
    n_wave_r400_4,
)

wave_b600_1 = np.linspace(
    wave_b600_1_start,
    wave_b600_1_start + wave_b600_1_bin * n_wave_b600_1,
    n_wave_b600_1,
)
wave_b600_2 = np.linspace(
    wave_b600_2_start,
    wave_b600_2_start + wave_b600_2_bin * n_wave_b600_2,
    n_wave_b600_2,
)
wave_b600_3 = np.linspace(
    wave_b600_3_start,
    wave_b600_3_start + wave_b600_3_bin * n_wave_b600_3,
    n_wave_b600_3,
)
wave_b600_4 = np.linspace(
    wave_b600_4_start,
    wave_b600_4_start + wave_b600_4_bin * n_wave_b600_4,
    n_wave_b600_4,
)


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


flux_r400 = (
    np.mean((flux_r400_1, flux_r400_2, flux_r400_3, flux_r400_4), axis=0) / 5.0
)
flux_b600 = (
    np.mean((flux_b600_1, flux_b600_2, flux_b600_3, flux_b600_4), axis=0)
    / 10.0
)

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

flux_averaged = (flux_r400_resampled + flux_b600[b600_mask]) / 2.0


gmos_flux_combined = np.concatenate(
    (
        flux_b600[wave_b600 < r400_limit],
        flux_averaged,
        flux_r400[wave_r400 > b600_limit],
    )
)
gmos_wave_combined = np.concatenate(
    (wave_b600[wave_b600 <= b600_limit], wave_r400[wave_r400 > b600_limit])
)


gmos_public_data_2 = np.genfromtxt(
    "gemini_gmos-gw170817/AT2017gfo_2017-08-21_00-15-00_Gemini-S_GMOS-S_None.dat"
)


gmos_flux_combined = spectres(
    gmos_public_data_2[:, 0], gmos_wave_combined, gmos_flux_combined
)

#   ___      ___    ___     ___    _____
#  / __|    | _ \  | _ \   /   \  |_   _|
#  \__ \    |  _/  |   /   | - |    | |
#  |___/   _|_|_   |_|_\   |_|_|   _|_|_
# _|"""""|_| """ |_|"""""|_|"""""|_|"""""|
# "`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'
#


sprat_fits = fits.open("lt-sprat-doaql/reduced_science_0.fits")[1]
sprat_data = sprat_fits.data
sprat_wave = np.linspace(
    sprat_fits.header["CRVAL1"],
    sprat_fits.header["CRVAL1"]
    + sprat_fits.header["CDELT1"] * sprat_fits.header["NAXIS1"],
    sprat_fits.header["NAXIS1"],
)

sprat_public_fits_1 = fits.open(
    "lt-sprat-doaql/v_e_20170618_20_1_0_2.fits.gz"
)[-1]
sprat_public_fits_2 = fits.open(
    "lt-sprat-doaql/v_e_20170618_20_2_0_2.fits.gz"
)[-1]
sprat_public_fits_3 = fits.open(
    "lt-sprat-doaql/v_e_20170618_20_3_0_2.fits.gz"
)[-1]
sprat_public_data = np.nanmean(
    (
        sprat_public_fits_1.data[0],
        sprat_public_fits_2.data[0],
        sprat_public_fits_3.data[0],
    ),
    axis=0,
)
sprat_public_wave = np.linspace(
    sprat_public_fits_1.header["CRVAL1"],
    sprat_public_fits_1.header["CRVAL1"]
    + sprat_public_fits_1.header["CDELT1"]
    * sprat_public_fits_1.header["NAXIS1"],
    sprat_public_fits_1.header["NAXIS1"],
)

sprat_data_resampled = spectres(sprat_public_wave, sprat_wave, sprat_data)


#    ___    _       ___   __   __   ___     ___
#   | __|  | |     / _ \  \ \ / /  |   \   / __|
#   | _|   | |__  | (_) |  \ V /   | |) |  \__ \
#  _|_|_   |____|  \___/   _|_|_   |___/   |___/
# _| """ |_|"""""|_|"""""|_| """ |_|"""""|_|"""""|
# "`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'
#

floyds_fits = fits.open("lco-floyds-iptf14hls/iptf14hls/IPTF14HLS_None.fits")[
    0
]
floyds_data = floyds_fits.data
floyds_wave = np.linspace(
    floyds_fits.header["CRVAL1"],
    floyds_fits.header["CRVAL1"]
    + floyds_fits.header["CDELT1"] * floyds_fits.header["NAXIS1"],
    floyds_fits.header["NAXIS1"],
)

floyds_public_ascii = np.genfromtxt(
    "lco-floyds-iptf14hls/iPTF14hls_2016-02-09_08-15-33_FTN_FLOYDS-N_iPTF.ascii"
)
floyds_public_wave = floyds_public_ascii[:, 0]
floyds_public_data = floyds_public_ascii[:, 1]

floyds_data_resampled = spectres(floyds_public_wave, floyds_wave, floyds_data)


#    ___    ___     ___     ___
#   | __|  / _ \   | _ \   / __|
#   | _|  | (_) |  |   /   \__ \
#  _|_|_   \___/   |_|_\   |___/
# _| """ |_|"""""|_|"""""|_|"""""|
# "`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'
#

fors_ascii = np.genfromtxt("vlt-fors-v418ser/v418ser_chip1_average.dat")
fors_wave = fors_ascii[:, 0]
fors_data = fors_ascii[:, 1] * 2.99792458e-6 / fors_wave**2.0

fors_fits = fits.open("vlt-fors-v418ser/reduced_science_0.fits")[
    "Flux resampled atm ext corrected"
]
fors_wave_aspired = np.linspace(
    fors_fits.header["CRVAL1"],
    fors_fits.header["CRVAL1"]
    + fors_fits.header["CDELT1"] * fors_fits.header["NAXIS1"],
    fors_fits.header["NAXIS1"],
)
fors_data_aspired = spectres(fors_wave, fors_wave_aspired, fors_fits.data)


# big plot here
plt.figure(1, figsize=(8, 8))
plt.clf()

# FLOYDS
plt.plot(
    floyds_public_wave, floyds_data_resampled / 8 + 13.5e-16, color="royalblue"
)
plt.plot(floyds_public_wave, floyds_public_data / 4 + 12e-16, color="crimson")
plt.text(8100, 1.65e-15, "LCO/FLOYDS - iPTF14hls")

# SPRAT
plt.plot(
    sprat_public_wave[10:],
    sprat_data_resampled[10:] * 3.0 + 7.5e-16,
    color="royalblue",
    label="ASPIRED reduction",
)
plt.plot(
    sprat_public_wave[10:],
    sprat_public_data[10:] + 6e-16,
    color="crimson",
    label="iraf-based reduction as per publication",
)
plt.text(8100, 8.5e-16, "LT/SPRAT - DO Aql")

# GMOS
chip_gap_mask = (gmos_public_data_2[:, 0] < 7956.0) | (
    gmos_public_data_2[:, 0] > 8009.0
)
plt.plot(
    gmos_public_data_2[:, 0][chip_gap_mask],
    8 * gmos_flux_combined[chip_gap_mask] + 1.5e-16,
    color="royalblue",
)
plt.plot(
    gmos_public_data_2[:, 0], 8 * gmos_public_data_2[:, 1], color="crimson"
)
plt.text(8200, 1.5e-16, "Gemini/GMOS - AT 2017gfo")

plt.xlim(3800, 9500)
plt.ylim(0, 2.25e-15)
plt.legend()

# FORS
plt.plot(fors_wave, 15 * fors_data_aspired * 21 + 1.5e-15, color="royalblue")
plt.plot(
    fors_wave,
    fors_data / 8.0 + 1.4e-15,
    color="olivedrab",
    label="starlink-based reduction as per publication",
)
plt.text(5100, 2.0e-15, "VLT/FORS2 - V418 Ser")


plt.xlim(3800, 9500)
plt.ylim(0, 2.25e-15)
plt.legend()


ax = plt.gca()
ax.set_yticks([])
plt.xlabel(r"Wavelength ($\mathrm{\AA}$)")
plt.ylabel("Arbitrary Flux (per $\mathrm{\AA}$)")
plt.tight_layout()
plt.savefig("fig_08_reduction_compared.pdf")
plt.savefig("fig_08_reduction_compared.png")
