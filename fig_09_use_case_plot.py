from astropy.io import fits
from matplotlib import pyplot as plt
from spectres import spectres
import numpy as np

plt.ion()

# DOLORES dmwd
dolores_fits_1 = fits.open("tng-dolores-dmwd/dmwd_blue_1_science_1.fits")[1]
dolores_fits_2 = fits.open("tng-dolores-dmwd/dmwd_blue_2_science_1.fits")[1]
dolores_fits_3 = fits.open("tng-dolores-dmwd/dmwd_blue_3_science_1.fits")[1]
dolores_fits_4 = fits.open("tng-dolores-dmwd/dmwd_blue_4_science_1.fits")[1]
dolores_wave_1 = np.linspace(
    dolores_fits_1.header["CRVAL1"],
    dolores_fits_1.header["CRVAL1"]
    + dolores_fits_1.header["CDELT1"] * dolores_fits_1.header["NAXIS1"],
    dolores_fits_1.header["NAXIS1"],
)
dolores_wave_2 = np.linspace(
    dolores_fits_2.header["CRVAL1"],
    dolores_fits_2.header["CRVAL1"]
    + dolores_fits_2.header["CDELT1"] * dolores_fits_2.header["NAXIS1"],
    dolores_fits_2.header["NAXIS1"],
)
dolores_wave_3 = np.linspace(
    dolores_fits_3.header["CRVAL1"],
    dolores_fits_3.header["CRVAL1"]
    + dolores_fits_3.header["CDELT1"] * dolores_fits_3.header["NAXIS1"],
    dolores_fits_3.header["NAXIS1"],
)
dolores_wave_4 = np.linspace(
    dolores_fits_4.header["CRVAL1"],
    dolores_fits_4.header["CRVAL1"]
    + dolores_fits_4.header["CDELT1"] * dolores_fits_4.header["NAXIS1"],
    dolores_fits_4.header["NAXIS1"],
)
dolores_wave = dolores_wave_1
dolores_flux_1 = dolores_fits_1.data
dolores_flux_2 = spectres(dolores_wave, dolores_wave_2, dolores_fits_2.data)
dolores_flux_3 = spectres(dolores_wave, dolores_wave_3, dolores_fits_3.data)
dolores_flux_4 = spectres(dolores_wave, dolores_wave_4, dolores_fits_4.data)

dolores_flux = np.mean(
    (dolores_flux_1, dolores_flux_2, dolores_flux_3, dolores_flux_4), axis=0
)

# DOLORES dm
dolores_fits_dm_1 = fits.open("tng-dolores-dmwd/dmwd_blue_1_science_2.fits")[1]
dolores_fits_dm_2 = fits.open("tng-dolores-dmwd/dmwd_blue_2_science_2.fits")[1]
dolores_fits_dm_3 = fits.open("tng-dolores-dmwd/dmwd_blue_3_science_2.fits")[1]
dolores_fits_dm_4 = fits.open("tng-dolores-dmwd/dmwd_blue_4_science_2.fits")[1]
dolores_wave_dm_1 = np.linspace(
    dolores_fits_dm_1.header["CRVAL1"],
    dolores_fits_dm_1.header["CRVAL1"]
    + dolores_fits_dm_1.header["CDELT1"] * dolores_fits_dm_1.header["NAXIS1"],
    dolores_fits_dm_1.header["NAXIS1"],
)
dolores_wave_dm_2 = np.linspace(
    dolores_fits_dm_2.header["CRVAL1"],
    dolores_fits_dm_2.header["CRVAL1"]
    + dolores_fits_dm_2.header["CDELT1"] * dolores_fits_dm_2.header["NAXIS1"],
    dolores_fits_dm_2.header["NAXIS1"],
)
dolores_wave_dm_3 = np.linspace(
    dolores_fits_dm_3.header["CRVAL1"],
    dolores_fits_dm_3.header["CRVAL1"]
    + dolores_fits_dm_3.header["CDELT1"] * dolores_fits_dm_3.header["NAXIS1"],
    dolores_fits_dm_3.header["NAXIS1"],
)
dolores_wave_dm_4 = np.linspace(
    dolores_fits_dm_4.header["CRVAL1"],
    dolores_fits_dm_4.header["CRVAL1"]
    + dolores_fits_dm_4.header["CDELT1"] * dolores_fits_dm_4.header["NAXIS1"],
    dolores_fits_dm_4.header["NAXIS1"],
)
dolores_dm_wave = dolores_wave_dm_1
dolores_flux_dm_1 = dolores_fits_dm_1.data
dolores_flux_dm_2 = spectres(
    dolores_dm_wave, dolores_wave_dm_2, dolores_fits_dm_2.data
)
dolores_flux_dm_3 = spectres(
    dolores_dm_wave, dolores_wave_dm_3, dolores_fits_dm_3.data
)
dolores_flux_dm_4 = spectres(
    dolores_dm_wave, dolores_wave_dm_4, dolores_fits_dm_4.data
)

dolores_dm_flux = np.mean(
    (
        dolores_flux_dm_1,
        dolores_flux_dm_2,
        dolores_flux_dm_3,
        dolores_flux_dm_4,
    ),
    axis=0,
)


# OSIRIS
r1000b_1 = np.load("gtc-osiris-zgpblap09/r1000b_1.npy")
r1000b_1_wave = r1000b_1[:, 0]
r1000b_1_flux = r1000b_1[:, 1]

r1000b_2 = np.load("gtc-osiris-zgpblap09/r1000b_2.npy")
r1000b_2_wave = r1000b_2[:, 0]
r1000b_2_flux = r1000b_2[:, 1]

r2500u_1 = np.load("gtc-osiris-zgpblap09/r2500u_1.npy")
r2500u_1_wave = r2500u_1[:, 0]
r2500u_1_flux = r2500u_1[:, 1]

r2500u_2 = np.load("gtc-osiris-zgpblap09/r2500u_2.npy")
r2500u_2_wave = r2500u_2[:, 0]
r2500u_2_flux = r2500u_2[:, 1]

r1000b_2_resampled = spectres(r1000b_1_wave, r1000b_2_wave, r1000b_2_flux)
r2500u_2_resampled = spectres(r2500u_1_wave, r2500u_2_wave, r2500u_2_flux)

ymax = np.nanmax((r2500u_1_flux, r2500u_2_flux))

plt.rcParams.update({"font.size": 16})
plt.rcParams.update({"legend.fontsize": 12})

plt.figure(1, figsize=(8, 10))
plt.clf()
# OSIRIS
plt.plot(
    r2500u_1_wave,
    np.nanmean((r2500u_1_flux, r2500u_2_resampled), axis=0) * 5 + 4.0e-14,
    color="royalblue",
    label="GTC/OSIRIS R2500U",
)
plt.plot(
    r1000b_1_wave[r1000b_1_wave > 4500.0],
    np.nanmean((r1000b_1_flux, r1000b_2_resampled), axis=0)[
        r1000b_1_wave > 4500.0
    ]
    * 1.35
    * 5.0
    + 4.0e-14,
    color="royalblue",
    alpha=0.6,
    label="GTC/OSIRIS R1000B",
)

# DOLORES
plt.plot(
    dolores_wave,
    dolores_flux + 1.5e-14,
    color="crimson",
    label="TNG/DOLORES dMWD",
)
plt.plot(
    dolores_dm_wave,
    dolores_dm_flux + 1.5e-14,
    color="crimson",
    alpha=0.6,
    label="TNG/DOLORES dM",
)
plt.text(3750, 2.3e-14, "Simultaneous Extraction", color="crimson", alpha=0.8)

plt.xlim(3500, 8000)
plt.ylim(0.0, 7e-14)
plt.xlabel(r"Wavelength / $\mathrm{\AA}$")
plt.ylabel(r"Flux / ( erg / cm / cm / s / $\mathrm{\AA}$)")
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig("fig_09_use_case_plots.pdf")
plt.savefig("fig_09_use_case_plots.png")
