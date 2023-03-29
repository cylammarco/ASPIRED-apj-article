import glob
import os

from aspired import spectral_reduction
from astropy.io import fits
from matplotlib import pyplot as plt
import numpy as np
from spectres import spectres

plt.ion()

#   ___      ___    ___     ___    _____
#  / __|    | _ \  | _ \   /   \  |_   _|
#  \__ \    |  _/  |   /   | - |    | |
#  |___/   _|_|_   |_|_\   |_|_|   _|_|_
# _|"""""|_| """ |_|"""""|_|"""""|_|"""""|
# "`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'"`-0-0-'
#

filelist = glob.glob(
    "ASPIRED-apj-article-data/lt-sprat-standards/*blue_science.fits"
)
sprat_data_list = []
sprat_wave_list = []

for f in filelist:
    sprat_fits = fits.open(f)[1]
    sprat_data_list.append(sprat_fits.data)
    sprat_wave_list.append(
        np.linspace(
            sprat_fits.header["CRVAL1"],
            sprat_fits.header["CRVAL1"]
            + sprat_fits.header["CDELT1"] * sprat_fits.header["NAXIS1"],
            sprat_fits.header["NAXIS1"],
        )
    )


sprat_wave = np.linspace(
    sprat_fits.header["CRVAL1"],
    sprat_fits.header["CRVAL1"]
    + sprat_fits.header["CDELT1"] * sprat_fits.header["NAXIS1"],
    sprat_fits.header["NAXIS1"],
)

sprat_data_resampled = []
for w, f in zip(sprat_wave_list, sprat_data_list):
    sprat_data_resampled.append(spectres(sprat_wave, w, f))

sprat_median = np.nanmedian(sprat_data_resampled[4:], axis=0)

# big plot here
plt.figure(1, figsize=(8, 6))
plt.clf()

sprat_sd = np.nanstd(
    (
        np.array(sprat_data_resampled[4:]).T
        / np.nanmedian(sprat_data_resampled[4:], axis=1)
    ).T
    * np.nanmedian(sprat_median)
    / sprat_median,
    axis=0,
)
plt.errorbar(
    sprat_wave,
    np.ones_like(sprat_median),
    yerr=sprat_sd,
    color="C0",
    alpha=0.3,
)

# SPRAT
for i, (j, name) in enumerate(zip(sprat_data_resampled, filelist)):
    if i < 4:
        continue
    date = name.split("_")[-3].split(os.sep)[-1]
    plt.plot(
        sprat_wave,
        j / np.nanmedian(j) * np.nanmedian(sprat_median) / sprat_median,
        color="black",
        alpha=0.2
        #        label=date
    )

onedspec = spectral_reduction.OneDSpec(
    log_level="INFO", log_file_name="sprat_doaql"
)
onedspec.load_standard(
    target="bd332642", library="irafirscal", cutoff=0.4, ftype="flux"
)

flux_resampled = (
    spectres(onedspec.fluxcal.standard_wave_true, sprat_wave, sprat_median)
    / np.nanmedian(sprat_median[(sprat_wave > 4500) & (sprat_wave < 7500)])
    * np.nanmedian(
        onedspec.fluxcal.standard_fluxmag_true[
            (onedspec.fluxcal.standard_wave_true > 4500)
            & (onedspec.fluxcal.standard_wave_true < 7500)
        ]
    )
)
"""
plt.plot(
    onedspec.fluxcal.standard_wave_true,
    onedspec.fluxcal.standard_fluxmag_true / flux_resampled,
    lw=2.5,
    label="Difference from the median"
)
"""

plt.hlines(1.0, 4000, 8000, ls="dashed", color="black")

plt.xlim(4250, 7750)
plt.ylim(0.85, 1.15)

ax = plt.gca()
plt.xlabel(r"Wavelength ($\mathrm{\AA}$)")
plt.ylabel("Fractional Residual Flux (per $\mathrm{\AA}$)")
plt.tight_layout()
plt.savefig("fig_08_sprat_standards_compared.pdf")
plt.savefig("fig_08_sprat_standards_compared.png")
