import sys
import numpy as np
from astropy.io import fits
from aspired import image_reduction
from aspired import spectral_reduction
from matplotlib.pyplot import *

ion()


# Line list
atlas = [
    5460.735,
    5769.598,
    5790.663,
    5852.488,
    5881.895,
    5944.834,
    5974.6273,
    6029.9969,
    6074.338,
    6096.163,
    6143.063,
    6163.594,
    6217.2812,
    6266.495,
    6304.7889,
    6334.4278,
    6382.9917,
    6402.248,
    6506.5281,
    6532.8822,
    6598.9529,
    6678.2762,
    6717.043,
    6929.4673,
    6965.4310,
    7032.4131,
    7067.218,
    7147.042,
    7173.9381,
    7245.1666,
    7272.9360,
    7383.98,
    7435.368,
    7503.869,
    7514.652,
    7635.106,
    7724.207,
    7948.176,
    8006.157,
    8014.786,
    8103.693,
    8115.311,
    8264.522,
    8377.608,
    8408.21,
    8424.648,
    8495.3598,
    8521.442,
    8591.2584,
    8634.6470,
    8654.3831,
    8783.7533,
    8853.8668,
    8919.5006,
]
element = ["HgArNe"] * len(atlas)

"""
atlas_Hg = [
    3650.153,
    3983.931,
    4046.563,
    4358.328,
    5128.442,
    5425.253,
    5460.735,
    5677.105,
    5888.939,
    6149.475,
    7944.555,
]


atlas_Ne = [
    3027.016,
    3028.864,
    3198.586,
    3323.745,
    3378.216,
    3520.4711,
    3766.259,
    3777.133,
    5400.5618,
    5852.4879,
    6029.9969,
    6074.3377,
    6143.0626,
    6163.5939,
    6217.2812,
    6266.4950,
    6382.9917,
    6402.248,
    6506.5281,
    6598.9529,
    6929.4673,
    7032.4131,
    7173.9381,
    7245.1666,
    8377.6080,
    8654.3831,
    8780.6226,
    8783.7533,
    11143.0200,
    11177.5240,
]

atlas_Ar = [
    4277.528,
    4348.064,
    4609.567,
    4726.868,
    4764.865,
    4806.020,
    4879.864,
    6965.431,
    7067.218,
    7503.869,
    7635.106,
    7948.176,
    8006.157,
    8014.786,
    8103.693,
    8115.311,
    9122.967,
    9657.786,
]

element_Hg = ["Hg"] * len(atlas_Hg)
element_Ne = ["Ne"] * len(atlas_Ne)
element_Ar = ["Ar"] * len(atlas_Ar)
"""

# ztf22aapmawo
sci_1_fits = fits.open("0303_ztf22aapmawo_04-07-2022.fits")[0]
sci_1_arc_fits = fits.open("0304_ztf22aapmawo_04-07-2022_comp.fits")[0]

std_fits = fits.open("0290_Calibration-Star_03-07-2022.fits")[0]
std_arc_fits = fits.open("0289_Calibration-Star_03-07-2022_comp.fits")[0]

flat_fits = fits.open("0305_ztf22aapmawo_04-07-2022.fits")[0]

sci_1_image_reduction = image_reduction.ImageReduction()
std_image_reduction = image_reduction.ImageReduction()

sci_1_image_reduction.add_light(
    sci_1_fits.data, sci_1_fits.header, sci_1_fits.header["EXPTIME"]
)
sci_1_image_reduction.add_flat(
    flat_fits.data, flat_fits.header, flat_fits.header["EXPTIME"]
)
sci_1_image_reduction.add_arc(sci_1_arc_fits.data, sci_1_arc_fits.header)
sci_1_image_reduction.reduce()

std_image_reduction.add_light(
    std_fits.data, std_fits.header, std_fits.header["EXPTIME"]
)
std_image_reduction.add_flat(
    flat_fits.data, flat_fits.header, flat_fits.header["EXPTIME"]
)
std_image_reduction.add_arc(std_arc_fits.data, std_arc_fits.header)
std_image_reduction.reduce()

# initialise the two aspired.TwoDSpec()
spatial_mask = np.arange(400, 600)
spec_mask = np.arange(25, 2068)

sci_1_twodspec = spectral_reduction.TwoDSpec(
    sci_1_image_reduction,
    spatial_mask=spatial_mask,
    spec_mask=spec_mask,
    cosmicray=True,
    fsmode="convolve",
    psfsize=15,
    gain=1.48,
    readnoise=3.89,
    sigclip=4.0,
)
std_twodspec = spectral_reduction.TwoDSpec(
    std_image_reduction,
    spatial_mask=spatial_mask,
    spec_mask=spec_mask,
    cosmicray=True,
    fsmode="convolve",
    psfsize=15,
    gain=1.48,
    readnoise=3.89,
    sigclip=4.0,
)

# Trace the spectrs
std_twodspec.ap_trace(
    nspec=1,
    ap_faint=10,
    fit_deg=2,
    resample_factor=5,
    save_fig=True,
    fig_type="jpg",
    filename="std_trace",
)
sci_1_twodspec.ap_trace(
    nspec=1,
    ap_faint=10,
    fit_deg=2,
    resample_factor=5,
    save_fig=True,
    fig_type="jpg",
    filename="sci_1_trace",
)

# Extract the spectra
sci_1_twodspec.ap_extract(
    apwidth=7,
    skysep=1,
    skywidth=5,
    skydeg=1,
    optimal=True,
    model="gauss",
    save_fig=True,
    fig_type="jpg",
    filename="sci_1_extract",
)
std_twodspec.ap_extract(
    apwidth=10,
    skywidth=10,
    skydeg=1,
    optimal=True,
    lowess_frac=0.1,
    save_fig=True,
    model="gauss",
    fig_type="jpg",
    filename="std_extract",
)

sci_1_twodspec.extract_arc_spec(display=False)
std_twodspec.extract_arc_spec(display=False)

# Now work in one dimension
sci_1_onedspec = spectral_reduction.OneDSpec()

sci_1_onedspec.from_twodspec(sci_1_twodspec, stype="science")

sci_1_onedspec.from_twodspec(std_twodspec, stype="standard")

sci_1_onedspec.find_arc_lines(
    prominence=1.0,
    refine=True,
    refine_window_width=3,
    top_n_peaks=100,
    display=True,
    save_fig=True,
    stype="science+standard",
)
sci_1_onedspec.initialise_calibrator(stype="science+standard")
sci_1_onedspec.set_hough_properties(
    num_slopes=2000,
    xbins=200,
    ybins=200,
    min_wavelength=4900,
    max_wavelength=9000,
    range_tolerance=500,
    stype="science+standard",
)
sci_1_onedspec.set_ransac_properties(
    filter_close=True,
    minimum_matches=20,
    ransac_tolerance=5.0,
    stype="science+standard",
)
sci_1_onedspec.add_user_atlas(
    elements=element,
    wavelengths=atlas,
    stype="science+standard",
)
sci_1_onedspec.do_hough_transform()
sci_1_onedspec.fit(
    max_tries=1000,
    fit_deg=5,
    fit_tolerance=5.0,
    candidate_tolerance=5.0,
    stype="science+standard",
    display=True,
    save_fig=True,
)

sci_1_onedspec.apply_wavelength_calibration()
sci_1_onedspec.load_standard(target="LTT4364", cutoff=0.4, ftype="flux")

# Get the sensitivity curves
sci_1_onedspec.get_sensitivity(method="polynomial", sens_deg=15)
sci_1_onedspec.apply_flux_calibration()

sci_1_onedspec.apply_atmospheric_extinction_correction()

sci_1_onedspec.get_telluric_profile()
sci_1_onedspec.get_telluric_strength()
sci_1_onedspec.apply_telluric_correction()

sci_1_onedspec.inspect_reduced_spectrum(
    wave_min=4970.0, wave_max=9000.0, save_fig=True, fig_type="png"
)

sci_1_onedspec.resample(wave_start=4970.0, wave_end=9000.0, wave_bin=2.0)
sci_1_onedspec.save_fits(
    output="flux_resampled_atm_ext_telluric_corrected", overwrite=True
)

flux_1 = sci_1_onedspec.science_spectrum_list[
    0
].flux_resampled_atm_ext_telluric_corrected
wavelength = sci_1_onedspec.science_spectrum_list[0].wave_resampled

figure(1, figsize=(12, 8))
clf()
plot(wavelength, flux_1)
xlim(4970.0, 9000.0)
ylim(np.nanpercentile(flux_1, 0), np.nanpercentile(flux_1, 99))
grid()
xlabel("Wavelength (A)")
ylabel("Flux (erg / s / cm / cm / A)")
tight_layout()
savefig("ztf22aapmawo_20220704.png")


flux_std = sci_1_onedspec.standard_spectrum_list[
    0
].flux_resampled_atm_ext_corrected
wavelength_std = sci_1_onedspec.standard_spectrum_list[0].wave_resampled

figure(2, figsize=(12, 8))
clf()
plot(wavelength_std, flux_std)
xlim(4970.0, 9000.0)
ylim(np.nanpercentile(flux_std, 0), np.nanpercentile(flux_std, 99))
grid()
xlabel("Wavelength (A)")
ylabel("Flux (erg / s / cm / cm / A)")
tight_layout()
savefig("ztf22aapmawo_standard_20220704.png")
