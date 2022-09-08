from astropy.io import fits
import numpy as np
from scipy import ndimage

from aspired import image_reduction
from aspired import spectral_reduction

# Line list
atlas = [
    4193.5,
    4385.77,
    4500.98,
    4524.68,
    4582.75,
    4624.28,
    4671.23,
    4697.02,
    4734.15,
    4807.02,
    4921.48,
    5028.28,
    5618.88,
    5823.89,
    5893.29,
    5934.17,
    6182.42,
    6318.06,
    6472.841,
    6595.56,
    6668.92,
    6728.01,
    6827.32,
    6976.18,
    7119.60,
    7257.9,
    7393.8,
    7584.68,
    7642.02,
    7740.31,
    7802.65,
    7887.40,
    7967.34,
    8057.258,
]
element = ["Xe"] * len(atlas)


science_light_1 = fits.open("v_e_20170618_20_1_0_2.fits.gz")
science_light_2 = fits.open("v_e_20170618_20_2_0_2.fits.gz")
science_light_3 = fits.open("v_e_20170618_20_3_0_2.fits.gz")
standard_light = fits.open("v_s_20170618_13_1_0_2.fits.gz")
arc = fits.open("v_a_20170618_21_1_0_1.fits.gz")
dark = fits.open("v_dark_1.fits.gz")

science_frame = image_reduction.ImageReduction(log_file_name=None)
science_frame.add_light(
    science_light_1[0].data,
    science_light_1[0].header,
    science_light_1[0].header["EXPTIME"],
)
science_frame.add_light(
    science_light_2[0].data,
    science_light_2[0].header,
    science_light_2[0].header["EXPTIME"],
)
science_frame.add_light(
    science_light_3[0].data,
    science_light_3[0].header,
    science_light_3[0].header["EXPTIME"],
)
science_frame.add_dark(dark[0].data, dark[0].header, dark[0].header["EXPTIME"])
science_frame.add_arc(arc[0].data, arc[0].header)
science_frame.reduce()
science_frame.inspect(
    filename="doaql",
    save_fig=True,
    fig_type="jpg",
    display=False,
)

standard_frame = image_reduction.ImageReduction(log_file_name=None)
standard_frame.add_light(
    standard_light[0].data,
    standard_light[0].header,
    standard_light[0].header["EXPTIME"],
)
standard_frame.add_dark(
    dark[0].data, dark[0].header, dark[0].header["EXPTIME"]
)
standard_frame.add_arc(arc[0].data, arc[0].header)
standard_frame.reduce()
standard_frame.inspect(
    save_fig=True,
    fig_type="jpg",
    display=False,
)


# Set the spectral direction is defaulted to saxis=1
# saxis = 1

# initialise the two aspired.TwoDSpec()
doaql_twodspec = spectral_reduction.TwoDSpec(
    science_frame, cosmicray=True, readnoise=5.7, log_file_name=None
)

hilt102_twodspec = spectral_reduction.TwoDSpec(
    standard_frame, cosmicray=True, readnoise=5.7, log_file_name=None
)

# automatically trace the spectrum
doaql_twodspec.ap_trace(nspec=1, display=False, save_fig=True, fig_type="jpg")
hilt102_twodspec.ap_trace(
    nspec=1, display=False, save_fig=True, fig_type="jpg", filename="ap_trace_standard"
)


# Optimal extracting spectrum by summing over the aperture along the trace
doaql_twodspec.ap_extract(
    spec_id=0,
    apwidth=12,
    skysep=12,
    skywidth=5,
    skydeg=1,
    optimal=False,
    display=False,
    save_fig=True,
    fig_type="jpg",
    filename="ap_extract"
)

hilt102_twodspec.ap_extract(
    apwidth=15,
    skysep=3,
    skywidth=5,
    skydeg=1,
    display=False,
    optimal=True,
    save_fig="jpg",
    filename="ap_extract_standard"
)

# Extract the 1D arc by aperture sum of the traces provided
doaql_twodspec.extract_arc_spec(
    spec_id=0, display=False, save_fig=True, renderer="jpg", fig_type="jpg"
)
hilt102_twodspec.extract_arc_spec(
    display=False, save_fig=True, renderer="jpg", fig_type="jpg"
)

doaql_onedspec = spectral_reduction.OneDSpec(
    log_level="INFO", log_file_name="sprat_doaql"
)
doaql_onedspec.from_twodspec(doaql_twodspec, spec_id=0, stype="science")
doaql_onedspec.from_twodspec(hilt102_twodspec, stype="standard")

# Find the peaks of the arc
doaql_onedspec.find_arc_lines(
    spec_id=0,
    top_n_peaks=35,
    prominence=1,
    stype="science",
    save_fig=True,
    fig_type="jpg",
    display=False,
)
doaql_onedspec.find_arc_lines(
    top_n_peaks=35,
    prominence=1,
    stype="standard",
    save_fig=True,
    fig_type="jpg",
    display=False,
)

# Configure the wavelength calibrator
doaql_onedspec.initialise_calibrator(spec_id=0, stype="science")
doaql_onedspec.initialise_calibrator(stype="standard")
doaql_onedspec.set_hough_properties(
    spec_id=0,
    num_slopes=500,
    xbins=100,
    ybins=100,
    min_wavelength=3500,
    max_wavelength=8000,
    range_tolerance=250,
    stype="science",
)
doaql_onedspec.set_hough_properties(
    num_slopes=500,
    xbins=100,
    ybins=100,
    min_wavelength=3500,
    max_wavelength=8000,
    range_tolerance=250,
    stype="standard",
)

doaql_onedspec.add_user_atlas(
    spec_id=0, elements=element, wavelengths=atlas, stype="science"
)
doaql_onedspec.add_user_atlas(
    elements=element, wavelengths=atlas, stype="standard"
)

doaql_onedspec.do_hough_transform(spec_id=0, stype="science")
doaql_onedspec.do_hough_transform(stype="standard")

doaql_onedspec.set_ransac_properties(
    spec_id=0, minimum_matches=18, stype="science"
)
doaql_onedspec.set_ransac_properties(minimum_matches=18, stype="standard")

doaql_onedspec.fit(
    spec_id=0,
    max_tries=1000,
    stype="science",
    save_fig=True,
    fig_type="jpg",
    display=False,
)
doaql_onedspec.fit(
    max_tries=1000,
    stype="standard",
    save_fig=True,
    fig_type="jpg",
    display=False,
)

doaql_onedspec.apply_wavelength_calibration(spec_id=0, stype="science")
doaql_onedspec.apply_wavelength_calibration(stype="standard")

doaql_onedspec.load_standard(
    target="BDp332642", cutoff=0.4, ftype="flux"
)

# Get the sensitivity curves
doaql_onedspec.get_sensitivity()
doaql_onedspec.apply_flux_calibration(spec_id=0, stype="science")
doaql_onedspec.apply_flux_calibration(stype="standard")


# atmospheric extinction correction
doaql_onedspec.set_atmospheric_extinction(location="orm")
doaql_onedspec.apply_atmospheric_extinction_correction(spec_id=0)

# telluric absorption correction
doaql_onedspec.get_telluric_profile()
doaql_onedspec.get_telluric_strength(spec_id=0)
doaql_onedspec.apply_telluric_correction(spec_id=0, stype="science")

doaql_onedspec.inspect_reduced_spectrum(
    spec_id=0,
    stype="science",
    save_fig=True,
    fig_type="jpg",
    display=False,
)

doaql_onedspec.create_fits(
    spec_id=0,
    output="flux_resampled_atm_ext_telluric_corrected",
)
doaql_onedspec.save_fits(
    spec_id=0,
    output="flux_resampled_atm_ext_telluric_corrected",
    overwrite=True,
)
