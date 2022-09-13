from astropy.io import fits
import glob
import numpy as np
from scipy import ndimage

from aspired import image_reduction
from aspired import spectral_reduction

# Line list
atlas = [
    3650.15,
    3888.700,
    4046.560,
    4077.810,
    4358.340,
    4471.500,
    4678.10,
    4713.200,
    4799.900,
    4921.929,
    5015.700,
    5085.800,
    5460.750,
    5769.590,
    5790.690,
    5875.620,
]
element = ["HeHgCd"] * len(atlas)


# Science frames
science_filelist = glob.glob("science/*.gz")
science_light_fits_list = []
for f in science_filelist:
    science_light_fits_list.append(fits.open(f))

arc1 = fits.open("hehgcd_arc/FORS2.2015-04-09T11_52_07.383.fits.gz")
arc2 = fits.open("hehgcd_arc/FORS2.2015-04-09T12_14_43.884.fits.gz")

flat_filelist = glob.glob("flat/*.gz")
flat_fits_list = []
for f in flat_filelist:
    flat_fits_list.append(fits.open(f))

bias_filelist = glob.glob("bias/*.gz")
bias_fits_list = []
for f in bias_filelist:
    bias_fits_list.append(fits.open(f))


# Add all the frames to the right places
science_frame = image_reduction.ImageReduction(log_file_name=None)
for f in science_light_fits_list:
    science_frame.add_light(
        f[0].data,
        f[0].header,
        f[0].header["EXPTIME"],
    )

science_frame.add_arc(
    arc1[0].data,
    arc1[0].header,
)

science_frame.add_arc(
    arc2[0].data,
    arc2[0].header,
)

for f in flat_fits_list:
    science_frame.add_flat(
        f[0].data,
        f[0].header,
        f[0].header["EXPTIME"],
    )

for f in bias_fits_list:
    science_frame.add_bias(
        f[0].data,
        f[0].header,
    )

science_frame.reduce()
science_frame.inspect(
    save_fig=True,
    fig_type="jpg",
    display=False,
    filename="v418ser_redcued_spectral_image",
)


# Standard frame
standard_light = fits.open("standard/FORS2.2015-04-09T08_56_00.429.fits.gz")

standard_frame = image_reduction.ImageReduction(log_file_name=None)
standard_frame.add_light(
    standard_light[0].data,
    standard_light[0].header,
    standard_light[0].header["EXPTIME"],
)

standard_frame.add_arc(
    arc1[0].data,
    arc1[0].header,
)

standard_frame.add_arc(
    arc2[0].data,
    arc2[0].header,
)

for f in flat_fits_list:
    standard_frame.add_flat(
        f[0].data,
        f[0].header,
        f[0].header["EXPTIME"],
    )

for f in bias_fits_list:
    standard_frame.add_bias(
        f[0].data,
        f[0].header,
    )

standard_frame.reduce()
standard_frame.inspect(
    save_fig=True,
    fig_type="jpg",
    display=False,
    filename="lt7379_redcued_spectral_image",
)

spatial_mask = np.arange(50, 300)

v418_twodspec = spectral_reduction.TwoDSpec(
    science_frame,
    readnoise=2.9,
    gain=1.43,
    log_file_name=None,
    spatial_mask=spatial_mask,
)

ltt7379_twodspec = spectral_reduction.TwoDSpec(
    standard_frame,
    readnoise=2.9,
    gain=1.43,
    log_file_name=None,
    spatial_mask=spatial_mask,
)

v418_twodspec.ap_trace(
    nspec=1,
    display=False,
    save_fig=True,
    fig_type="jpg",
    filename="ap_trace_v418",
)
ltt7379_twodspec.ap_trace(
    nspec=1,
    display=False,
    save_fig=True,
    fig_type="jpg",
    filename="ap_trace_ltt7379",
)


v418_twodspec.ap_extract(
    apwidth=8,
    skysep=3,
    skywidth=6,
    skydeg=1,
    optimal=True,
    display=False,
    save_fig=True,
    fig_type="jpg",
    filename="ap_extract_v418",
)

ltt7379_twodspec.ap_extract(
    apwidth=12,
    skysep=3,
    skywidth=6,
    skydeg=1,
    display=False,
    optimal=True,
    save_fig="jpg",
    filename="ap_extract_ltt7379",
)


# Extract the 1D arc by aperture sum of the traces provided
v418_twodspec.extract_arc_spec(
    display=False, save_fig=True, renderer="jpg", fig_type="jpg"
)
ltt7379_twodspec.extract_arc_spec(
    display=False, save_fig=True, renderer="jpg", fig_type="jpg"
)

v418_onedspec = spectral_reduction.OneDSpec(
    log_level="INFO", log_file_name="vlt_fors_v418ser"
)
v418_onedspec.from_twodspec(v418_twodspec, stype="science")
v418_onedspec.from_twodspec(ltt7379_twodspec, stype="standard")


# Find the peaks of the arc
v418_onedspec.find_arc_lines(
    top_n_peaks=15,
    prominence=0.01,
    save_fig=True,
    fig_type="jpg",
    display=False,
)

# Configure the wavelength calibrator
v418_onedspec.initialise_calibrator()
v418_onedspec.set_hough_properties(
    num_slopes=500,
    xbins=100,
    ybins=100,
    min_wavelength=3500,
    max_wavelength=6000,
    range_tolerance=250,
)

v418_onedspec.add_user_atlas(elements=element, wavelengths=atlas, vacuum=False)

v418_onedspec.do_hough_transform()

v418_onedspec.set_ransac_properties(minimum_matches=5)

v418_onedspec.fit(
    fit_deg=2,
    candidate_tolerance=5.,
    fit_tolerance=10.,
    max_tries=1000,
    save_fig=True,
    fig_type="jpg",
    display=False,
)

v418_onedspec.apply_wavelength_calibration()


v418_onedspec.load_standard(target="ltt7379", cutoff=0.4, ftype="flux")

# Get the sensitivity curves
v418_onedspec.get_sensitivity()
v418_onedspec.apply_flux_calibration()


v418_onedspec.set_atmospheric_extinction(location='cp')
v418_onedspec.apply_atmospheric_extinction_correction(science_airmass=1.207, standard_airmass=1.101)

v418_onedspec.inspect_reduced_spectrum(
    save_fig=True,
    fig_type="jpg",
    display=False,
)

v418_onedspec.create_fits(
    output="*",
)
v418_onedspec.save_fits(
    output="*",
    overwrite=True,
)
