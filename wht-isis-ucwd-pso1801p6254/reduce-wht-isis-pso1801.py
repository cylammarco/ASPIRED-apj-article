from aspired import image_reduction
from aspired import spectral_reduction
from astropy.io import fits
import numpy as np
import os


# Set the spectral and spatial direction
Saxis = 0


# RED ARM HERE
# RED ARM HERE
# RED ARM HERE
# RED ARM HERE
# RED ARM HERE

# spec mask
spatial_mask = np.arange(475, 575)
spec_mask = np.arange(200, 4000)


science_frame = image_reduction.ImageReduction()
science_frame.add_filelist("isis_pso1801p6254_red.list")
science_frame.set_properties(cosmicray=True, fsmode="median", saxis=Saxis)
science_frame.load_data()
science_frame.reduce()
science_frame.inspect(
    display=False, save_fig=True, filename="reduced_image_pso1801p6254_red"
)
# science_frame.savefits(overwrite=True)

standard_frame = image_reduction.ImageReduction()
standard_frame.add_filelist("isis_g93m48_red.list")
standard_frame.set_properties(saxis=Saxis)
standard_frame.load_data()
standard_frame.reduce()
standard_frame.inspect(
    display=False, save_fig=True, filename="reduced_image_g93m48_red"
)


# initialise the two spectral_reduction.TwoDSpec()
pso = spectral_reduction.TwoDSpec(
    science_frame,
    spatial_mask=spatial_mask,
    spec_mask=spec_mask,
    readnoise=4.5,
    cosmicray=False,
    gain=0.98,
    airmass=1.498530,
    seeing=1.1,
    verbose=False,
)
g93 = spectral_reduction.TwoDSpec(
    standard_frame,
    spatial_mask=spatial_mask,
    spec_mask=spec_mask,
    readnoise=4.5,
    cosmicray=False,
    gain=0.98,
    airmass=1.130636,
    seeing=1.1,
    verbose=False,
)

# automatically trace the spectrum
pso.ap_trace(
    display=False,
    nwindow=15,
    percentile=25,
    shift_tol=10.0,
    ap_faint=20,
    fit_deg=1,
    filename="pso_trace_red",
)

g93.ap_trace(display=False, fig_type="png", filename="g93_trace_red")

pso.spectrum_list[0].trace_sigma = g93.spectrum_list[0].trace_sigma

# Optimal extracting spectrum by summing over the aperture along the trace
pso.ap_extract(
    apwidth=15,
    optimal=True,
    display=False,
    fig_type="png",
    filename="pso_extract_red",
)
g93.ap_extract(
    apwidth=20,
    skysep=3,
    skywidth=5,
    skydeg=1,
    optimal=True,
    display=False,
    fig_type="png",
    filename="g93_extract_red",
)

pso.extract_arc_spec(
    display=False, fig_type="png", filename="arc_spec_science_red"
)
g93.extract_arc_spec(
    display=False, fig_type="png", filename="arc_spec_standard_red"
)


pso_reduced = spectral_reduction.OneDSpec()
pso_reduced.from_twodspec(pso, stype="science")
pso_reduced.from_twodspec(g93, stype="standard")


pso_reduced.find_arc_lines(
    prominence=2.0,
    display=False,
    stype="science+standard",
    fig_type="png",
    filename="arc_lines_red",
)


pso_reduced.initialise_calibrator(stype="science+standard")
pso_reduced.set_hough_properties(
    num_slopes=10000,
    xbins=500,
    ybins=500,
    min_wavelength=7000.0,
    max_wavelength=10500.0,
    range_tolerance=200,
    stype="science+standard",
)
pso_reduced.set_ransac_properties(filter_close=True, minimum_matches=13)
pso_reduced.add_atlas(
    elements=["Cu", "Ne", "Ar"],
    min_atlas_wavelength=7000.0,
    max_atlas_wavelength=10500.0,
    min_intensity=5,
    pressure=80000.0,
    temperature=286.25,
    relative_humidity=30.0,
    stype="science+standard",
)
pso_reduced.do_hough_transform(brute_force=False)


pso_reduced.fit(
    max_tries=1000,
    fit_tolerance=3.0,
    display=False,
    fig_type="png",
    stype="science+standard",
)


pso_reduced.apply_wavelength_calibration(stype="science+standard")


pso_reduced.load_standard(
    target="g93_48",
    library="esowdstan",
    cutoff=0.4,
)
pso_reduced.inspect_standard(
    fig_type="png", filename="literature_standard_red"
)

pso_reduced.get_sensitivity()
pso_reduced.inspect_sensitivity(fig_type="png", filename="sensitivity_red")

pso_reduced.apply_flux_calibration(stype="science+standard")


pso_reduced.set_atmospheric_extinction(location="orm")
pso_reduced.apply_atmospheric_extinction_correction()


pso_reduced.get_telluric_profile()
pso_reduced.get_telluric_strength()
pso_reduced.apply_telluric_correction()


pso_reduced.inspect_reduced_spectrum(
    wave_min=6200.0,
    wave_max=10500.0,
    stype="science",
    save_fig=True,
    fig_type="png",
    filename="pso_reduced_spectrum_red",
)

pso_reduced.inspect_reduced_spectrum(
    wave_min=6200.0,
    wave_max=10500.0,
    stype="standard",
    save_fig=True,
    fig_type="png",
    filename="g93_reduced_spectrum_red",
)

pso_reduced.create_fits(
    output="flux_resampled_atm_ext_telluric_corrected",
)

pso_reduced.save_fits(
    output="flux_resampled_atm_ext_telluric_corrected",
    filename="pso_reduced_red",
    overwrite=True,
)

# BLUE ARM HERE
# BLUE ARM HERE
# BLUE ARM HERE
# BLUE ARM HERE
# BLUE ARM HERE


atlas = [
    3247.54,
    3273.96,
    3337.84,
    3350.92,
    3520.472,
    3561.03,
    3568.53,
    3582.35,
    3588.44,
    3593.64,
    3694.20,
    3713.08,
    3727.08,
    3766.12,
    3850.58,
    3868.53,
    3925.72,
    3946.10,
    4052.92,
    4072.00,
    4103.91,
    4131.72,
    4158.59,
    4181.88,
    4190.87,
    4200.68,
    4229.87,
    4259.36,
    4266.29,
    4272.17,
    4277.53,
    4300.10,
    4309.24,
    4331.20,
    4335.56,
    4426.00,
    4474.76,
    4481.81,
    4510.73,
    4545.05,
    4579.35,
    4589.9,
    4609.57,
    4657.90,
    4726.87,
    4764.86,
    4806.02,
    4847.81,
    4879.66,
    4965.08,
    5017.18,
    5105.54,
    5153.23,
    5218.20,
    5330.78,
    5341.09,
    5400.56,
]
elements = ["CuNeAr"] * len(atlas)

spatial_mask = np.arange(450, 550)
spec_mask = np.arange(500, 3200)


if os.path.exists("reduced_image_g93m48_blue.fits"):

    standard_frame = fits.open("reduced_image_g93m48_blue.fits")

else:

    standard_frame = image_reduction.ImageReduction()
    standard_frame.add_filelist("isis_g93m48_blue.list")
    standard_frame.set_properties(saxis=Saxis)
    standard_frame.load_data()
    standard_frame.reduce()
    standard_frame.inspect(
        display=False,
        save_fig=True,
        fig_type="png",
        filename="reduced_image_g93m48_blue",
    )
    standard_frame.save_fits(filename="reduced_image_g93m48_blue")


g93 = spectral_reduction.TwoDSpec(
    standard_frame,
    spatial_mask=spatial_mask,
    spec_mask=spec_mask,
    readnoise=4.5,
    cosmicray=False,
    gain=0.98,
    airmass=1.132391,
    seeing=1.1,
    verbose=False,
)


if g93.arc is None:
    g93_arc = fits.open(
        "isis_pso1801p6254_g93m48_raw/wht20180801_02701260.fits.fz"
    )[1]
    g93.add_arc(g93_arc.data.T, g93_arc.header)
    g93.apply_mask_to_arc()

# Update for only the science frame
# spec mask
spatial_mask = np.arange(410, 460)
spec_mask = np.arange(800, 3200)


if os.path.exists("reduced_image_pso1801p6254_blue.fits"):

    science_frame = fits.open("reduced_image_pso1801p6254_blue.fits")

else:

    science_frame = image_reduction.ImageReduction()
    science_frame.add_filelist("isis_pso1801p6254_blue.list")
    science_frame.set_properties(cosmicray=True, fsmode="median", saxis=Saxis)
    science_frame.load_data()
    science_frame.reduce()
    science_frame.inspect(
        display=False,
        save_fig=True,
        fig_type="png",
        filename="reduced_image_pso1801p6254_blue",
    )
    science_frame.save_fits(filename="reduced_image_pso1801p6254_blue")


# initialise the two spectral_reduction.TwoDSpec()
pso = spectral_reduction.TwoDSpec(
    science_frame,
    spatial_mask=spatial_mask,
    spec_mask=spec_mask,
    readnoise=4.5,
    cosmicray=True,
    psfsize=11,
    gain=0.98,
    airmass=1.498595,
    seeing=1.1,
    verbose=False,
)
if pso.arc is None:
    pso_arc = fits.open(
        "isis_pso1801p6254_g93m48_raw/wht20180801_02701246.fits.fz"
    )[1]
    pso.add_arc(pso_arc.data.T, pso_arc.header)
    pso.apply_mask_to_arc()

# automatically trace the spectrum
pso.ap_trace(
    display=False,
    nwindow=8,
    percentile=10,
    shift_tol=10.0,
    ap_faint=10,
    fit_deg=1,
    save_fig=True,
    fig_type="png",
    filename="pso_trace_blue",
)


g93.ap_trace(
    display=False,
    fit_deg=1,
    save_fig=True,
    fig_type="png+browser",
    filename="g93_trace_blue",
)


# Optimal extracting spectrum by summing over the aperture along the trace
pso.ap_extract(
    apwidth=5,
    skysep=3,
    skywidth=5,
    optimal=True,
    display=False,
    save_fig=True,
    fig_type="png",
    filename="pso_extract_blue",
)
g93.ap_extract(
    apwidth=20,
    skysep=3,
    skywidth=5,
    skydeg=1,
    optimal=True,
    display=False,
    save_fig=True,
    fig_type="png",
    filename="g93_extract_blue",
)

pso.extract_arc_spec(
    display=False, fig_type="png", filename="arc_spec_science_blue"
)
g93.extract_arc_spec(
    display=False, fig_type="png", filename="arc_spec_standard_blue"
)


pso_reduced = spectral_reduction.OneDSpec()
pso_reduced.from_twodspec(pso, stype="science")
pso_reduced.from_twodspec(g93, stype="standard")


pso_reduced.find_arc_lines(
    prominence=6.0,
    display=False,
    save_fig=True,
    stype="science+standard",
    fig_type="png",
    filename="arc_lines_blue",
)


pso_reduced.initialise_calibrator(stype="science+standard")
pso_reduced.set_hough_properties(
    num_slopes=5000,
    xbins=500,
    ybins=500,
    min_wavelength=3200.0,
    max_wavelength=5500.0,
    range_tolerance=250,
    stype="science+standard",
)
pso_reduced.set_ransac_properties(minimum_matches=38)
pso_reduced.add_user_atlas(
    elements=elements,
    wavelengths=atlas,
    vacuum=True,
    stype="science+standard",
)
pso_reduced.do_hough_transform(brute_force=False)

pso_reduced.fit(
    fit_deg=4,
    max_tries=1000,
    fit_tolerance=10.0,
    candidate_tolerance=3.0,
    display=False,
    save_fig=True,
    fig_type="png",
    stype="science+standard",
)


pso_reduced.apply_wavelength_calibration(stype="science+standard")


pso_reduced.load_standard(
    target="g93_48",
    library="esowdstan",
    cutoff=0.4,
)
pso_reduced.inspect_standard(
    fig_type="png", filename="literature_standard_blue"
)


pso_reduced.get_sensitivity(
    lowess_frac=0.2, method="polynomial", sens_deg=15, smooth=True
)
pso_reduced.inspect_sensitivity(
    fig_type="png", save_fig=True, filename="sensitivity_blue"
)

pso_reduced.apply_flux_calibration(stype="science+standard")


pso_reduced.set_atmospheric_extinction(location="orm")
pso_reduced.apply_atmospheric_extinction_correction()

pso_reduced.inspect_reduced_spectrum(
    wave_min=3000.0,
    wave_max=5500.0,
    stype="science",
    save_fig=True,
    fig_type="png",
    filename="pso_reduced_spectrum_blue",
)

pso_reduced.inspect_reduced_spectrum(
    wave_min=3000.0,
    wave_max=5500.0,
    stype="standard",
    save_fig=True,
    fig_type="png",
    filename="g93_reduced_spectrum_blue",
)

pso_reduced.create_fits(
    output="flux_resampled_atm_ext_corrected",
)

pso_reduced.save_fits(
    output="flux_resampled_atm_ext_corrected",
    filename="pso_reduced_blue",
    overwrite=True,
)
