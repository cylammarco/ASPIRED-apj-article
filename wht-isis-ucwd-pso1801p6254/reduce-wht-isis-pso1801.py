import numpy as np
from astropy.io import fits
from aspired import image_reduction
from aspired import spectral_reduction

# Set the spectral and spatial direction
Saxis = 0

science_frame = image_reduction.ImageReduction()
science_frame.add_filelist("isis_pso1801p6254.list")
science_frame.set_properties(cosmicray=True, fsmode="median", saxis=Saxis)
science_frame.load_data()
science_frame.reduce()
science_frame.inspect(filename="reduced_image_pso1801p6254")
# science_frame.savefits(overwrite=True)

standard_frame = image_reduction.ImageReduction()
standard_frame.add_filelist("isis_g93m48.list")
standard_frame.set_properties(saxis=Saxis)
standard_frame.load_data()
standard_frame.reduce()
standard_frame.inspect(filename="reduced_image_g93m48")


# spec mask
spatial_mask = np.arange(475, 575)
spec_mask = np.arange(200, 4000)

# initialise the two spectral_reduction.TwoDSpec()
pso = spectral_reduction.TwoDSpec(
    science_frame,
    spatial_mask=spatial_mask,
    spec_mask=spec_mask,
    readnoise=4.5,
    cosmicray=False,
    gain=0.98,
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
    seeing=1.1,
    verbose=False,
)

# automatically trace the spectrum
pso.ap_trace(
    display=True,
    nwindow=15,
    percentile=25,
    shift_tol=10.0,
    ap_faint=20,
    fit_deg=1,
    
    filename="pso_trace",
)
print(np.nanmean(pso.spectrum_list[0].trace))
print(np.nanmean(pso.spectrum_list[0].trace_sigma))

g93.ap_trace(display=True, renderer="jpg+notebook", filename="g93_trace")
print(np.nanmean(g93.spectrum_list[0].trace))
print(np.nanmean(g93.spectrum_list[0].trace_sigma))


pso.spectrum_list[0].trace_sigma = g93.spectrum_list[0].trace_sigma

# Optimal extracting spectrum by summing over the aperture along the trace
pso.ap_extract(
    apwidth=15,
    optimal=True,
    display=True,
    renderer="jpg+notebook",
    filename="pso_extract",
)
g93.ap_extract(
    apwidth=20,
    skysep=3,
    skywidth=5,
    skydeg=1,
    optimal=True,
    display=True,
    renderer="jpg+notebook",
    filename="g93_extract",
)

pso.extract_arc_spec(
    display=False, renderer="jpg+notebook", filename="arc_spec_science"
)
g93.extract_arc_spec(
    display=False, renderer="jpg+notebook", filename="arc_spec_standard"
)


pso_reduced = spectral_reduction.OneDSpec()
pso_reduced.from_twodspec(pso, stype="science")
pso_reduced.from_twodspec(g93, stype="standard")


pso_reduced.find_arc_lines(
    prominence=2.0,
    display=True,
    stype="science+standard",
    renderer="jpg+notebook",
    filename="arc_lines",
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
    display=True,
    fig_type="iframe+png",
    stype="science+standard",
)


pso_reduced.apply_wavelength_calibration(stype="science+standard")


pso_reduced.load_standard(
    target="g93_48",
    library="esohststan",
    cutoff=0.4,
)
pso_reduced.inspect_standard(
    renderer="jpg+notebook", filename="literature_standard"
)


pso_reduced.get_sensitivity()
pso_reduced.inspect_sensitivity(
    renderer="jpg+notebook", filename="sensitivity"
)

pso_reduced.apply_flux_calibration(stype="science+standard")


pso_reduced.set_atmospheric_extinction(location="orm")
pso_reduced.apply_atmospheric_extinction_correction(
    science_airmass=1.498530, standard_airmass=1.13
)


pso_reduced.get_telluric_profile()
pso_reduced.get_telluric_strength()
pso_reduced.apply_telluric_correction()


pso_reduced.inspect_reduced_spectrum(
    wave_min=7000.0,
    wave_max=10500.0,
    stype="science",
    save_fig=True,
    renderer="jpg+notebook",
    filename="pso_reduced_spectrum",
)

pso_reduced.inspect_reduced_spectrum(
    wave_min=7000.0,
    wave_max=10500.0,
    stype="standard",
    save_fig=True,
    renderer="jpg+notebook",
    filename="g93_reduced_spectrum",
)

pso_reduced.create_fits(
    output="flux_resampled_atm_ext_telluric_corrected",
)

pso_reduced.save_fits(
    output="flux_resampled_atm_ext_telluric_corrected",
    filename="pso_reduced_red",
    overwrite=True,
)
