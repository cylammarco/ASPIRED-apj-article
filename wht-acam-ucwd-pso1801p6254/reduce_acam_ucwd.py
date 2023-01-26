from aspired import image_reduction
from aspired import spectral_reduction
from astropy.io import fits
import numpy as np
import os


# Set the spectral and spatial direction
Saxis = 0

# spec mask
spatial_mask = np.arange(500, 1500)
spec_mask = np.arange(800, 2300)


science_frame = image_reduction.ImageReduction()
science_frame.add_filelist("wd.list")
science_frame.set_properties(cosmicray=True, fsmode="median", saxis=Saxis)
science_frame.load_data()
science_frame.reduce()
science_frame.inspect(
    display=False, save_fig=True, filename="wd"
)
# science_frame.savefits(overwrite=True)

standard_frame = image_reduction.ImageReduction()
standard_frame.add_filelist("standard.list")
standard_frame.set_properties(saxis=Saxis)
standard_frame.load_data()
standard_frame.reduce()
standard_frame.inspect(
    display=False, save_fig=True, filename="standard"
)

# initialise the two spectral_reduction.TwoDSpec()
pso = spectral_reduction.TwoDSpec(
    science_frame,
    spatial_mask=spatial_mask,
    spec_mask=spec_mask,
    readnoise=4.5,
    cosmicray=False,
    verbose=False,
)
ross640 = spectral_reduction.TwoDSpec(
    standard_frame,
    spatial_mask=spatial_mask,
    spec_mask=spec_mask,
    cosmicray=False,
    verbose=False,
)

# automatically trace the spectrum
pso.ap_trace(
    display=True,
    ap_faint=20,
    fit_deg=1,
    filename="pso_trace",
)

ross640.ap_trace(display=True, fig_type="png", filename="ross640_trace")

pso.spectrum_list[0].trace_sigma = ross640.spectrum_list[0].trace_sigma

# Optimal extracting spectrum by summing over the aperture along the trace
pso.ap_extract(
    optimal=True,
    display=True,
    fig_type="png",
    filename="pso_extract",
)
ross640.ap_extract(
    optimal=True,
    display=True,
    fig_type="png",
    filename="ross640_extract",
)

pso.extract_arc_spec(
    display=False, fig_type="png", filename="arc_spec_science"
)
ross640.extract_arc_spec(
    display=False, fig_type="png", filename="arc_spec_standard"
)


pso_reduced = spectral_reduction.OneDSpec()
pso_reduced.from_twodspec(pso, stype="science")
pso_reduced.from_twodspec(ross640, stype="standard")


pso_reduced.find_arc_lines(
    prominence=2.0,
    display=True,
    stype="science+standard",
    fig_type="png",
    filename="arc_lines",
)


atlas = [
    5852.50,
    5881.90,
    5944.83,
    5975.53,
    6030.00,
    6096.16,
    6143.06,6217.28,
    6266.50,
    630.79,
    6334.43,
    6402.25,
    6538.11,
    6604.85,
    6678.28,
    6717.04,
    6929.47,
    7032.41,
    7173.94,
    7245.17,
    7439.00,
    7488.87,
    7535.77,
    8136.41,
    8300.33,
    8377.61,
    8418.43,
    8495.36,
    8591.28,
    8654.38,
    8681.92,
    8780.62
]
elements = ["CuNe"] * len(atlas)

pso_reduced.initialise_calibrator(stype="science+standard")
pso_reduced.set_hough_properties(
    num_slopes=2000,
    xbins=200,
    ybins=200,
    min_wavelength=4000.0,
    max_wavelength=9000.0,
    range_tolerance=500,
    stype="science+standard",
)
pso_reduced.set_ransac_properties(minimum_matches=20)
pso_reduced.add_user_atlas(
    elements=elements,
    wavelengths=atlas,
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
    save_fig=True,
    fig_type="png",
    stype="science+standard",
)


pso_reduced.apply_wavelength_calibration(stype="science+standard")


pso_reduced.load_standard(
    target="r640",
    library="ing_oke",
    cutoff=0.4,
)
pso_reduced.inspect_standard(
    fig_type="png", filename="literature_standard"
)

pso_reduced.get_sensitivity()
pso_reduced.inspect_sensitivity(fig_type="png", filename="sensitivity")

pso_reduced.apply_flux_calibration(stype="science+standard")


pso_reduced.set_atmospheric_extinction(location="orm")
pso_reduced.apply_atmospheric_extinction_correction()


pso_reduced.get_telluric_profile()
pso_reduced.get_telluric_strength()
pso_reduced.apply_telluric_correction()


pso_reduced.inspect_reduced_spectrum(
    wave_min=4000.0,
    wave_max=9000.0,
    stype="science",
    save_fig=True,
    fig_type="png",
    filename="pso_reduced_spectrum",
)

pso_reduced.inspect_reduced_spectrum(
    wave_min=4000.0,
    wave_max=9000.0,
    stype="standard",
    save_fig=True,
    fig_type="png",
    filename="ross640_reduced_spectrum",
)

pso_reduced.create_fits(
    output="flux_resampled_atm_ext_telluric_corrected",
)

pso_reduced.save_fits(
    output="flux_resampled_atm_ext_telluric_corrected",
    filename="pso_reduced",
    overwrite=True,
)

np.savetxt(
    "acam-pixel-wavelength-solution-pairs.txt",
    np.column_stack(
        [
            pso_reduced.science_spectrum_list[0].calibrator.matched_peaks,
            pso_reduced.science_spectrum_list[0].calibrator.matched_atlas,
        ]
    ),
    delimiter=",",
)

np.savetxt(
    "acam-effective-pixel-spectrum.txt",
    np.column_stack(
        [
            pso_reduced.science_spectrum_list[0].pixel_list,
            pso_reduced.science_spectrum_list[0].arc_spec,
        ]
    ),
    delimiter=",",
)
