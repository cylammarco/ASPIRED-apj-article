import glob
import os
import sys
import warnings

from aspired import image_reduction
from aspired import spectral_reduction
from astropy.io import fits
import numpy as np
from scipy import signal

from gmos_fieldflattening import create_pixel_array


try:
    base_path = os.path.abspath(os.path.dirname(__file__))
except:
    base_path = os.path.abspath(os.path.dirname(__name__))


atlas_lines = [
    4609.567,
    4764.865,
    4806.02,
    4879.864,
    4965.080,
    5017.163,
    5062.037,
    5105.54,
    5153.23,
    5187.746,
    5218.20,
    5782.13,
    6677.282,
    6752.834,
    6871.289,
    6965.431,
    7030.25,
    7067.218,
    7147.04,
    7206.98,
    7272.936,
    7353.29,
    7383.981,
    7503.87,
    7514.65,
    7635.105,
    7723.98,
    7948.176,
    8006.16,
    8014.79,
    8103.69,
    8115.311,
    8264.522,
    8408.21,
    8424.648,
    8521.442,
    8667.944,
    9122.967,
    9194.64,
    9224.499,
    9291.53,
    9354.220,
    9657.786,
    9784.503,
]
elements = ["CuAr"] * len(atlas_lines)

# Get the filepaths

# R400 first
science_filepath = glob.glob("gemini_gmos_light_r400/*.fits")
science_arc_filepath = glob.glob("gemini_gmos_arc_r400/*.fits")
standard_filepath = glob.glob("gemini_gmos_standard_light_r400/*.fits")
standard_arc_filepath = glob.glob("gemini_gmos_standard_arc_r400/*.fits")

spatial_mask = np.arange(50, 600)

n_science = len(science_filepath)

science_twodspec = np.array([None] * n_science, dtype="object")
onedspec = np.array([None] * n_science, dtype="object")

standard_fits = fits.open(standard_filepath[0])[0]
standard_flattened = standard_fits.data
standard_header = standard_fits.header
standard_northsouth = standard_header["NS"].lower()
standard_binx = standard_header["BINX"]
standard_pixel_list = create_pixel_array(standard_northsouth, standard_binx)
standard_atlas_temperature = float(standard_header["TAMBIENT"]) + 273.15
standard_atlas_pressure = float(standard_header["PRESSUR2"])
standard_atlas_relative_humidity = float(standard_header["HUMIDITY"])

standard_arc_fits = fits.open(standard_arc_filepath[0])[0]
standard_arc_flattened = standard_arc_fits.data
standard_arc_header = standard_arc_fits.header
standard_arc_northsouth = standard_arc_header["NS"].lower()
standard_arc_binx = standard_arc_header["BINX"]


standard_twodspec = spectral_reduction.TwoDSpec(
    standard_fits.data,
    spatial_mask=spatial_mask,
    saxis=1,
    flip=True,
    cosmicray=True,
    readnoise=0.0,
    gain=1.64,
    psfsize=11,
)
standard_twodspec.ap_trace(
    fit_deg=2,
    save_fig=True,
    filename="gw170817_r400_standard_trace",
    display=False,
)
standard_twodspec.ap_extract(
    apwidth=25,
    save_fig=True,
    algorithm="horne86",
    filename="gw170817_r400_standard_extraction",
    display=False,
)

standard_twodspec.add_arc(np.flip(standard_arc_flattened))
standard_twodspec.extract_arc_spec()


standard_onedspec = spectral_reduction.OneDSpec()
standard_onedspec.from_twodspec(standard_twodspec, stype="standard")

standard_onedspec.find_arc_lines(prominence=0.25, stype="standard")

standard_onedspec.initialise_calibrator(
    peaks=standard_pixel_list[
        standard_onedspec.standard_spectrum_list[0].peaks
    ],
    stype="standard",
)
standard_onedspec.set_calibrator_properties(
    pixel_list=standard_pixel_list, stype="standard"
)
standard_onedspec.set_hough_properties(
    num_slopes=10000.0,
    xbins=500,
    ybins=500,
    min_wavelength=4900.0,
    max_wavelength=9500.0,
    range_tolerance=500,
    stype="standard",
)
standard_onedspec.add_user_atlas(
    elements=elements,
    wavelengths=atlas_lines,
    pressure=standard_atlas_pressure,
    temperature=standard_atlas_temperature,
    relative_humidity=standard_atlas_relative_humidity,
    vacuum=False,
    stype="standard",
)
standard_onedspec.do_hough_transform(stype="standard")
standard_onedspec.set_ransac_properties(
    sample_size=6, minimum_matches=20, stype="standard"
)
# wavelength solution for target
standard_onedspec.fit(
    max_tries=1000,
    fit_deg=4,
    candidate_tolerance=3.0,
    fit_tolerance=10.0,
    stype="standard",
)
standard_onedspec.apply_wavelength_calibration(stype="standard")

standard_onedspec.load_standard(target="EG274", library="esoxshooter")
standard_onedspec.get_sensitivity(lowess_frac=0.01)
standard_onedspec.inspect_sensitivity(
    save_fig=True, filename="gw170817_r400_sensitivity", display=False
)
standard_onedspec.apply_wavelength_calibration(stype="standard")
standard_onedspec.apply_flux_calibration(stype="standard")
standard_onedspec.get_telluric_profile()
standard_onedspec.apply_telluric_correction(stype="standard")
standard_onedspec.resample(stype="standard")
standard_onedspec.inspect_reduced_spectrum(
    save_fig=True,
    filename="gw170817_r400_standard",
    stype="standard",
    display=False,
)


for i in range(n_science):
    # Science frames must exist
    science_fits = fits.open(science_filepath[i])[0]
    science_flattened = science_fits.data
    science_header = science_fits.header
    science_northsouth = science_header["NS"].lower()
    science_binx = science_header["BINX"]
    # Initialise the 1D spec
    onedspec[i] = spectral_reduction.OneDSpec()
    science_pixel_list = create_pixel_array(science_northsouth, science_binx)
    science_atlas_temperature = float(science_header["TAMBIENT"]) + 273.15
    science_atlas_pressure = float(science_header["PRESSUR2"])
    science_atlas_relative_humidity = float(science_header["HUMIDITY"])
    # Make sure standard_noise, gain, seeing, exptime are None for the
    # above keywords to be used.
    science_twodspec[i] = spectral_reduction.TwoDSpec(
        science_flattened.copy(),
        spatial_mask=spatial_mask,
        saxis=1,
        flip=True,
        cosmicray=True,
        readnoise=0.0,
        gain=1.64,
        psfsize=11,
    )
    science_twodspec[i].ap_trace(
        fit_deg=2,
        save_fig=True,
        filename="gw170817_r400_trace_{}".format(i),
    )
    science_twodspec[i].ap_extract(
        save_fig=True,
        algorithm="horne86",
        filename="gw170817_r400_extraction_{}".format(i),
    )
    #################################################################
    # # # ##### #   # ##### #     ##### #   # ##### ##### #   #
    # # # #   # #   # #     #     #     ##  # #       #   #   #
    # # # ##### #   # ##### #     ##### # # # # ###   #   #####
    # # # #   #  # #  #     #     #     #  ## #   #   #   #   #
    ##### #   #   #   ##### ##### ##### #   # #####   #   #   #
    ##### ##### #     ##### ####  ####  ##### ##### ##### ##### #   #
    #     #   # #       #   #   # #   # #   #   #     #   #   # ##  #
    #     ##### #       #   ####  ####  #####   #     #   #   # # # #
    #     #   # #       #   #   # #  #  #   #   #     #   #   # #  ##
    ##### #   # ##### ##### ####  #   # #   #   #   ##### ##### #   #
    #################################################################
    science_arc_fits = fits.open(science_arc_filepath[0])[0]
    science_arc_flattened = science_arc_fits.data
    science_arc_header = science_arc_fits.header
    science_arc_northsouth = science_arc_header["NS"].lower()
    science_arc_binx = science_arc_header["BINX"]
    if science_arc_binx == science_binx:
        science_twodspec[i].add_arc(science_arc_flattened)
        # apply_mask handles the flip
        science_twodspec[i].apply_mask_to_arc()
        science_twodspec[i].extract_arc_spec()
        onedspec[i].from_twodspec(science_twodspec[i], stype="science")
    # Check arc binning matches other frames
    # GMOS-South didn't take any CuAr R400 in 1x1 binning since 2018
    else:
        onedspec[i].from_twodspec(science_twodspec[i], stype="science")
        binx_ratio = science_arc_binx / science_binx
        spec = onedspec.science_twodspec[0]
        len_trace = len(spec.trace)
        trace = np.nanmean(spec.trace) / binx_ratio
        trace_width = np.nanmean(spec.trace_sigma) * 3.0 / binx_ratio
        arc_trace = np.flip(science_arc_flattened)[
            max(0, int(trace - trace_width - 1)) : min(
                int(trace + trace_width), len_trace
            ),
            :,
        ]
        arc_spec = np.nanmedian(arc_trace, axis=0)
        arc_spec_resampled = signal.resample(arc_spec, len_trace)
        onedspec[i].add_arc_spec(arc_spec=arc_spec_resampled)
    onedspec[i].find_arc_lines(prominence=0.25, stype="science")
    onedspec[i].initialise_calibrator(
        peaks=science_pixel_list[onedspec[i].science_spectrum_list[0].peaks],
        stype="science",
    )
    onedspec[i].set_calibrator_properties(
        pixel_list=science_pixel_list, stype="science"
    )
    onedspec[i].set_hough_properties(
        num_slopes=10000.0,
        xbins=500,
        ybins=500,
        min_wavelength=4900.0,
        max_wavelength=9500.0,
        range_tolerance=500,
        stype="science",
    )
    onedspec[i].add_user_atlas(
        elements=elements,
        wavelengths=atlas_lines,
        pressure=science_atlas_pressure,
        temperature=science_atlas_temperature,
        relative_humidity=science_atlas_relative_humidity,
        vacuum=False,
        stype="science",
    )
    onedspec[i].do_hough_transform(stype="science")
    onedspec[i].set_ransac_properties(
        sample_size=6, minimum_matches=15, stype="science"
    )
    # wavelength solution for target
    onedspec[i].fit(
        max_tries=1000,
        fit_deg=4,
        candidate_tolerance=3.0,
        fit_tolerance=10.0,
        stype="science",
        save_fig=True,
        filename="gw170817_r400_wavelength_solution_{}".format(i),
    )
    onedspec[i].apply_wavelength_calibration(stype="science")
    onedspec[i].standard_spectrum_list[
        0
    ] = standard_onedspec.standard_spectrum_list[0]
    onedspec[i].fluxcal = standard_onedspec.fluxcal
    onedspec[i].sensitivity_curve_available = True
    onedspec[i].apply_flux_calibration(stype="science")
    onedspec[i].get_continuum()
    onedspec[i].get_telluric_profile()
    onedspec[i].get_telluric_strength()
    onedspec[i].apply_telluric_correction(stype="science")
    onedspec[i].set_atmospheric_extinction(location="ls")
    onedspec[i].apply_atmospheric_extinction_correction()
    onedspec[i].resample()

    onedspec[i].inspect_reduced_spectrum(
        save_fig=True,
        filename="gw170817_r400_{}".format(i),
        stype="science",
        display=False,
    )
    onedspec[i].create_fits(
        output="flux_resampled_atm_ext_telluric_corrected",
    )
    onedspec[i].save_fits(
        output="flux_resampled_atm_ext_telluric_corrected",
        filename="gw170817_r400_{}".format(i),
        stype="science",
        overwrite=True,
    )


# Now the B600
# Now the B600
# Now the B600
# Now the B600
# Now the B600
# Now the B600
# Now the B600
# Now the B600
# Now the B600
# Now the B600
# Now the B600
# Now the B600
# Now the B600
# Now the B600
# Now the B600
# Now the B600
# Now the B600

atlas_lines = [
    3925.72,
    3946.10,
    4042.89,
    4072.11,
    4103.91,
    4131.72,
    4158.59,
    4181.88,
    4190.87,
    4200.67,
    4259.36,
    4300.1,
    4309.24,
    4331.20,
    4348.06,
    4370.75,
    4379.67,
    4400.76,
    4426.0,
    4481.81,
    4510.73,
    4545.05,
    4579.35,
    4589.90,
    4609.567,
    4726.87,
    4735.91,
    4764.865,
    4806.02,
    4847.81,
    4879.864,
    4965.080,
    5017.163,
    5495.87,
    5558.70,
    5606.73,
    5650.70,
    5782.13,
    6032.13,
    6043.22,
    6059.37,
    6114.92,
    6145.44,
    6172.28,
    6212.50,
    6243.12,
    6296.87,
    6307.66,
    6384.72,
    6416.31,
    6677.282,
    6752.834,
    6871.289,
    6965.431,
]
elements = ["CuAr"] * len(atlas_lines)

spatial_mask = np.arange(50, 600)

science_filepath = glob.glob("gemini_gmos_light_b600/*.fits")
science_arc_filepath = glob.glob("gemini_gmos_arc_b600/*.fits")
standard_filepath = glob.glob("gemini_gmos_standard_light_b600/*.fits")
standard_arc_filepath = glob.glob("gemini_gmos_standard_arc_b600/*.fits")

n_science = len(science_filepath)

science_twodspec = np.array([None] * n_science, dtype="object")
onedspec = np.array([None] * n_science, dtype="object")

standard_fits = fits.open(standard_filepath[0])[0]
standard_flattened = standard_fits.data
standard_header = standard_fits.header
standard_northsouth = standard_header["NS"].lower()
standard_binx = standard_header["BINX"]
standard_pixel_list = create_pixel_array(standard_northsouth, standard_binx)
standard_atlas_temperature = float(standard_header["TAMBIENT"]) + 273.15
standard_atlas_pressure = float(standard_header["PRESSUR2"])
standard_atlas_relative_humidity = float(standard_header["HUMIDITY"])

standard_arc_fits = fits.open(standard_arc_filepath[0])[0]
standard_arc_flattened = standard_arc_fits.data
standard_arc_header = standard_arc_fits.header
standard_arc_northsouth = standard_arc_header["NS"].lower()
standard_arc_binx = standard_arc_header["BINX"]


standard_twodspec = spectral_reduction.TwoDSpec(
    standard_fits.data,
    spatial_mask=spatial_mask,
    saxis=1,
    flip=True,
    cosmicray=True,
    readnoise=0.0,
    gain=1.64,
    psfsize=11,
)

standard_twodspec.ap_trace(fit_deg=2)
standard_twodspec.ap_extract(algorithm="horne86")

standard_twodspec.add_arc(np.flip(standard_arc_flattened))
standard_twodspec.extract_arc_spec()


standard_onedspec = spectral_reduction.OneDSpec()
standard_onedspec.from_twodspec(standard_twodspec, stype="standard")

standard_onedspec.find_arc_lines(prominence=0.5, stype="standard")

standard_onedspec.initialise_calibrator(
    peaks=standard_pixel_list[
        standard_onedspec.standard_spectrum_list[0].peaks
    ],
    stype="standard",
)
standard_onedspec.set_calibrator_properties(
    pixel_list=standard_pixel_list, stype="standard"
)
standard_onedspec.set_hough_properties(
    num_slopes=5000,
    xbins=250,
    ybins=250,
    min_wavelength=3800.0,
    max_wavelength=7000.0,
    range_tolerance=200,
    stype="standard",
)
standard_onedspec.add_user_atlas(
    elements=elements,
    wavelengths=atlas_lines,
    pressure=standard_atlas_pressure,
    temperature=standard_atlas_temperature,
    relative_humidity=standard_atlas_relative_humidity,
    vacuum=False,
    stype="standard",
)
standard_onedspec.do_hough_transform(stype="standard")
standard_onedspec.set_ransac_properties(
    sample_size=5, minimum_matches=16, stype="standard"
)
# wavelength solution for target
standard_onedspec.fit(
    max_tries=1000,
    fit_deg=2,
    candidate_tolerance=10.0,
    fit_tolerance=25.0,
    stype="standard",
    display=False,
)
standard_onedspec.apply_wavelength_calibration(stype="standard")

standard_onedspec.load_standard(target="EG274", library="esoxshooter")
standard_onedspec.get_sensitivity(lowess_frac=0)
standard_onedspec.inspect_sensitivity(
    save_fig=True, filename="gw170817_b600_sensitivity", display=False
)

for i in range(n_science):
    # Science frames must exist
    science_fits = fits.open(science_filepath[i])[0]
    science_flattened = science_fits.data
    science_header = science_fits.header
    science_northsouth = science_header["NS"].lower()
    science_binx = science_header["BINX"]
    # Initialise the 1D spec
    onedspec[i] = spectral_reduction.OneDSpec()
    science_pixel_list = create_pixel_array(science_northsouth, science_binx)
    science_atlas_temperature = float(science_header["TAMBIENT"]) + 273.15
    science_atlas_pressure = float(science_header["PRESSUR2"])
    science_atlas_relative_humidity = float(science_header["HUMIDITY"])
    # Make sure standard_noise, gain, seeing, exptime are None for the
    # above keywords to be used.
    science_twodspec[i] = spectral_reduction.TwoDSpec(
        science_flattened.copy(),
        spatial_mask=spatial_mask,
        saxis=1,
        flip=True,
        cosmicray=True,
        readnoise=0.0,
        gain=1.64,
        psfsize=11,
    )
    science_twodspec[i].ap_trace(
        fit_deg=2,
        save_fig=True,
        filename="gw170817_b600_trace_{}".format(i),
    )
    science_twodspec[i].ap_extract(
        save_fig=True,
        algorithm="horne86",
        filename="gw170817_b600_extraction_{}".format(i),
    )
    #################################################################
    # # # ##### #   # ##### #     ##### #   # ##### ##### #   #
    # # # #   # #   # #     #     #     ##  # #       #   #   #
    # # # ##### #   # ##### #     ##### # # # # ###   #   #####
    # # # #   #  # #  #     #     #     #  ## #   #   #   #   #
    ##### #   #   #   ##### ##### ##### #   # #####   #   #   #
    ##### ##### #     ##### ####  ####  ##### ##### ##### ##### #   #
    #     #   # #       #   #   # #   # #   #   #     #   #   # ##  #
    #     ##### #       #   ####  ####  #####   #     #   #   # # # #
    #     #   # #       #   #   # #  #  #   #   #     #   #   # #  ##
    ##### #   # ##### ##### ####  #   # #   #   #   ##### ##### #   #
    #################################################################
    science_arc_fits = fits.open(science_arc_filepath[0])[0]
    science_arc_flattened = science_arc_fits.data
    science_arc_header = science_arc_fits.header
    science_arc_northsouth = science_arc_header["NS"].lower()
    science_arc_binx = science_arc_header["BINX"]
    if science_arc_binx == science_binx:
        science_twodspec[i].add_arc(science_arc_flattened)
        # apply_mask handles the flip
        science_twodspec[i].apply_mask_to_arc()
        science_twodspec[i].extract_arc_spec()
        onedspec[i].from_twodspec(science_twodspec[i], stype="science")
    # Check arc binning matches other frames
    # GMOS-South didn't take any CuAr R400 in 1x1 binning since 2018
    else:
        onedspec[i].from_twodspec(science_twodspec[i], stype="science")
        binx_ratio = science_arc_binx / science_binx
        spec = onedspec.science_twodspec[0]
        len_trace = len(spec.trace)
        trace = np.nanmean(spec.trace) / binx_ratio
        trace_width = np.nanmean(spec.trace_sigma) * 3.0 / binx_ratio
        arc_trace = np.flip(science_arc_flattened)[
            max(0, int(trace - trace_width - 1)) : min(
                int(trace + trace_width), len_trace
            ),
            :,
        ]
        arc_spec = np.nanmedian(arc_trace, axis=0)
        arc_spec_resampled = signal.resample(arc_spec, len_trace)
        onedspec[i].add_arc_spec(arc_spec=arc_spec_resampled)
    onedspec[i].find_arc_lines(prominence=0.25, stype="science")
    onedspec[i].initialise_calibrator(
        peaks=science_pixel_list[onedspec[i].science_spectrum_list[0].peaks],
        stype="science",
    )
    onedspec[i].set_calibrator_properties(
        pixel_list=science_pixel_list, stype="science"
    )
    onedspec[i].set_hough_properties(
        num_slopes=5000,
        xbins=250,
        ybins=250,
        min_wavelength=3800.0,
        max_wavelength=7000.0,
        range_tolerance=200,
        stype="science",
    )
    onedspec[i].add_user_atlas(
        elements=elements,
        wavelengths=atlas_lines,
        pressure=science_atlas_pressure,
        temperature=science_atlas_temperature,
        relative_humidity=science_atlas_relative_humidity,
        vacuum=False,
        stype="science",
    )
    onedspec[i].do_hough_transform(stype="science")
    onedspec[i].set_ransac_properties(
        sample_size=5, minimum_matches=16, stype="science"
    )
    # wavelength solution for target
    onedspec[i].fit(
        max_tries=1000,
        fit_deg=2,
        candidate_tolerance=10.0,
        fit_tolerance=25.0,
        stype="science",
        display=False,
        save_fig=True,
        filename="gw170817_b600_wavelength_solution_{}".format(i),
    )
    onedspec[i].apply_wavelength_calibration(stype="science")
    onedspec[i].standard_spectrum_list[
        0
    ] = standard_onedspec.standard_spectrum_list[0]
    onedspec[i].fluxcal = standard_onedspec.fluxcal
    onedspec[i].sensitivity_curve_available = True
    onedspec[i].apply_flux_calibration()
    onedspec[i].set_atmospheric_extinction(location="ls")
    onedspec[i].apply_atmospheric_extinction_correction()
    onedspec[i].resample()
    onedspec[i].inspect_reduced_spectrum(
        save_fig=True,
        filename="gw170817_b600_{}".format(i),
        stype="science",
        display=False,
    )
    onedspec[i].create_fits(
        output="flux_resampled_atm_ext_corrected",
        stype="science",
    )
    onedspec[i].save_fits(
        output="flux_resampled_atm_ext_corrected",
        filename="gw170817_b600_{}".format(i),
        stype="science",
        overwrite=True,
    )

np.savetxt(
    "gmos-pixel-wavelength-solution-pairs.txt",
    np.column_stack(
        [
            onedspec[i].science_spectrum_list[0].calibrator.matched_peaks,
            onedspec[i].science_spectrum_list[0].calibrator.matched_atlas,
        ]
    ),
    delimiter=",",
)

np.savetxt(
    "gmos-effective-pixel-spectrum.txt",
    np.column_stack(
        [
            onedspec[i].science_spectrum_list[0].pixel_list,
            onedspec[i].science_spectrum_list[0].arc_spec,
        ]
    ),
    delimiter=",",
)
