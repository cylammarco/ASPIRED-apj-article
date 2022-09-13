from astropy.nddata.decorators import SUPPORTED_PROPERTIES
import copy
import numpy as np
from aspired import spectral_reduction
from aspired import image_reduction
from spectres import spectres
from matplotlib import pyplot as plt

spatial_mask = np.arange(500, 1500)
spec_mask_red = np.arange(50, 2100)
spec_mask_blue = np.arange(250, 2100)

# Science frame
dmwd_red_1_image_reduction = image_reduction.ImageReduction(
    log_level="INFO", log_file_name=None
)
dmwd_red_1_image_reduction.add_filelist("dmwd_red_1.list")
dmwd_red_1_image_reduction.set_detector_properties(gain="E_GAIN")
dmwd_red_1_image_reduction.load_data()
dmwd_red_1_image_reduction.reduce()

dmwd_red_2_image_reduction = image_reduction.ImageReduction(
    log_level="INFO", log_file_name=None
)
dmwd_red_2_image_reduction.add_filelist("dmwd_red_2.list")
dmwd_red_2_image_reduction.set_detector_properties(gain="E_GAIN")
dmwd_red_2_image_reduction.load_data()
dmwd_red_2_image_reduction.reduce()

dmwd_red_3_image_reduction = image_reduction.ImageReduction(
    log_level="INFO", log_file_name=None
)
dmwd_red_3_image_reduction.add_filelist("dmwd_red_3.list")
dmwd_red_3_image_reduction.set_detector_properties(gain="E_GAIN")
dmwd_red_3_image_reduction.load_data()
dmwd_red_3_image_reduction.reduce()

dmwd_blue_1_image_reduction = image_reduction.ImageReduction(
    log_level="INFO", log_file_name=None
)
dmwd_blue_1_image_reduction.add_filelist("dmwd_blue_1.list")
dmwd_blue_1_image_reduction.set_detector_properties(gain="E_GAIN")
dmwd_blue_1_image_reduction.load_data()
dmwd_blue_1_image_reduction.reduce()

dmwd_blue_2_image_reduction = image_reduction.ImageReduction(
    log_level="INFO", log_file_name=None
)
dmwd_blue_2_image_reduction.add_filelist("dmwd_blue_2.list")
dmwd_blue_2_image_reduction.set_detector_properties(gain="E_GAIN")
dmwd_blue_2_image_reduction.load_data()
dmwd_blue_2_image_reduction.reduce()

dmwd_blue_3_image_reduction = image_reduction.ImageReduction(
    log_level="INFO", log_file_name=None
)
dmwd_blue_3_image_reduction.add_filelist("dmwd_blue_3.list")
dmwd_blue_3_image_reduction.set_detector_properties(gain="E_GAIN")
dmwd_blue_3_image_reduction.load_data()
dmwd_blue_3_image_reduction.reduce()

dmwd_blue_4_image_reduction = image_reduction.ImageReduction(
    log_level="INFO", log_file_name=None
)
dmwd_blue_4_image_reduction.add_filelist("dmwd_blue_4.list")
dmwd_blue_4_image_reduction.set_detector_properties(gain="E_GAIN")
dmwd_blue_4_image_reduction.load_data()
dmwd_blue_4_image_reduction.reduce()

dmwd_red_1_twodspec = spectral_reduction.TwoDSpec(
    dmwd_red_1_image_reduction,
    spatial_mask=spatial_mask,
    spec_mask=spec_mask_red,
    cosmic=True,
    log_level="INFO",
    log_file_name=None,
)

dmwd_red_2_twodspec = spectral_reduction.TwoDSpec(
    dmwd_red_2_image_reduction,
    spatial_mask=spatial_mask,
    spec_mask=spec_mask_red,
    cosmic=True,
    log_level="INFO",
    log_file_name=None,
)

dmwd_red_3_twodspec = spectral_reduction.TwoDSpec(
    dmwd_red_3_image_reduction,
    spatial_mask=spatial_mask,
    spec_mask=spec_mask_red,
    cosmic=True,
    log_level="INFO",
    log_file_name=None,
)

dmwd_blue_1_twodspec = spectral_reduction.TwoDSpec(
    dmwd_blue_1_image_reduction,
    spatial_mask=spatial_mask,
    spec_mask=spec_mask_blue,
    cosmic=True,
    log_level="INFO",
    log_file_name=None,
)

dmwd_blue_2_twodspec = spectral_reduction.TwoDSpec(
    dmwd_blue_2_image_reduction,
    spatial_mask=spatial_mask,
    spec_mask=spec_mask_blue,
    cosmic=True,
    log_level="INFO",
    log_file_name=None,
)

dmwd_blue_3_twodspec = spectral_reduction.TwoDSpec(
    dmwd_blue_3_image_reduction,
    spatial_mask=spatial_mask,
    spec_mask=spec_mask_blue,
    cosmic=True,
    log_level="INFO",
    log_file_name=None,
)

dmwd_blue_4_twodspec = spectral_reduction.TwoDSpec(
    dmwd_blue_4_image_reduction,
    spatial_mask=spatial_mask,
    spec_mask=spec_mask_blue,
    cosmic=True,
    log_level="INFO",
    log_file_name=None,
)

dmwd_red_1_twodspec.ap_trace(nspec=3, fit_deg=2, display=False)
dmwd_red_1_twodspec.ap_extract(spec_id=[0], apwidth=8, display=False)
dmwd_red_1_twodspec.ap_extract(
    spec_id=[1], apwidth=6, skysep=[3, 15], display=False
)
dmwd_red_1_twodspec.ap_extract(
    spec_id=[2], apwidth=6, skysep=[15, 3], display=False
)

dmwd_red_2_twodspec.ap_trace(nspec=3, fit_deg=2, display=False)
dmwd_red_2_twodspec.ap_extract(spec_id=[0], apwidth=8, display=False)
dmwd_red_2_twodspec.ap_extract(
    spec_id=[1], apwidth=6, skysep=[3, 15], display=False
)
dmwd_red_2_twodspec.ap_extract(
    spec_id=[2], apwidth=6, skysep=[15, 3], display=False
)

dmwd_red_3_twodspec.ap_trace(nspec=3, fit_deg=2, display=False)
dmwd_red_3_twodspec.ap_extract(spec_id=[0], apwidth=8, display=False)
dmwd_red_3_twodspec.ap_extract(
    spec_id=[1], apwidth=6, skysep=[3, 15], display=False
)
dmwd_red_3_twodspec.ap_extract(
    spec_id=[2], apwidth=6, skysep=[15, 3], display=False
)

dmwd_blue_1_twodspec.ap_trace(nspec=3, fit_deg=2, display=False)
dmwd_blue_1_twodspec.ap_extract(spec_id=[0], apwidth=8, display=False)
dmwd_blue_1_twodspec.ap_extract(
    spec_id=[1], apwidth=6, skysep=[3, 15], display=False
)
dmwd_blue_1_twodspec.ap_extract(
    spec_id=[2], apwidth=6, skysep=[15, 3], display=False
)

dmwd_blue_2_twodspec.ap_trace(nspec=3, fit_deg=2, display=False)
dmwd_blue_2_twodspec.ap_extract(spec_id=[0], apwidth=8, display=False)
dmwd_blue_2_twodspec.ap_extract(
    spec_id=[1], apwidth=6, skysep=[3, 15], display=False
)
dmwd_blue_2_twodspec.ap_extract(
    spec_id=[2], apwidth=6, skysep=[15, 3], display=False
)

dmwd_blue_3_twodspec.ap_trace(nspec=3, fit_deg=2, display=False)
dmwd_blue_3_twodspec.ap_extract(spec_id=[0], apwidth=8, display=False)
dmwd_blue_3_twodspec.ap_extract(
    spec_id=[1], apwidth=6, skysep=[3, 15], display=False
)
dmwd_blue_3_twodspec.ap_extract(
    spec_id=[2], apwidth=6, skysep=[15, 3], display=False
)

dmwd_blue_4_twodspec.ap_trace(nspec=3, fit_deg=2, display=False)
dmwd_blue_4_twodspec.ap_extract(spec_id=[0], apwidth=8, display=False)
dmwd_blue_4_twodspec.ap_extract(
    spec_id=[1], apwidth=6, skysep=[3, 15], display=False
)
dmwd_blue_4_twodspec.ap_extract(
    spec_id=[2], apwidth=6, skysep=[15, 3], display=False
)

dmwd_red_1_twodspec.extract_arc_spec(display=False)
dmwd_red_2_twodspec.extract_arc_spec(display=False)
dmwd_red_3_twodspec.extract_arc_spec(display=False)
dmwd_blue_1_twodspec.extract_arc_spec(display=False)
dmwd_blue_2_twodspec.extract_arc_spec(display=False)
dmwd_blue_3_twodspec.extract_arc_spec(display=False)
dmwd_blue_4_twodspec.extract_arc_spec(display=False)

# Standard frames
BDp253941_red_image_reduction = image_reduction.ImageReduction(
    log_level="INFO", log_file_name=None
)
BDp253941_red_image_reduction.add_filelist("BDp33d2642_red.list")
BDp253941_red_image_reduction.set_properties(
    cosmicray=False, psfmodel="gauss", gain="E_GAIN"
)
BDp253941_red_image_reduction.load_data()
BDp253941_red_image_reduction.reduce()

BDp253941_red_twodspec = spectral_reduction.TwoDSpec(
    BDp253941_red_image_reduction,
    spatial_mask=spatial_mask,
    spec_mask=spec_mask_red,
    cosmicray=True,
    log_level="INFO",
    log_file_name=None,
)
BDp253941_red_twodspec.ap_trace(nspec=1, display=False)
BDp253941_red_twodspec.ap_extract(apwidth=15, display=False)
BDp253941_red_twodspec.extract_arc_spec(display=False)

BDp253941_blue_image_reduction = image_reduction.ImageReduction(
    log_level="INFO", log_file_name=None
)
BDp253941_blue_image_reduction.add_filelist("BDp253941_blue.list")
BDp253941_blue_image_reduction.set_properties(
    cosmicray=False, psfmodel="gauss", gain="E_GAIN"
)
BDp253941_blue_image_reduction.load_data()
BDp253941_blue_image_reduction.reduce()

BDp253941_blue_twodspec = spectral_reduction.TwoDSpec(
    BDp253941_blue_image_reduction,
    spatial_mask=spatial_mask,
    spec_mask=spec_mask_blue,
    cosmicray=True,
    log_level="INFO",
    log_file_name=None,
)
BDp253941_blue_twodspec.ap_trace(nspec=1, display=False)
BDp253941_blue_twodspec.ap_extract(apwidth=15, display=False)
BDp253941_blue_twodspec.extract_arc_spec(display=False)

# Handle 1D Science spectrum
dmwd_red_1_onedspec = spectral_reduction.OneDSpec(
    log_level="INFO", log_file_name=None
)
dmwd_red_1_onedspec.from_twodspec(dmwd_red_1_twodspec, stype="science")
dmwd_red_1_onedspec.from_twodspec(
    copy.copy(BDp253941_red_twodspec), stype="standard"
)

dmwd_red_2_onedspec = spectral_reduction.OneDSpec(
    log_level="INFO", log_file_name=None
)
dmwd_red_2_onedspec.from_twodspec(dmwd_red_2_twodspec, stype="science")
dmwd_red_2_onedspec.from_twodspec(
    copy.copy(BDp253941_red_twodspec), stype="standard"
)

dmwd_red_3_onedspec = spectral_reduction.OneDSpec(
    log_level="INFO", log_file_name=None
)
dmwd_red_3_onedspec.from_twodspec(dmwd_red_3_twodspec, stype="science")
dmwd_red_3_onedspec.from_twodspec(
    copy.copy(BDp253941_red_twodspec), stype="standard"
)

dmwd_blue_1_onedspec = spectral_reduction.OneDSpec(
    log_level="INFO", log_file_name=None
)
dmwd_blue_1_onedspec.from_twodspec(dmwd_blue_1_twodspec, stype="science")
dmwd_blue_1_onedspec.from_twodspec(
    copy.copy(BDp253941_blue_twodspec), stype="standard"
)

dmwd_blue_2_onedspec = spectral_reduction.OneDSpec(
    log_level="INFO", log_file_name=None
)
dmwd_blue_2_onedspec.from_twodspec(dmwd_blue_2_twodspec, stype="science")
dmwd_blue_2_onedspec.from_twodspec(
    copy.copy(BDp253941_blue_twodspec), stype="standard"
)

dmwd_blue_3_onedspec = spectral_reduction.OneDSpec(
    log_level="INFO", log_file_name=None
)
dmwd_blue_3_onedspec.from_twodspec(dmwd_blue_3_twodspec, stype="science")
dmwd_blue_3_onedspec.from_twodspec(
    copy.copy(BDp253941_blue_twodspec), stype="standard"
)

dmwd_blue_4_onedspec = spectral_reduction.OneDSpec(
    log_level="INFO", log_file_name=None
)
dmwd_blue_4_onedspec.from_twodspec(dmwd_blue_4_twodspec, stype="science")
dmwd_blue_4_onedspec.from_twodspec(
    copy.copy(BDp253941_blue_twodspec), stype="standard"
)

# Extract arc spectrum
dmwd_red_1_onedspec.find_arc_lines(
    prominence=1.0, display=False, stype="science+standard"
)
dmwd_red_2_onedspec.find_arc_lines(
    prominence=1.0, display=False, stype="science+standard"
)
dmwd_red_3_onedspec.find_arc_lines(
    prominence=1.0, display=False, stype="science+standard"
)
dmwd_blue_1_onedspec.find_arc_lines(
    prominence=1.0, display=False, stype="science+standard"
)
dmwd_blue_2_onedspec.find_arc_lines(
    prominence=1.0, display=False, stype="science+standard"
)
dmwd_blue_3_onedspec.find_arc_lines(
    prominence=1.0, display=False, stype="science+standard"
)
dmwd_blue_4_onedspec.find_arc_lines(
    prominence=1.0, display=False, stype="science+standard"
)

# Configure the wavelength calibrator
dmwd_red_1_onedspec.initialise_calibrator(stype="science+standard")
dmwd_red_2_onedspec.initialise_calibrator(stype="science+standard")
dmwd_red_3_onedspec.initialise_calibrator(stype="science+standard")
dmwd_blue_1_onedspec.initialise_calibrator(stype="science+standard")
dmwd_blue_2_onedspec.initialise_calibrator(stype="science+standard")
dmwd_blue_3_onedspec.initialise_calibrator(stype="science+standard")
dmwd_blue_4_onedspec.initialise_calibrator(stype="science+standard")

dmwd_red_1_onedspec.set_hough_properties(
    num_slopes=5000,
    xbins=500,
    ybins=500,
    min_wavelength=5000,
    max_wavelength=10500,
    stype="science+standard",
)
dmwd_red_2_onedspec.set_hough_properties(
    num_slopes=5000,
    xbins=500,
    ybins=500,
    min_wavelength=5000,
    max_wavelength=10500,
    stype="science+standard",
)
dmwd_red_3_onedspec.set_hough_properties(
    num_slopes=5000,
    xbins=500,
    ybins=500,
    min_wavelength=5000,
    max_wavelength=10500,
    stype="science+standard",
)

dmwd_blue_1_onedspec.set_hough_properties(
    num_slopes=1000,
    xbins=100,
    ybins=100,
    min_wavelength=3000,
    max_wavelength=8500,
    stype="science+standard",
)
dmwd_blue_2_onedspec.set_hough_properties(
    num_slopes=1000,
    xbins=100,
    ybins=100,
    min_wavelength=3000,
    max_wavelength=8500,
    stype="science+standard",
)
dmwd_blue_3_onedspec.set_hough_properties(
    num_slopes=1000,
    xbins=100,
    ybins=100,
    min_wavelength=3000,
    max_wavelength=8500,
    stype="science+standard",
)
dmwd_blue_4_onedspec.set_hough_properties(
    num_slopes=1000,
    xbins=100,
    ybins=100,
    min_wavelength=3000,
    max_wavelength=8500,
    stype="science+standard",
)

dmwd_red_1_onedspec.set_ransac_properties(
    filter_close=True, stype="science+standard"
)
dmwd_red_2_onedspec.set_ransac_properties(
    filter_close=True, stype="science+standard"
)
dmwd_red_3_onedspec.set_ransac_properties(
    filter_close=True, stype="science+standard"
)
dmwd_blue_1_onedspec.set_ransac_properties(
    filter_close=True, stype="science+standard"
)
dmwd_blue_2_onedspec.set_ransac_properties(
    filter_close=True, stype="science+standard"
)
dmwd_blue_3_onedspec.set_ransac_properties(
    filter_close=True, stype="science+standard"
)
dmwd_blue_4_onedspec.set_ransac_properties(
    filter_close=True, stype="science+standard"
)

atlas = [
    5769.598,
    5820.1558,
    5852.4879,
    5944.8342,
    5975.5340,
    6029.9969,
    6074.3377,
    6096.1631,
    6143.0626,
    6163.5939,
    6217.2812,
    6266.4950,
    6334.4278,
    6382.9917,
    6506.5281,
    6532.8822,
    6598.9529,
    6678.2762,
    6717.0430,
    6929.4673,
    7173.9381,
    7245.1666,
    7272.9360,
    7383.9800,
    7435.7800,
    7488.8712,
    7587.4136,
    7601.5457,
    7948.176,
    8059.5048,
    8112.9012,
    8136.4054,
    8263.2426,
    8298.1099,
    8495.3598,
    8521.4420,
    8654.3831,
    8704.1116,
    8776.7505,
    8928.6934,
    9122.9670,
    9657.7860,
]
elements = ["ArKrNeHg"] * len(atlas)

dmwd_red_1_onedspec.add_user_atlas(
    elements=elements, wavelengths=atlas, stype="science+standard"
)
dmwd_red_2_onedspec.add_user_atlas(
    elements=elements, wavelengths=atlas, stype="science+standard"
)
dmwd_red_3_onedspec.add_user_atlas(
    elements=elements, wavelengths=atlas, stype="science+standard"
)
dmwd_blue_1_onedspec.add_atlas(elements=["He"], stype="science+standard")
dmwd_blue_2_onedspec.add_atlas(elements=["He"], stype="science+standard")
dmwd_blue_3_onedspec.add_atlas(elements=["He"], stype="science+standard")
dmwd_blue_4_onedspec.add_atlas(elements=["He"], stype="science+standard")

dmwd_red_1_onedspec.do_hough_transform()
dmwd_red_2_onedspec.do_hough_transform()
dmwd_red_3_onedspec.do_hough_transform()
dmwd_blue_1_onedspec.do_hough_transform()
dmwd_blue_2_onedspec.do_hough_transform()
dmwd_blue_3_onedspec.do_hough_transform()
dmwd_blue_4_onedspec.do_hough_transform()

# Solve for the pixel-to-wavelength solution
dmwd_red_1_onedspec.fit(
    max_tries=2000, stype="science+standard", display=False
)
dmwd_red_2_onedspec.fit(
    max_tries=2000, stype="science+standard", display=False
)
dmwd_red_3_onedspec.fit(
    max_tries=2000, stype="science+standard", display=False
)
dmwd_blue_1_onedspec.fit(
    max_tries=5000, stype="science+standard", display=False
)
dmwd_blue_2_onedspec.fit(
    max_tries=5000, stype="science+standard", display=False
)
dmwd_blue_3_onedspec.fit(
    max_tries=5000, stype="science+standard", display=False
)
dmwd_blue_4_onedspec.fit(
    max_tries=5000, stype="science+standard", display=False
)

# Apply the wavelength calibration and display it
dmwd_red_1_onedspec.apply_wavelength_calibration(stype="science+standard")
dmwd_red_2_onedspec.apply_wavelength_calibration(stype="science+standard")
dmwd_red_3_onedspec.apply_wavelength_calibration(stype="science+standard")
dmwd_blue_1_onedspec.apply_wavelength_calibration(stype="science+standard")
dmwd_blue_2_onedspec.apply_wavelength_calibration(stype="science+standard")
dmwd_blue_3_onedspec.apply_wavelength_calibration(stype="science+standard")
dmwd_blue_4_onedspec.apply_wavelength_calibration(stype="science+standard")

# Get the standard from the library
dmwd_red_1_onedspec.load_standard(target="bd33d2642", library="esohststan")
dmwd_red_2_onedspec.load_standard(target="bd33d2642", library="esohststan")
dmwd_red_3_onedspec.load_standard(target="bd33d2642", library="esohststan")
dmwd_blue_1_onedspec.load_standard(target="bd253941", library="irafirscal")
dmwd_blue_2_onedspec.load_standard(target="bd253941", library="irafirscal")
dmwd_blue_3_onedspec.load_standard(target="bd253941", library="irafirscal")
dmwd_blue_4_onedspec.load_standard(target="bd253941", library="irafirscal")

dmwd_red_1_onedspec.get_sensitivity(method="polynomial", frac=0.02)
dmwd_red_2_onedspec.get_sensitivity(method="polynomial", frac=0.02)
dmwd_red_3_onedspec.get_sensitivity(method="polynomial", frac=0.02)
dmwd_blue_1_onedspec.get_sensitivity(method="polynomial", frac=0.025)
dmwd_blue_2_onedspec.get_sensitivity(method="polynomial", frac=0.025)
dmwd_blue_3_onedspec.get_sensitivity(method="polynomial", frac=0.025)
dmwd_blue_4_onedspec.get_sensitivity(method="polynomial", frac=0.025)

dmwd_red_1_onedspec.apply_flux_calibration(stype="science+standard")
dmwd_red_2_onedspec.apply_flux_calibration(stype="science+standard")
dmwd_red_3_onedspec.apply_flux_calibration(stype="science+standard")
dmwd_blue_1_onedspec.apply_flux_calibration(stype="science+standard")
dmwd_blue_2_onedspec.apply_flux_calibration(stype="science+standard")
dmwd_blue_3_onedspec.apply_flux_calibration(stype="science+standard")
dmwd_blue_4_onedspec.apply_flux_calibration(stype="science+standard")

# Apply Telluric Correction
dmwd_red_1_onedspec.get_telluric_profile()
dmwd_red_2_onedspec.get_telluric_profile()
dmwd_red_3_onedspec.get_telluric_profile()
dmwd_blue_1_onedspec.get_telluric_profile()
dmwd_blue_2_onedspec.get_telluric_profile()
dmwd_blue_3_onedspec.get_telluric_profile()
dmwd_blue_4_onedspec.get_telluric_profile()

dmwd_red_1_onedspec.get_telluric_strength()
dmwd_red_2_onedspec.get_telluric_strength()
dmwd_red_3_onedspec.get_telluric_strength()
dmwd_blue_1_onedspec.get_telluric_strength()
dmwd_blue_2_onedspec.get_telluric_strength()
dmwd_blue_3_onedspec.get_telluric_strength()
dmwd_blue_4_onedspec.get_telluric_strength()

dmwd_red_1_onedspec.apply_telluric_correction()
dmwd_red_2_onedspec.apply_telluric_correction()
dmwd_red_3_onedspec.apply_telluric_correction()
dmwd_blue_1_onedspec.apply_telluric_correction()
dmwd_blue_2_onedspec.apply_telluric_correction()
dmwd_blue_3_onedspec.apply_telluric_correction()
dmwd_blue_4_onedspec.apply_telluric_correction()

# Apply atmospheric extinction correction
dmwd_red_1_onedspec.set_atmospheric_extinction(location="orm")
dmwd_red_2_onedspec.set_atmospheric_extinction(location="orm")
dmwd_red_3_onedspec.set_atmospheric_extinction(location="orm")
dmwd_blue_1_onedspec.set_atmospheric_extinction(location="orm")
dmwd_blue_2_onedspec.set_atmospheric_extinction(location="orm")
dmwd_blue_3_onedspec.set_atmospheric_extinction(location="orm")
dmwd_blue_4_onedspec.set_atmospheric_extinction(location="orm")

dmwd_red_1_onedspec.apply_atmospheric_extinction_correction()
dmwd_red_2_onedspec.apply_atmospheric_extinction_correction()
dmwd_red_3_onedspec.apply_atmospheric_extinction_correction()
dmwd_blue_1_onedspec.apply_atmospheric_extinction_correction()
dmwd_blue_2_onedspec.apply_atmospheric_extinction_correction()
dmwd_blue_3_onedspec.apply_atmospheric_extinction_correction()
dmwd_blue_4_onedspec.apply_atmospheric_extinction_correction()

dmwd_red_1_onedspec.inspect_reduced_spectrum(
    wave_min=5500, wave_max=10000, display=False
)
dmwd_red_2_onedspec.inspect_reduced_spectrum(
    wave_min=5500, wave_max=10000, display=False
)
dmwd_red_3_onedspec.inspect_reduced_spectrum(
    wave_min=5500, wave_max=10000, display=False
)
dmwd_blue_1_onedspec.inspect_reduced_spectrum(
    wave_min=3000, wave_max=8500, display=False
)
dmwd_blue_2_onedspec.inspect_reduced_spectrum(
    wave_min=3000, wave_max=8500, display=False
)
dmwd_blue_3_onedspec.inspect_reduced_spectrum(
    wave_min=3000, wave_max=8500, display=False
)
dmwd_blue_4_onedspec.inspect_reduced_spectrum(
    wave_min=3000, wave_max=8500, display=False
)

dmwd_red_1_onedspec.create_fits(
    output="flux_resampled_atm_ext_telluric_corrected",
)
dmwd_red_2_onedspec.create_fits(
    output="flux_resampled_atm_ext_telluric_corrected",
)
dmwd_red_3_onedspec.create_fits(
    output="flux_resampled_atm_ext_telluric_corrected",
)
dmwd_blue_1_onedspec.create_fits(
    output="flux_resampled_atm_ext_telluric_corrected",
)
dmwd_blue_2_onedspec.create_fits(
    output="flux_resampled_atm_ext_telluric_corrected",
)
dmwd_blue_3_onedspec.create_fits(
    output="flux_resampled_atm_ext_telluric_corrected",
)
dmwd_blue_4_onedspec.create_fits(
    output="flux_resampled_atm_ext_telluric_corrected",
)

dmwd_red_1_onedspec.save_fits(
    output="flux_resampled_atm_ext_telluric_corrected",
    filename="dmwd_red_1",
    overwrite=True,
)
dmwd_red_2_onedspec.save_fits(
    output="flux_resampled_atm_ext_telluric_corrected",
    filename="dmwd_red_2",
    overwrite=True,
)
dmwd_red_3_onedspec.save_fits(
    output="flux_resampled_atm_ext_telluric_corrected",
    filename="dmwd_red_3",
    overwrite=True,
)
dmwd_blue_1_onedspec.save_fits(
    output="flux_resampled_atm_ext_telluric_corrected",
    filename="dmwd_blue_1",
    overwrite=True,
)
dmwd_blue_2_onedspec.save_fits(
    output="flux_resampled_atm_ext_telluric_corrected",
    filename="dmwd_blue_2",
    overwrite=True,
)
dmwd_blue_3_onedspec.save_fits(
    output="flux_resampled_atm_ext_telluric_corrected",
    filename="dmwd_blue_3",
    overwrite=True,
)
dmwd_blue_4_onedspec.save_fits(
    output="flux_resampled_atm_ext_telluric_corrected",
    filename="dmwd_blue_4",
    overwrite=True,
)

dmwd_red_1_onedspec.resample(wave_start=5500.0, wave_end=10400.0, wave_bin=2.5)
dmwd_red_2_onedspec.resample(wave_start=5500.0, wave_end=10400.0, wave_bin=2.5)
dmwd_red_3_onedspec.resample(wave_start=5500.0, wave_end=10400.0, wave_bin=2.5)
dmwd_blue_1_onedspec.resample(wave_start=3500.0, wave_end=8000.0, wave_bin=2.5)
dmwd_blue_2_onedspec.resample(wave_start=3500.0, wave_end=8000.0, wave_bin=2.5)
dmwd_blue_3_onedspec.resample(wave_start=3500.0, wave_end=8000.0, wave_bin=2.5)
dmwd_blue_4_onedspec.resample(wave_start=3500.0, wave_end=8000.0, wave_bin=2.5)

wave_red = dmwd_red_1_onedspec.science_spectrum_list[0].wave_resampled
wave_blue = dmwd_blue_1_onedspec.science_spectrum_list[0].wave_resampled

flux_red_1 = dmwd_red_1_onedspec.science_spectrum_list[
    0
].flux_resampled_atm_ext_telluric_corrected
flux_red_2 = dmwd_red_2_onedspec.science_spectrum_list[
    0
].flux_resampled_atm_ext_telluric_corrected
flux_red_3 = dmwd_red_3_onedspec.science_spectrum_list[
    0
].flux_resampled_atm_ext_telluric_corrected
flux_blue_1 = dmwd_blue_1_onedspec.science_spectrum_list[
    0
].flux_resampled_atm_ext_telluric_corrected
flux_blue_2 = dmwd_blue_2_onedspec.science_spectrum_list[
    0
].flux_resampled_atm_ext_telluric_corrected
flux_blue_3 = dmwd_blue_3_onedspec.science_spectrum_list[
    0
].flux_resampled_atm_ext_telluric_corrected
flux_blue_4 = dmwd_blue_4_onedspec.science_spectrum_list[
    0
].flux_resampled_atm_ext_telluric_corrected

plt.figure(1, figsize=(12, 6))
plt.clf()
plt.plot(wave_red, flux_red_1 * 18)
plt.plot(wave_red, flux_red_2 * 22.5)
plt.plot(wave_red, flux_red_3 * 27.5)
plt.plot(wave_blue, flux_blue_1)
plt.plot(wave_blue, flux_blue_2)
plt.plot(wave_blue, flux_blue_3)
plt.plot(wave_blue, flux_blue_4)
plt.ylim(0, 6.5e-14)
plt.xlim(3500, 10000)
plt.grid()
plt.xlabel("Wavelength / A")
plt.ylabel("Flux")
plt.tight_layout()
plt.savefig("kstar_20210727.png")

flux_red_1 = dmwd_red_1_onedspec.science_spectrum_list[
    1
].flux_resampled_atm_ext_telluric_corrected
flux_red_2 = dmwd_red_2_onedspec.science_spectrum_list[
    1
].flux_resampled_atm_ext_telluric_corrected
flux_red_3 = dmwd_red_3_onedspec.science_spectrum_list[
    1
].flux_resampled_atm_ext_telluric_corrected
flux_blue_1 = dmwd_blue_1_onedspec.science_spectrum_list[
    1
].flux_resampled_atm_ext_telluric_corrected
flux_blue_2 = dmwd_blue_2_onedspec.science_spectrum_list[
    1
].flux_resampled_atm_ext_telluric_corrected
flux_blue_3 = dmwd_blue_3_onedspec.science_spectrum_list[
    1
].flux_resampled_atm_ext_telluric_corrected
flux_blue_4 = dmwd_blue_4_onedspec.science_spectrum_list[
    1
].flux_resampled_atm_ext_telluric_corrected

plt.figure(2, figsize=(12, 6))
plt.clf()
plt.plot(wave_red, flux_red_1 * 18)
plt.plot(wave_red, flux_red_2 * 22.5)
plt.plot(wave_red, flux_red_3 * 27.5)
plt.plot(wave_blue, flux_blue_1)
plt.plot(wave_blue, flux_blue_2)
plt.plot(wave_blue, flux_blue_3)
plt.plot(wave_blue, flux_blue_4)
plt.ylim(0, 7.0e-14)
plt.xlim(3500, 10000)
plt.grid()
plt.xlabel("Wavelength / A")
plt.ylabel("Flux")
plt.tight_layout()
plt.savefig("dmwd_20210727.png")

flux_red_1 = dmwd_red_1_onedspec.science_spectrum_list[
    2
].flux_resampled_atm_ext_telluric_corrected
flux_red_2 = dmwd_red_2_onedspec.science_spectrum_list[
    2
].flux_resampled_atm_ext_telluric_corrected
flux_red_3 = dmwd_red_3_onedspec.science_spectrum_list[
    2
].flux_resampled_atm_ext_telluric_corrected
flux_blue_1 = dmwd_blue_1_onedspec.science_spectrum_list[
    2
].flux_resampled_atm_ext_telluric_corrected
flux_blue_2 = dmwd_blue_2_onedspec.science_spectrum_list[
    2
].flux_resampled_atm_ext_telluric_corrected
flux_blue_3 = dmwd_blue_3_onedspec.science_spectrum_list[
    2
].flux_resampled_atm_ext_telluric_corrected
flux_blue_4 = dmwd_blue_4_onedspec.science_spectrum_list[
    2
].flux_resampled_atm_ext_telluric_corrected

plt.figure(3, figsize=(12, 6))
plt.clf()
plt.plot(wave_red, flux_red_1 * 18)
plt.plot(wave_red, flux_red_2 * 22.5)
plt.plot(wave_red, flux_red_3 * 27.5)
plt.plot(wave_blue, flux_blue_1)
plt.plot(wave_blue, flux_blue_2)
plt.plot(wave_blue, flux_blue_3)
plt.plot(wave_blue, flux_blue_4)
plt.ylim(0, 2e-14)
plt.xlim(3500, 10000)
plt.grid()
plt.xlabel("Wavelength / A")
plt.ylabel("Flux")
plt.tight_layout()
plt.savefig("dm_20210727.png")
