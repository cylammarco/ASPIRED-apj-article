import os

from astropy.io import fits
from scipy import ndimage

from gmos_fieldflattening import *

try:
    base_path = os.path.abspath(os.path.dirname(__file__))
except:
    base_path = os.path.abspath(os.path.dirname(__name__))

# Get the working directory and output directory
output_folder_path = base_path


# Get the folder paths R400 first
light_folder = "gemini_gmos_light_r400"
standard_light_folder = "gemini_gmos_standard_light_r400"
flat_folder = "gemini_gmos_flat"
bias_folder = "gemini_gmos_bias"

arc_folder = "gemini_gmos_arc_r400"
standard_arc_folder = "gemini_gmos_standard_arc_r400"


# Get only the files with the right extension
light_path = []
for i, path_i in enumerate(glob.glob(os.path.join(light_folder, "*.bz2"))):
    light_path.append(path_i)

standard_light_path = []
for i, path_i in enumerate(
    glob.glob(os.path.join(standard_light_folder, "*.bz2"))
):
    standard_light_path.append(path_i)

# For all cases, flat frame must exist
flat_path = []
for i, path_i in enumerate(glob.glob(os.path.join(flat_folder, "*.bz2"))):
    flat_path.append(path_i)


if os.path.exists("gemini_gmos_bias/bias_master.npy"):
    bias_master = np.load("gemini_gmos_bias/bias_master.npy")
    bias_binx = int(2048 / np.shape(bias_master)[2])
    bias_biny = int(4176 / np.shape(bias_master)[1])
else:
    # Get bias data
    bias_fits, bias_binx, bias_biny, bias_northsouth = make_bias_master(
        bias_folder
    )
    bias_master = bias_fits.data
    bias_header = bias_fits.header

# If arc is available
arc_path = []
for i, path_i in enumerate(glob.glob(os.path.join(arc_folder, "*.bz2"))):
    arc_path.append(path_i)

standard_arc_path = []
for i, path_i in enumerate(
    glob.glob(os.path.join(standard_arc_folder, "*.bz2"))
):
    standard_arc_path.append(path_i)

# Reduction stars here

# Reconstruct the light frames
light_rc_frames = []
light_header = []
light_binx = []
light_biny = []
light_northsouth = []

for i, path_i in enumerate(light_path):
    light_temp = gmos_hamamatsu(fits.open(path_i))
    light_rc_frames.append(light_temp[0])
    light_binx.append(light_temp[2])
    light_biny.append(light_temp[3])
    light_northsouth.append(light_temp[4])
    light_header.append(light_temp[5])

light_binx = light_binx[0]
light_biny = light_biny[0]

# Reconstruct the standard light frames
standard_light_rc_frames = []
standard_light_header = []
standard_light_binx = []
standard_light_biny = []
standard_light_northsouth = []

for i, path_i in enumerate(standard_light_path):
    standard_light_temp = gmos_hamamatsu(fits.open(path_i))
    standard_light_rc_frames.append(standard_light_temp[0])
    standard_light_binx.append(standard_light_temp[2])
    standard_light_biny.append(standard_light_temp[3])
    standard_light_northsouth.append(standard_light_temp[4])
    standard_light_header.append(standard_light_temp[5])

standard_light_binx = standard_light_binx[0]
standard_light_biny = standard_light_biny[0]

# Reconstruct the flat frames
flat_rc_frames = []

for i, path_i in enumerate(flat_path):
    flat_temp = gmos_hamamatsu(fits.open(path_i))
    flat_rc_frames.append(flat_temp[0])
    flat_binx = flat_temp[2]
    flat_biny = flat_temp[3]
    flat_northsouth = flat_temp[4]


# Reconstruct the arc frames
arc_rc_frames = []
arc_header = []
arc_binx = []
arc_biny = []

for i, path_i in enumerate(arc_path):
    arc_temp = gmos_hamamatsu(fits.open(path_i))
    arc_rc_frames.append(arc_temp[0])
    arc_binx.append(arc_temp[2])
    arc_biny.append(arc_temp[3])
    arc_northsouth = arc_temp[4]
    arc_header.append(arc_temp[5])

arc_binx = arc_binx[0]
arc_biny = arc_biny[0]


# Reconstruct the arc frames
standard_arc_rc_frames = []
standard_arc_header = []
standard_arc_binx = []
standard_arc_biny = []

for i, path_i in enumerate(standard_arc_path):
    standard_arc_temp = gmos_hamamatsu(fits.open(path_i))
    standard_arc_rc_frames.append(standard_arc_temp[0])
    standard_arc_binx.append(standard_arc_temp[2])
    standard_arc_biny.append(standard_arc_temp[3])
    standard_arc_northsouth = standard_arc_temp[4]
    standard_arc_header.append(standard_arc_temp[5])

standard_arc_binx = standard_arc_binx[0]
standard_arc_biny = standard_arc_biny[0]


# Bias subtraction
for frame_i in light_rc_frames:
    frame_i -= bias_master

for frame_i in standard_light_rc_frames:
    frame_i -= bias_master

for frame_i in flat_rc_frames:
    frame_i -= bias_master

# Mean combine the flat frames
# Get the relative detector response from the flats, correcting for both
# CCD response and vignetting across the focal plane
flat = np.nanmean(flat_rc_frames, axis=0)
flat_normed, flat_normalisation, _ = normalise_flat(
    flat=flat,
    binx=flat_binx,
    biny=flat_biny,
    northsouth="south",
    return_imgdata=True,
)

# Flatten the light frames
light_flattened = []
for frame_i in light_rc_frames:
    light_temp = flatten_image(
        image=frame_i,
        flat=flat_normed,
        binx=flat_binx,
        biny=flat_biny,
        normalisation=flat_normalisation,
        return_imgdata=False,
    )
    light_flattened.append(light_temp)


standard_light_flattened = []
for frame_i in standard_light_rc_frames:
    standard_light_temp = flatten_image(
        image=frame_i,
        flat=flat_normed,
        binx=flat_binx,
        biny=flat_biny,
        normalisation=flat_normalisation,
        return_imgdata=False,
    )
    standard_light_flattened.append(standard_light_temp)


# Flatten the arc frames
arc_flattened = []
for frame_i in arc_rc_frames:
    arc_temp = flatten_image(
        image=frame_i,
        flat=flat_normed,
        binx=flat_binx,
        biny=flat_biny,
        normalisation=flat_normalisation,
        return_imgdata=False,
    )
    arc_flattened.append(arc_temp)

standard_arc_flattened = []
for frame_i in standard_arc_rc_frames:
    standard_arc_temp = flatten_image(
        image=frame_i,
        flat=flat_normed,
        binx=flat_binx,
        biny=flat_biny,
        normalisation=flat_normalisation,
        return_imgdata=False,
    )
    standard_arc_flattened.append(standard_arc_temp)


frame_path = light_path + arc_path + standard_light_path + standard_arc_path
frame_header = (
    light_header + arc_header + standard_light_header + standard_arc_header
)
frame_flattened = (
    light_flattened
    + arc_flattened
    + standard_light_flattened
    + standard_arc_flattened
)

# Save as npy or FITS
for i, (frame_i, frame_header_i, name_i) in enumerate(zip(frame_flattened, frame_header, frame_path)):
    flattened_image = fits.PrimaryHDU(frame_i, header=fits.Header())
    # Add the names of the frames to header
    # flattened_image.header['COMMENT'] = "The frames."
    for i in range(len(flat_path)):
        flattened_image.header.set(
            keyword="FLAT" + str(i + 1),
            value=flat_path[i].split("/")[-1].split(".")[0],
            comment="Flat frame " + str(i + 1),
        )
    flattened_image.header.set(
        keyword="SCLIP",
        value="3",
        comment="Number of sigma used for outlier clipping.",
    )
    flattened_image.header.set(
        keyword="NS", value="south", comment="GMOS North or South."
    )
    flattened_image.header.set(
        keyword="BINX",
        value=light_binx,
        comment="Binning in the spectral direction.",
    )
    flattened_image.header.set(
        keyword="BINY",
        value=light_biny,
        comment="Binning in the spatial direction.",
    )
    flattened_image.header.set(
        keyword="PRESSUR2",
        value=frame_header_i["PRESSUR2"],
        comment="Pressure (Pa).",
    )
    flattened_image.header.set(
        keyword="TAMBIENT",
        value=frame_header_i["TAMBIENT"],
        comment="Ambient temperature (C).",
    )
    flattened_image.header.set(
        keyword="HUMIDITY",
        value=frame_header_i["HUMIDITY"],
        comment="Relative humidity (0-100%).",
    )
    flattened_image.header.set(
        keyword="EXPTIME",
        value=frame_header_i["EXPTIME"],
        comment="Exposure time",
    )
    flattened_image.header.set(
        keyword="AIRMASS",
        value=frame_header_i["AIRMASS"],
        comment="Airmass",
    )
    flattened_image.writeto(
        os.path.join(
            output_folder_path,
            name_i.split("/")[-1].split(".")[0] + "_flattened.fits",
        ),
        overwrite=True,
    )


# Repeat for B600
# Repeat for B600
# Repeat for B600
# Repeat for B600
# Repeat for B600
# Repeat for B600
# Repeat for B600
# Repeat for B600
# Repeat for B600
# Repeat for B600
# Repeat for B600
# Repeat for B600
# Repeat for B600
# Repeat for B600
# Repeat for B600
# Repeat for B600
# Repeat for B600

# Get the folder paths R400 first
light_folder = "gemini_gmos_light_b600"
standard_light_folder = "gemini_gmos_standard_light_b600"
arc_folder = "gemini_gmos_arc_b600"
standard_arc_folder = "gemini_gmos_standard_arc_b600"

flat_folder = "gemini_gmos_flat"

# Get only the files with the right extension
light_path = []
for i, path_i in enumerate(glob.glob(os.path.join(light_folder, "*.bz2"))):
    light_path.append(path_i)

# Get only the files with the right extension
standard_light_path = []
for i, path_i in enumerate(
    glob.glob(os.path.join(standard_light_folder, "*.bz2"))
):
    standard_light_path.append(path_i)

# If arc is available
arc_path = []
for i, path_i in enumerate(glob.glob(os.path.join(arc_folder, "*.bz2"))):
    arc_path.append(path_i)


standard_arc_path = []
for i, path_i in enumerate(
    glob.glob(os.path.join(standard_arc_folder, "*.bz2"))
):
    standard_arc_path.append(path_i)


# Reduction stars here

# Reconstruct the light frames
light_rc_frames = []
light_header = []
light_binx = []
light_biny = []
light_northsouth = []

for i, path_i in enumerate(light_path):
    light_temp = gmos_hamamatsu(fits.open(path_i))
    light_rc_frames.append(light_temp[0])
    light_binx.append(light_temp[2])
    light_biny.append(light_temp[3])
    light_northsouth.append(light_temp[4])
    light_header.append(light_temp[5])

light_binx = light_binx[0]
light_biny = light_biny[0]


# Reconstruct the standard light frames
standard_light_rc_frames = []
standard_light_header = []
standard_light_binx = []
standard_light_biny = []
standard_light_northsouth = []

for i, path_i in enumerate(standard_light_path):
    standard_light_temp = gmos_hamamatsu(fits.open(path_i))
    standard_light_rc_frames.append(standard_light_temp[0])
    standard_light_binx.append(standard_light_temp[2])
    standard_light_biny.append(standard_light_temp[3])
    standard_light_northsouth.append(standard_light_temp[4])
    standard_light_header.append(standard_light_temp[5])

standard_light_binx = standard_light_binx[0]
standard_light_biny = standard_light_biny[0]


# Reconstruct the arc frames
arc_rc_frames = []
arc_header = []
arc_binx = []
arc_biny = []

for i, path_i in enumerate(arc_path):
    arc_temp = gmos_hamamatsu(fits.open(path_i))
    arc_rc_frames.append(arc_temp[0])
    arc_binx.append(arc_temp[2])
    arc_biny.append(arc_temp[3])
    arc_northsouth = arc_temp[4]
    arc_header.append(arc_temp[5])

arc_binx = arc_binx[0]
arc_biny = arc_biny[0]


# Reconstruct the arc frames
standard_arc_rc_frames = []
standard_arc_header = []
standard_arc_binx = []
standard_arc_biny = []

for i, path_i in enumerate(standard_arc_path):
    standard_arc_temp = gmos_hamamatsu(fits.open(path_i))
    standard_arc_rc_frames.append(standard_arc_temp[0])
    standard_arc_binx.append(standard_arc_temp[2])
    standard_arc_biny.append(standard_arc_temp[3])
    standard_arc_northsouth = standard_arc_temp[4]
    standard_arc_header.append(standard_arc_temp[5])

standard_arc_binx = standard_arc_binx[0]
standard_arc_biny = standard_arc_biny[0]


# Bias subtraction
for frame_i in light_rc_frames:
    frame_i -= bias_master

frame_i = standard_light_rc_frames[0]
zone1 = ndimage.zoom(
    frame_i[0],
    (standard_light_biny / bias_biny, standard_light_binx / light_binx),
)
zone2 = ndimage.zoom(
    frame_i[1],
    (standard_light_biny / bias_biny, standard_light_binx / light_binx),
)
zone3 = ndimage.zoom(
    frame_i[2],
    (standard_light_biny / bias_biny, standard_light_binx / light_binx),
)
_frame_i_zone1 = zone1 - bias_master[0]
_frame_i_zone2 = zone2 - bias_master[1]
_frame_i_zone3 = zone3 - bias_master[2]
standard_light_rc_frames[0] = (_frame_i_zone1, _frame_i_zone2, _frame_i_zone3)

for frame_i in flat_rc_frames:
    frame_i -= bias_master

frame_i = standard_arc_rc_frames[0]
zone1 = ndimage.zoom(
    frame_i[0],
    (standard_arc_biny / bias_biny, standard_arc_binx / bias_binx),
)
zone2 = ndimage.zoom(
    frame_i[1],
    (standard_arc_biny / bias_biny, standard_arc_binx / bias_binx),
)
zone3 = ndimage.zoom(
    frame_i[2],
    (standard_arc_biny / bias_biny, standard_arc_binx / bias_binx),
)
_frame_i_zone1 = zone1 - bias_master[0]
_frame_i_zone2 = zone2 - bias_master[1]
_frame_i_zone3 = zone3 - bias_master[2]
standard_arc_rc_frames[0] = (_frame_i_zone1, _frame_i_zone2, _frame_i_zone3)

# Flatten the light frames
light_flattened = []
for frame_i in light_rc_frames:
    light_temp = flatten_image(
        image=frame_i,
        flat=flat_normed,
        binx=flat_binx,
        biny=flat_biny,
        normalisation=flat_normalisation,
        return_imgdata=False,
    )
    light_flattened.append(light_temp)


standard_light_flattened = []
for frame_i in standard_light_rc_frames:
    standard_light_temp = flatten_image(
        image=frame_i,
        flat=flat_normed,
        binx=flat_binx,
        biny=flat_biny,
        normalisation=flat_normalisation,
        return_imgdata=False,
    )
    standard_light_flattened.append(standard_light_temp)


# Flatten the arc frames
arc_flattened = []
for frame_i in arc_rc_frames:
    arc_temp = flatten_image(
        image=frame_i,
        flat=flat_normed,
        binx=flat_binx,
        biny=flat_biny,
        normalisation=flat_normalisation,
        return_imgdata=False,
    )
    arc_flattened.append(arc_temp)

standard_arc_flattened = []
for frame_i in standard_arc_rc_frames:
    standard_arc_temp = flatten_image(
        image=frame_i,
        flat=flat_normed,
        binx=flat_binx,
        biny=flat_biny,
        normalisation=flat_normalisation,
        return_imgdata=False,
    )
    standard_arc_flattened.append(standard_arc_temp)


frame_path = light_path + arc_path + standard_light_path + standard_arc_path
frame_header = (
    light_header + arc_header + standard_light_header + standard_arc_header
)
frame_flattened = (
    light_flattened
    + arc_flattened
    + standard_light_flattened
    + standard_arc_flattened
)

# Save as npy or FITS
for i, (frame_i, frame_header_i, name_i) in enumerate(zip(frame_flattened, frame_header, frame_path)):
    flattened_image = fits.PrimaryHDU(frame_i, header=fits.Header())
    # Add the names of the frames to header
    # flattened_image.header['COMMENT'] = "The frames."
    for i in range(len(flat_path)):
        flattened_image.header.set(
            keyword="FLAT" + str(i + 1),
            value=flat_path[i].split("/")[-1].split(".")[0],
            comment="Flat frame " + str(i + 1),
        )
    flattened_image.header.set(
        keyword="SCLIP",
        value="3",
        comment="Number of sigma used for outlier clipping.",
    )
    flattened_image.header.set(
        keyword="NS", value="south", comment="GMOS North or South."
    )
    flattened_image.header.set(
        keyword="BINX",
        value=light_binx,
        comment="Binning in the spectral direction.",
    )
    flattened_image.header.set(
        keyword="BINY",
        value=light_biny,
        comment="Binning in the spatial direction.",
    )
    flattened_image.header.set(
        keyword="PRESSUR2",
        value=frame_header_i["PRESSUR2"],
        comment="Pressure (Pa).",
    )
    flattened_image.header.set(
        keyword="TAMBIENT",
        value=frame_header_i["TAMBIENT"],
        comment="Ambient temperature (C).",
    )
    flattened_image.header.set(
        keyword="HUMIDITY",
        value=frame_header_i["HUMIDITY"],
        comment="Relative humidity (0-100%).",
    )
    flattened_image.header.set(
        keyword="EXPTIME",
        value=frame_header_i["EXPTIME"],
        comment="Exposure time",
    )
    flattened_image.header.set(
        keyword="AIRMASS",
        value=frame_header_i["AIRMASS"],
        comment="Airmass",
    )
    flattened_image.writeto(
        os.path.join(
            output_folder_path,
            name_i.split("/")[-1].split(".")[0] + "_flattened.fits",
        ),
        overwrite=True,
    )
