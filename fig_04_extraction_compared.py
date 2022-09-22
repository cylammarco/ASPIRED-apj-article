from aspired import image_reduction
from aspired import spectral_reduction
from astropy.io import fits
import copy
import glob
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np


# Standard frame
standard_light = fits.open(
    "vlt-fors-v418ser/standard/FORS2.2015-04-09T08_56_00.429.fits.gz"
)

flat_filelist = glob.glob("vlt-fors-v418ser/flat/*.gz")
flat_fits_list = []
for f in flat_filelist:
    flat_fits_list.append(fits.open(f))

bias_filelist = glob.glob("vlt-fors-v418ser/bias/*.gz")
bias_fits_list = []
for f in bias_filelist:
    bias_fits_list.append(fits.open(f))


standard_frame = image_reduction.ImageReduction(log_file_name=None)
standard_frame.add_light(
    standard_light[0].data,
    standard_light[0].header,
    standard_light[0].header["EXPTIME"],
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
    fig_type="jpg",
    display=False,
)

spatial_mask = np.arange(50, 300)


ltt7379_twodspec = spectral_reduction.TwoDSpec(
    standard_frame,
    readnoise=2.9,
    gain=1.43,
    log_file_name=None,
    spatial_mask=spatial_mask,
)


ltt7379_twodspec.ap_trace(
    nspec=1,
    display=False,
)

# Tophat
ltt7379_twodspec.ap_extract(
    spec_id=0, apwidth=10, skywidth=5, skysep=0, skydeg=1, optimal=False
)
count_tophat = copy.deepcopy(ltt7379_twodspec.spectrum_list[0].count)
count_tophat_err = copy.deepcopy(ltt7379_twodspec.spectrum_list[0].count_err)
img_residual_tophat = copy.deepcopy(ltt7379_twodspec.img_residual)

# Horne86
ltt7379_twodspec.ap_extract(
    spec_id=0,
    apwidth=10,
    skywidth=5,
    skysep=0,
    skydeg=1,
    optimal=True,
    algorithm="horne86",
)
count_horne86 = copy.deepcopy(ltt7379_twodspec.spectrum_list[0].count)
count_horne86_err = copy.deepcopy(ltt7379_twodspec.spectrum_list[0].count_err)
img_residual_horne86 = copy.deepcopy(ltt7379_twodspec.img_residual)

# Marsh89
ltt7379_twodspec.ap_extract(
    spec_id=0,
    apwidth=10,
    skywidth=5,
    skysep=0,
    skydeg=1,
    optimal=True,
    algorithm="marsh89",
    pord=4,
    nreject=0,
    qmode="fast-linear",
)
count_marsh89 = copy.deepcopy(ltt7379_twodspec.spectrum_list[0].count)
count_marsh89_err = copy.deepcopy(ltt7379_twodspec.spectrum_list[0].count_err)
img_residual_marsh89 = copy.deepcopy(ltt7379_twodspec.img_residual)

fig = plt.figure(1, figsize=(6, 6))
fig.clf()

gs = GridSpec(5, 1, figure=fig)

ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1:3])
ax3 = fig.add_subplot(gs[3:])

data_logged = np.log10(ltt7379_twodspec.img)
data_logged[np.isnan(data_logged)] = np.nanmin(data_logged)

ax1.imshow(data_logged, origin="lower", aspect="auto", vmax=5.0)
ax1.plot(
    np.array(ltt7379_twodspec.spectrum_list[0].trace) - 0.5,
    lw=1,
    ls=":",
    color="black",
)
ax1.set_xticks([])
ax1.set_ylim(35, 115)
ax1.set_ylabel("Pixel (Spatial)")

ax2.plot(count_tophat, color="red", label="Tophat", lw=1)
ax2.plot(count_horne86, color="green", label="Horne86", lw=1)
ax2.plot(count_marsh89, color="blue", label="Marsh89", lw=1)


ax2.set_xticks([])
ax2.set_ylabel(r"e$^{-}$ count")

ax2.legend(framealpha=0.9)

ax3.plot(count_tophat_err, color="red", label="Tophat", lw=1)
ax3.plot(count_horne86_err, color="green", label="Horne86", lw=1)
ax3.plot(count_marsh89_err, color="blue", label="Marsh89", lw=1)

ax3.set_xlabel("Pixel (Dispersion)")
ax3.set_ylabel(r"e$^-$ count uncertainty")

ax1.set_xlim(0, 2048)
ax2.set_xlim(0, 2048)
ax3.set_xlim(0, 2048)

ax2.set_ylim(10000, 2e5)
#ax2.set_yscale('log')

ax3.set_ylim(0, 3e3)
# ax3.set_yscale('log')

fig.tight_layout()
fig.subplots_adjust(hspace=0, left=0.15)
fig.savefig("fig_04_extraction_compared.jpg")
fig.savefig("fig_04_extraction_compared.pdf")


fig2 = plt.figure(2, figsize=(6, 6))
fig2.clf()

ax1 = fig2.add_subplot(3, 1, 1)
ax2 = fig2.add_subplot(3, 1, 2)
ax3 = fig2.add_subplot(3, 1, 3)

ax1.imshow(
    np.log10(img_residual_tophat), origin="lower", aspect="auto", vmax=5.0
)
ax2.imshow(
    np.log10(img_residual_horne86), origin="lower", aspect="auto", vmax=5.0
)
ax3.imshow(
    np.log10(img_residual_marsh89), origin="lower", aspect="auto", vmax=5.0
)

ax1.set_ylim(35, 115)
ax2.set_ylim(35, 115)
ax3.set_ylim(35, 115)

ax1.set_xticks([])
ax2.set_xticks([])

fig.tight_layout()
fig.subplots_adjust(hspace=0)

fig.savefig("fig_04b_extraction_residual_compared.jpg")
fig.savefig("fig_04b_extraction_residual_compared.pdf")
