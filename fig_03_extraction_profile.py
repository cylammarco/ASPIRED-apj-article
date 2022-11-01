from astropy.io import fits
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon, FancyArrowPatch, Rectangle
from aspired import spectral_reduction
import numpy as np
from scipy.ndimage import rotate
from statsmodels.nonparametric.smoothers_lowess import lowess

fits_file = fits.open('ogg2m001-en06-20160111-0005-e00.fits.fz')[1]
data = fits_file.data
header = fits_file.header

red_spatial_mask = np.arange(0, 330)
red_spec_mask = np.arange(15, 1500)

twodspec = spectral_reduction.TwoDSpec(data,
                                       header=header,
                                       spatial_mask=red_spatial_mask,
                                       spec_mask=red_spec_mask,
                                       cosmicray=True,
                                       sigclip=3.,
                                       readnoise=3.5,
                                       gain=2.3,
                                       log_level='CRITICAL',
                                       log_file_name=None)

twodspec.ap_trace(nspec=1,
                  nwindow=30,
                  ap_faint=5,
                  trace_width=25,
                  resample_factor=5,
                  shift_tol=120,
                  fit_deg=7,
                  display=False)

trace = np.array(twodspec.spectrum_list[0].trace)

# trace is at (spatial) pix 70 for (dispersion) pix 600
extraction_slice = Polygon(
    [i for i in zip(np.arange(len(trace)), trace - 30)] +
    [i for i in zip(np.arange(len(trace))[::-1], trace[::-1] + 30)],
    hatch='\\/...',
    edgecolor='C1',
    alpha=0.3,
    facecolor='none',
    label="Region for extraction")
# arrows
sky_arrow_1 = FancyArrowPatch((750, 72), (750, 93),
                              color='C2',
                              arrowstyle='<|-|>',
                              mutation_scale=8)
sky_arrow_2 = FancyArrowPatch((750, 122), (750, 143),
                              color='C2',
                              arrowstyle='<|-|>',
                              mutation_scale=8)
source_arrow = FancyArrowPatch((750, 92), (750, 123),
                               color='C3',
                               arrowstyle='<|-|>',
                               mutation_scale=8)

# boxes
box_1 = Rectangle((230, 52), 20, 70, edgecolor='black', facecolor='none', lw=1)
box_2 = Rectangle((1150, 128),
                  20,
                  70,
                  edgecolor='black',
                  facecolor='none',
                  lw=1)

slice_1 = twodspec.img[50:120, 240]
slice_2 = twodspec.img[125:195, 1160]

lowess_fit_1 = lowess(slice_1[25:51], np.arange(74, 100), frac=0.1)[:, 1]
lowess_fit_2 = lowess(slice_2[25:51], np.arange(149, 175), frac=0.1)[:, 1]

fig = plt.figure(1, figsize=(6, 6))
fig.clf()

ax1 = fig.add_subplot(2, 1, 1)
ax2 = fig.add_subplot(2, 2, 3)
ax3 = fig.add_subplot(2, 2, 4)

ax1.imshow(np.log10(twodspec.img),
           origin='lower',
           aspect='auto',
           vmin=np.nanpercentile(np.log10(twodspec.img), 10),
           vmax=np.nanpercentile(np.log10(twodspec.img), 97),
           cmap='gray_r')

# Extraction slice
ax1.add_patch(extraction_slice)
# Sky arrows
ax1.add_patch(sky_arrow_1)
ax1.add_patch(sky_arrow_2)
# Source arrow
ax1.add_patch(source_arrow)
# Boxes
ax1.add_patch(box_1)
ax1.add_patch(box_2)

# Sky
ax1.plot(trace + 30, c='C2', label="Region for sky estimation")
ax1.plot(trace + 20, c='C2')
ax1.plot(trace - 20, c='C2')
ax1.plot(trace - 30, c='C2')

# Source
ax1.plot(trace + 10, c='C3', label="Region for the source count")
ax1.plot(trace - 10, c='C3')
ax1.legend()
ax1.set_xlabel('Pixel (Dispersion Direction)')
ax1.set_ylabel('Pixel (Spatial Direction)')
ax1.xaxis.tick_top()
ax1.xaxis.set_label_position('top')
ax1.set_ylim(20, 260)
ax1.set_xlim(0, 1450)

ax1.text(300, 40, 'Slice 1')
ax1.text(1200, 120, 'Slice 2')

vmin = min(slice_1) - 10
vmax = max(slice_1) * 1.15

ax2.plot(np.arange(49, 119), slice_1/max(slice_1), color='black')
ax2.plot(np.arange(74, 100), lowess_fit_1/max(lowess_fit_1) * 1.15, color='C0')
ax2.vlines(trace[240] - 1, 0.0, 1.2, ls=":", color='black')
ax2.vlines(trace[240] + 10 - 1, 0.0, 1.2, ls=":", color='C3')
ax2.vlines(trace[240] - 10 - 1, 0.0, 1.2, ls=":", color='C3')
ax2.vlines(trace[240] + 20 - 1, 0.0, 1.2, ls=":", color='C2')
ax2.vlines(trace[240] + 30 - 1, 0.0, 1.2, ls=":", color='C2')
ax2.vlines(trace[240] - 20 - 1, 0.0, 1.2, ls=":", color='C2')
ax2.vlines(trace[240] - 30 - 1, 0.0, 1.2, ls=":", color='C2')

ax3.plot(np.arange(125, 195), slice_2/max(slice_2), color='black', label='Centroid (Trace)')
ax3.plot(np.arange(150, 176),
         lowess_fit_2/max(lowess_fit_2) * 1.15,
         color='C0',
         label='LOWESS profile')
ax3.vlines(trace[1160], 0, vmax, ls=":", color='black')
ax3.vlines(trace[1160] + 10,
           0.0,
           1.2,
           ls=":",
           color='C3',
           label='Source region')
ax3.vlines(trace[1160] - 10, 0, 1.2, ls=":", color='C3')
ax3.vlines(trace[1160] + 20,
           0.0,
           1.2,
           ls=":",
           color='C2',
           label='Sky region')
ax3.vlines(trace[1160] + 30, 0, 1.2, ls=":", color='C2')
ax3.vlines(trace[1160] - 20, 0, 1.2, ls=":", color='C2')
ax3.vlines(trace[1160] - 30, 0, 1.2, ls=":", color='C2')

ax2.set_xlim(45, 125)
ax3.set_xlim(128, 198)

ax2.set_ylim(0.45, 1.2)
ax3.set_ylim(0.45, 1.2)
ax3.set_yticks([])

ax2.text(50, 1.1, 'Slice 1')
ax3.text(135, 1.1, 'Slice 2')

ax2.set_ylabel('Normalised Electron count')
ax2.set_xlabel('Pixel (Spatial Direction)')
ax3.set_xlabel('Pixel (Spatial Direction)')

ax3.legend(loc='lower right')

plt.tight_layout()
plt.subplots_adjust(wspace=0, hspace=0.05)
plt.savefig('fig_03_extraction_profile.jpg')
plt.savefig('fig_03_extraction_profile.pdf')
