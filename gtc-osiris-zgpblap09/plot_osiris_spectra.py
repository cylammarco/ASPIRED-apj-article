from matplotlib import pyplot as plt
from spectres import spectres
import numpy as np

plt.ion()

r1000b_1 = np.load('r1000b_1.npy')
r1000b_1_wave = r1000b_1[:, 0]
r1000b_1_flux = r1000b_1[:, 1]

r1000b_2 = np.load('r1000b_2.npy')
r1000b_2_wave = r1000b_2[:, 0]
r1000b_2_flux = r1000b_2[:, 1]

r2500u_1 = np.load('r2500u_1.npy')
r2500u_1_wave = r2500u_1[:, 0]
r2500u_1_flux = r2500u_1[:, 1]

r2500u_2 = np.load('r2500u_2.npy')
r2500u_2_wave = r2500u_2[:, 0]
r2500u_2_flux = r2500u_2[:, 1]

r1000b_2_resampled = spectres(r1000b_1_wave, r1000b_2_wave, r1000b_2_flux)
r2500u_2_resampled = spectres(r2500u_1_wave, r2500u_2_wave, r2500u_2_flux)


lines_H = [6562.79, 4861.35, 4340.472, 4101.734, 3970.075, 3889.064, 3835.397]
lines_HeI = [10830, 7065, 6678, 5876, 5016, 4922, 4713, 4472, 4026, 3965, 3889]
lines_HeI_highres = [
    3819.607, 3819.76, 3888.6489, 3964.729, 4009.27, 4026.191, 4120.82,
    4143.76, 4387.929, 4437.55, 4471.479
]
lines_HeII = [5411.52, 4685.7, 4542.8, 4340, 4200, 3203.10]
lines_CII = [
    3918.978, 3920.693, 4267.003, 4267.258, 5145.16, 5151.09, 5889.77, 6578.05,
    6582.88, 7231.32, 7236.42
]
lines_NII = [3995.00, 4630.54, 5005.15, 5679.56, 6482.05, 6610.56]
lines_OI = [9263, 8446, 7775, 7774, 7772, 6158]
lines_OI_forbidden = [6363, 6300, 5577]
lines_OII = [6721, 6641, 4649, 4416, 4317, 4076, 3982, 3973, 3911, 3749, 3713]
lines_OII_forbidden = [3729, 3726]
lines_OIII_forbidden = [5007, 4959, 4363]
lines_NaI = [8195, 8183, 5896, 5890]
lines_MgI = [
    8807, 5528, 5184, 5173, 5167, 4703, 4571, 3838, 3832, 3829, 2852, 2780
]
lines_MgII = [9632, 9244, 9218, 8235, 8214, 7896, 7877, 4481, 2803, 2796, 2791]
lines_SiII = [7850, 6371, 6347, 5979, 5958, 5670, 5056, 5041, 4131, 4128, 3856]
lines_SII = [6715, 6397, 6313, 6305, 6287, 5647, 5640, 5606, 5454, 5433, 4163]
lines_CaII = [8662, 8542, 8498, 3969, 3934, 3737, 3706, 3180, 3159]
lines_CaII_forbidden = [7324, 7292]
lines_FeII = [5363, 5235, 5198, 5169, 5018, 4924, 4549, 4515, 4352, 4303]
lines_FeIII = [5158, 5129, 4432, 4421, 4397]

lines_OI_atm = [5577, 6300.3, 7774.2]
lines_OH_atm = [6863, 7340, 7523, 7662, 7750, 7821, 7913, 8025]
lines_O2_atm = [6277.7, 7605, 6869]
lines_H2O_atm = [
    6940.7, 7185.9, 7234.1, 7253.2, 7276.2, 7292.0, 7303.2, 7318.4
]

ymax = np.nanmax((r2500u_1_flux, r2500u_2_flux))

model = np.genfromtxt('Husfeld/t35g50e25.dat.txt')
w = model[:, 0]
f = model[:, 1] / 3.5e23

plt.rcParams.update({'font.size': 16})
plt.rcParams.update({'legend.fontsize': 12})

fig, (ax1, ax2) = plt.subplots(2, figsize=(15, 9))

ax1.plot(r2500u_1_wave,
         np.nanmean((r2500u_1_flux, r2500u_2_resampled), axis=0),
         color='0.0',
         label='R2500U average')
ax1.plot(r2500u_1_wave, r2500u_1_flux, lw=0.5, label='Epoch 3')
ax1.plot(r2500u_2_wave, r2500u_2_flux, lw=0.5, label='Epoch 4')
ax1.vlines(lines_H, 0, ymax, color='C0', label='H')
ax1.vlines(lines_HeI_highres, 0, ymax, color='C1', label='He I')
ax1.vlines(lines_HeII, 0, ymax, color='C2', label='He II')
ax1.vlines(lines_CII, 0, ymax, color='C3', label='C II')
ax1.vlines(lines_NII, 0, ymax, color='C4', label='N II')
ax1.vlines(lines_OII, 0, ymax, color='C5', label='O II')
ax1.vlines(lines_MgII, 0, ymax, color='C6', label='Mg II')
# C7 is grey, reserved for atmosphere...
ax1.vlines(lines_SiII, 0, ymax, color='C8', label='Si II')
ax1.vlines(lines_CaII, 0, ymax, color='C9', label='Ca II')

ax1.plot(w,
         f,
         lw=3,
         color='red',
         alpha=0.5,
         label='Husfeld 35000K log(g)=5.0 Y=0.25')
ax1.plot(w, f * 1.5, lw=3, color='green', alpha=0.5)

ax1.set_xlim(3800, 4600)
ax1.set_ylim(1.5e-15, 7.0e-15)
ax1.set_xlabel('Wavelength / A')
ax1.set_ylabel('Flux / ( erg / cm / cm / s / A)')
ax1.grid()
ax1.legend(loc='upper right', ncol=2)

ax2.plot(r1000b_1_wave,
         np.nanmean((r1000b_1_flux, r1000b_2_resampled), axis=0) * 1.35,
         color='0.2',
         label='R1000B average')
ax2.plot(r1000b_1_wave, r1000b_1_flux * 1.35, lw=0.5, label='Epoch 1')
ax2.plot(r1000b_2_wave, r1000b_2_flux * 1.35, lw=0.5, label='Epoch 2')
ax2.vlines(lines_H, 0, ymax, color='C0', label='H')
ax2.vlines(lines_HeI, 0, ymax, color='C1', label='He I')
ax2.vlines(lines_HeI_highres, 0, ymax, color='C1')
ax2.vlines(lines_HeII, 0, ymax, color='C2', label='He II')
ax2.vlines(lines_CII, 0, ymax, color='C3', label='C II')
ax2.vlines(lines_NII, 0, ymax, color='C4', label='N II')
ax2.vlines(lines_OII, 0, ymax, color='C5', label='O II')
ax2.vlines(lines_MgII, 0, ymax, color='C6', label='Mg II')
# C7 is grey, reserved for atmosphere...
ax2.vlines(lines_SiII, 0, ymax, color='C8', label='Si II')

ax2.plot(w,
         f,
         lw=3,
         color='red',
         alpha=0.5,
         label='Husfeld 35000K log(g)=5.0 Y=0.25')
ax2.plot(w, f * 1.5, lw=3, color='green', alpha=0.5)

ax2.vlines(lines_OI_atm,
           0,
           ymax,
           color='grey',
           alpha=0.25,
           lw=8,
           label='Atmosphere')
ax2.vlines(lines_O2_atm, 0, ymax, color='grey', alpha=0.25, lw=8)
ax2.vlines(lines_OH_atm, 0, ymax, color='grey', alpha=0.25, lw=8)
ax2.vlines(lines_H2O_atm, 0, ymax, color='grey', alpha=0.25, lw=8)

ax2.set_xlim(4600, 7800)
ax2.set_ylim(0.0, 3.5e-15)
ax2.set_xlabel('Wavelength / A')
ax2.set_ylabel('Flux / ( erg / cm / cm / s / A)')
ax2.grid()
ax2.legend(ncol=2)

plt.subplots_adjust(left=0.05, bottom=0.06, right=0.98, top=0.95, hspace=0.15)
plt.suptitle('ZGP-BLAP-09 (GTC/OSIRIS)')

plt.savefig('ZTF_BLAP_09_GTC.png')
plt.savefig('ZTF_BLAP_09_GTC.pdf')
