import numpy as np

from preamble import *
import math
from scipy import integrate
from check_beta import *
import matplotlib.colors as mcolors
from matplotlib.font_manager import FontProperties
import matplotlib.lines as mlines

font = FontProperties()
font.set_weight('bold')

# 2 - Constants & functions

c = 2.99792458e+10  # [cm / c]
h = 6.582119569e-22  # [Mev * с]
k = 8.617333262e-11  # [Mev / K]
T_nu = 1.945e+0  # [K]
T_gamma = 2.7277e+0  # [K]
GG = 6.67e-8  # [cm^3 * g^-1 * s^-2]
kt = h * c ** 3 / (8 * np.pi * GG * 10 ** 8)
print(kt / k)
print(T_gamma * 10 ** 11)


# ### TEST FUNCTIONS ###
# def n_pbh_f(f_pbh, m_pbh, sigma):
#     if sigma == 0:
#         return (f_pbh / 4.52386e+44) * (1.e+15 / m_pbh)
#     else:
#         return 2.24e-45 * f_pbh * (1.e+15 / m_pbh) * np.exp(- sigma ** 2 / 2) / \
#             (0.5 * (1 + math.erf((sigma ** 2 - np.log(5.e+14 / m_pbh)) / (math.sqrt(2) * sigma))))
# def n_to_f(n, m, sigma):
#     return 4.52386e+44 * (m / 1.e+15) * n * np.exp(sigma ** 2 / 2) * \
#         (0.5 * (1 + math.erf((sigma ** 2 - np.log10(8.e+14 / m)) / (math.sqrt(2) * sigma))))


def n_pbh_beta(beta, m_pbh, sigma):
    if sigma == 0:
        return 1.19e-27 * beta * (1.e+15 / m_pbh) ** (3 / 2)
    else:
        return 1.19e-27 * beta * (1.e+15 / m_pbh) ** (3 / 2) * np.exp(-9 * sigma ** 2 / 8)


def data_extract_nu(path_pri: str, path_e: str, path_mu: str, path_tau: str):
    data_nu_pri = (np.genfromtxt(data_cosmology_new + path_pri, skip_header=1))
    data_sec_e = (np.genfromtxt(data_cosmology_new + path_e, skip_header=1))
    data_sec_mu = (np.genfromtxt(data_cosmology_new + path_mu, skip_header=1))
    data_sec_tau = (np.genfromtxt(data_cosmology_new + path_tau, skip_header=1))
    data_final = (np.array(data_nu_pri[:, 1]) + np.array(data_sec_e[:, 1]) + np.array(data_sec_mu[:, 1]) +
                  np.array(data_sec_tau[:, 1]))
    return np.array(data_nu_pri[:, 0]), data_final


def data_extract_nu_new(path_pri: str, path_e: str, path_mu: str, path_tau: str):
    cosmodata = r"D:\Downloads\BlackHawk\scripts\cosmology_scripts\CosmologyData"
    data_nu_pri = (np.genfromtxt(cosmodata + path_pri, skip_header=1))
    data_sec_e = (np.genfromtxt(cosmodata + path_e, skip_header=1))
    data_sec_mu = (np.genfromtxt(cosmodata + path_mu, skip_header=1))
    data_sec_tau = (np.genfromtxt(cosmodata + path_tau, skip_header=1))
    data_final = (np.array(data_nu_pri[:, 1]) + np.array(data_sec_e[:, 1]) + np.array(data_sec_mu[:, 1]) +
                  np.array(data_sec_tau[:, 1]))
    return np.array(data_nu_pri[:, 0]), data_final


def data_extract_ph(path_pri: str, path_sec: str):
    cosmodata = r"D:\Downloads\BlackHawk\scripts\cosmology_scripts\CosmologyData"
    data_ph_pri = (np.genfromtxt(cosmodata + path_pri, skip_header=1))
    data_ph_sec = (np.genfromtxt(cosmodata + path_sec, skip_header=1))
    data_final = np.array(data_ph_pri[:, 1]) + np.array(data_ph_sec[:, 1])
    return np.array(data_ph_pri[:, 0]), data_final


def relic_nu(points):
    return (24 * np.pi * points ** 2) / ((c ** 2 * (2 * np.pi * h) ** 3) * (np.exp(points / (k * T_nu)) + 1))


#################################
# Основной результат - нейтрино #
#################################

f, ax = plt.subplots(1, 1, figsize=(fig_width * 0.9, fig_width * 0.75))
f.subplots_adjust(left=0.1, right=0.95)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$E\,\,[$МэВ$]$')
ax.set_ylabel(r'$\mathrm{d}\Phi/\mathrm{d}E\,\, [$МэВ$^{-1}\cdot$см$^{-2}\cdot$с$^{-1}$$]$')
ax.set_ylim(1.e-22, 1.e+24)
ax.set_xlim(1.e-13, 1.e+11)
ax.set_xticks(np.logspace(-10, 10, num=5), minor=False)
ax.set_xticks(np.logspace(-13, 11, num=25), labels=[], minor=True)
ax.set_yticks(np.logspace(-20, 20, num=5), minor=False)
ax.set_yticks(np.logspace(-22, 24, num=24), labels=[], minor=True)

# ##### Nu-background ##### #

# Relic Neutrinos
x = np.logspace(-12, -8, 50)
ax.plot(x, relic_nu(x), linewidth=4, color='cyan')
ax.plot(x, relic_nu(x), linewidth=2, color='darkblue')

ax.axvline(x=0.1, color='black', linewidth=1, linestyle='dashed')

# Thermal Sun Neutrinos
data_TNF = np.genfromtxt(data_here + r'\Sun-thermal.dat.txt', skip_header=4)
ax.plot(np.array(data_TNF[:, 0]) * 1.e-6, np.array(data_TNF[:, 1]) * 1.e+6, linestyle='dashdot', linewidth=2,
        color='orange')

# Solar Neutrinos
# pp
x = np.logspace(-6.5, -0.40, 50)
ax.plot(x, 832.7e+10 * x ** 2 * (1 - 2.4 * x), linewidth=4, color='yellow')
ax.plot((10 ** (-0.40), 10 ** (-0.40)), (10 ** 10.6, 10 ** 4.2), linewidth=4, color='yellow')
# B8
data_SNF_B8 = np.genfromtxt(data_here + r'\Sun-nuclear-B8.dat.txt', skip_header=4)
ax.plot(np.array(data_SNF_B8[20:, 0]) * 1.e-6, np.array(data_SNF_B8[20:, 1]) * 1.e+6, linewidth=4, color='yellow')

ax.plot(x, 832.7e+10 * x ** 2 * (1 - 2.4 * x), linewidth=2, color='indianred')
ax.plot((10 ** (-0.40), 10 ** (-0.40)), (10 ** 10.6, 10 ** 4.2), linewidth=2, color='indianred')
ax.plot(np.array(data_SNF_B8[20:, 0]) * 1.e-6, np.array(data_SNF_B8[20:, 1]) * 1.e+6, linewidth=2, color='indianred')

# Atmospheric Neutrinos
data_ANF = np.genfromtxt(data_here + r'\Atmospheric.dat.txt', skip_header=4)
ax.plot(np.array(data_ANF[:, 0]) * 1.e-6, (np.array(data_ANF[:, 1]) + np.array(data_ANF[:, 2])) * 1.e+6,
        linewidth=4, color='lightcyan')
ax.plot(np.array(data_ANF[:, 0]) * 1.e-6, (np.array(data_ANF[:, 1]) + np.array(data_ANF[:, 2])) * 1.e+6,
        linewidth=2, color='cornflowerblue')

# Data from old supernova
osn_x, osn_y = np.loadtxt(data_here + r'\dataoldsn.csv', delimiter=' ', unpack=True)
ax.plot(osn_x, osn_y, linewidth=2, color='sienna')

# # Antineutrinos from BBN (Ivanchick)
# Antinu_n
data_antinu_n = np.genfromtxt(data_here + r'\antinu_n.txt')  # , skip_header=15)
ax.plot(np.array(data_antinu_n[:, 0]), np.array(data_antinu_n[:, 1]), linewidth=2, color='red', linestyle='solid')
# Antinu_tau
data_antinu_t = np.genfromtxt(data_here + r'\antinu_t.txt')
ax.plot(np.array(data_antinu_t[:, 0]), np.array(data_antinu_t[:, 1]), linewidth=2, color='magenta')

# #### For main results nu_1 #### #

# 10^8 (log) sigma = 2
energies_8, nu_all_8 = data_extract_nu_new(r"\ra_1.0e+08_2.0_pri.txt",
                                           r"\ra_1.0e+08_2.0_e_sec.txt",
                                           r"\ra_1.0e+08_2.0_mu_sec.txt",
                                           r"\ra_1.0e+08_2.0_tau_sec.txt")
ax.plot(energies_8 * 1.e+3, n_pbh_beta(beta_shtrix(1.e+8), 1.e+8, 2.0) * (c / (4 * np.pi)) * nu_all_8 / 1.e+3,
        label=r"$\sigma = 2.0,\,\, M_{PBH}=2\cdot10^{8}$ г", linewidth=3, linestyle='dashed', color='lime')

# 2 * 10^9 (mono) (Pythia)
energies, nu_all = data_extract_nu_new(r"\ra_2.0e+09_0.0_pri.txt",
                                       r"\ra_2.0e+09_0.0_e_sec.txt",
                                       r"\ra_2.0e+09_0.0_mu_sec.txt",
                                       r"\ra_2.0e+09_0.0_tau_sec.txt")
ax.plot(energies * 1.e+3, n_pbh_beta(beta_shtrix(1.e+9), 1.e+9, 0) * (c / (4 * np.pi)) * nu_all / 1.e+3,
        label=r"$\sigma = 0.01,\,\, M_{PBH}=2\cdot10^{9}$ г", linewidth=3, linestyle='dashed', color='indigo')

# 10^11 (mono) NEW
energies_11_new, nu_all_11_new = data_extract_nu_new(r"\ra_1.0e+11_0.0_pri.txt",
                                                     r"\ra_1.0e+11_0.0_e_sec.txt",
                                                     r"\ra_1.0e+11_0.0_mu_sec.txt",
                                                     r"\ra_1.0e+11_0.0_tau_sec.txt")
ax.plot(energies_11_new * 1.e+3,
        n_pbh_beta(beta_shtrix(1.e+11), 1.e+11, 0.0) * (c / (4 * np.pi)) * nu_all_11_new / 1.e+3,
        label=r"$\sigma = 0.01,\,\, M_{PBH}=10^{11}$ г", linewidth=5, linestyle='dashed', color='olive')

# #### For second results nu_2 #### #

# # M = 10^11 (sigma = 1.0)
energies_11_1, nu_all_11_1 = data_extract_nu_new(r"\ra_1.0e+11_1.0_pri.txt",
                                                 r"\ra_1.0e+11_1.0_e_sec.txt",
                                                 r"\ra_1.0e+11_1.0_mu_sec.txt",
                                                 r"\ra_1.0e+11_1.0_tau_sec.txt")
# ax.plot(energies_11_1 * 1.e+3,
#         n_pbh_beta(beta_shtrix(1.e+11), 1.e+11, 1.0) * (c / (4 * np.pi)) * nu_all_11_1 / 1.e+3,
#         label=r"$\sigma = 1.0,\,\, M_{PBH}=10^{11}$ г", linewidth=3, linestyle='dashed', color='indianred')
#
# # M = 10^13 (sigma = 2.0)
energies_13_2, nu_all_13_2 = data_extract_nu_new(r"\ra_1.0e+13_2.0_pri.txt",
                                                 r"\ra_1.0e+13_2.0_e_sec.txt",
                                                 r"\ra_1.0e+13_2.0_mu_sec.txt",
                                                 r"\ra_1.0e+13_2.0_tau_sec.txt")
# ax.plot(energies_13_2 * 1.e+3,
#         n_pbh_beta(beta_shtrix(1.e+13), 1.e+13, 2.0) * (c / (4 * np.pi)) * nu_all_13_2 / 1.e+3,
#         label=r"$\sigma = 2.0,\,\, M_{PBH}=10^{13}$ г", linewidth=3, linestyle='dashed', color='lawngreen')
#
# # M = 10^15 (sigma = 3.0)
energies_15_3, nu_all_15_3 = data_extract_nu_new(r"\ra_1.0e+15_3.0_pri.txt",
                                                 r"\ra_1.0e+15_3.0_e_sec.txt",
                                                 r"\ra_1.0e+15_3.0_mu_sec.txt",
                                                 r"\ra_1.0e+15_3.0_tau_sec.txt")
# ax.plot(energies_15_3 * 1.e+3,
#         n_pbh_beta(1.e-26, 1.e+15, 3.0) * (c / (4 * np.pi)) * nu_all_15_3 / 1.e+3,
#         label=r"$\sigma = 3.0,\,\, M_{PBH}=10^{15}$ г", linewidth=3, linestyle='dashed', color='royalblue')

# ## Подписи и заштрихованная область  ## #

yn11 = n_pbh_beta(beta_shtrix(1.e+11), 1.e+11, 0.0) * (c / (4 * np.pi)) * nu_all_11_new / 1.e+3
yn9 = n_pbh_beta(beta_shtrix(1.e+9), 1.e+9, 0.0) * (c / (4 * np.pi)) * nu_all / 1.e+3
xn = energies * 1.e+3
ax.fill_between(xn, yn11, yn9, where=(yn9 > yn11), color='grey', alpha=0.6)
xr = np.logspace(-12, -8, 50)
ax.fill_between(xr, relic_nu(xr), 1e-12, color='white', alpha=1)
xs = np.logspace(-6.5, -0.40, 50)
ys = 832.7e+10 * xs ** 2 * (1 - 2.4 * xs)
ax.fill_between(xs, ys, 1, color='white', alpha=1)
ax.fill_between(np.array(data_TNF[:, 0]) * 1.e-6, np.array(data_TNF[:, 1]) * 1.e+6, 1, color='white', alpha=1)
xs_new = np.logspace(-4, 6, 100)
ax.fill_between(xs_new, 1, 1.e-25, color='white', alpha=1)

ax.text(3.e-9, 1.e+19, r'C$\nu$B', fontsize=18, color='cyan')
ax.text(3.e-9, 1.e+19, r'C$\nu$B', fontsize=18, color='darkblue')

ax.text(1.e-11, 5.e+11, r'$\widetilde{\nu}_{BBN}(n)$', fontsize=18, color='red')
ax.text(2.e-12, 3.e+3, r'$\widetilde{\nu}_{BBN}(t)$', fontsize=18, color='magenta')

ax.text(1.e-3, 4.e+9, r'$\nu_{\mathrm{sun}}$', fontsize=22, color='yellow', rotation=37)
ax.text(1.e-3, 4.e+9, r'$\nu_{\mathrm{sun}}$', fontsize=22, color='indianred', rotation=37)

ax.text(1.e+4, 3.e-8, r'$\nu_{\mathrm{atm}}$', fontsize=22, color='lightcyan', rotation=-50)
ax.text(1.e+4, 3.e-8, r'$\nu_{\mathrm{atm}}$', fontsize=22, color='cornflowerblue', rotation=-50)

ax.text(4.e+1, 2.e-1, r'$\nu_{\mathrm{SN}}$', fontsize=20, color='saddlebrown')

ax.legend(loc='upper right')

fig_name = figure_folder + r"\nu_results_1.pdf"
plt.savefig(fig_name)

f.clf()

########
# GUNS #
########

f, ax = plt.subplots(1, 1, figsize=(fig_width * 0.9, fig_width * 0.75))
f.subplots_adjust(left=0.1, right=0.95)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$E\,\,[$МэВ$]$')
ax.set_ylabel(r'$\mathrm{d}\Phi/\mathrm{d}E\,\, [$МэВ$^{-1}\cdot$см$^{-2}\cdot$с$^{-1}$$]$')
ax.set_ylim(1.e-22, 1.e+24)
ax.set_xlim(1.e-13, 1.e+11)
ax.set_xticks(np.logspace(-10, 10, num=5), minor=False)
ax.set_xticks(np.logspace(-13, 11, num=25), labels=[], minor=True)
ax.set_yticks(np.logspace(-20, 20, num=5), minor=False)
ax.set_yticks(np.logspace(-22, 24, num=24), labels=[], minor=True)

# ##### Nu-background ##### #

# Relic Neutrinos
x = np.logspace(-12, -8, 50)
ax.plot(x, relic_nu(x), linewidth=4, color='cyan')
ax.plot(x, relic_nu(x), linewidth=2, color='darkblue')

ax.axvline(x=0.1, color='black', linewidth=1, linestyle='dashed')

# Thermal Sun Neutrinos
data_TNF = np.genfromtxt(data_here + r'\Sun-thermal.dat.txt', skip_header=4)
ax.plot(np.array(data_TNF[:, 0]) * 1.e-6, np.array(data_TNF[:, 1]) * 1.e+6, linestyle='dashdot', linewidth=2,
        color='orange')

# Solar Neutrinos
# pp
x = np.logspace(-6.5, -0.40, 50)
ax.plot(x, 832.7e+10 * x ** 2 * (1 - 2.4 * x), linewidth=4, color='yellow')
ax.plot((10 ** (-0.40), 10 ** (-0.40)), (10 ** 10.6, 10 ** 4.2), linewidth=4, color='yellow')
# B8
data_SNF_B8 = np.genfromtxt(data_here + r'\Sun-nuclear-B8.dat.txt', skip_header=4)
ax.plot(np.array(data_SNF_B8[20:, 0]) * 1.e-6, np.array(data_SNF_B8[20:, 1]) * 1.e+6, linewidth=4, color='yellow')

ax.plot(x, 832.7e+10 * x ** 2 * (1 - 2.4 * x), linewidth=2, color='indianred')
ax.plot((10 ** (-0.40), 10 ** (-0.40)), (10 ** 10.6, 10 ** 4.2), linewidth=2, color='indianred')
ax.plot(np.array(data_SNF_B8[20:, 0]) * 1.e-6, np.array(data_SNF_B8[20:, 1]) * 1.e+6, linewidth=2, color='indianred')

# Atmospheric Neutrinos
data_ANF = np.genfromtxt(data_here + r'\Atmospheric.dat.txt', skip_header=4)
ax.plot(np.array(data_ANF[:, 0]) * 1.e-6, (np.array(data_ANF[:, 1]) + np.array(data_ANF[:, 2])) * 1.e+6,
        linewidth=4, color='lightcyan')
ax.plot(np.array(data_ANF[:, 0]) * 1.e-6, (np.array(data_ANF[:, 1]) + np.array(data_ANF[:, 2])) * 1.e+6,
        linewidth=2, color='cornflowerblue')

# Data from old supernova
osn_x, osn_y = np.loadtxt(data_here + r'\dataoldsn.csv', delimiter=' ', unpack=True)
ax.plot(osn_x, osn_y, linewidth=2, color='sienna')

# # Antineutrinos from BBN (Ivanchick)
# Antinu_n
data_antinu_n = np.genfromtxt(data_here + r'\antinu_n.txt')  # , skip_header=15)
ax.plot(np.array(data_antinu_n[:, 0]), np.array(data_antinu_n[:, 1]), linewidth=2, color='red', linestyle='solid')
# Antinu_tau
data_antinu_t = np.genfromtxt(data_here + r'\antinu_t.txt')
ax.plot(np.array(data_antinu_t[:, 0]), np.array(data_antinu_t[:, 1]), linewidth=2, color='magenta')

# ## Подписи и заштрихованная область  ## #
#
# yn11 = n_pbh_beta(beta_shtrix(1.e+11), 1.e+11, 0.0) * (c / (4 * np.pi)) * nu_all_11_new / 1.e+3
# yn9 = n_pbh_beta(beta_shtrix(1.e+9), 1.e+9, 0.0) * (c / (4 * np.pi)) * nu_all / 1.e+3
# xn = energies * 1.e+3
# ax.fill_between(xn, yn11, yn9, where=(yn9 > yn11), color='grey', alpha=0.6)
# xr = np.logspace(-12, -8, 50)
# ax.fill_between(xr, relic_nu(xr), 1e-12, color='white', alpha=1)
# xs = np.logspace(-6.5, -0.40, 50)
# ys = 832.7e+10 * xs ** 2 * (1 - 2.4 * xs)
# ax.fill_between(xs, ys, 1, color='white', alpha=1)
# ax.fill_between(np.array(data_TNF[:, 0]) * 1.e-6, np.array(data_TNF[:, 1]) * 1.e+6, 1, color='white', alpha=1)
# xs_new = np.logspace(-4, 6, 100)
# ax.fill_between(xs_new, 1, 1.e-25, color='white', alpha=1)

ax.text(3.e-9, 1.e+19, r'C$\nu$B', fontsize=18, color='cyan')
ax.text(3.e-9, 1.e+19, r'C$\nu$B', fontsize=18, color='darkblue')

ax.text(1.e-11, 5.e+11, r'$\widetilde{\nu}_{BBN}(n)$', fontsize=18, color='red')
ax.text(2.e-12, 3.e+3, r'$\widetilde{\nu}_{BBN}(t)$', fontsize=18, color='magenta')

ax.text(1.e-3, 4.e+9, r'$\nu_{\mathrm{sun}}$', fontsize=22, color='yellow', rotation=37)
ax.text(1.e-3, 4.e+9, r'$\nu_{\mathrm{sun}}$', fontsize=22, color='indianred', rotation=37)

ax.text(1.e+4, 3.e-8, r'$\nu_{\mathrm{atm}}$', fontsize=22, color='lightcyan', rotation=-50)
ax.text(1.e+4, 3.e-8, r'$\nu_{\mathrm{atm}}$', fontsize=22, color='cornflowerblue', rotation=-50)

ax.text(4.e+1, 2.e-1, r'$\nu_{\mathrm{SN}}$', fontsize=20, color='saddlebrown')

# ax.legend(loc='upper right')

fig_name = figure_folder + r"\guns.pdf"
plt.savefig(fig_name)

f.clf()

#######################################################
# ## Photons + Neutrinos (for 1 maxE distribution) ## #
#######################################################

f, ax = plt.subplots(1, 1, figsize=(fig_width * 0.85, fig_width * 0.65))
f.subplots_adjust(left=0.1, right=0.95)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$E,\,\, $МэВ')
ax.set_ylabel(r'$\mathrm{d}\Phi/\mathrm{d}E,\,\, $МэВ$^{-1}\cdot$см$^{-2}\cdot$с$^{-1}$')
ax.set_ylim(1.e-22, 1.e+24)
ax.set_xlim(1.e-13, 1.e+11)
ax.set_xticks(np.logspace(-10, 10, num=5), minor=False)
ax.set_xticks(np.logspace(-13, 14, num=28), labels=[], minor=True)
ax.set_yticks(np.logspace(-20, 20, num=5), minor=False)
ax.set_yticks(np.logspace(-22, 24, num=24), labels=[], minor=True)

# # CMB
# x1 = np.logspace(-13, -8.3, 100)
# ax.plot(x1, (8 * np.pi * x1 ** 2) / ((c ** 2 * (2 * np.pi * h) ** 3) * (np.exp(x1 / (k * T_gamma)) - 1)),
#         label='CMB', linewidth=4, color='deepskyblue')
#
# # CIB + COB
# file_CIB = data_here + r'\ebl_dominguez11.out.txt'
# data_CIB = np.genfromtxt(file_CIB, skip_header=10, usecols=(0, 1))
# data_x_CIB = 2 * np.pi * h * c / (1.e-4 * np.array(data_CIB[:, 0]))
# ax.plot(data_x_CIB, 1.e-13 * 6.24150964e+12 * np.array(data_CIB[:, 1] / np.array(data_x_CIB) ** 2),
#         label='CIB + COB', linewidth=4, color='green')
#
# # CUB (upper limit) - just bad approximately
# xc = np.logspace(-5, -3, 50)
# ax.plot(xc, 10.15e+8 / ((xc / 3.999e-5) ** 0.82 + (xc / 3.999e-5) ** 3.88),
#         label='CUB', linewidth=4, color='gold')
#
# # CXB
# x2 = np.logspace(-3, 2, 100)
# ax.plot(x2, 10.15e+1 / ((x2 / 2.999e-2) ** 1.32 + (x2 / 2.999e-2) ** 2.88),
#         label='CXB', linewidth=4, color='crimson')
#
# # CGB
# x3 = np.logspace(2, 6.3, 100)
# i_100 = 1.66e-7
# gamma = 2.28
# E_cut = 2.67e+5
# ax.plot(x3, i_100 * (x3 / 100) ** (-gamma) * np.exp(-x3 / E_cut), linewidth=4, color='darkblue')  #, label='CGB' )

x_ph, y_ph = np.loadtxt(data_here + r'\ph_background.csv', delimiter=' ', unpack=True)
ax.plot(x_ph, y_ph, linewidth=4, color='indianred')
ax.plot(x_ph, y_ph, linewidth=2, color='yellow')

x_nu, y_nu = np.loadtxt(data_here + r'\nu_background.csv', delimiter=' ', unpack=True)
ax.plot(x_nu, y_nu, linewidth=4, color='darkblue')
ax.plot(x_nu, y_nu, linewidth=2, color='cyan')

# ax.plot(x_nu[:34], y_nu[:34], linewidth=4, color='cyan')
# ax.plot(x_nu[:34], y_nu[:34], linewidth=2, color='darkblue')

energies_9_3_ph, ph_all_9_3_ph = data_extract_ph(r"\ra_2.4e+09_2.8_ph_pri.txt",
                                                 r"\ra_2.4e+09_2.8_ph_sec.txt")
ax.plot(energies_9_3_ph[110:] * 1.e+3,
        n_pbh_beta(beta_shtrix(2.4e+9), 2.4e+9, 2.8) * 3.e-4 * (c / (4 * np.pi)) * ph_all_9_3_ph[110:] / 1.e+3,
        label=r"$\sigma = 2.78,\,\, M_{PBH}=2.4\times10^{9}$ г", linewidth=2, linestyle='dashed', color='orange')

energies_9_3_nu, nu_all_9_3_nu = data_extract_nu_new(r"\ra_2.4e+09_2.8_pri.txt",
                                                     r"\ra_2.4e+09_2.8_e_sec.txt",
                                                     r"\ra_2.4e+09_2.8_mu_sec.txt",
                                                     r"\ra_2.4e+09_2.8_tau_sec.txt")
ax.plot(energies_9_3_nu * 1.e+3,
        n_pbh_beta(beta_shtrix(2.4e+9), 2.4e+9, 2.8) * (c / (4 * np.pi)) * nu_all_9_3_nu / 1.e+3,
        label=r"$\sigma = 2.78,\,\, M_{PBH}=2.4\times10^{9}$ г", linewidth=2, linestyle='dashed', color='darkblue')

# Создаем легенду
yel_line = mlines.Line2D([], [], color='orange', linestyle='dashed', label=r'$\gamma$')
blue_line = mlines.Line2D([], [], color='darkblue', linestyle='dashed',
                          label=r'$\nu_{e, \mu, \tau},\,\widetilde{\nu}_{e, \mu, \tau}$')
dashed_line = mlines.Line2D([], [], color='black', linestyle='--',
                            label=r'$\sigma = 2.78, M_{c} = 2.4 \times 10^{9}$ г')

ax.text(1.e-7, 1.0e+15, r'$\gamma$-фон', fontsize=18, color='yellow')
ax.text(1.e-7, 1.0e+15, r'$\gamma$-фон', fontsize=18, color='indianred')
ax.text(1.e+0, 1.e+8, r'$\nu$-фон', fontsize=18, color='cyan')
ax.text(1.e+0, 1.e+8, r'$\nu$-фон', fontsize=18, color='darkblue')

plt.legend(handles=[dashed_line, blue_line, yel_line], loc='upper right')

fig_name = figure_folder + r"\sum_results_maxE.pdf"
plt.savefig(fig_name)

f.clf()

#######################################################################
# ## SAME Photons + Neutrinos (for 1 maxE distribution) BUT in ERG ## #
#######################################################################

f, ax = plt.subplots(1, 1, figsize=(fig_width * 0.85, fig_width * 0.65))
f.subplots_adjust(left=0.1, right=0.95)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$E,\,\, $МэВ')
ax.set_ylabel(r'$E\:\mathrm{d}\Phi/\mathrm{d}E,\,\, $см$^{-2}\cdot$с$^{-1}$')
ax.set_ylim(1.e-22, 1.e+24)
ax.set_xlim(1.e-13, 1.e+11)
ax.set_xticks(np.logspace(-10, 10, num=5), minor=False)
ax.set_xticks(np.logspace(-13, 14, num=28), labels=[], minor=True)
ax.set_yticks(np.logspace(-20, 20, num=5), minor=False)
ax.set_yticks(np.logspace(-22, 24, num=24), labels=[], minor=True)

ax.plot(x_ph, x_ph * y_ph, linewidth=4, color='indianred')
ax.plot(x_ph, x_ph * y_ph, linewidth=2, color='yellow')

ax.plot(x_nu, x_nu * y_nu, linewidth=4, color='darkblue')
ax.plot(x_nu, x_nu * y_nu, linewidth=2, color='cyan')

ax.plot(energies_9_3_ph[110:] * 1.e+3,
        n_pbh_beta(beta_shtrix(2.4e+9), 2.4e+9, 2.8) * 3.e-4 * energies_9_3_ph[110:] * 1.e+3 * (
                c / (4 * np.pi)) * ph_all_9_3_ph[110:] / 1.e+3,
        label=r"$\sigma = 2.78,\,\, M_{PBH}=2.4\times10^{9}$ г", linewidth=2, linestyle='dashed', color='orange')

ax.plot(energies_9_3_nu * 1.e+3,
        n_pbh_beta(beta_shtrix(2.4e+9), 2.4e+9, 2.8) * energies_9_3_nu * 1.e+3 * (
                c / (4 * np.pi)) * nu_all_9_3_nu / 1.e+3,
        label=r"$\sigma = 2.78,\,\, M_{PBH}=2.4\times10^{9}$ г", linewidth=2, linestyle='dashed', color='darkblue')

# Создаем легенду
yel_line = mlines.Line2D([], [], color='orange', linestyle='dashed', label=r'$\gamma$')
blue_line = mlines.Line2D([], [], color='darkblue', linestyle='dashed',
                          label=r'$\nu_{e, \mu, \tau},\,\widetilde{\nu}_{e, \mu, \tau}$')
dashed_line = mlines.Line2D([], [], color='black', linestyle='--',
                            label=r'$\sigma = 2.78, M_{c} = 2.4 \times 10^{9}$ г')

ax.text(1.e-7, 1.0e+8, r'$\gamma$-фон', fontsize=18, color='yellow')
ax.text(1.e-7, 1.0e+8, r'$\gamma$-фон', fontsize=18, color='indianred')
ax.text(1.e+0, 1.e+8, r'$\nu$-фон', fontsize=18, color='cyan')
ax.text(1.e+0, 1.e+8, r'$\nu$-фон', fontsize=18, color='darkblue')

ax.text(1.e-11, 1.e+6, r'C$\nu$B', fontsize=18, color='cyan')
ax.text(1.e-11, 1.e+6, r'C$\nu$B', fontsize=18, color='darkblue')

plt.legend(handles=[dashed_line, blue_line, yel_line], loc='upper right')

fig_name = figure_folder + r"\sum_results_max_in_erg.pdf"
plt.savefig(fig_name)

f.clf()

####################################################################################
# ## 2 GRAPHICS for Photons + Neutrinos (for 1 maxE distribution and 1 max_erg) ## #
####################################################################################

f, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(fig_width * 0.85, fig_width * 0.85),
                             gridspec_kw={'height_ratios': [3, 2]}, sharex='all')
f.subplots_adjust(wspace=0.0, hspace=0.0, left=0.1, right=0.95)

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylabel(r'$\mathrm{d}\Phi/\mathrm{d}E,\,\, $МэВ$^{-1}\cdot$см$^{-2}\cdot$с$^{-1}$')
ax1.set_ylim(1.e-22, 1.e+24)
ax1.set_xlim(1.e-13, 1.e+11)
ax1.tick_params(axis='both', which='major', top=True, right=True, labelright=False, labeltop=False)
ax1.set_xticks(np.logspace(-10, 10, num=5), minor=False)
ax1.set_xticks(np.logspace(-13, 11, num=25), labels=[], minor=True)
ax1.set_yticks(np.logspace(-20, 20, num=5), minor=False)
ax1.set_yticks(np.logspace(-22, 24, num=24), labels=[], minor=True)

# ax2.tick_params(axis='both', which='major', top=True, right=True, labelright=False, labeltop=False)
# ax2.set_xticks(np.logspace(-10, 10, num=5), minor=False)
# ax2.set_xticks(np.logspace(-13, 14, num=28), labels=[], minor=True)
# ax2.set_yticks(np.logspace(-20, 20, num=5), minor=False)
# ax2.set_yticks(np.logspace(-22, 24, num=24), labels=[], minor=True)

ax2.set_ylabel(r'$E\:\mathrm{d}\Phi/\mathrm{d}E,\,\, $см$^{-2}\cdot$с$^{-1}$')
ax2.set_xlabel(r'$E,\,\, $МэВ')
ax2.tick_params(axis='both', which='major', top=False, right=True, labelright=False, labeltop=False)
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_ylim(1.e-14, 1.e+14)
ax2.set_xlim(1.e-13, 1.e+11)
ax2.set_yticks(np.logspace(-10, 10, num=3), minor=False)
ax2.set_yticks(np.logspace(-14, 15, num=30), labels=[], minor=True)
ax2.set_xticks(np.logspace(-10, 10, num=5), minor=False)
ax2.set_xticks(np.logspace(-13, 11, num=25), labels=[], minor=True)

plt.setp(ax1.get_xticklabels(), visible=False)

ax1.plot(x_ph, y_ph, linewidth=4, color='indianred')
ax1.plot(x_ph, y_ph, linewidth=2, color='yellow')

ax1.plot(x_nu, y_nu, linewidth=4, color='darkblue')
ax1.plot(x_nu, y_nu, linewidth=2, color='cyan')

ax1.plot(energies_9_3_ph[110:] * 1.e+3,
         n_pbh_beta(beta_shtrix(2.4e+9), 2.4e+9, 2.8) * 3.e-4 * (
                 c / (4 * np.pi)) * ph_all_9_3_ph[110:] / 1.e+3,
         label=r"$\sigma = 2.78,\,\, M_{PBH}=2.4\times10^{9}$ г", linewidth=2, linestyle='dashed', color='orange')

ax1.plot(energies_9_3_nu * 1.e+3,
         n_pbh_beta(beta_shtrix(2.4e+9), 2.4e+9, 2.8) * (
                 c / (4 * np.pi)) * nu_all_9_3_nu / 1.e+3,
         label=r"$\sigma = 2.78,\,\, M_{PBH}=2.4\times10^{9}$ г", linewidth=2, linestyle='dashed', color='darkblue')

# Создаем легенду
yel_line = mlines.Line2D([], [], color='orange', linestyle='dashed', label=r'$\gamma$')
blue_line = mlines.Line2D([], [], color='darkblue', linestyle='dashed',
                          label=r'$\nu_{e, \mu, \tau},\,\widetilde{\nu}_{e, \mu, \tau}$')
dashed_line = mlines.Line2D([], [], color='black', linestyle='--',
                            label=r'$\sigma = 2.78, M_{c} = 2.4 \times 10^{9}$ г')

ax1.legend(handles=[dashed_line, blue_line, yel_line], loc='upper right')

ax2.plot(x_ph, x_ph * y_ph, linewidth=4, color='indianred')
ax2.plot(x_ph, x_ph * y_ph, linewidth=2, color='yellow')

ax2.plot(x_nu, x_nu * y_nu, linewidth=4, color='darkblue')
ax2.plot(x_nu, x_nu * y_nu, linewidth=2, color='cyan')

ax2.plot(energies_9_3_ph[110:] * 1.e+3,
         n_pbh_beta(beta_shtrix(2.4e+9), 2.4e+9, 2.8) * 3.e-4 * energies_9_3_ph[110:] * 1.e+3 * (
                 c / (4 * np.pi)) * ph_all_9_3_ph[110:] / 1.e+3,
         label=r"$\sigma = 2.78,\,\, M_{PBH}=2.4\times10^{9}$ г", linewidth=2, linestyle='dashed', color='orange')

ax2.plot(energies_9_3_nu * 1.e+3,
         n_pbh_beta(beta_shtrix(2.4e+9), 2.4e+9, 2.8) * energies_9_3_nu * 1.e+3 * (
                 c / (4 * np.pi)) * nu_all_9_3_nu / 1.e+3,
         label=r"$\sigma = 2.78,\,\, M_{PBH}=2.4\times10^{9}$ г", linewidth=2, linestyle='dashed', color='darkblue')

fig_name = figure_folder + r"\two_graph_results.pdf"
plt.savefig(fig_name)

f.clf()

#######################################################
# ## Photons + Neutrinos (for 3 maxE distribution) ## #
#######################################################

f, ax = plt.subplots(1, 1, figsize=(fig_width * 0.85, fig_width * 0.65))
f.subplots_adjust(left=0.1, right=0.95)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$E,\,\, $МэВ')
ax.set_ylabel(r'$\mathrm{d}\Phi/\mathrm{d}E,\,\, $МэВ$^{-1}\cdot$см$^{-2}\cdot$с$^{-1}$')
ax.set_ylim(1.e-22, 1.e+24)
ax.set_xlim(1.e-13, 1.e+11)
ax.set_xticks(np.logspace(-10, 10, num=5), minor=False)
ax.set_xticks(np.logspace(-13, 14, num=28), labels=[], minor=True)
ax.set_yticks(np.logspace(-20, 20, num=5), minor=False)
ax.set_yticks(np.logspace(-22, 24, num=24), labels=[], minor=True)

ax.plot(x_ph, y_ph, linewidth=4, color='indianred')
ax.plot(x_ph, y_ph, linewidth=2, color='yellow')

ax.plot(x_nu, y_nu, linewidth=4, color='darkblue')
ax.plot(x_nu, y_nu, linewidth=2, color='cyan')

# M = 10^11 (sigma = 1.0)
energies_11_1_ph, nu_all_11_1_ph = data_extract_ph(r"\ra_1.0e+11_1.0_ph_pri.txt",
                                                   r"\ra_1.0e+11_1.0_ph_sec.txt")
ax.plot(energies_11_1 * 1.e+3,
        n_pbh_beta(beta_shtrix(1.e+11), 1.e+11, 1.0) * (c / (4 * np.pi)) * nu_all_11_1 / 1.e+3,
        label=r"$\sigma = 1.0,\,\, M_{PBH}=10^{11}$ г", linewidth=2, linestyle='dashed', color='indianred')
ax.plot(energies_11_1_ph[105:] * 1.e+3,
        n_pbh_beta(beta_shtrix(1.e+11), 1.e+11, 1.0) * (c / (4 * np.pi)) * nu_all_11_1_ph[105:] / 1.e+3,
        label=r"$\sigma = 1.0,\,\, M_{PBH}=10^{11}$ г", linewidth=2, linestyle='dotted', color='indianred')

# M = 10^13 (sigma = 2.0)
energies_13_2_ph, nu_all_13_2_ph = data_extract_ph(r"\ra_1.0e+13_2.0_ph_pri.txt",
                                                   r"\ra_1.0e+13_2.0_ph_sec.txt")
ax.plot(energies_13_2 * 1.e+3,
        n_pbh_beta(beta_shtrix(1.e+13), 1.e+13, 2.0) * (c / (4 * np.pi)) * nu_all_13_2 / 1.e+3,
        label=r"$\sigma = 2.0,\,\, M_{PBH}=10^{13}$ г", linewidth=2, linestyle='dashed', color='lawngreen')
ax.plot(energies_13_2_ph[105:] * 1.e+3,
        n_pbh_beta(beta_shtrix(1.e+13) * 0.03, 1.e+13, 2.0) * (c / (4 * np.pi)) * nu_all_13_2_ph[105:] / 1.e+3,
        label=r"$\sigma = 2.0,\,\, M_{PBH}=10^{13}$ г", linewidth=2, linestyle='dotted', color='lawngreen')

# M = 10^15 (sigma = 3.0)
energies_15_3_ph, nu_all_15_3_ph = data_extract_ph(r"\ra_1.0e+15_3.0_ph_pri.txt",
                                                   r"\ra_1.0e+15_3.0_ph_sec.txt")
ax.plot(energies_15_3 * 1.e+3,
        n_pbh_beta(1.e-26, 1.e+15, 3.0) * (c / (4 * np.pi)) * nu_all_15_3 / 1.e+3,
        label=r"$\sigma = 3.0,\,\, M_{PBH}=10^{15}$ г", linewidth=2, linestyle='dashed', color='royalblue')
ax.plot(energies_15_3_ph[105:] * 1.e+3,
        n_pbh_beta(1.e-26, 1.e+15, 3.0) * (c / (4 * np.pi)) * nu_all_15_3_ph[105:] / 1.e+3,
        label=r"$\sigma = 3.0,\,\, M_{PBH}=10^{15}$ г", linewidth=2, linestyle='dotted', color='royalblue')

# custom legend
red_line = mlines.Line2D([], [], color='indianred', label=r'$M = 10^{11}\:\text{г}, \sigma = 1.0$')
blue_line = mlines.Line2D([], [], color='royalblue', label=r'$M = 10^{15}\:\text{г}, \sigma = 3.0$')
green_line = mlines.Line2D([], [], color='lawngreen', label=r'$M = 10^{13}\:\text{г}, \sigma = 2.0$')
dashed_line = mlines.Line2D([], [], color='black', linestyle='--',
                            label=r'$\nu_{e, \mu, \tau},\,\widetilde{\nu}_{e, \mu, \tau}$')
dotted_line = mlines.Line2D([], [], color='black', linestyle=':', label=r'$\gamma$')

ax.text(1.e-7, 1.0e+15, r'$\gamma$', fontsize=18, color='yellow')
ax.text(1.e-7, 1.0e+15, r'$\gamma$', fontsize=18, color='indianred')
ax.text(1.e-8, 1.e+6, r'$\nu$', fontsize=18, color='cyan')
ax.text(1.e-8, 1.e+6, r'$\nu$', fontsize=18, color='darkblue')

plt.legend(handles=[red_line, green_line, blue_line, dashed_line, dotted_line], loc='upper right')

fig_name = figure_folder + r"\sum_results_three_distr.pdf"
plt.savefig(fig_name)

f.clf()

####################################################
# Comparison for E-multiplied #
####################################################

# f, ax = plt.subplots(1, 1, figsize=(fig_width * 0.85, fig_width * 0.65))
# f.subplots_adjust(left=0.1, right=0.95)
# ax.set_xscale('log')
# ax.set_yscale('log')
# ax.set_xlabel(r'$E\,\,[$МэВ$]$')
# ax.set_ylabel(r'$E\,\mathrm{d}\Phi/\mathrm{d}E\,\, [$см$^{2}\cdot$с$^{-1}$$]$')
# ax.set_ylim(1.e-22, 1.e+24)
# ax.set_xlim(1.e-13, 1.e+11)
# ax.set_xticks(np.logspace(-10, 10, num=5), minor=False)
# ax.set_xticks(np.logspace(-13, 14, num=28), labels=[], minor=True)
# ax.set_yticks(np.logspace(-20, 20, num=5), minor=False)
# ax.set_yticks(np.logspace(-22, 24, num=24), labels=[], minor=True)
#
# # Plotting Relic Neutrinos
# x = np.logspace(-12, -8, 100)
# ax.plot(x, x * relic_nu(x), linewidth=4, color='cyan')
# ax.plot(x, x * relic_nu(x), label='Реликтовые', linewidth=1, color='darkblue')
#
# ax.axvline(x=1.e+0, color='black', linewidth=1, linestyle='dashed')
#
# # Energy density for: 10^9 (log) sigma=3 (Pythia) NU
# energies_e9_3, nu_all_e9_3 = data_extract_nu(r"\results_alternate_nu_pri_log_9_3.txt",
#                                              r"\results_alternate_nu_e_sec_log_9_3.txt",
#                                              r"\results_alternate_nu_mu_sec_log_9_3.txt",
#                                              r"\results_alternate_nu_tau_sec_log_9_3.txt")
#
# E_nu = energies_e9_3 * 1.e+3
# n_PBH = n_pbh_beta(beta_shtrix(1.e+9), 1.e+9, 3.0)
# y_calc = n_PBH * (c / (4 * np.pi)) * E_nu * nu_all_e9_3 / 1.e+3
# ax.plot(E_nu, y_calc, label=r"$\sigma = 3.0,\,\, M_{PBH}=10^{9}$ г", linewidth=3, linestyle='dashed', color='indigo')
#
# area_relic = integrate.trapz(x * relic_nu(x), x)
# area_calc = integrate.trapz(y_calc, E_nu)
#
# print("relic:", "{:e}".format(area_relic))
# print("calc:", "{:e}".format(area_calc))
# print('calc / relic:', area_calc / area_relic)
#
# ax.legend(loc='upper right')
#
# fig_name = figure_folder + r"\energy_comparison_E.pdf"
# plt.savefig(fig_name)
#
# f.clf()
