import numpy as np
from preamble import *

# 2 - Constants & functions

c = 2.99792458e+10  # [cm / c]
h = 6.582119569e-22  # [Mev * с]
k = 8.617333262e-11  # [Mev / K]
T_nu = 1.945e+0  # [K]
T_gamma = 2.7277e+0  # [K]
GG = 6.67e-8  # [cm^3 * g^-1 * s^-2]


def Tpbh(m: float) -> float:
    return h * c ** 3 / (8 * np.pi * GG * k * m)  # K


def T_universe(t: float) -> float:
    return 1.52e+10 * t ** (-1/2)  # K


def Tpbh_2(t: float) -> float:
    return 3.04e-13 / t


def t_init(m):
    return 2 * GG * m / c ** 3


def t_evap(m):
    return 2.0e+18 * (m * 10 ** (-15)) ** 3


Mpbhs = [1.e+8, 1.e+9, 1.e+10, 1.e+11, 1.e+12, 1.e+13, 1.e+14, 1.e+15, 1.e+16, 1.e+17, 1e+18, 1.e+19, 1.e+20]
print(T_universe(1))
for M in Mpbhs:
    print('For M = {:.1e}:'.format(M))
    print('t_init = {:.2e}'.format(t_init(M)))
    print('T_U = {:.2e}'.format(T_universe(t_init(M))))
    print('Tpbh = {:.2e}'.format(Tpbh(M)))
    print('Tpbh_2 = {:.2e}'.format(Tpbh_2(t_init(M))))
    print()

data = r"D:\Downloads\BlackHawk\results"

Mt_15 = np.genfromtxt(data + r'\1.0e+15_0.0' + r'\life_evolutions.txt', skip_header=4)
ts_15, Ms_15 = Mt_15[:, 0], Mt_15[:, 1]

Mt_11 = np.genfromtxt(data + r'\1.0e+11_0.0' + r'\life_evolutions.txt', skip_header=4)
ts_11, Ms_11 = Mt_11[:, 0], Mt_11[:, 1]

Mt_13 = np.genfromtxt(data + r'\1.0e+13_0.0' + r'\life_evolutions.txt', skip_header=4)
ts_13, Ms_13 = Mt_13[:, 0], Mt_13[:, 1]

M2 = np.logspace(9, 16, 100)

f, ax = plt.subplots(1, 1, figsize=(fig_width*0.85, fig_width * 0.65))
f.subplots_adjust(left=0.1, right=0.95)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$M_{\mathrm{PBH}},$ г')
ax.set_ylabel(r'$t,$ с')
ax.set_ylim(1.e+0, 1.e+24)
ax.set_xlim(1.e9, 1.e+16)

ax.plot(Ms_15, ts_15, linewidth=2, color='royalblue', linestyle='dotted')
ax.plot(Ms_13, ts_13, linewidth=2, color='lawngreen', linestyle='dotted')
ax.plot(Ms_11, ts_11, linewidth=2, color='indianred', linestyle='dotted')
ax.plot(M2, t_init(M2), linewidth=2, color='indianred', linestyle='dashed')
ax.plot(M2, t_evap(M2), linewidth=2, color='violet', linestyle='dashed', label=r'$\tau_{evap}(M_i)$')
ax.axhline(y=13.8e+9 * 365. * 24. * 3600., color='black', linewidth=1, linestyle='dashed')
ax.axvline(x=5.5e+14, color='black', linewidth=1, linestyle='dashed')

ax.text(1.e+12, 3.0e+18, r'$M=10^{15}\,$ г', fontsize=15, color='royalblue')
ax.text(3.e+10, 1.e+12, r'$M=10^{13}\,$ г', fontsize=15, color='lawngreen')
ax.text(3.1e+9, 1.0e+6, r'$M=10^{11}\,$ г', fontsize=15, color='indianred')

f.text(0.2, 0.65, r'$t_U = 4.35 \times 10^{17}\,$ с', ha='center', fontsize=15)
ax.text(3.4e+14, 3.0e+7, r'$M_*=5.1 \times 10^{14}\,$ г', fontsize=15, rotation=90)

ax.legend(loc='upper right')
fig_name = figure_folder + r"\mass_via_time.pdf"
plt.savefig(fig_name)
f.clf()

# ########## CHECK ########## #
f, ax = plt.subplots(1, 1, figsize=(fig_width*0.85, fig_width * 0.65))
f.subplots_adjust(left=0.1, right=0.95)
ax.set_xscale('linear')
ax.set_yscale('linear')
ax.set_xlabel(r'$\sigma$')
ax.set_ylabel(r'$n$')
# ax.set_ylim(1.e+0, 1.e+24)
ax.set_xlim(0.0, 4.0)

sigmas = np.linspace(0.0, 5.0, 100)

n = np.exp(-9 * sigmas ** 2 / 8)

ax.plot(sigmas, n)

fig_name = figure_folder + r"\check.pdf"
plt.savefig(fig_name)
f.clf()

