from preamble import *
import matplotlib.lines as mlines


def lognormal_mass_spec(mu, sigma, a_1, m):
    return a_1 * np.exp(-((np.log(m / mu)) ** 2) / (2 * sigma ** 2)) / (sigma * np.sqrt(2 * np.pi))


def format_ticks(axis):
    axis.set_yscale('log')
    axis.set_xscale('log')
    axis.set_ylim(1.e-8, 1.e+2)
    axis.set_yticks(np.logspace(-7, 1, 5), ['$10^{-7}$', '$10^{-5}$', '$10^{-3}$', '$10^{-1}$', '$10^{1}$'])
    if axis != ax4:
        axis.set_xlim(1.e+9, 1.e+15)
        axis.set_xticks(np.logspace(10, 14, 5), ['$10^{10}$', '', '$10^{12}$', '', '$10^{14}$'])


def format_axes(axis):
    # Убираем обводку графика
    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)
    axis.spines['bottom'].set_visible(False)
    axis.spines['left'].set_visible(False)

    # Убираем метки (тики) на осях вместе с надписями
    axis.tick_params(axis='both', which='both', bottom=False, top=False, left=False, right=False,
                     labelbottom=False, labelleft=False, labelright=False, labeltop=False)


c = 2.99792458e+10  # [cm / c]
h = 6.582119569e-22  # [Mev * с]
k = 8.617333262e-11  # [Mev / K]
T_nu = 1.945e+0  # [K]
T_gamma = 2.7277e+0  # [K]
GG = 6.67e-8  # [cm^3 * g^-1 * s^-2]
kt = h * c ** 3 / (8 * np.pi * GG * 1.e+26)
print(kt / k)
print(T_gamma)

# 3 - Plotting mass spectrum via lognormal

# f, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(fig_width, fig_width / 3.0))
# f.subplots_adjust(wspace=0.0, bottom=0.15, left=0.1, right=0.95)
# f.text(0.52, 0.03, r'$M,\,\text{г}$', ha='center')
# f.text(0.03, 0.5, r'$d\mathcal{N}/dM$', va='center', rotation='vertical')
# for ax in (ax1, ax2, ax3, ax4):
#     format_ticks(ax)
#
# m_log = np.logspace(1, 15, 1000)
#
# # # Dirac # #
# ax1.axvline(x=1.e+12, ymax=0.8, color='black', linewidth='1', linestyle='solid')
#
# ax1.plot(m_log, lognormal_mass_spec(1.e+12, 0.03, 0.1, m_log), linewidth=1, color='darkgrey', linestyle='solid')
#
# # # power-law
# ax2.yaxis.set_major_formatter(ticker.NullFormatter())
#
# gamma = -1. / 2.
# a_2 = (gamma - 1) / (10 ** (14 * (gamma - 1)) - 10 ** (11 * (gamma - 1)))  #
# m_plaw = np.logspace(11, 14, 100)
# dNdM_plaw = a_2 * m_plaw ** (gamma - 1.)
#
# ax2.plot(m_plaw, dNdM_plaw, linewidth=1, color='black')
# ax2.axvline(x=1.e+11, ymax=0.8, color='black', linewidth='1', linestyle='dashed')
# ax2.axvline(x=1.e+14, ymax=0.4, color='black', linewidth='1', linestyle='dashed')
#
# ax2.plot(m_log, lognormal_mass_spec(1.e+12, 0.9, 6., m_log), linewidth=1, color='darkgrey', linestyle='solid')
#
# # critical collapse
# ax3.yaxis.set_major_formatter(ticker.NullFormatter())
#
# mf, k = 1.e+12, 2.85
# a_3 = (k / mf ** k) * (1 / (np.exp(-10 ** (k * 10) / mf ** k) - np.exp(-1000 ** (k * 5) / mf ** k)))  # 1.79e-34
# m_cc = np.logspace(9, 14, 100)
# dNdM_cc = a_3 * (m_cc ** k) * np.exp(-(m_cc / mf) ** k)
#
# ax3.plot(m_cc, dNdM_cc, linewidth=1, color='black', linestyle='solid')
#
# ax3.plot(m_log, lognormal_mass_spec(8.e+11, 0.8, 1., m_log), linewidth=1, color='grey', linestyle='solid')
#
# # peak theory
# ax4.set_xlim(1.e+5, 1.e+13)
# ax4.set_xticks(np.logspace(6, 12, 4), ['$10^{6}$', '$10^{8}$', '$10^{10}$', '$10^{12}$'])
# data_spectrum = np.genfromtxt(r"C:\cygwin64\home\BlackHawk\results\ms_peak_theory\BH_spectrum.txt", skip_header=3)
# ax4.yaxis.set_major_formatter(ticker.NullFormatter())
# ax4.plot(data_spectrum[:, 0], 1.e+30 * data_spectrum[:, 1], linewidth=1, color='black', linestyle='solid')
# ax4.plot(m_log, lognormal_mass_spec(1.e+9, 1., 10., m_log), linewidth=1, color='grey', linestyle='solid')
#
# fig_name = figure_folder + r"\lognormal_approx.pdf"
# plt.savefig(fig_name)
#
# f.clf()

# # ---------------------------------------------- #
# # ------------- 3 mass spectrum ---------------- #
# # ---------------------------------------------- #

f, ax = plt.subplots(1, 1, figsize=(fig_width * 0.85, fig_width * 0.65))
f.subplots_adjust(left=0.1, right=0.95)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$M, \,$г')
ax.set_ylabel(r'$d\mathcal{N}/dM, \text{с}^{-1} \text{ГэВ}^{-1}$')
ax.set_ylim(1.e-10, 1.e+2)
ax.set_xlim(1.e+9, 1.e+16)

# -------------------- log-normal ----------------- #

mu, sigma, a_1 = 1.e+12, 0.5, 1.  # mean and standard deviation
m_log = np.logspace(8, 16, 1000)
label_log = r"$\psi_{\mathrm{ln}} \propto \frac{1}{M} \exp{\left [-\frac{\ln^2{(M/M_{\mathrm{PBH}})}}{2\sigma^2}\right]}$"


def dNdM_log(mu, sigma, a_1, m_log):
    return a_1 * np.exp(-((np.log(m_log / mu)) ** 2) / (2 * sigma ** 2)) / (sigma * np.sqrt(2 * np.pi))


# -------------------- power-law ------------------ #

gamma = -1. / 2.
a_2 = (gamma - 1) / (10 ** (15 * (gamma - 1)) - 10 ** (10 * (gamma - 1)))  # 1.5e+15
m_plaw = np.logspace(10, 15, 1000)
dNdM_plaw = a_2 * m_plaw ** (gamma - 1.)
label_plaw = r"$\psi_{\mathrm{pl}} \propto M^{\gamma-2}$"

# -------------------- critical collapse ------------------ #

mf, k = 1.e+12, 2.85
a_3 = (k / mf ** k) * (1 / (np.exp(-10 ** (k * 10) / mf ** k) - np.exp(-1000 ** (k * 5) / mf ** k)))  # 1.79e-34
m_cc = np.logspace(8, 20, 1000)
dNdM_cc = a_3 * (m_cc ** k) * np.exp(-(m_cc / mf) ** k)
label_cc = r"$\psi_{\mathrm{cc}} \propto M^{1.85}\exp{\left [-\left(\frac{M}{M_{f}}\right)^{2.85}\right]}$"

# --------------------------------------------------#

ax.plot(m_log, dNdM_log(mu, sigma, a_1, m_log), label=label_log, linewidth=2, color='blue', linestyle='solid')

ax.plot(m_plaw, dNdM_plaw, label=label_plaw, linewidth=2, color='lime', linestyle='dashed')
ax.axvline(x=1.e+10, ymax=0.85, color='black', linewidth='1', linestyle='dashed')
ax.axvline(x=1.e+15, ymax=0.22, color='black', linewidth='1', linestyle='dashed')

ax.plot(m_cc, dNdM_cc, label=label_cc, linewidth=2, color='red', linestyle='dashed')
# ax.grid(which='major', axis='both', linewidth=1, color='grey')
ax.legend(loc='upper right')

fig_name = figure_folder + r"\three_mass_spectrum.pdf"

plt.savefig(fig_name)
f.clf()

# # ----------------------------------------------- #
# # ------------- nu-inst-spectrum ---------------- #
# # ----------------------------------------------- #

f, ax = plt.subplots(1, 1, figsize=(fig_width * 0.9, fig_width * 0.65))
f.subplots_adjust(left=0.1, right=0.95)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$E,\,$ МэВ')
ax.set_ylabel(r'$\mathrm{d}^{2}N/\mathrm{d}t\mathrm{d}E,\, \,$МэВ$^{-1}\cdot$с$^{-1}$')
ax.set_ylim(1.e+15, 1.e+27)
ax.set_xlim(1.e-2, 3.e+6)

# ax.set_yticks(np.logspace(15, 27, 12),
#               ['$10^{15}$', '$10^{16}$', '$10^{17}$', '$10^{18}$', '$10^{19}$', '$10^{20}$', '$10^{21}$', '$10^{22}$',
#                '$10^{23}$', '$10^{24}$', '$10^{25}$', '$10^{27}$'])

data_primary_e15 = np.genfromtxt(data_folder + r"\e15\instantaneous_primary_spectra.txt", skip_header=2)
data_secondary_e15 = np.genfromtxt(data_folder + r"\e15\instantaneous_secondary_spectra.txt", skip_header=2)

data_primary_e13 = np.genfromtxt(data_folder + r"\e13\instantaneous_primary_spectra.txt", skip_header=2)
data_secondary_e13 = np.genfromtxt(data_folder + r"\e13\instantaneous_secondary_spectra.txt", skip_header=2)

data_primary_e11 = np.genfromtxt(data_folder + r"\e11\instantaneous_primary_spectra.txt", skip_header=2)
data_secondary_e11 = np.genfromtxt(data_folder + r"\e11\instantaneous_secondary_spectra.txt", skip_header=2)

max15nu = np.max(data_primary_e15[:, 6])
max15ph = np.max(data_primary_e15[:, 1])
rel = max15nu / max15ph
print(rel)

# ax.plot(1.e+3 * data_primary_e13[:, 0], 1.e-3 * data_primary_e13[:, 1] * rel,
#         linewidth=3, color='orange', linestyle='dotted')

ax.plot(1.e+3 * data_primary_e15[:, 0], 1.e-3 * data_primary_e15[:, 6],
        linewidth=3, color='royalblue', linestyle='dotted', label=r'$\nu_{e, \mu, \tau}^{\mathrm{pri}}$')
ax.plot(1.e+3 * data_primary_e15[:, 0],
        1.e-3 * (data_secondary_e15[:, 3] + data_secondary_e15[:, 4] + data_secondary_e15[:, 5] - data_primary_e15[:,
                                                                                                  6]),
        linewidth=3, color='royalblue', linestyle='dashed', label=r'$\nu_{e, \mu, \tau}^{\mathrm{sec}}$')
ax.plot(1.e+3 * data_primary_e15[:, 0],
        1.e-3 * (data_secondary_e15[:, 3] + data_secondary_e15[:, 4] + data_secondary_e15[:, 5]),
        linewidth=3, color='royalblue', linestyle='solid')

ax.plot(1.e+3 * data_primary_e13[:, 0], 1.e-3 * data_primary_e13[:, 6],
        linewidth=3, color='lawngreen', linestyle='dotted')
ax.plot(1.e+3 * data_primary_e13[:, 0],
        1.e-3 * (data_secondary_e13[:, 3] + data_secondary_e13[:, 4] + data_secondary_e13[:, 5] - data_primary_e13[:,
                                                                                                  6]),
        linewidth=3, color='lawngreen', linestyle='dashed')
ax.plot(1.e+3 * data_primary_e13[:, 0],
        1.e-3 * (data_secondary_e13[:, 3] + data_secondary_e13[:, 4] + data_secondary_e13[:, 5]),
        linewidth=3, color='lawngreen', linestyle='solid')

ax.plot(1.e+3 * data_primary_e11[:, 0], 1.e-3 * data_primary_e11[:, 6],
        linewidth=3, color='indianred', linestyle='dotted')
ax.plot(1.e+3 * data_primary_e11[:, 0],
        1.e-3 * (data_secondary_e11[:, 3] + data_secondary_e11[:, 4] + data_secondary_e11[:, 5] - data_primary_e11[:,
                                                                                                  6]),
        linewidth=3, color='indianred', linestyle='dashed')
ax.plot(1.e+3 * data_primary_e11[:, 0],
        1.e-3 * (data_secondary_e11[:, 3] + data_secondary_e11[:, 4] + data_secondary_e11[:, 5]),
        linewidth=3, color='indianred', linestyle='solid')

# Создаем легенду
red_line = mlines.Line2D([], [], color='indianred', label=r'$M = 10^{11}$ г')
blue_line = mlines.Line2D([], [], color='royalblue', label=r'$M = 10^{15}$ г')
green_line = mlines.Line2D([], [], color='lawngreen', label=r'$M = 10^{13}$ г')
dashed_line = mlines.Line2D([], [], color='black', linestyle='--',
                            label=r'$\nu_{e, \mu, \tau}^{\mathrm{sec}}$')
dotted_line = mlines.Line2D([], [], color='black', linestyle=':',
                            label=r'$\nu_{e, \mu, \tau}^{\mathrm{pri}}$')
solid_line = mlines.Line2D([], [], color='black', linestyle='solid',
                           label=r'$\nu_{e, \mu, \tau}^{\mathrm{pri}} + \nu_{e, \mu, \tau}^{\mathrm{sec}}$')

plt.legend(handles=[red_line, green_line, blue_line, solid_line, dashed_line, dotted_line], loc='upper right')

fig_name = figure_folder + r"\nu_inst_spec.pdf"
plt.savefig(fig_name)
f.clf()

# # ---------------------------------------------------- #
# # ------------- gamma- VS nu-spectrum ---------------- #
# # ---------------------------------------------------- #

font = {  # 'family': 'serif',
    # 'color': 'red',
    # 'weight': 'light',
    'size': 18
}

box1 = {'facecolor': 'whitesmoke',
        'edgecolor': 'steelblue',
        'boxstyle': 'round'
        }

box2 = {'facecolor': 'whitesmoke',
        'edgecolor': 'darkorange',
        'boxstyle': 'round'
        }

f, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(fig_width, fig_width / 2.0), sharex='all',
                             sharey='all')  # , figsize=(20*.66, 10*.66))
f.subplots_adjust(wspace=0.0, bottom=0.15, left=0.1, right=0.95)
f.text(0.52, 0.06, r'$E,\,\text{МэВ}$', ha='center', fontsize=20)
ax1.set_xscale('log')
ax1.set_yscale('log')
# ax1.set_xlabel(r'$E,\,$ МэВ$')
ax1.set_ylabel(r'$\frac{\mathrm{d}^{2}N}{\mathrm{d}t\mathrm{d}E},\, \,[$МэВ$^{-1}\cdot$с$^{-1}]$')
ax1.set_ylim(1.e+15, 1.e+27)
ax1.set_xlim(1.e-2, 3.e+6)
ax1.tick_params(axis='both', which='major', top=True, right=True, labelright=False, labeltop=False)
ax1.grid(which='major', axis='both', linestyle='--', linewidth=0.5)
ax1.set_yticks(np.logspace(16, 26, 6),
               ['$10^{16}$', '$10^{18}$', '$10^{20}$', '$10^{22}$', '$10^{24}$', '$10^{26}$'])
ax1.set_xticks(np.logspace(-1, 5, 4),
               ['$10^{-1}$', '$10^{1}$', '$10^{3}$', '$10^{5}$'])
ax1.text(3.e+4, 1.7e+26, r'$\nu_{e,\nu,\tau}(\widetilde{\nu}_{e,\nu,\tau})$', color='steelblue', fontdict=font,
         backgroundcolor='white', bbox=box1)

ax1.text(3.e-1, 2.e+19, r'$M=10^{15}\,$г', rotation=-30)
ax1.text(1.e+2, 1.e+21, r'$M=10^{13}\,$г', rotation=-42)
ax1.text(1.e+3, 1.e+23, r'$M=10^{11}\,$г', rotation=-40)

ax2.text(3.e-1, 2.e+19, r'$M=10^{15}\,$г', rotation=-30)
ax2.text(1.e+2, 1.e+21, r'$M=10^{13}\,$г', rotation=-42)
ax2.text(1.e+3, 1.e+23, r'$M=10^{11}\,$г', rotation=-40)

ax1.plot(1.e+3 * data_primary_e15[:, 0], 1.e-3 * data_primary_e15[:, 6],
         linewidth=2, color='cyan', linestyle='dotted')
ax1.plot(1.e+3 * data_primary_e15[:, 0],
         1.e-3 * (data_secondary_e15[:, 3] + data_secondary_e15[:, 4] + data_secondary_e15[:, 5] - data_primary_e15[:,
                                                                                                   6]),
         linewidth=2, color='cyan', linestyle='dashed')
ax1.plot(1.e+3 * data_primary_e15[:, 0],
         1.e-3 * (data_secondary_e15[:, 3] + data_secondary_e15[:, 4] + data_secondary_e15[:, 5]),
         linewidth=2, color='cyan', linestyle='solid')

ax1.plot(1.e+3 * data_primary_e13[:, 0], 1.e-3 * data_primary_e13[:, 6],
         linewidth=2, color='cyan', linestyle='dotted')
ax1.plot(1.e+3 * data_primary_e13[:, 0],
         1.e-3 * (data_secondary_e13[:, 3] + data_secondary_e13[:, 4] + data_secondary_e13[:, 5] - data_primary_e13[:,
                                                                                                   6]),
         linewidth=2, color='cyan', linestyle='dashed')
ax1.plot(1.e+3 * data_primary_e13[:, 0],
         1.e-3 * (data_secondary_e13[:, 3] + data_secondary_e13[:, 4] + data_secondary_e13[:, 5]),
         linewidth=2, color='cyan', linestyle='solid')

ax1.plot(1.e+3 * data_primary_e11[:, 0], 1.e-3 * data_primary_e11[:, 6],
         linewidth=2, color='cyan', linestyle='dotted')
ax1.plot(1.e+3 * data_primary_e11[:, 0],
         1.e-3 * (data_secondary_e11[:, 3] + data_secondary_e11[:, 4] + data_secondary_e11[:, 5] - data_primary_e11[:,
                                                                                                   6]),
         linewidth=2, color='cyan', linestyle='dashed')
ax1.plot(1.e+3 * data_primary_e11[:, 0],
         1.e-3 * (data_secondary_e11[:, 3] + data_secondary_e11[:, 4] + data_secondary_e11[:, 5]),
         linewidth=2, color='cyan', linestyle='solid')

ax1.plot(1.e+3 * data_primary_e15[:, 0], 1.e-3 * data_primary_e15[:, 6],
         linewidth=1, color='darkblue', linestyle='dotted')
ax1.plot(1.e+3 * data_primary_e15[:, 0],
         1.e-3 * (data_secondary_e15[:, 3] + data_secondary_e15[:, 4] + data_secondary_e15[:, 5] - data_primary_e15[:,
                                                                                                   6]),
         linewidth=1, color='darkblue', linestyle='dashed')
ax1.plot(1.e+3 * data_primary_e15[:, 0],
         1.e-3 * (data_secondary_e15[:, 3] + data_secondary_e15[:, 4] + data_secondary_e15[:, 5]),
         linewidth=1, color='darkblue', linestyle='solid')

ax1.plot(1.e+3 * data_primary_e13[:, 0], 1.e-3 * data_primary_e13[:, 6],
         linewidth=1, color='darkblue', linestyle='dotted')
ax1.plot(1.e+3 * data_primary_e13[:, 0],
         1.e-3 * (data_secondary_e13[:, 3] + data_secondary_e13[:, 4] + data_secondary_e13[:, 5] - data_primary_e13[:,
                                                                                                   6]),
         linewidth=1, color='darkblue', linestyle='dashed')
ax1.plot(1.e+3 * data_primary_e13[:, 0],
         1.e-3 * (data_secondary_e13[:, 3] + data_secondary_e13[:, 4] + data_secondary_e13[:, 5]),
         linewidth=1, color='darkblue', linestyle='solid')

ax1.plot(1.e+3 * data_primary_e11[:, 0], 1.e-3 * data_primary_e11[:, 6],
         linewidth=1, color='darkblue', linestyle='dotted')
ax1.plot(1.e+3 * data_primary_e11[:, 0],
         1.e-3 * (data_secondary_e11[:, 3] + data_secondary_e11[:, 4] + data_secondary_e11[:, 5] - data_primary_e11[:,
                                                                                                   6]),
         linewidth=1, color='darkblue', linestyle='dashed')
ax1.plot(1.e+3 * data_primary_e11[:, 0],
         1.e-3 * (data_secondary_e11[:, 3] + data_secondary_e11[:, 4] + data_secondary_e11[:, 5]),
         linewidth=1, color='darkblue', linestyle='solid')

ax2.grid(which='major', axis='both', linestyle='--', linewidth=0.5)
# ax2.text(1.4e+5, 1.7e+26, r'$\gamma$', color='yellow', fontdict=font)
ax2.text(4.e+5, 1.7e+26, r'\,$\gamma$\,', color='darkorange', fontdict=font, bbox=box2)

ax2.plot(1.e+3 * data_primary_e15[:, 0], 1.e-3 * data_primary_e15[:, 1],
         linewidth=2, color='yellow', linestyle='dotted')
ax2.plot(1.e+3 * data_primary_e15[:, 0], 1.e-3 * (data_secondary_e15[:, 1] - data_primary_e15[:, 1]),
         linewidth=2, color='yellow', linestyle='dashed')
ax2.plot(1.e+3 * data_primary_e15[:, 0], 1.e-3 * data_secondary_e15[:, 1],
         linewidth=2, color='yellow', linestyle='solid')

ax2.plot(1.e+3 * data_primary_e13[:, 0], 1.e-3 * data_primary_e13[:, 1],
         linewidth=2, color='yellow', linestyle='dotted')
ax2.plot(1.e+3 * data_primary_e13[137:, 0], 1.e-3 * (data_secondary_e13[137:, 1] - data_primary_e13[137:, 1]),
         linewidth=2, color='yellow', linestyle='dashed')
ax2.plot(1.e+3 * data_primary_e13[137:, 0], 1.e-3 * data_secondary_e13[137:, 1],
         linewidth=2, color='yellow', linestyle='solid')

ax2.plot(1.e+3 * data_primary_e11[:, 0], 1.e-3 * data_primary_e11[:, 1],
         linewidth=2, color='yellow', linestyle='dotted')
ax2.plot(1.e+3 * data_primary_e11[:, 0], 1.e-3 * (data_secondary_e11[:, 1] - data_primary_e11[:, 1]),
         linewidth=2, color='yellow', linestyle='dashed')
ax2.plot(1.e+3 * data_primary_e11[:, 0], 1.e-3 * data_secondary_e11[:, 1],
         linewidth=2, color='yellow', linestyle='solid')

ax2.plot(1.e+3 * data_primary_e15[:, 0], 1.e-3 * data_primary_e15[:, 1],
         linewidth=1, color='indianred', linestyle='dotted')
ax2.plot(1.e+3 * data_primary_e15[:, 0], 1.e-3 * (data_secondary_e15[:, 1] - data_primary_e15[:, 1]),
         linewidth=1, color='indianred', linestyle='dashed')
ax2.plot(1.e+3 * data_primary_e15[:, 0], 1.e-3 * data_secondary_e15[:, 1],
         linewidth=1, color='indianred', linestyle='solid')

ax2.plot(1.e+3 * data_primary_e13[:, 0], 1.e-3 * data_primary_e13[:, 1],
         linewidth=1, color='indianred', linestyle='dotted')
ax2.plot(1.e+3 * data_primary_e13[137:, 0], 1.e-3 * (data_secondary_e13[137:, 1] - data_primary_e13[137:, 1]),
         linewidth=1, color='indianred', linestyle='dashed')
ax2.plot(1.e+3 * data_primary_e13[137:, 0], 1.e-3 * data_secondary_e13[137:, 1],
         linewidth=1, color='indianred', linestyle='solid')

ax2.plot(1.e+3 * data_primary_e11[:, 0], 1.e-3 * data_primary_e11[:, 1],
         linewidth=1, color='indianred', linestyle='dotted')
ax2.plot(1.e+3 * data_primary_e11[:, 0], 1.e-3 * (data_secondary_e11[:, 1] - data_primary_e11[:, 1]),
         linewidth=1, color='indianred', linestyle='dashed')
ax2.plot(1.e+3 * data_primary_e11[:, 0], 1.e-3 * data_secondary_e11[:, 1],
         linewidth=1, color='indianred', linestyle='solid')

# Строим подогнанный график в собственных осях для ax2
x = np.logspace(-2, 3.25, 1000)
ax3 = ax2.twinx()
ax4 = ax2.twiny()
ax3.set_yscale('log')
ax4.set_xscale('log')
ax3.set_ylim(1.e+15, 1.e+27)
ax4.set_xlim(1.7e+2, 3.e+6)
ax4.plot(x, x ** 6, linewidth=2, color='yellow', linestyle='solid')
ax4.plot(x, x ** 6, linewidth=1, color='indianred', linestyle='solid')
for ax in (ax3, ax4):
    format_axes(ax)
for ax in (ax2, ax3):
    ax.tick_params(axis='both', which='major', top=True, right=True, labelright=False, labeltop=False)

fig_name = figure_folder + r"\inst_spec_e15_e13_e11.pdf"
plt.savefig(fig_name)
f.clf()

# ### M dN/dM ### #
f, ax = plt.subplots(1, 1, figsize=(fig_width * 0.65, fig_width * 0.5))
f.subplots_adjust(left=0.1, right=0.95)
ax.set_xscale('log')
# ax.set_yscale('log')
ax.set_xlabel(r'$M,\,$ г')
ax.set_ylabel(r'$\mathrm{d}N/\mathrm{d}M,\:\text{г}^{-1}\:\text{см}^{-3}$')
# ax.set_ylim(1.e+15, 1.e+27)
ax.set_xlim(1.e+7, 1.e+11)

ms = np.logspace(6, 12, 1000)
ax.plot(ms, dNdM_log(1.e+9, 0.5, 1., ms), color='indigo', label=r'$\sigma = 0.5$', linewidth=3)
ax.plot(ms, dNdM_log(1.e+9, 1.0, 1., ms), color='olive', label=r'$\sigma = 1.0$', linewidth=3)
ax.plot(ms, dNdM_log(1.e+9, 2.0, 1., ms), color='teal', label=r'$\sigma = 2.0$', linewidth=3)

ax.legend(loc='upper right')

fig_name = figure_folder + r"\logspec.pdf"
plt.savefig(fig_name)

f.clf()


