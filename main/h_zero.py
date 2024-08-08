# from preamble import *
from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt

data_here = r'C:\Users\Luffy\PycharmProjects\visual_script\data'
figure_folder = r"C:\Article_PBH\Figures"


def rho_crit(h: float) -> float:
    return 1.878 * 1.e-29 * h ** 2  # г * cm^-3


def H(z, h: float, is_pbh: bool) -> float:
    if is_pbh:
        return ((8 * np.pi * GG) * (omega_m * rho_crit(h) * z ** 3 + rho_l + rho_r * z ** 4 +
                                    rho_nu_pbh * 2) / 3) ** (1 / 2)
    else:
        return ((8 * np.pi * GG) * (omega_m * rho_crit(h) * z ** 3 + rho_l + rho_r * z ** 4) / 3) ** (1 / 2)


def t(x):
    H0 = H(0.0, h0, False)
    return 1 / (H0 * x * (omega_m * x ** (-3) + omega_l + omega_r * x ** (-4)) ** (1 / 2))


# data for H(z) formula
GG = 6.67e-8  # [cm^3 * g^-1 * s^-2]
h0_min = 67.40
h0 = 0.674
h0_max = 74.03
zs = np.logspace(0, 10, 1000)  # (1+z)

g_to_mev = 5.60959e+26
mev_to_g = 1 / 5.60959e+26
rho_r = 7.8e-34  # г * см^-3
rho_l = 60.3e-31  # г * см^-3
rho_nu_rel = 3.167e-34  # г * см^-3

# For M_c = 2.4e+9 and sigma = 2.8
max_relation = 2.34e-2
# rho_nu_pbh = max_relation * rho_nu_rel  # г * см^-3

omega_m = 0.315
omega_l = rho_l / rho_crit(h0)
omega_r = rho_r / rho_crit(h0)

ts = np.logspace(1, 17, 100)
data = np.genfromtxt(data_here + r'\rhos_2.0.txt')
rhos = data[:, 1]
rho_nu_pbh = np.average(rhos)

H0 = H(10 ** 7, h0, False)
print('H0_std =', H0)
H0_pbh = H(10 ** 7, h0, True)
print('H0_pbh =', H0_pbh)
print((H0_pbh - H0) / H0)


z_max = 1.e+5
a, b = 0, 1 / (1 + z_max)
result = integrate.quad(t, a, b)
t_today = 13.8e+9 * 365. * 24. * 3600.
# print('t_today in sec:', "{:e}".format(t_today))
print("t({:e}) = {:e}".format(z_max, result[0]))
print()


if __name__ == '__main__':
    f, ax = plt.subplots(1, 1, figsize=(6, 4))
    f.subplots_adjust(left=0.1, right=0.95)
    ax.set_xlabel(r'$\log_{10}(1+z)$', fontsize=13)
    ax.set_ylabel(r'$\delta H(z)$', fontsize=13)

    ax.set_ylim(-0.3e-5, 2e-5)
    ax.set_xlim(0, 2)

    H_s = H(zs, h0, False)
    H_p = H(zs, h0, True)
    y_p = (H_p - H_s) / H_s
    y_0 = (H_s - H_s) / H_s

    zs_new = np.log10(zs)

    ax.plot(zs_new, y_p, linewidth=3, color='royalblue', linestyle='solid', label=r'$H(z)$')
    ax.plot(zs_new, y_0, linewidth=3, color='black', linestyle='dashed', label=r'$H_s(z)$')

    ax.axvline(x=np.log10(0.4), color='grey', linewidth=2, linestyle='dashed')
    ax.text(1.25, 0.05, r'$z_* = 1090$', fontsize=18, rotation=90)

    plt.gca().invert_xaxis()

    ax.legend(loc='upper right')
    fig_name = figure_folder + r"\h_results.pdf"
    plt.show()

    f.clf()

    # ###### RHOS ###### #
    f, ax = plt.subplots(1, 1, figsize=(10, 6))
    f.subplots_adjust(left=0.15, right=0.95, bottom=0.2)
    ax.set_ylim(1.e-36, 1e-33)
    # ax.set_xlim(0, 2)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$t, $ с', fontsize=13)
    ax.set_ylabel(r'$\rho_{\nu}^{PBH}, $ г см$^{-3}$', fontsize=13)
    ax.tick_params(axis='x', labelsize=13)  # Размер подписей к тикам на оси X
    ax.tick_params(axis='y', labelsize=13)

    # Пример данных
    x = np.array([1, 100, 1e5, 1e10, 1e15])  # X данные
    y = np.array([1e-36, 1e-35, 1e-34, 5e-34, 1e-33])  # Y данные

    # Логарифмическое преобразование данных
    log_x = np.log10(ts[1:])
    log_y = np.log10(rhos[1:])

    # Аппроксимация логарифмированных данных полиномом 1-й степени (линейная регрессия)
    coefficients = np.polyfit(log_x, log_y, 3)
    polynomial = np.poly1d(coefficients)

    # Генерация аппроксимационных данных
    log_x_fit = np.linspace(min(log_x), max(log_x), 500)
    log_y_fit = polynomial(log_x_fit)

    ax.scatter(ts[1:], rhos[1:], marker='x', s=30, color='olive', linestyle='solid')
    ax.plot(10 ** log_x_fit, 10 ** log_y_fit, color='indigo', linewidth=3)

    # fig_name = figure_folder + r"\rhos_results.pdf"
    plt.show()

    f.clf()

    # #### PLot N_eff #### #
    f, ax = plt.subplots(1, 1, figsize=(10, 6))
    f.subplots_adjust(left=0.15, right=0.95, bottom=0.2)
    ax.set_ylim(2, 5)
    ax.set_xlim(0, 3)
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    ax.set_xlabel(r'$\log_{10}(1+z)$', fontsize=13)
    ax.set_ylabel(r'$N_{\mathrm{eff}}$', fontsize=13)
    ax.tick_params(axis='x', labelsize=13)  # Размер подписей к тикам на оси X
    ax.tick_params(axis='y', labelsize=13)

    rho_gamma = 4.64e-34
    fac = 7/8 * (4/11) ** (4/3)


    def Neff(z, ispbh: bool):
        if ispbh:
            return (rho_nu_rel * z ** 4 + rho_nu_pbh) / (fac * rho_gamma * z ** 4)
        else:
            return (rho_nu_rel * z ** 4) / (fac * rho_gamma * z ** 4)


    ax.plot(zs_new, Neff(zs, True), color='royalblue', linewidth=3, label=r'$N_{\mathrm{eff}}(z)$')
    ax.plot(zs_new, Neff(zs, False), color='gray', linestyle='dashed', linewidth=3, label=r'$N^s_{\mathrm{eff}}(z)$')
    plt.gca().invert_xaxis()
    ax.legend(loc='upper right')
    plt.show()

    f.clf()






