import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm


def poly_expression(xs, ys, deg=10):
    log_x, log_y = np.log10(xs), np.log10(ys)  # Преобразование исходных данных в линейную шкалу
    coeffs = np.polyfit(log_x, log_y, deg=deg)  # Аппроксимация полиномом на линейной шкале
    return np.poly1d(coeffs)


def beta_omega_constraint(m):
    return 0.282 * 7.06e-18 * (m / 1.e+15) ** 0.5


# Define the function using your piecewise conditions
def beta_shtrix(m_pbh):
    if m_pbh < bbn_x[29]:  # 8.977355e+09
        return beta_omega_constraint(m_pbh)
    elif bbn_x[29] <= m_pbh <= bbn_x[-1]:
        return 10 ** poly_bbn(np.log10(m_pbh))
    elif m_pbh > bbn_x[-1]:
        return 10 ** poly_cmb(np.log10(m_pbh))


# comoving number density today
def n_pbh_beta(m_pbh, sigma):
    return 1.19e-27 * beta_shtrix(m_pbh) * (1.e+15 / m_pbh) ** (3 / 2) * np.exp(-9 * (sigma ** 2) / 8)


data_here = r'C:\Users\Luffy\PycharmProjects\visual_script\data'

bbn_x, bbn_y = np.loadtxt(data_here + r'\beta_bbn_data.csv', delimiter=' ', unpack=True)
cmb_x, cmb_y = np.loadtxt(data_here + r'\beta_cmb_data.csv', delimiter=' ', unpack=True)

bbn_new_x = np.linspace(min(np.log10(bbn_x[20:])), max(np.log10(bbn_x[20:])), 1000)
poly_bbn = poly_expression(bbn_x[20:], bbn_y[20:], deg=36)

cmb_new_x = np.linspace(min(np.log10(cmb_x)), max(np.log10(cmb_x)), 1000)
poly_cmb = poly_expression(cmb_x, cmb_y, deg=36)

f, ax = plt.subplots(1, 1)
f.subplots_adjust(left=0.1, right=0.95)
ax.set_xscale('linear')
ax.set_yscale('linear')
# ax.set_xlabel(r'$M\,\,[$г$]$')
ax.set_xlabel(r'$\log_{10}(M\,/\,$г$)$')
ax.set_ylabel(r'$\sigma$')
# ax.set_xlim([1.e+7, 1.e+12])
ax.set_xlim([7.0, 11.5])
ax.set_ylim([0.0, 4.0])

# ####################################################### #
# ### Строи общую Z и закрасшиваем области под маской ### # (imshow)
# ####################################################### #

bins = 1000

# Возьмем первую кривую (M, sigma), ограничивающую сверху и проаппроксимируем
x_all, y_all = np.loadtxt(data_here + r'\msigma_best.csv', delimiter=' ', unpack=True)
x_all_lin = np.log10(x_all)
cfs = np.polyfit(x_all_lin, y_all, 9)
pol = np.poly1d(cfs)

x_lin = np.linspace(min(x_all_lin), max(x_all_lin), bins)
y_pol = pol(x_lin)

x_all_log = 10 ** x_lin
ax.plot(x_lin, y_pol, color='blue')  ###

# Задаем двумерные сетки на одномерных осях
x_log = np.logspace(np.log10(min(x_all)), np.log10(max(x_all)), bins)
y_lin = np.linspace(0.0, 4.0, bins)
X_all, Y_all = np.meshgrid(x_log, y_lin)  # двумерные сетки

# Вторая кривая (M, sigma), задающая ограничение снизу
x_line, y_line = np.loadtxt(data_here + r'\msigma_lineup.csv', delimiter=' ', unpack=True)
cfs_line = np.polyfit(x_line, y_line, 4)
pol_line = np.poly1d(cfs_line)
x_line = np.linspace(min(x_line), max(x_line), bins)
y_line = pol_line(x_line)

# ### ПРОБЛЕМА С ЛОГАРИФМИЧЕСКИМ МАСШТАБОМ ДЛЯ ВТОРОЙ ЛИНИИ (ДЕЛАЮ ОБРЕЗАНИЕ)

ax.plot(x_line, y_line, color='blue')

# Вычисляем цвет
Z = np.zeros((bins, bins))
for i in range(bins):
    for j in range(bins):
        # z[i][j] = n_pbh_beta(X_all[i][j], Y_all[i][j])
        Z[i][j] = n_pbh_beta(x_log[i], y_lin[j])

Y_poly = np.tile(y_pol, (bins, 1))
Y_line = np.tile(y_line, (bins, 1))
mask = (Y_all > Y_poly) | (Y_all < np.log10(Y_line))  # Маска для значений Z
Z_masked = np.copy(Z)  # Применяем маску: значениям Z, удовлетворяющим маске, назначаем значение `np.nan`
Z_masked[mask] = np.nan  # np.nan значения будут показываться как прозрачные на графике
vmax = np.nanmax(Z_masked)

# Создаем сплошной цветной фон
cmap = plt.get_cmap('Blues_r')
color_mesh = ax.imshow(Z_masked, extent=(x_lin.min(), x_lin.max(), y_lin.min(), y_lin.max()),
                       origin='lower', cmap=cmap, norm=LogNorm())  # vmax=vmax)

f.colorbar(color_mesh, ax=ax, label='$n_{PBH}(\\sigma, M)$')

ax.fill_between(x_line, y_line, min(y_line), color='white', alpha=1)

plt.show()
