from preamble import *


def poly_expression(x, y, deg=10):
    log_x, log_y = np.log10(x), np.log10(y)  # Converting source data to linear scale
    coeffs = np.polyfit(log_x, log_y, deg=deg)  # Polynomial approximation on a linear scale
    return np.poly1d(coeffs)


def beta_omega_constraint(m):
    return 0.282 * 7.06e-18 * (m / 1.e+15) ** 0.5


# Define the function using your piecewise conditions
def beta_shtrix(m_pbh):
    if m_pbh < bbn_x[29]:  # 8.977355e+09
        return beta_omega_constraint(m_pbh)
    elif bbn_x[29] <= m_pbh <= bbn_x[-1]:
        return 10 ** poly_bbn(np.log10(m_pbh))
    elif bbn_x[-1] < m_pbh <= egb_x[20]:
        return 10 ** poly_cmb(np.log10(m_pbh))
    elif egb_x[20] < m_pbh <= egb_x[31]:
        return 10 ** poly_egb_1(np.log10(m_pbh))
    elif m_pbh > egb_x[31]:
        return 10 ** poly_egb_2(np.log10(m_pbh))


bbn_x, bbn_y = np.loadtxt(data_here + r'\beta_bbn_data.csv', delimiter=' ', unpack=True)
cmb_x, cmb_y = np.loadtxt(data_here + r'\beta_cmb_data.csv', delimiter=' ', unpack=True)
egb_x, egb_y = np.loadtxt(data_here + r'\egb_data.csv', delimiter=' ', unpack=True)

bbn_new_x = np.linspace(min(np.log10(bbn_x[20:])), max(np.log10(bbn_x[20:])), 1000)
poly_bbn = poly_expression(bbn_x[20:], bbn_y[20:], deg=36)

cmb_new_x = np.linspace(min(np.log10(cmb_x[:31])), max(np.log10(cmb_x[:31])), 1000)
poly_cmb = poly_expression(cmb_x, cmb_y, deg=36)

egb_new_x = np.linspace(min(np.log10(egb_x[20:])), max(np.log10(egb_x[20:])), 1000)

egb_new_x_1 = np.linspace(min(np.log10(egb_x[20:32])), max(np.log10(egb_x[20:32])), 1000)
poly_egb_1 = poly_expression(egb_x[20:32], egb_y[20:32], deg=36)

egb_new_x_2 = np.linspace(min(np.log10(egb_x[31:])), max(np.log10(egb_x[31:])), 1000)
poly_egb_2 = poly_expression(egb_x[31:], egb_y[31:], deg=36)

if __name__ == '__main__':
    f, ax = plt.subplots(1, 1, figsize=(fig_width * 0.85, fig_width * 0.65))
    f.subplots_adjust(left=0.1, right=0.95)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r"$M$, г")
    ax.set_ylabel(r"$\beta'$")
    ax.set_ylim(1.e-33, 1.e-17)
    ax.set_xlim(1.e+9, 1.e+17)
    ax.set_xticks(np.logspace(10, 16, num=7), minor=False)
    ax.set_xticks(np.logspace(9, 17, num=81), labels=[], minor=True)
    ax.set_yticks(np.logspace(-33, -17, num=9), minor=False)
    ax.yaxis.set_minor_locator(plt.FixedLocator(np.logspace(10, 16)))

    ax.fill_between(bbn_x[29:], bbn_y[29:], 1.e-16, color='grey', alpha=0.3)
    ax.fill_between(10 ** cmb_new_x, 10 ** poly_cmb(cmb_new_x), 1.e-16, color='grey', alpha=0.3)
    m_s = np.logspace(9, 17, 1000)
    mask = m_s <= bbn_x[29]
    ax.fill_between(m_s, beta_omega_constraint(m_s), 1.e-16, where=mask, color='grey', alpha=0.3)
    ax.fill_between(egb_x[19:], egb_y[19:], 1.e-16, color='grey', alpha=0.3)

    ax.plot(bbn_x, bbn_y, label='BBN', linewidth=4, color='magenta')
    ax.plot(10 ** bbn_new_x, 10 ** poly_bbn(bbn_new_x), linewidth=1, color='black')
    ax.axvline(x=bbn_x[29], color='black', linewidth=1, linestyle='dashed')

    m_s = np.linspace(1.e+9, 1.e+17, 100)
    ax.plot(m_s, beta_omega_constraint(m_s), linewidth=2, color='darkblue', linestyle='dashed',
            label=r"$\Omega_{\mathrm{PBH}}$")
    # ax.scatter(1.e+10, 10 ** polynomial(np.log10(1.e+10)), marker='*', s=100, edgecolors='r')

    ax.plot(cmb_x[:32], cmb_y[:32], label='CMB', linewidth=4, color='orange')
    ax.plot(10 ** cmb_new_x, 10 ** poly_cmb(cmb_new_x), linewidth=1, color='black')
    ax.axvline(x=cmb_x[0], color='black', linewidth=1, linestyle='dashed')

    ax.plot(egb_x[19:], egb_y[19:], label='EGB', linewidth=4, color='red')
    ax.plot(10 ** egb_new_x_1, 10 ** poly_egb_1(egb_new_x_1), linewidth=1, color='black')
    ax.plot(10 ** egb_new_x_2, 10 ** poly_egb_2(egb_new_x_2), linewidth=1, color='black')
    ax.axvline(x=egb_x[20], color='black', linewidth=1, linestyle='dashed')

    ax.legend(loc='upper right', framealpha=1.0)
    fig_name = figure_folder + r"\constraints.pdf"
    plt.savefig(fig_name)
    f.clf()
