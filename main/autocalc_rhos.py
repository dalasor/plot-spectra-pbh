import os, time
import subprocess
from scipy import integrate
from check_beta import *

time_ini = time.time()
path = os.getcwd()
print(path)

if not os.path.exists("CosmologyData"):
    os.system("mkdir CosmologyData")

dct = {'./stack1.x': 'pri',
       './stack2.x': 'e_sec',
       './stack3.x': 'mu_sec',
       './stack4.x': 'tau_sec'}

cosmodata = r"D:\Downloads\BlackHawk\scripts\cosmology_scripts\CosmologyData"
c = 2.99792458e+10  # [cm / c]
mev_to_g = 1 / 5.60959e+26
rho_nu_rel = 3.16e-34  # [г * см^-3]

rhos = []


def n_pbh_beta(beta: float, m_pbh: float, sigma: float) -> float:
    if sigma == 0:
        return 1.19e-27 * beta * (1.e+15 / m_pbh) ** (3 / 2)
    else:
        return 1.19e-27 * beta * (1.e+15 / m_pbh) ** (3 / 2) * np.exp(-9 * sigma ** 2 / 8)


def data_extract_nu(path_pri: str, path_e: str, path_mu: str, path_tau: str):
    data_nu_pri = (np.genfromtxt(cosmodata + path_pri, skip_header=1))
    data_sec_e = (np.genfromtxt(cosmodata + path_e, skip_header=1))
    data_sec_mu = (np.genfromtxt(cosmodata + path_mu, skip_header=1))
    data_sec_tau = (np.genfromtxt(cosmodata + path_tau, skip_header=1))
    data_final = (np.array(data_nu_pri[:, 1]) + np.array(data_sec_e[:, 1]) + np.array(data_sec_mu[:, 1]) +
                  np.array(data_sec_tau[:, 1]))
    return np.array(data_nu_pri[:, 0]), data_final


# Функция для запуска одной программы на C и сохранение результатов
def run_c_program(fname: str, fpath: str, m, s) -> None:
    subprocess.run([fname, fpath, "500", "2", "1", "1"], stdout=subprocess.DEVNULL)
    os.system("cp ./ra_{}.txt ./CosmologyData/ra_{:.1e}_{:.1f}_{}.txt".format(dct[fname], m, s, dct[fname]))


# Функция для запуска 4 программ на C в цикле и нахождение плотности энергии для текущего распределения
def get_all_files_distribution(m: float, s: float) -> None:
    fpath = "/cygdrive/d/Downloads/BlackHawk/results/{:.1e}_{:.1f}/".format(m, s)  # ????
    print('######## M = {:.1e}, sigma = {:.1f} #######'.format(m, s))
    run_c_program("./stack1.x", fpath, m, s)
    print('Done 1 - pri')
    run_c_program("./stack2.x", fpath, m, s)
    print('Done 2 - e_sec')
    run_c_program("./stack3.x", fpath, m, s)
    print('Done 3 - mu_sec')
    run_c_program("./stack4.x", fpath, m, s)
    print('Done 4 - tau_sec')

    energies, nu_all = data_extract_nu(r"\ra_{:.1e}_{:.1f}_pri.txt".format(m, s),
                                       r"\ra_{:.1e}_{:.1f}_e_sec.txt".format(m, s),
                                       r"\ra_{:.1e}_{:.1f}_mu_sec.txt".format(m, s),
                                       r"\ra_{:.1e}_{:.1f}_tau_sec.txt".format(m, s))

    E_nu = energies * 1.e+3
    n_PBH = n_pbh_beta(beta_shtrix(m), m, s)
    y_calc = n_PBH * (c / (4 * np.pi)) * E_nu * nu_all / 1.e+3
    area_calc = integrate.trapz(y_calc, E_nu)
    rho_nu_pbh = area_calc * mev_to_g / c  # г * см^-3
    rho_ratio = rho_nu_pbh / rho_nu_rel

    rhos.append(rho_ratio)

    print("For {:.1e}_{:.1f}:".format(m, s), rho_ratio)

    with open('rhos.txt', 'a') as file:
        file.write('{:.2e}\t{:.2e}\n'.format(rho_nu_pbh, rho_ratio))


if __name__ == "__main__":

    num = 42
    Mpbhs = [1.e8, 1.e9, 1.e10, 1.e11, 1.e12, 1.e13, 1.e14]
    sigmas = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
    count = 1

    for M in Mpbhs:
        for sigma in sigmas:
            get_all_files_distribution(M, sigma)
            print("Done {}/{}".format(count, num))
            print()
            count += 1

    print("Minutes elapsed:", (time.time() - time_ini) / 60.)
