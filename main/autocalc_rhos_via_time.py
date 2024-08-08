import os, time
import fileinput, sys
import subprocess
from scipy import integrate
from check_beta import *

time_ini = time.time()
path = os.getcwd()

if not os.path.exists("DataRhos"):
    os.system("mkdir DataRhos")

dct = {'./stack1.x': 'pri',
       './stack2.x': 'e_sec',
       './stack3.x': 'mu_sec',
       './stack4.x': 'tau_sec'}

c_files = {'./stack1.x': './stack1.c',
           './stack2.x': './stack2.c',
           './stack3.x': './stack3.c',
           './stack4.x': './stack4.c'}

cosmodata = r"D:\Downloads\BlackHawk\scripts\cosmology_scripts\DataRhos"
c = 2.99792458e+10  # [cm / c]
mev_to_g = 1 / 5.60959e+26
rho_nu_rel = 3.16e-34  # [г * см^-3]
rho_ph_rel = 7.8e-34  # [г * см^-3]

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


def data_extract_ph(path_pri: str, path_sec: str):
    data_ph_pri = (np.genfromtxt(cosmodata + path_pri, skip_header=1))
    data_ph_sec = (np.genfromtxt(cosmodata + path_sec, skip_header=1))
    data_final = np.array(data_ph_pri[:, 1]) + np.array(data_ph_sec[:, 1])
    return np.array(data_ph_pri[:, 0]), data_final


def replaceAll(file, searchExp, replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = replaceExp
        sys.stdout.write(line)


# Функция для запуска одной программы на C и сохранение результатов
def run_c_program(fname: str, fpath: str, m, s, t) -> None:
    replaceAll(c_files[fname], "double t_today =", "double t_today = {:.1e};\n".format(t))
    subprocess.run(['make', '-f', 'Makefile{}'.format(fname[-3])], stdout=subprocess.DEVNULL)
    subprocess.run([fname, fpath, "500", "2", "1", "0"], stdout=subprocess.DEVNULL)
    os.system("cp ./ra_{}.txt ./DataRhos/ra_{:.1e}_{:.1f}_{}_{:.1e}_new.txt".format(dct[fname], m, s, dct[fname], t))


# Функция для запуска 4 программ на C в цикле и нахождение плотности энергии для текущего распределения
def get_all_files_distribution(m: float, s: float, t: float) -> None:
    fpath = "/cygdrive/d/Downloads/BlackHawk/results/{:.1e}_{:.1f}_new/".format(m, s)  # ????
    print('######## Nuetrinos: M = {:.1e}, sigma = {:.1f}, t = {:.1e} #######'.format(m, s, t))
    run_c_program("./stack1.x", fpath, m, s, t)
    print('Done 1 - pri')
    run_c_program("./stack2.x", fpath, m, s, t)
    print('Done 2 - e_sec')
    run_c_program("./stack3.x", fpath, m, s, t)
    print('Done 3 - mu_sec')
    run_c_program("./stack4.x", fpath, m, s, t)
    print('Done 4 - tau_sec')

    energies, nu_all = data_extract_nu(r"\ra_{:.1e}_{:.1f}_pri_{:.1e}_new.txt".format(m, s, t),
                                       r"\ra_{:.1e}_{:.1f}_e_sec_{:.1e}_new.txt".format(m, s, t),
                                       r"\ra_{:.1e}_{:.1f}_mu_sec_{:.1e}_new.txt".format(m, s, t),
                                       r"\ra_{:.1e}_{:.1f}_tau_sec_{:.1e}_new.txt".format(m, s, t))

    E_nu = energies * 1.e+3
    n_PBH = n_pbh_beta(beta_shtrix(m), m, s)
    y_calc = n_PBH * (c / (4 * np.pi)) * E_nu * nu_all / 1.e+3
    area_calc = integrate.trapz(y_calc, E_nu)
    rho_nu_pbh = area_calc * mev_to_g / c  # г * см^-3

    rhos.append(rho_nu_pbh)

    print("For {:.1e}_{:.1f}_{:.1e}:".format(m, s, t), rho_nu_pbh)

    with open('rhos_2.0.txt', 'a') as file:
        file.write(r'{:.2e}\t{:.2e}\n'.format(t, rho_nu_pbh))


# def get_ph_files_distribution(m: float, s: float) -> None:
#     fpath = "/cygdrive/d/Downloads/BlackHawk/results/{:.1e}_{:.1f}/".format(m, s)  # ????
#     print('######## Photons: M = {:.1e}, sigma = {:.1f} #######'.format(m, s))
#     run_c_program("./stack5.x", fpath, m, s)
#     print('Done 1 - ph_pri')
#     run_c_program("./stack6.x", fpath, m, s)
#     print('Done 2 - ph_sec')
#
#     energies, ph_all = data_extract_ph(r"\ra_{:.1e}_{:.1f}_ph_pri.txt".format(m, s),
#                                        r"\ra_{:.1e}_{:.1f}_ph_sec.txt".format(m, s))
#
#     E_ph = energies * 1.e+3
#     n_PBH = n_pbh_beta(beta_shtrix(m), m, s)
#     y_calc = n_PBH * (c / (4 * np.pi)) * E_ph * ph_all / 1.e+3
#     area_calc = integrate.trapz(y_calc, E_ph)
#     rho_ph_pbh = area_calc * mev_to_g / c  # г * см^-3
#     rho_ratio = rho_ph_pbh / rho_ph_rel
#
#     print("For ph {:.1e}_{:.1f}:".format(m, s), rho_ratio)


if __name__ == "__main__":

    num = 100  # 42

    ts = np.logspace(0, 17.638, num)

    Mpbhs = [1.0e+11]
    sigmas = [0.0]
    count = 1

    for M in Mpbhs:
        for sigma in sigmas:
            for t in ts:
                get_all_files_distribution(M, sigma, t)  # nu
                # get_ph_files_distribution(M, sigma, t)  # ph
                print("Done {}/{}".format(count, num))
                print()
                count += 1

    print("Minutes elapsed:", (time.time() - time_ini) / 60.)
