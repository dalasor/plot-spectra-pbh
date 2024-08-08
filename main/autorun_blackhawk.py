import os, fileinput, sys, time
import numpy as np
import subprocess

time_ini = time.time()
path = os.getcwd()
print(path)


def replaceAll(file, searchExp, replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = replaceExp
        sys.stdout.write(line)


Mpbhs = [1.e8, 1.e9, 1.e10, 1.e11, 1.e12, 1.e13, 1.e14]

# Choose mass spectrum: 0 = monochromatic, 11 = lognormal
spectrum_choice = 11

# For lognormal:
# sigmas = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
sigmas = [3.5, 4.0]
BHnumber = 50  # number of masses

param_file = path + '/parameters.txt'

count = 1
num = 14  # 42

# --- MAIN --- #

for M in Mpbhs:
    for sigma in sigmas:

        Emin = 1.e-15
        Emax = 1.e+15
        Enumber = 500

        replaceAll(param_file, "destination_folder =",
                   "destination_folder = {:.1e}_{:.1f} \t\t\t\t\t\t\t\t\t # name of the output folder in results/ \n".format(M, sigma))
        replaceAll(param_file, "spectrum_choice = ",
                   "spectrum_choice = {:} \t\t\t\t\t\t\t\t\t # form of the BH distribution: 0=Dirac, 11=log-normal for the number, 2=power-law, 3=critical collapse, 4=peak theory, 5=uniform -1=user-defined \n".format(
                       spectrum_choice))

        if spectrum_choice == 0:
            replaceAll(param_file, "Mmin = ",
                       "Mmin = {:.1e} \t\t\t\t\t\t\t\t\t # lowest BH mass in g (larger than the Planck mass) \n".format(M))
            replaceAll(param_file, "BH_number = ",
                       "BH_number = {:} \t\t\t\t\t\t\t\t\t # number of BH masses (should be the number of tabulated masses if spectrum_choice=5) \n".format(
                           1))

        if spectrum_choice == 11:
            replaceAll(param_file, "BH_number = ",
                       "BH_number = {:} \t\t\t\t\t\t\t\t\t # number of BH masses (should be the number of tabulated masses if spectrum_choice=5) \n".format(
                           BHnumber))
            replaceAll(param_file, "Mmin = ",
                       "Mmin = {:.1e} \t\t\t\t\t\t\t\t\t # lowest BH mass in g (larger than the Planck mass) \n".format(1.e+8))
            replaceAll(param_file, "Mmax = ",
                       "Mmax = {:.1e} \t\t\t\t\t\t\t\t\t # highest BH mass in g (larger than the Planck mass) \n".format(1.e+15))
            replaceAll(param_file, "stand_dev_lognormal = ",
                       "stand_dev_lognormal = {:.1f} \t\t\t\t\t\t\t\t\t # dimensionless variance of the log-normal distribution \n".format(
                           sigma))
            replaceAll(param_file, "crit_mass_lognormal = ",
                       "crit_mass_lognormal = {:.1e} \t\t\t\t\t\t\t\t\t # # characteristic mass of the log-normal distribution in g \n".format(
                           M))

        replaceAll(param_file, "Emin = ",
                   "Emin = {:.1e} \t\t\t\t\t\t\t\t\t # minimal energy in GeV of the primary particles \n".format(Emin))
        replaceAll(param_file, "Emax = ",
                   "Emax = {:.1e} \t\t\t\t\t\t\t\t\t # maximal energy in GeV of the primary particles \n".format(Emax))
        replaceAll(param_file, "E_number = ",
                   "E_number = {:} \t\t\t\t\t\t\t\t\t # number of primary particles energies to be simulated \n".format(
                       Enumber))

        # Строка ответов, содержащая 42 "y", разделённых символами новой строки
        input_data = ("y\n" * 42).encode()

        # Запуск программы C с автоматическим ответом "y" 42 раза на входные запросы
        process = subprocess.Popen(["./BlackHawk_tot.x", "parameters.txt"], stdin=subprocess.PIPE)

        # Передача всех ответов в stdin
        process.stdin.write(input_data)
        process.stdin.flush()

        # Закрытие stdin и ожидание завершения процесса
        process.stdin.close()
        process.wait()

        print("Done {}/{}".format(count, num))
        count += 1

print("Minutes elapsed:", (time.time() - time_ini) / 60.)
