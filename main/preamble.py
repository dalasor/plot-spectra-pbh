#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 13, 2024

@author: Julian

Data received for BlackHawk v. 2.1
"""

# 1 - Necessary importations

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import rcParams
import matplotlib as mpl
import numpy as np

mpl.use("pgf")

fig_width_pt = 1000  # Get this from LaTeX using \showthe\columnwidth
inches_per_pt = 1.0 / 72.27  # Convert pt to inch
golden_mean = (np.sqrt(5) - 1.0) / 2.0  # Aesthetic ratio
fig_width = fig_width_pt * inches_per_pt  # width in inches
fig_height = fig_width * golden_mean  # height in inches
fig_size = [fig_width, fig_height]

params = {'backend': 'ps',
          'axes.labelsize': 20,
          'axes.linewidth': 1,
          'font.size': 15,
          'font.family': 'serif',
          "font.weight": "ultralight",
          'ytick.direction': 'in',
          'xtick.direction': 'in',
          'figure.titlesize': 25,
          'legend.fontsize': 18,
          'xtick.labelsize': 20,
          'ytick.labelsize': 20,
          'xtick.major.pad': 5,
          'text.usetex': True,
          'figure.figsize': fig_size,
          "pgf.preamble": "\n".join([
              r"\usepackage{url}",
              r"\usepackage{unicode-math}",
              r"\usepackage{amsmath}",
              r"\usepackage[russian]{babel}"
          ])}
rcParams.update(params)

# 2 - Folder definition

# Here put the path
figure_folder = r"C:\Article_PBH\Figures"  # куда хочешь сохранять графики
data_folder = r"D:\Downloads\BlackHawk_v2.1\version_finale\results"  # используется в "plot_distributions"
data_folder_v2 = r"C:\cygwin64\home\BlackHawk_v2.2\results"  # не нужен
data_cosmology = r"D:\Downloads\BlackHawk_v2.1\version_finale\scripts\cosmology_scripts"  # местонажождение stack-нутых спектров № 1
data_cosmology_new = r"C:\cygwin64\home\BlackHawk\scripts\cosmology_scripts"  # местонажождение stack-нутых спектров № 2
data_here = r'C:\Users\Luffy\PycharmProjects\visual_script\data'  # местонахождение данных для построения известных кривых
