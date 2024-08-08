# Plotting time-integrated neutrino spectra from Primordial Black Holes (PBH) with log-normal mass distribution

## Overview

To simulate the life of evaporating primordial black holes and numerically calculate energy spectra, the public code [BlackHawk v.2.1](https://blackhawk.hepforge.org/) is used. 

This code on **Python 3.10** can be used as a template to plot your own plots using **matplotlib 3.8.4**.

## Data

The "Data" folder contains files that are used in plotting which are not from **BlackHawk**. They can be divided into 3 groups:

- [Constraints](http://dx.doi.org/10.1088/1361-6633/ac1e31) on abundance of PBH - dimensionless parameter $\beta'$ - the fraction of the PBH energy density in the total energy density of the Universe at the time of their formation:
  1. `beta_bbn_data.csv` - obtained from Big Bang Nucleosynthesis
  2. `beta_cmb_data.csv` - obtained from Cosmic Microwave Background
  3. `egb_data.csv` - obtained from ExtraGalactic Background

- [$\nu$](http://dx.doi.org/10.1016/j.astropartphys.2020.102537)- and [$\gamma$](http://dx.doi.org/10.1177/0003702818767133)-background data:
  1. `antinu_n.txt` -
  2. `antinu_t.txt` -
  3. `Sun-thermal.dat` -
  4. `Sun-nuclear-pp.dat` -
  5. `Sun-nuclear-hep.dat` -
  6. `Sun-nuclear-B8.dat` -
  7. `dataoldsn.csv` -
  8. `Atmospheric.dat` -
  9. `nu_background.txt` -
  10. `ph_background.txt` -

- Additional / Intermediate
