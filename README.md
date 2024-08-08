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

- $\nu$-background ([link](http://dx.doi.org/10.1016/j.astropartphys.2020.102537)) and $\gamma$-background ([link](http://dx.doi.org/10.1177/0003702818767133)) data:
  1. `antinu_n.txt` and `antinu_t.txt` - antineutrinos from the decay of neutrons and tritons during primary nucleosynthesis
  2. `Sun-thermal.dat` - predicted thermal solar neutrino
  3. `Sun-nuclear-pp.dat`, `Sun-nuclear-hep.dat` and `Sun-nuclear-B8.dat` - neutrinno from nuclear fusion reactions on Sun
  4. `dataoldsn.csv` - diffuse neutrino background from distant supernovae (DSBN)
  5. `Atmospheric.dat` - atmospheric neutrino 
  6. `nu_background.txt` and `ph_background.txt` - combined full background for neutrinos and photons

- Additional / Intermediate

## Code
