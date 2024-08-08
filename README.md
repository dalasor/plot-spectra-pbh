# Plotting time-integrated neutrino spectra from Primordial Black Holes (PBH) with log-normal mass distribution

## Overview

To simulate the life of evaporating primordial black holes and numerically calculate energy spectra, the public code [BlackHawk v.2.1](https://blackhawk.hepforge.org/) is used. 

This code on **Python 3.10** can be used as a template to plot your own plots using **matplotlib 3.8.4**.

## Data

The "Data" folder contains files that are used in plotting. They can be divided into 3 groups:

- Constraints on abundance of PBH - dimensionless parameter $\beta'$ - the fraction of the PBH energy density in the total energy density of the Universe at the time of their formation:
  1. `beta_bbn_data.csv` - obtained from Big Bang Nucleosynthesis
  2. `beta_cmb_data.csv` - obtained from Cosmic Microwave Background
  3. `egb_data.csv` - obtained from ExtraGalactic Background

- Background $\nu$ and $\gamma$ data
