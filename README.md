# Plotting time-integrated neutrino spectra from Primordial Black Holes (PBH) with log-normal mass distribution

## Overview

To simulate the life of evaporating primordial black holes and numerically calculate energy spectra, the public code [BlackHawk v.2.1](https://blackhawk.hepforge.org/) is used. 

This code on **Python 3.10** can be used as a template to plot your own plots using **matplotlib 3.8.4**.

Goal: 
- calculate and visualize the neutrino flux densities generated by PBH throughout their existence in the Universe, taking into account strict constraints on the abundance of PBH and cosmological redshift $z$ for various parameters;
- show an assessment of the influence of new taken into account neutrinos on the evolution of the universe through the calculation of energy densities $\rho$ and the relative deviation of the Hubble parameter $\delta H(z)$.

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

- Additional / Intermediate / Optional

The rest of the input data for directly calculating the flux density that appears in the code is obtained as a result of BlackHawk simulations.

## Code

Here is a short description of the scripts included in the "main" where the corresponding calculations are defined:

- `autorun_blackhawk.py` - Automatically runs long-running BlackHawk simulations for various distributions by modifying the parameter file and controlling I/O flow;
- `preamble.py` - Pre-settings for final graphs using _rcParams_ such as sizes, latex enabled, etc.;
- `check_beta.py` - Describes a function that glues piecewise conditions of constraints used in other scripts, builds a plot of existing constraints (`constraints.pdf`);
- `plot_distributions.py` - Comparison instantaneous $\nu-$ and $\gamma-$spectra for 3 mass distributions (`three_mass_spectrum.pdf`, `nu_inst_spec.pdf`, `inst_spec_e15_e13_e11.pdf`); 
- `mass_via_time.py` - Builds a graph of the dependence of mass PBH on the their lifetime (`mass_via_time.pdf`);
- `plot_results.py` - Main results. Builds a graph of time-integrated neutrino flux densities for various distributions. Flux densities calculate according to:

$$ \frac{\mathrm{d}\Phi_{\nu}}{\mathrm{d}E_{\nu}} = \frac{c}{4\pi} \frac{\mathrm{d}N_{\nu}}{\mathrm{d}E_{\nu}} = n_{\mathrm{PBH}}(t_{0}) \int_{t_{\mathrm{min}}}^{t_{0}} \mathrm{d}t (1+z) \frac{\mathrm{d}^{2}N_{\nu}^{\mathrm{tot}}(E_{\nu}(1+z), M(t))}{\mathrm{d}t\mathrm{d}E_{\nu}} ,$$

где $t_{\mathrm{min}} = 1$ s -- time of neutrino splitting off from matter in early Universe, $z(t)$ -- redshift, $n_{\mathrm{PBH}}(t_{0})$ -- current PBH number density with $t_0 = 13.8 \text{billion years}=4.35\times 10^{17}$ s, $c$ - speed of light.

- `allowed_window.py` - 
- `h_zero.py` -
- `rhos_3d.py` -
- `autocalc_rhos.py` -
- `autocalc_rhos_via_time` -

## Some figures

- From `plot_results.py`:
  
  <img src="/Figures/jpg/nu_results_1_cropped.jpg" width="700">
  
- From `plot_distributions.py`:
  
  <img src="/Figures/jpg/inst_spec_e15_e13_e11_cropped.jpg" width="700">

- From `allowed_window.py`:

  <img src="/Figures/jpg/window_cropped.jpg" width="700">

- From `rhos_3d.py`:

  <img src="/Figures/jpg/2dplot_best_cropped.jpg" width="700">
