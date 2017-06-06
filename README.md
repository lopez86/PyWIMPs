# WimpSim

This package contains various python tools for simulation &amp; analysis
of dark matter direct detection experiments.

Everything is still under construction. 
The goal is for the code here to all be in Python 3 with standard Python packages such as NumPy and SciPy. I might make a ROOT interface but the core code should be independent of ROOT. This ensures that the code will be easy to install and use by non-experts.

## Requirements

So far, the following dependencies are needed:

 * Python 3 (likely 3.4 or above) - primary language
 * NumPy - numerical calculations
 * AstroPy - astrophysics libraries (coordinate transforms)

For the examples, you will also need (depending on the specific example):
 * Matplotlib
 * PyROOT (ROOT with Python bindings)

It is very likely that the following will eventually be needed:
 
 * SciPy - more numerical tools
 * Basemap - extra map plotting tools for Matplotlib


## Features

 * Standard dark matter-nucleus interaction model:
   * Standard Halo Model: Truncated Maxwellian velocity distribution
   * Isotropic cross section
   * Several basic form factors
 * Monte Carlo simulation of recoils using the standard halo and cross section assumptions
   * Weighted sampling for building histograms and distributions, calculating weights, etc. (one throwing uniformly over a region and another drawing from a Maxwell-Boltzmann distribution)
   * Un-weighted event-by-event sampling using (1) a basic rejection sampling method and (2) a Markov Chain Monte Carlo (MCMC) using the Metropolis-Hastings algorithm.

## Future Features
 * Nucleus to nucleon normalization
 * Data for common nuclei
 * Implementation of some basic limit-setting
   * Background-free counting
   * Feldman-Cousins procedure
   * Maximum Gap (Yellen)
   * Annual modulation limits
 * Detector effects
   * Efficiency curves
   * Reconstruction effects (smearing of reconstructed properties)
 * Examples of various plots and calculations
   * Recoil distribution skymaps
   * Sidereal modulation skymaps
   * Recoil energy distributions
   * Annual modulation rate
 * References and readings on dark matter
 * Maybe/Might be fun
   * Simplified parameterized simulation of a LUX or XENON type detector
   * Inelastic dark matter
   * Q^2-dependent cross sections
   * Coherent neutrino elastic scattering


