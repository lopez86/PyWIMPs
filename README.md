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

It is very likely that the following will eventually be needed:
 
 * SciPy - more numerical tools
 * Matplotlib - Matlab-like plotting in python
 * Basemap - extra map plotting tools for Matplotlib


## Features

 * Standard dark matter-nucleus interaction model:
   * Standard Halo Model: Truncated Maxwellian velocity distribution
   * Isotropic cross section
   * Several basic form factors
 * Monte Carlo simulation of recoils given some basic properties
   * Weighted sampling for building histograms and distributions, calculating weights, etc.
   * Un-weighted event-by-event sampling

## Future Features
 * Nucleus to nucleon normalization
 * Data for common nuclei
 * Implementation of some basic limit-setting
   * Background-free counting
   * Feldman-Cousins procedure
   * Maximum Gap (Yellen)
   * Annual modulation?
 * Detector effects
   * Efficiency curves
   * Reconstruction effects (smearing of reconstructed properties)
 * Plotting tools
   * Recoil distribution skymaps
   * Sidereal modulation skymaps
   * Recoil energy distributions
   * Annual modulation rate
 * References and readings on dark matter
 * Maybe/Might be fun
   * More sampling methods (MCMC?)
   * Simplified parameterized simulation of a LUX or XENON type detector
   * Inelastic dark matter
   * Q^2-dependent cross sections
   * Coherent neutrino elastic scattering


