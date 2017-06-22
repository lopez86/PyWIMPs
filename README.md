# PyWIMPs

This package contains various python tools for simulation &amp; analysis
of dark matter direct detection experiments. If you're a member of the 
HEP or astro community and might want to contribute something, let me know.

## Requirements

So far, the following dependencies are needed:

 * Python 3 (likely 3.4 or above) - primary language
 * NumPy - numerical calculations
 * AstroPy - astrophysics libraries (coordinate transforms)
 * SciPy - used for various statistics things

For the examples, you will also need (depending on the specific example):
 * Matplotlib
 * PyROOT (ROOT with Python bindings)
 * Basemap - extra map plotting tools for Matplotlib

## Features

 * Standard dark matter-nucleus interaction model:
   * Standard Halo Model: Truncated Maxwellian velocity distribution
   * Isotropic cross section
   * Various form factors
   * Nucleus to nucleon normalization
 * Monte Carlo simulation of recoils using the standard halo and cross section assumptions
   * Weighted sampling for building histograms and distributions, calculating weights, etc. (one throwing uniformly over a region and another drawing from a Maxwell-Boltzmann distribution)
   * Un-weighted event-by-event sampling using (1) a basic rejection sampling method and (2) a Markov Chain Monte Carlo (MCMC) using the Metropolis-Hastings algorithm.
 * Some basic limit setting for a simple counting analysis:
   * Background-free counting
   * Feldman-Cousins confidence intervals
   * CLs limits
 * Detector effects: very basic classes for:
   * Efficiency curves
   * Reconstruction effects
   * Realistically, the user will need to make custon classes for their experiment
 * Examples:
   * Running threaded processes
   * Annual modulation curves
   * Limit plot generation 
   * Comparison of sampling methods
   * MCMC tuning

## Future Features
 * Data for common nuclei
 * More limit setting stuff
   * Maximum Gap (Yellen)
   * Annual modulation limits
   * Bayesian limits
   * Parameter fitting for positive results
   * Detector/model systematics treatment (easier in Bayesian case?)
 * Examples of various plots and calculations
   * Recoil distribution skymaps
   * Sidereal modulation skymaps
 * References and readings on dark matter
 * Maybe/Might be fun
   * Simplified parameterized simulation of a LUX or XENON type detector
   * Inelastic dark matter
   * Q^2-dependent cross sections
   * Coherent neutrino elastic scattering
