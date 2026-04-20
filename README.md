# TiH Thermodynamic Fit for NASA/CEA

This repository contains a reconstructed and cleaned workflow for deriving TiH thermodynamic coefficients from partition-function data and fitting them to the NASA polynomial form used in equilibrium chemistry solvers (e.g., CEA-style codes).

## Overview

TiH is not included in many standard thermochemical databases. This workflow was developed to construct TiH thermodynamic inputs for use in stellar atmosphere and chemical equilibrium modeling.

The steps are:

1. Start from a tabulated partition function \( Q(T) \)
2. Compute thermodynamic quantities:
   - \( C_p/R \)
   - \( H/RT \)
   - \( S/R \)
3. Fit these to the NASA polynomial form used by CEA
4. Export coefficients for use in equilibrium chemistry calculations

These coefficients can then be used in chemistry solvers to compute TiH abundances.

## Context (SPHINX)

This work was developed as part of the **SPHINX stellar atmosphere modeling framework**, which focuses on self-consistent modeling of M-dwarf atmospheres and spectra. Accurate thermochemistry is critical for modeling molecular abundances and understanding stellar contamination in exoplanet observations.

## What this repository does NOT include

- TiH opacity or cross-section calculations  
- Radiative transfer or spectral synthesis  
- Direct abundance plotting (e.g., TiH mixing fraction vs temperature)

Those steps occur downstream after these coefficients are used in a chemistry solver.

## References

- Barklem, P. S., & Collet, R. (2016)  
  *Partition functions and equilibrium constants for diatomic molecules and atoms of astrophysical interest*,  
  Astronomy & Astrophysics, 588, A96  
  https://doi.org/10.1051/0004-6361/201526961  

- Burrows, A., Dulick, M., Bauschlicher, C. W., Jr., Bernath, P. F., Ram, R. S., Sharp, C. M., & Milsom, J. A. (2005)  
  *Spectroscopic Constants, Abundances, and Opacities of the TiH Molecule*,  
  The Astrophysical Journal, 624, 988–1002  
  https://doi.org/10.1086/429366  

- McBride, B. J., Gordon, S., & Reno, M. A. (1993)  
  *Coefficients for Calculating Thermodynamic and Transport Properties of Individual Species*,  
  NASA Reference Publication 1311  

- Gordon, S., & McBride, B. J. (1994)  
  *Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications*,  
  NASA Reference Publication 1311  

- Iyer, A. R., Line, M. R., Muirhead, P. S., Fortney, J. J., & Gharib-Nezhad, E. (2023)  
  *The SPHINX M-dwarf Spectral Grid. I. Benchmarking New Model Atmospheres to Derive Fundamental M-dwarf Properties*,  
  The Astrophysical Journal, 944, 41  
  https://doi.org/10.3847/1538-4357/acabc2  

- Iyer, A. R., Line, M. R., Muirhead, P. S., Fortney, J. J., & Faherty, J. K. (2026)  
  *The SPHINX M Dwarf Spectral Grid. II. New Model Atmospheres and Spectra to Derive the Fundamental Properties of Mid-to-late Type M Dwarfs*,  
  The Astrophysical Journal, 998, 88  
  https://doi.org/10.3847/1538-4357/ae285d  

## Notes

This repository preserves the original fitting approach used historically. It is intended as a transparent and documented reconstruction of the method rather than an optimized modern implementation.
