# TiH thermo fit for CEA-style chemistry input

This repository contains a cleaned version of an older TiH thermochemistry script used to derive NASA/CEA-style polynomial coefficients for TiH from tabulated partition-function data.

## What the script does

The workflow is:

1. Start from a tabulated TiH partition function \(Q(T)\)
2. Compute thermodynamic quantities from the partition function:
   - \(C_p/R\)
   - \(H/RT\)
   - \(S/R\)
3. Fit those curves with the standard NASA polynomial form used by CEA-style equilibrium chemistry solvers
4. Export a CEA-style coefficient block for TiH

These coefficients can then be used in an equilibrium chemistry code to include TiH as a species.

## What the script does not do

This script does **not**:
- compute TiH line opacities or cross sections
- generate the final TiH abundance or mixing-fraction plot directly

Those steps happen downstream, after the fitted thermo coefficients are passed into a chemistry solver.

## Scientific context

This script was originally used for Iyer et al. 2023 and Iyer et al. 2026 to add TiH thermodynamics to a chemistry workflow because TiH was not available in the default thermo database being used. The partition-function input was based on Barklem & Collet (2016), and the resulting TiH behavior was checked against Burrows et al. (2005).
https://iopscience.iop.org/article/10.3847/1538-4357/acabc2/meta
https://iopscience.iop.org/article/10.3847/1538-4357/ae285d
## References

- Barklem, P. S., & Collet, R. 2016, *Partition functions and equilibrium constants for diatomic molecules and atoms of astrophysical interest*, A&A, 588, A96. https://doi.org/10.1051/0004-6361/201526961
- Burrows, A., Dulick, M., Bauschlicher, C. W., Jr., Bernath, P. F., Ram, R. S., Sharp, C. M., & Milsom, J. A. 2005, *Spectroscopic Constants, Abundances, and Opacities of the TiH Molecule*, ApJ, 624, 988–1002. https://doi.org/10.1086/429366
- McBride, B. J., Gordon, S., & Reno, M. A. 1993, *Coefficients for Calculating Thermodynamic and Transport Properties of Individual Species*, NASA RP-1311
- Gordon, S., & McBride, B. J. 1994, *Computer Program for Calculation of Complex Chemical Equilibrium Compositions and Applications*, NASA RP-1311

## Notes

This cleaned version preserves the legacy fitting logic closely so the historical workflow is easy to follow. It is intended primarily as a documented reconstruction of the original method.
