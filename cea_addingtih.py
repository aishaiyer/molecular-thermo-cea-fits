"""
TiH thermodynamic coefficient reconstruction for NASA CEA-style chemistry input.

This script reproduces the legacy workflow used to derive TiH thermodynamic
coefficients from a tabulated partition function, then fit those quantities to
the NASA polynomial form used by CEA-style equilibrium chemistry solvers.

Workflow
--------
1. Read a tabulated TiH partition function Q(T).
2. Compute Cp/R, H/RT, and S/R from partition-function relations.
3. Fit those thermodynamic curves with NASA/CEA-style polynomial coefficients.
4. Export a CEA-style coefficient block for use in equilibrium chemistry.

What this script does NOT do
----------------------------
- It does not compute TiH line opacities or cross sections.
- It does not generate the final TiH abundance/mixing-fraction plot directly.
  That plot would be produced downstream after these coefficients are passed
  into a chemistry solver such as NASA CEA or a CEA-style thermo module.

References
----------
Barklem, P. S., & Collet, R. 2016,
"Partition functions and equilibrium constants for diatomic molecules and atoms
of astrophysical interest," Astronomy & Astrophysics, 588, A96.
https://doi.org/10.1051/0004-6361/201526961

Burrows, A., Dulick, M., Bauschlicher, C. W., Jr., Bernath, P., Ram, R. S.,
Sharp, C. M., & Milsom, J. A. 2005,
"Spectroscopic Constants, Abundances, and Opacities of the TiH Molecule,"
The Astrophysical Journal, 624, 988-1002.
https://doi.org/10.1086/429366

McBride, B. J., Gordon, S., & Reno, M. A. 1993,
"Coefficients for Calculating Thermodynamic and Transport Properties of
Individual Species," NASA Reference Publication 1311.

Gordon, S., & McBride, B. J. 1994,
"Computer Program for Calculation of Complex Chemical Equilibrium Compositions
and Applications," NASA Reference Publication 1311.

Iyer, A. R., Line, M. R., Muirhead, P. S., Fortney, J. J.,
& Gharib-Nezhad, E. 2023,
"The SPHINX M-dwarf Spectral Grid. I. Benchmarking New Model Atmospheres
to Derive Fundamental M-dwarf Properties,"
The Astrophysical Journal, 944, 41.
https://doi.org/10.3847/1538-4357/acabc2

Iyer, A. R., Line, M. R., Muirhead, P. S., Fortney, J. J.,
& Faherty, J. K. 2026,
"The SPHINX M Dwarf Spectral Grid. II. New Model Atmospheres and Spectra
to Derive the Fundamental Properties of Mid-to-late Type M Dwarfs,"
The Astrophysical Journal, 998, 88.
https://doi.org/10.3847/1538-4357/ae285d
"""

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit


R_GAS = 8.31446261815324  # J K^-1 mol^-1


# Barklem & Collet (2016) partition-function table for TiH.
# This is the thermodynamic input used to derive NASA/CEA-style coefficients.
BARKLEM_T = np.array([
    1.00000e-05, 1.00000e-04, 1.00000e-03, 1.00000e-02, 1.00000e-01,
    1.50000e-01, 2.00000e-01, 3.00000e-01, 5.00000e-01, 7.00000e-01,
    1.00000e+00, 1.30000e+00, 1.70000e+00, 2.00000e+00, 3.00000e+00,
    5.00000e+00, 7.00000e+00, 1.00000e+01, 1.50000e+01, 2.00000e+01,
    3.00000e+01, 5.00000e+01, 7.00000e+01, 1.00000e+02, 1.30000e+02,
    1.70000e+02, 2.00000e+02, 2.50000e+02, 3.00000e+02, 5.00000e+02,
    7.00000e+02, 1.00000e+03, 1.50000e+03, 2.00000e+03, 3.00000e+03,
    4.00000e+03, 5.00000e+03, 6.00000e+03, 7.00000e+03, 8.00000e+03,
    9.00000e+03, 1.00000e+04,
])

BARKLEM_Q_TIH = np.array([
    1.40000e+01, 1.40000e+01, 1.40000e+01, 1.40000e+01, 1.40000e+01,
    1.40000e+01, 1.40000e+01, 1.40000e+01, 1.40000e+01, 1.40000e+01,
    1.40000e+01, 1.40000e+01, 1.40000e+01, 1.40000e+01, 1.40000e+01,
    1.40001e+01, 1.40027e+01, 1.40383e+01, 1.43008e+01, 1.48634e+01,
    1.66856e+01, 2.21996e+01, 2.95947e+01, 4.36225e+01, 6.06378e+01,
    8.70653e+01, 1.09163e+02, 1.49564e+02, 1.93935e+02, 4.12370e+02,
    7.14423e+02, 1.37894e+03, 3.22730e+03, 6.28151e+03, 1.74678e+04,
    3.83544e+04, 7.26767e+04, 1.23609e+05, 1.93324e+05, 2.83017e+05,
    3.93090e+05, 5.23344e+05,
])


@dataclass
class ThermoCurves:
    temperature: np.ndarray
    cp_over_r: np.ndarray
    h_over_rt: np.ndarray
    s_over_r: np.ndarray


@dataclass
class CEACoefficients:
    a1: float
    a2: float
    a3: float
    a4: float
    a5: float
    a6: float
    a7: float
    a8: float
    a9: float

    @property
    def as_array(self) -> np.ndarray:
        return np.array([
            self.a1, self.a2, self.a3, self.a4, self.a5,
            self.a6, self.a7, self.a8, self.a9,
        ])

    def format_cea_block(
        self,
        species: str = "TiH",
        t_low: float = 298.15,
        t_high: float = 1500.0,
        t_mid: float = 1500.0,
        phase: str = "G",
    ) -> str:
        coeffs = self.as_array
        return (
            f"{species:<16s}{phase:>2s}   {t_low:8.3f}  {t_high:8.3f}  {t_mid:8.3f}\n"
            f" {coeffs[0]: .8E} {coeffs[1]: .8E} {coeffs[2]: .8E} {coeffs[3]: .8E} {coeffs[4]: .8E}\n"
            f" {coeffs[5]: .8E} {coeffs[6]: .8E} {coeffs[7]: .8E} {coeffs[8]: .8E}"
        )


def nasa_cp_over_r(
    temperature: np.ndarray,
    a1: float, a2: float, a3: float, a4: float, a5: float, a6: float, a7: float
) -> np.ndarray:
    return (
        a1 * temperature ** -2
        + a2 * temperature ** -1
        + a3
        + a4 * temperature
        + a5 * temperature ** 2
        + a6 * temperature ** 3
        + a7 * temperature ** 4
    )


def nasa_h_over_rt(
    temperature: np.ndarray,
    a1: float, a2: float, a3: float, a4: float,
    a5: float, a6: float, a7: float, a8: float
) -> np.ndarray:
    return (
        -a1 * temperature ** -2
        + a2 * np.log(temperature) * temperature ** -1
        + a3
        + 0.5 * a4 * temperature
        + (a5 * temperature ** 2) / 3.0
        + (a6 * temperature ** 3) / 4.0
        + (a7 * temperature ** 4) / 5.0
        + a8 / temperature
    )


def nasa_s_over_r(
    temperature: np.ndarray,
    a1: float, a2: float, a3: float, a4: float,
    a5: float, a6: float, a7: float, a9: float
) -> np.ndarray:
    return (
        -0.5 * a1 * temperature ** -2
        - a2 * temperature ** -1
        + a3 * np.log(temperature)
        + a4 * temperature
        + 0.5 * a5 * temperature ** 2
        + (a6 * temperature ** 3) / 3.0
        + 0.25 * a7 * temperature ** 4
        + a9
    )


def compute_thermo_from_partition_function(
    temperature: np.ndarray,
    partition_function: np.ndarray,
    fit_temperature: np.ndarray,
) -> ThermoCurves:
    """
    Compute Cp/R, H/RT, and S/R from a partition function table.

    Notes
    -----
    This follows the same logic as the original script:
    1. spline-fit Q(T)
    2. derive Q', Q''
    3. compute thermodynamic functions from partition-function relations
    """
    spline_q = UnivariateSpline(temperature, partition_function, s=0)

    q = spline_q(fit_temperature)
    q_prime = spline_q.derivative(1)(fit_temperature) * fit_temperature
    q_double_prime = (
        spline_q.derivative(2)(fit_temperature) * fit_temperature ** 2
        + 2.0 * q_prime
    )

    cp_over_r = (q_double_prime / q) - (q_prime / q) ** 2

    h0 = R_GAS * 298.15 * (
        spline_q.derivative(1)(298.15) * 298.15 / spline_q(298.15)
    )
    h_over_rt = np.log(q) + (h0 / (R_GAS * fit_temperature))
    s_over_r = (q_prime / q) + np.log(q)

    return ThermoCurves(
        temperature=fit_temperature,
        cp_over_r=cp_over_r,
        h_over_rt=h_over_rt,
        s_over_r=s_over_r,
    )


def fit_cea_coefficients(
    curves: ThermoCurves,
    initial_guess: np.ndarray,
) -> CEACoefficients:
    """
    Fit NASA/CEA-style coefficients.

    This preserves the spirit of the original script:
    - fit Cp first for a1..a7
    - fit H for a1..a8
    - fit S for a1..a7,a9
    - iterate once more to stabilize the solution
    """
    t = curves.temperature

    p_cp, _ = curve_fit(
        nasa_cp_over_r,
        t,
        curves.cp_over_r,
        p0=initial_guess[:7],
        maxfev=50000,
    )

    p_h, _ = curve_fit(
        nasa_h_over_rt,
        t,
        curves.h_over_rt,
        p0=np.append(p_cp, initial_guess[7]),
        maxfev=50000,
    )

    p_s, _ = curve_fit(
        nasa_s_over_r,
        t,
        curves.s_over_r,
        p0=np.append(p_h[:7], initial_guess[8]),
        maxfev=50000,
    )

    p_cp, _ = curve_fit(
        nasa_cp_over_r,
        t,
        curves.cp_over_r,
        p0=p_s[:7],
        maxfev=50000,
    )

    p_h, _ = curve_fit(
        nasa_h_over_rt,
        t,
        curves.h_over_rt,
        p0=np.append(p_cp, p_h[-1]),
        maxfev=50000,
    )

    p_s, _ = curve_fit(
        nasa_s_over_r,
        t,
        curves.s_over_r,
        p0=np.append(p_h[:7], p_s[-1]),
        maxfev=50000,
    )

    return CEACoefficients(
        a1=p_cp[0],
        a2=p_cp[1],
        a3=p_cp[2],
        a4=p_cp[3],
        a5=p_cp[4],
        a6=p_cp[5],
        a7=p_cp[6],
        a8=p_h[-1],
        a9=p_s[-1],
    )


def plot_fit_residuals(curves: ThermoCurves, coeffs: CEACoefficients) -> None:
    """Plot residuals between the derived thermo curves and the NASA polynomial fit."""
    t = curves.temperature
    c = coeffs.as_array

    cp_model = nasa_cp_over_r(t, *c[:7])
    h_model = nasa_h_over_rt(t, *c[:8])
    s_model = nasa_s_over_r(t, *np.append(c[:7], c[8]))

    plt.figure(figsize=(8, 5))
    plt.plot(t, curves.cp_over_r - cp_model, label="Cp/R residual")
    plt.plot(t, curves.h_over_rt - h_model, label="H/RT residual")
    plt.plot(t, curves.s_over_r - s_model, label="S/R residual")
    plt.axhline(0.0, linewidth=1)
    plt.xlabel("Temperature (K)")
    plt.ylabel("Residual")
    plt.legend()
    plt.tight_layout()
    plt.show()


def main() -> None:
    # Temperature range actually used in the original script for the fit.
    fit_temperature = np.arange(298.15, 1500.0, 0.1)

    # Original seed from the legacy script.
    initial_guess = np.array([
        -3.750684436e04,
         6.229478175e02,
        -9.360137531e-02,
         8.631692229e-03,
        -7.672598171e-06,
         3.428472624e-09,
        -6.107038964e-13,
         4.716527478e04,
         2.629609365e01,
    ])

    curves = compute_thermo_from_partition_function(
        temperature=BARKLEM_T,
        partition_function=BARKLEM_Q_TIH,
        fit_temperature=fit_temperature,
    )

    coeffs = fit_cea_coefficients(curves, initial_guess)

    print("CEA-style coefficients for TiH:")
    for index, value in enumerate(coeffs.as_array, start=1):
        print(f"a{index} = {value: .8E}")

    print("\nCEA block:\n")
    print(coeffs.format_cea_block())

    plot_fit_residuals(curves, coeffs)


if __name__ == "__main__":
    main()