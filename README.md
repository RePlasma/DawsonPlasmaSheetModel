# 1DPSM

1D plasma sheet model by John Dawson, implemented in Matlab.

---------------------------
###  Authors / Core Developers

  - [Filipe Cruz](https://github.com/filipe-dcruz)
  - [Óscar Amaro](https://github.com/OsAmaro)

---------------------------
### Outline

```/breakcold``` Wavebreaking in a cold plasma. Considers the electron sheet at the far left flowing to the right with constant velocity. This velocity is constantly fixed to be the same.

> breakcold_aux.m - function will return results of the system for a given parameters including the sheet velocity

> breakcold.m - uses the function in breakcold_aux.m to calculate the maximum value of the electric field, of the maximum velocity and number of collision and produce a plot.

```/conservation``` Conservation of energy study for warm plasmas with a waterbag distribution.

> conservation_aux.m - function will return error of energy of the system for a given parameters

> conservation.m - uses the function in conservation_aux.m to obtain a plot for error of the energy

```/drag``` Drag on a fast sheet. Conservation of energy study for warm plasmas with a waterbag distribution.

> drag_aux.m - calculates results for a a system where the left particle as a given initial velocity.

> drag.m - get plot of the sheet velocity with time for different initial velocities.

```/thermalize``` Evolution of a waterbag towards a maxwellian.
> thermalize.m - returns plot of velocity distribution for the initial and final time

```/wakefield``` Shape of the electric field from an external driver.

```/waterbag``` Generation of velocity distributions (waterbag and maxwellian).

```/wavelength``` Get wavelength from phasespace of the wake of an external driver plasma sheet.

---------------------------
### References:

[1] J. Dawson, The Physics of Fluids 5, 445 (1962), https://aip.scitation.org/doi/pdf/10.1063/1.1706638, URL https://aip.scitation.org/doi/abs/10.1063/1. 1706638.
