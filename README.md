# helmholtz_scattering problem
We consider both sound-soft and sound-hard diffraction problems by a radius one sphere.
We solve those problems using Bempp.
We compute the coefficients of the spherical harmonic functions.
We compare the numerical error of the Cauchy data.
Also, we validate the calculus of the integral operators in "validation_integral_operators.py"

# Notes:
Tested with Bempp 3.0.3.
To get the directory do: 

`git clone https://github.com/pescap/helmholtz_scattering.git`

Sound-soft problem tested with a number of precisions and frequencies.

Work in progress for sound-hard problem. (not stable).

A little documentation available in "scattering.tex" (or "scattering.pdf").

Calculus of volumic exact solutions not done yet.

# References:

http://www.bempp.org/

http://volta.byu.edu/winzip/scalar_sphere.pdf
