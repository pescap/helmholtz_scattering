# helmholtz_scattering problem
We consider both sound-soft and sound-hard diffraction problems by a radius one sphere.
We solve those problems using bempp.
We compute the coefficients of the spherical harmonic functions.
We compare the error of the Cauchy Data.
Also, we validate the calculus of the integral operators in "validation_integral_operators.py"

#Notes:
Tested with bempp 3.0.3.
To get the directory do: 

inline 'git clone https://github.com/pescap/helmholtz_scattering.git'

Sound-soft problem tested with a number of precision and frequencies.
Work in progress for Sound-hard problem.
Calculus of volumic exact solutions not done yet.

# References:

http://www.bempp.org/
http://volta.byu.edu/winzip/scalar_sphere.pdf
