# -*- coding: utf-8 -*-
# Paul Escapil-Inchauspé 

#######
# Exterior sound-hard scattering problem by a radius 1 sphere using Bempp 
# We use a direct formulation
#######

from __future__ import division
import bempp.api as bem
import numpy as np
import scipy.sparse.linalg

from exact_solution import Helmholtz_Exact_Solution


# Parameters: 
#   kappa: wavenumber
#   h : mesh width          

kappa = 1
precision = 10
h = 2.0 * np.pi / (precision * kappa)

def U_inc(x, n, domain_index, result):
    result[0] = np.exp(1j * kappa * x[0])

def gamma_1_U_inc(x, n, domain_index, result):
    result[0] = 1j * kappa * np.exp(1j * kappa * x[0]) * n[0]


grid = bem.shapes.sphere(h=h)
space = bem.function_space(grid, "P", 1) # Can be applied to any function space

# Assembly of the hypersingular

W = bem.operators.boundary.helmholtz.hypersingular(space, space, space, kappa)
Wop = W.weak_form()

# We work with the total field U

neumann_data = bem.GridFunction(space, fun=gamma_1_U_inc)
lhs = Wop
rhs_0 = neumann_data.projections()
dirichlet_coeff_U, info = scipy.sparse.linalg.gmres(lhs, rhs_0, tol=1E-5)

# Convert into bempp grid 

dirichlet_U = bem.GridFunction(space, coefficients=dirichlet_coeff_U)


dirichlet_U.plot()

sol = Helmholtz_Exact_Solution(kappa)

dirichlet_U_exact_fun = sol.uExactDirichletTrace
dirichlet_U_exact = bem.GridFunction(space, fun=sol.uExactBoundaryDirichletTrace)

#neumann_U_exact.plot()

# Error 

print "relative error direct: ", dirichlet_U.relative_error(dirichlet_U_exact_fun)

