# -*- coding: utf-8 -*-
# Paul Escapil-Inchauspé 

#######
# Exterior soft-sphere problem by using Bempp 
# We solve the problem using both direct and combined formulations
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

def combined_data(x, n, domain_index, result):
    result[0] = 1j * kappa * np.exp(1j * kappa * x[0]) * (n[0]-1) 


grid = bem.shapes.sphere(h=h)
space = bem.function_space(grid, "P", 2) # Can be applied to any function space

# Assembly of operators

V = bem.operators.boundary.helmholtz.single_layer(space, space, space, kappa)
Vop = V.weak_form()
ID = bem.operators.boundary.sparse.identity(space, space, space)
K = bem.operators.boundary.helmholtz.double_layer(space, space, space, kappa)
KP = bem.operators.boundary.helmholtz.adjoint_double_layer(space, space, space,kappa)

dirichlet_data = bem.GridFunction(space, fun=U_inc)
neumann_data = bem.GridFunction(space, fun=gamma_1_U_inc)

# We work with the total field U

# Direct formulation

lhs = Vop
rhs_0 = dirichlet_data.projections()
neumann_coeff_U, info = scipy.sparse.linalg.gmres(lhs, rhs_0, tol=1E-5)

# Combined formulation

lhs_combined =  .5*ID + KP - 1j * kappa * V
lhs_combined = lhs_combined.weak_form()

rhs_combined = bem.GridFunction(space, fun=combined_data)
rhs_combined = rhs_combined.projections()

neumann_coeff_U_combined, info = scipy.sparse.linalg.gmres(lhs_combined, rhs_combined, tol=1E-5)

# Convert into bempp grid 

neumann_U = bem.GridFunction(space, coefficients=neumann_coeff_U)
neumann_U_combined = bem.GridFunction(space, coefficients=neumann_coeff_U_combined)

neumann_U.plot()
#neumann_U_combined.plot()

sol = Helmholtz_Exact_Solution(kappa)

neumann_U_exact_fun = sol.uExactNeumannTrace
neumann_U_exact = bem.GridFunction(space, fun=sol.uExactBoundaryNeumannTrace)

#neumann_U_exact.plot()

# Error 

print "relative error direct: ", neumann_U.relative_error(neumann_U_exact_fun)
print "relative error combined: ", neumann_U_combined.relative_error(neumann_U_exact_fun)
