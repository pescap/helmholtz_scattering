# -*- coding: utf-8 -*-
# Paul Escapil-Inchauspé 

#######
# We verify the convergence of the natural norm of the boundary integral operators
# We use the diagonalization properties of the spherical harmonics
# Allows to a comparison to analytical values
#######

import numpy as np
import bempp.api as bem

h = 1

gr = sphere_bemtool(h)
gr.plot()


bem.global_parameters.assembly.boundary_operator_assembly_type = 'dense'

print gr.leaf_view.entity_count(2), "Number of dofs"

points = gr.leaf_view.vertices.T

def Get_theta(x):
	return np.arccos(x[2])

def Get_phi(x):
	res = np.arctan2(x[1],x[0])	
	if res < 0:
		res += 2 * np.pi
	return res

sp = bem.function_space(gr, 'P', 1)

import scipy
from scipy.special import sph_harm, sph_jn, sph_yn, lpn

p = 1
kappa = 1

def Y(x,n,d,r):
	phi = Get_phi(x)
	theta = Get_theta(x)
	r[0] = sph_harm(p,p,phi,theta)
	#r[0] = sph_harm(p,p,theta,phi)

Uinc = bem.GridFunction(sp, fun = Y)
#Uinc.plot()


V = bem.operators.boundary.helmholtz.single_layer(sp,sp,sp,kappa)
K = bem.operators.boundary.helmholtz.double_layer(sp,sp,sp,kappa)
W = bem.operators.boundary.helmholtz.hypersingular(sp,sp,sp,kappa)
M = bem.operators.boundary.sparse.identity(sp,sp,sp)

Vmat = bem.as_matrix(V.weak_form())
Kmat = bem.as_matrix(K.weak_form())
Wmat = bem.as_matrix(W.weak_form())
Mmat = bem.as_matrix(M.weak_form())

# Orthogonality

Utemp = Uinc.projections()
Err = np.vdot(Utemp,Uinc.coefficients)

print abs(np.real(Err)-1), 'Err Orthogonality'

# V

j_p, dj_p = sph_jn(p, kappa)
y_p, dy_p = sph_yn(p, kappa)
h_p = j_p + 1j * y_p
dh_p = dj_p + 1j * dy_p


#Temp = V * Uinc
#Res_V = np.vdot(Uinc.coefficients, Temp.projections())

# Error for Single Layer

Temp = np.dot(Vmat, Uinc.coefficients)

print Temp
Res_V = np.vdot(Uinc.coefficients, Temp)
Ref_V = 1j * kappa * j_p[p] * h_p[p]

print np.abs(Res_V-Ref_V)/np.abs(Ref_V), 'Error V'

# Error for Double Layer

Temp = np.dot(Kmat, Uinc.coefficients)
Res_K = np.vdot(Uinc.coefficients, Temp)
Ref_K= 1j /2. * kappa**2 * (dj_p[p] * h_p[p] + j_p[p]*dh_p[p])

print np.abs(Res_K-Ref_K)/np.abs(Ref_K), 'Error K'

# Error for Hypersingular

Temp = np.dot(Wmat, Uinc.coefficients)
Res_W = np.vdot(Uinc.coefficients, Temp)
Ref_W = -1j * kappa**3 * dj_p[p] * dh_p[p]

print np.abs(Res_W-Ref_W)/np.abs(Ref_W), 'Error W'

print Ref_K, "ref_K"




