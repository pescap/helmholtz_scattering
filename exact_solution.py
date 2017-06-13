# Compare the numerical and analytical solutions 

# Define functions that evaluate the Dirichlet and Neumann trace of the
# analytical solution, expressed as a series of spherical wave functions
# Done for an exterior sound-soft and soft-hard scattering problem by a radius 1 sphere. 
# Could be adapted to transmission problem

from __future__ import division
import numpy as np
from scipy.special import sph_jn, sph_yn, lpn

class Helmholtz_Exact_Solution():
    """Interface for getting exact solution of a soft sphere
    scattering problem with Helhmoltz equations

    Parameters:
    k : wavenumber
    """

    def __init__(self, k):
        self.k = k
        self.kExt = k

        l_max = 200
        l = np.arange(0, l_max + 1)

        jExt, djExt = sph_jn(l_max, self.kExt)
        yExt, dyExt = sph_yn(l_max, self.kExt)
        hExt = jExt + 1j * yExt


        aInc = (2 * l + 1) * 1j ** l

        np.seterr(divide='ignore',  invalid='ignore')

        cBound = 1/(1j * self.kExt) / hExt
        cDir = jExt / hExt

        for l in range(l_max + 1):
            if abs(cBound[l]) < 1e-16:
                # neglect all further terms
                l_max = l - 1
                aInc = aInc[:l]
                cBound = cBound[:l]
                cDir = cDir[:l]
                break


        self.cDir = cDir
        self.aInc = aInc
        self.cBound = cBound
        self.l_max = l_max

    def uExactDirichletTrace(self, point):
        x, y, z = point
        r = np.sqrt(x**2 + y**2 + z**2)
        hD, dhD = sph_jn(self.l_max, self.kExt) + 1j*sph_yn(self.l_max, self.kExt * r)
        Y, dY = lpn(self.l_max, x / r)
        return (self.cDir * hD * Y).sum()

    def uExactBoundaryDirichletTrace(self, point, normal, domain_index, result):
        x, y, z = point
        r = np.sqrt(x**2 + y**2 + z**2)
        hD, dhD = sph_jn(self.l_max, self.kExt) + 1j*sph_yn(self.l_max, self.kExt * r)
        Y, dY = lpn(self.l_max, x / r)
        result[0] = (self.cDir * hD * Y).sum()


    def uExactBoundaryNeumannTrace(self, point, normal, domain_index, result):
        x, y, z = point
        n_x, n_y, n_z = normal
        r = np.sqrt(x**2 + y**2 + z**2)
        Y, dY = lpn(self.l_max, x / r)
        result[0] = (self.aInc * self.cBound * Y).sum()
        #return result

    def uExactNeumannTrace(self, point):
        x, y, z = point
        r = np.sqrt(x**2 + y**2 + z**2)
        Y, dY = lpn(self.l_max, x / r)
        val = (self.aInc * self.cBound * Y).sum()
        return val



# Test for sound-soft scattering
# Parameters: 
#   kappa: wavenumber
#   h : mesh width          

import bempp.api as bem
k = 1
h = 0.5
grid = bem.shapes.sphere(h=h)
#if __name__ == "__main__":
 