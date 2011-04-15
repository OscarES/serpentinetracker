#!/usr/bin/python

#    Copyright 2009,  Stephen Molloy, Stewart Boogert
# 
#    This file is part of Serpentine.
#
#    Serpentine is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Serpentine is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Serpentine.  If not, see <http://www.gnu.org/licenses/>.
#
"""Module containing the standard beam representations for Serpentine sims."""

import numpy as np
from numpy.random import normal
import matplotlib.pylab as plt
from globals import electron_mass, proton_mass

# Superclass that defaults all new beams to being electrons
class Beam:
    """__init__(self, P=1, Q=1e-9, chargesign=-1, restmass=electron_mass)
    The primary beam class from which all others inherit.
    This will generate a beam with all the necessary methods.
    It will contain a single particle at [0,0,0,0,0,P] with the
    charge and restmass of an electron."""
    def __init__(self, P=1, Q=1e-9, chargesign=-1, restmass=electron_mass):
        self.chargesign = chargesign
        self.restmass   = restmass
        self.x          = np.array([0, 0, 0, 0, 0, P])
        self.x.shape    = (6, 1)
        self.Q          = Q

    def MakePositrons(self):
        "Convert the bunch species to positrons"
        self.chargesign = 1
        self.restmass   = electron_mass

    def MakeElectrons(self):
        "Convert the bunch species to electrons"
        self.chargesign = -1
        self.restmass   = electron_mass

    def MakeProtons(self):
        "Convert the bunch species to electrons"
        self.chargesign = 1
        self.restmass   = proton_mass

    def PlotXSpace(self, marker='bx'):
        "Plot the x phase-space of the bunch"
        plt.plot(self.x[0, :], self.x[1, :], marker)
        plt.xlabel('x Position / m')
        plt.ylabel('x Angle / rad')

    def PlotYSpace(self, marker='bx'):
        "Plot the y phase-space of the bunch"
        plt.plot(self.x[2, :], self.x[3, :], marker)
        plt.xlabel('y Position / m')
        plt.ylabel('y Angle / rad')

    def PlotLongSpace(self, marker='bx'):
        "Plot the longitudinal phase-space of the bunch"
        plt.plot(self.x[4, :], self.x[5, :], marker)
        plt.xlabel('Distance / m')
        plt.ylabel('Energy / GeV')
        
    def PlotXLongSpace(self, marker='bx'):
        "Plot the y Vs. z of the bunch"
        plt.plot(self.x[4, :], self.x[0, :], marker)
        plt.xlabel('Distance / m')
        plt.ylabel('x Position / m')
        
    def PlotYLongSpace(self, marker='bx'):
        "Plot the y Vs. z of the bunch"
        plt.plot(self.x[4, :], self.x[2, :], marker)
        plt.xlabel('Distance / m')
        plt.ylabel('y Position / m')

    def PlotPxLongSpace(self, marker='bx'):
        "Plot the Px Vs. z of the bunch"
        plt.plot(self.x[4, :], self.x[1, :], marker)
        plt.xlabel('Distance / m')
        plt.ylabel('x Angle / rad')
        
    def PlotPyLongSpace(self, marker='bx'):
        "Plot the Py Vs. z of the bunch"
        plt.plot(self.x[4, :], self.x[3, :], marker)
        plt.xlabel('Distance / m')
        plt.ylabel('y Angle / rad')

# =========================================================
# A beam with a gaussian spread in 6D
class GaussBeam(Beam):
    """__init__(self, N=1000, Q=1e-9, pos=np.array([0, 0, 0, 0, 0, 1]), 
    sig=np.array([1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 0.01]), chargesign=-1,
    restmass=electron_mass)
    Create a multi-particle beam with a Gaussian spread in each of the 
    6 dimensions"""
    def __init__(self, N=1000, Q=1e-9, pos=np.array([0, 0, 0, 0, 0, 1]), 
        sig=np.array([1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 0.01]), chargesign=-1, 
        restmass=electron_mass):
        Beam.__init__(self, chargesign, restmass)
        self.x = np.array([
            normal(loc=pos[0], scale=sig[0], size=N),
            normal(loc=pos[1], scale=sig[1], size=N),
            normal(loc=pos[2], scale=sig[2], size=N),
            normal(loc=pos[3], scale=sig[3], size=N),
            normal(loc=pos[4], scale=sig[4], size=N),
            normal(loc=pos[5], scale=sig[5], size=N)
            ])
        for dof in range(6):
            self.x[dof, :] -= np.mean(self.x[dof, :]) - pos[dof]
        self.Q = np.ones(N) * Q/N

# Similar to GaussBeam, but calculated from initial Twiss params
class TwissGaussBeam(GaussBeam):
    """__init__(self, twiss, N=1000, pos=np.array([0, 0, 0, 0, 0, 1]),
    Q=1e-9, chargesign=-1, restmass=electron_mass)
    Create a multi-particle beam with a Gaussian spread in each of the 
    6 dimensions.  The statistics of the beam distribution are calculated
    from the Twiss parameters in twiss."""
    def __init__(self, twiss, N=1000, pos=np.array([0, 0, 0, 0, 0, 1]), 
        Q=1e-9, chargesign=-1, restmass=electron_mass):
        P = pos[5]
        gammasquared = ( (P*1e9)**2 + restmass**2 ) / restmass**2
        betagamma = np.sqrt(gammasquared - 1)
        sig = np.zeros(6)
        sig[0] = np.sqrt((twiss.nemitx/betagamma) * twiss.betax)
        sig[1] = np.sqrt((twiss.nemitx/betagamma) / twiss.betax)
        sig[2] = np.sqrt((twiss.nemity/betagamma) * twiss.betay)
        sig[3] = np.sqrt((twiss.nemity/betagamma) / twiss.betay)
        sig[4] = twiss.sigz
        sig[5] = twiss.sigP
        
        rmat = np.eye(6, 6)
        rmat[1, 0] = -twiss.alphax / twiss.betax
        rmat[3, 2] = -twiss.alphay / twiss.betay

        rmat[5, 4] = twiss.pz_cor
        rmat[0, 5] = twiss.etax / P
        rmat[1, 5] = twiss.etaxp / P
        rmat[2, 5] = twiss.etay / P
        rmat[3, 5] = twiss.etayp / P

        GaussBeam.__init__(self, N, Q, pos, sig, chargesign, restmass)

        self.x = np.dot(rmat, self.x)
        self.x[5, :] += P

        self.x[0, :] *= sig[0]/np.std(self.x[0, :])
        self.x[1, :] *= sig[1]/np.std(self.x[1, :])
        self.x[2, :] *= sig[2]/np.std(self.x[2, :])
        self.x[3, :] *= sig[3]/np.std(self.x[3, :])
        self.x[4, :] *= sig[4]/np.std(self.x[4, :])
        self.x[5, :] *= sig[5]/np.std(self.x[5, :])

        self.x[0, :] -= np.mean(self.x[0, :]) - pos[0]
        self.x[1, :] -= np.mean(self.x[1, :]) - pos[1]
        self.x[2, :] -= np.mean(self.x[2, :]) - pos[2]
        self.x[3, :] -= np.mean(self.x[3, :]) - pos[3]
        self.x[4, :] -= np.mean(self.x[4, :]) - pos[4]
        self.x[5, :] -= np.mean(self.x[5, :]) - pos[5]

if __name__ == '__main__':
    mytwiss = {}
    mytwiss.betax = 6.85338806855804
    mytwiss.alphax = 1.11230788371885
    mytwiss.etax = 3.89188697330735e-012
    mytwiss.etaxp = 63.1945125619190e-015
    mytwiss.betay = 2.94129410712918
    mytwiss.alphay = -1.91105724003646
    mytwiss.etay = 0
    mytwiss.etayp = 0
    mytwiss.Nemitx = 5.08807339588144e-006
    mytwiss.Nemity = 50.8807339588144e-009
    mytwiss.sigz = 8.00000000000000e-003
    mytwiss.sigP = 1.03999991965541e-003
    mytwiss.pz_cor = 0
    mytwiss.x = 0
    mytwiss.xp = 0
    mytwiss.yp = 0
    mytwiss.y = 0
    mytwiss.s = 0
    mytwiss.P = 1.3
    mybeam = TwissGaussBeam(mytwiss, N=10000, Q=1.60217646200000e-009)
    plt.figure(1)
    mybeam.PlotXSpace()
    plt.figure(2)
    mybeam.PlotYSpace()
    plt.figure(3)
    mybeam.PlotLongSpace()

