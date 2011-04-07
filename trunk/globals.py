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
# Define all the global variables here to import once
#from numpy import sqrt
from scipy.special import cbrt

electron_mass = 0.51099892e6 # eV/c**2
positron_mass = electron_mass
proton_mass   = 938.272013e6 # eV/c**2
c_light = 299792458 # m/s
e_charge = 1.60217653e-19 # C
Brho1GeV = 1e9 / c_light # T.m / GeV/c


class lietrackparams:
    def __init__(self, dlengths=[], klengths=[]):
        self.dlengths = dlengths
        self.klengths = klengths

lietrackarray = []
lietrackarray.append(False)
lietrackarray.append(lietrackparams(dlengths=[1.0],klengths=[1.0]))
lietrackarray.append(lietrackparams(dlengths=[0.5,0.5],klengths=[1.0]))
lietrackarray.append(False)
lietrackarray.append(lietrackparams(dlengths=[(1/12)*(4+2*cbrt(2)+cbrt(4)),(1/2)-((1/12)*(4+2*cbrt(2)+cbrt(4))),(1/2)-((1/12)*(4+2*cbrt(2)+cbrt(4))),(1/12)*(4+2*cbrt(2)+cbrt(4))],klengths=[(1/6)*(4+2*cbrt(2)+cbrt(4)),1-((1/3)*(4+2*cbrt(2)+cbrt(4))),(1/6)*(4+2*cbrt(2)+cbrt(4))]))

