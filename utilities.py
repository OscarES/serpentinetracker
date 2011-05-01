#!/usr/bin/python

# Changes to BasicCav and Cavity by Gemmma Smith
#27/01/2011

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
"""A few useful functions"""
import numpy as np

def DriftRmat(L):
    """Return the R matrix for a drift"""
    R = np.array([
        [1, L, 0, 0, 0, 0],
        [0, 1, 0, 0, 0, 0],
        [0, 0, 1, L, 0, 0],
        [0, 0, 0, 1, 0, 0],
        [0, 0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0, 1],
        ])
    return R

def DriftTmat():
    """Return the T matrix for a drift"""
    T = np.zeros((6, 6, 6))
    return T

def SplitParams(param):
    """Split a parameter into a 1x2 np.array"""
    outparam = np.array([0.0, 0.0])
    if np.array(param).size == 1:
        outparam[0] = param
        outparam[1] = param
    elif np.array(param).size > 1:
        outparam[0] = param[0]
        outparam[1] = param[1]
    return outparam

def RotMats(alpha):
    """Calculate the matrix required to rotate the beam by a
    particular angle"""
    C = np.cos(-alpha)
    S = np.sin(-alpha)
    r_inrot = np.array([
        [C , 0 , S , 0 , 0 , 0],
        [0 , C , 0 , S , 0 , 0],
        [-S, 0 , C , 0 , 0 , 0],
        [0 , -S, 0 , C , 0 , 0],
        [0 , 0 , 0 , 0 , 1 , 0],
        [0 , 0 , 0 , 0 , 0 , 1],
        ])
    C = np.cos(alpha)
    S = np.sin(alpha)
    r_outrot = np.array([
        [C , 0 , S , 0 , 0 , 0],
        [0 , C , 0 , S , 0 , 0],
        [-S, 0 , C , 0 , 0 , 0],
        [0 , -S, 0 , C , 0 , 0],
        [0 , 0 , 0 , 0 , 1 , 0],
        [0 , 0 , 0 , 0 , 0 , 1],
        ])
    return r_inrot, r_outrot

