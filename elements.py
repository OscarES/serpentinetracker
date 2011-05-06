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
"""A module defining the Element superclass, and
individual subclasses for each element type."""

import numpy as np
from copy import deepcopy
import matplotlib.cbook as cb
import types
from numpy.random import normal
from globals import electron_mass, c_light, Brho1GeV, lietrackarray
from scipy.special import jn
import beamline as BL
import beamrep
from utilities import DriftRmat, DriftTmat, SplitParams, RotMats

# Twiss object class for addition to the Element superclass
class Twiss(object):
    """A class to hold Twiss parameters"""
    def __init__(self, 
                 betax  = None,  betay = None, 
                 alphax = None, alphay = None,
                 etax   = None,   etay = None,
                 etapx  = None,  etapy = None,
                 phix   = None,   phiy = None,
                 nemitx = None, nemity = None,
                 sigz   = None,   sigP = None,
                 pz_cor = None):

        self.betax = betax
        self.betay = betay

        self.alphax = alphax
        self.alphay = alphay

        self.etax = etax
        self.etay = etay

        self.etapx = etapx
        self.etapy = etapy

        self.phix = phix
        self.phiy = phiy

        self.nemitx = nemitx
        self.nemity = nemity

        self.sigz = sigz
        self.sigP = sigP
        self.pz_cor = pz_cor

    def __repr__(self) :
        pass

# 'Element' superclass, defining common qualities of all elements
# Also defines CalcRmat for a drift, which will then be overridden
# in classes where this is not true.
class Element:
    """The Element superclass.  All beamline element classes should
    inherit from this class"""
    def __init__(self, name, L=1, P=1, S=0, aper=0, apershape='ELLIPTICAL', 
        is_diag=False):
        self.name = name                                         # Element name
        self.L = L                                               # Element length
        self.P = P                                               # Design momentum 
        self.S = S                                               # S position of start of element
        self.aper = aper
        self.offset = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])   # (x,px,y,py,z,P) 
        self.is_diag = is_diag
        if self.is_diag:
            self.DiagOut = DiagOut()
        self.twiss = Twiss()
        self.R = None
        self.T = None
        self.Tin = None
        self.Tout = None
        self.CalcRmat()
        self.CalcTmat()

    def SetOffset(self, offset):
        """Set the offset parameter of the element"""
        if type(offset) == np.ndarray:
            self.offset = offset
        elif type(offset) == list:
            self.offset = np.array(offset)

    def GetOffset(self) :
        """Return the offset parameter"""
        return self.offset.copy()

    def SetL(self, L) :
        """Set the length, L, of the element"""
        self.L = L

    def GetL(self) :
        """Return the length, L, of the element"""
        return self.L 

    def SetP(self, P) :
        """Set the design momentum of the element"""
        self.P = P

    def GetP(self) :
        """Return the design momentum of the element"""
        return self.P

    def SetS(self, S) :
        """Set the S location of the element"""
        self.S = S

    def GetS(self) :
        """Return the S location of the element"""
        return self.S
 
    def __str__(self):
        retstr = ""
        methodwrappertype = type(object().__hash__)
        print "Element type = " + self.__class__.__name__
        for i in dir(self):
            iattr = getattr(self, i)
            if (isinstance(iattr, types.BuiltinMethodType) or 
                    isinstance(iattr, methodwrappertype) or
                    isinstance(iattr, types.MethodType) or
                    isinstance(iattr, types.FunctionType) or
                    isinstance(iattr, types.TypeType)):
                continue
            elif (cb.is_numlike(iattr) and 
                i != "R" and 
                i != "T" and 
                i != "Tin" and 
                i != "Tout" and 
                not i.startswith("_")):
                retstr += i+' = '+str(iattr)+"\n"
            elif cb.is_string_like(iattr) and not i.startswith("_"):
                retstr += i+' = '+iattr+"\n"
        return retstr

    def CalcRmat(self, P=None):
        """Calculate and set the R matrix of the element.
        Assumes a drift matrix for this default element."""
        self.R = DriftRmat(self.L)

    def CalcTmat(self, P=None):
        """Calculate and set the T matrix of the element.
        Assumes a drift matrix for this default element."""
        self.T = DriftTmat()

    def TrackThruEle(self, beam_in):
        """Tracks a beam, beam_in, through the element, and returns the
        output bea,."""
        beam_out = deepcopy(beam_in)
        beam_out.x = np.dot(self.R, beam_in.x)
        return beam_out

# ===============================================================
# 'Element' has a subclass for each of the main element types
#  Magnets, cavities (not BPMs), diagnostics, and drifts & collimators
class BasicMag(Element):
    """A class describing the operation of magnets.  Inherits from Element."""
    def __init__(self, name, L=1, P=1, S=0, aper=0, apershape='ELLIPTICAL',
        is_diag=False, B=0, tilt=0):
        self.B = B
        self.tilt = tilt
        Element.__init__(self, name, L, P, S, aper, apershape, is_diag)

    def SetB(self, B) :
        """Sets the B field of the element"""
        self.B = B

    def GetB(self) :
        """Returns the B field of the element"""
        return self.B

    def SetTilt(self, tilt) :
        """Sets the tilt of the element"""
        self.tilt = tilt

    def GetTilt(self) :
        """Returns the tilt of the element"""
        return self.tilt
    
    def TrackThruEle(self, beam_in):
        """Tracks a beam through the element"""
        [r_in, r_out] = RotMats(-self.tilt)
        beam_out = deepcopy(beam_in)
        if self.is_diag:
            self.DiagOut.centroid = np.mean(beam_out.x, 1)
        beam_out.x = np.dot(r_in, beam_out.x)
        doflist = range(6)
        zeromat = np.zeros((6, 6, 6))
        for partnum in range(beam_out.x.shape[1]):
            self.CalcRmat((beam_out.x[5, partnum] * self.P) + self.P)
            self.CalcTmat((beam_out.x[5, partnum] * self.P) + self.P)

            beam_out.x[:, partnum] = np.dot(self.R, beam_out.x[:, partnum])

            if not np.allclose(self.T, zeromat):
                for dof in doflist:
                    for Tdof1 in doflist:
                        for Tdof2 in doflist:
                            beam_out.x[dof, partnum] += \
                                self.T[dof, Tdof1, Tdof2] * \
                                beam_in.x[Tdof1, partnum] * \
                                beam_in.x[Tdof2, partnum]

        self.CalcRmat()
        self.CalcTmat()
        beam_out.x = np.dot(r_out, beam_out.x)
        return beam_out

class BasicCav(Element):
    """A class describing RF cavities.  Inherits from Element"""
    def __init__(self, name, L=1, P=1, S=0, aper=0, apershape='ELLIPTICAL', 
        is_diag=False, phi=0, freq=1e9, numdrift=2, egain=0, designQ=1e9, 
        slices=1, loading=0, restmass=electron_mass):

        self.P        = P
        self.phi      = phi
        self.restmass = restmass
        self.freq     = freq
        self.egain    = egain
        self.numdrift = numdrift
        self.energy0  = np.sqrt((P*1e9)**2 + self.restmass**2)
        self.gamma0   = self.energy0/self.restmass
        self.beta0    = np.sqrt(1 - 1/(self.gamma0**2))
        self.alpha    = self.egain / (P*1e9)
        self.k        = self.freq*2*np.pi/c_light
        self.designQ  = designQ
        self.slices   = slices
        self.loading  = loading
        Element.__init__(self, name, L, P, S, aper, apershape, is_diag)
    
class BasicDiag(Element):
    """A class describing diagnostics .  Inherits from Element"""
    def __init__(self, name, L=1, P=1, S=0, aper=0, apershape='ELLIPTICAL',
        is_diag=True, res=np.array([0, 0])):

        self.DiagOut = DiagOut()

        if type(res) == int or type(res) == float or type(res) == np.double:
            self.res = np.array([res, res])
        else:
            self.res = res
        Element.__init__(self, name, L, P, S, aper, apershape, is_diag)
    
    def Processor(self, beam_in):
        """The (empty) processor method for BasicDiag"""
        pass

class DiagOut:
    """A class to gather all the diagnostics data"""
    def __init__(self):
        self.x_centroid, self.xp_centroid = None, None
        self.y_centroid, self.yp_centroid = None, None
        self.centroid  = None
        self.beam_dist = None
        self.S_pos     = None
        self.charge    = None

class BasicDrift(Element):
    """A class describing drifts .  Inherits from Element.
    Will be used in the future to describe things like collimators"""
    pass

# ===============================================================
# Each of the magnet types is a subclass of BasicMag
class Xcor(BasicMag):
    """X corrector class"""
    def TrackThruEle(self, beam_in):
        """Track a beam through the element"""
        beam_out = deepcopy(beam_in)
        if self.is_diag:
            self.DiagOut.centroid = np.mean(beam_out.x, 1)
        beam_out.x = np.dot(DriftRmat(self.L/2.), beam_out.x)
        for partnum in range(beam_in.x.shape[1]):
            momentum = (beam_in.x[5, partnum] * self.P) + self.P
            Brho = Brho1GeV * momentum
            theta = self.B / Brho
            beam_out.x[1, partnum] += theta
        beam_out.x = np.dot(DriftRmat(self.L/2.), beam_out.x)
        return beam_out

class Ycor(BasicMag):
    """Y corrector class"""
    def TrackThruEle(self, beam_in):
        """Track a beam through the element"""
        beam_out = deepcopy(beam_in)
        if self.is_diag:
            self.DiagOut.centroid = np.mean(beam_out.x, 1)
        beam_out.x = np.dot(DriftRmat(self.L/2.), beam_out.x)
        for partnum in range(beam_in.x.shape[1]):
            momentum = (beam_in.x[5, partnum] * self.P) + self.P
            Brho = Brho1GeV * momentum
            theta = self.B / Brho
            beam_out.x[3, partnum] += theta
        beam_out.x = np.dot(DriftRmat(self.L/2.), beam_out.x)
        return beam_out

class XYcor(BasicMag):
    """XY corrector class"""
    pass

class Quad(BasicMag):
    """Quad class"""
    def CalcRmat(self, P=-1):
        """Calculate the R matrix of the element"""
        if P == -1:
            P = self.P
        L = self.L
        if self.B == 0:
            self.R = DriftRmat(L)
        elif self.B < 0:
            Brho = Brho1GeV * P
            K = np.sqrt(-(self.B/L) / Brho)
            CH = np.cosh(K*L)
            SH = np.sinh(K*L)
            C = np.cos(K*L)
            S = np.sin(K*L)
            if hasattr(self, 'R') and not self.R==None:
                self.R -= self.R
                self.R[0, 0] = CH
                self.R[0, 1] = SH/K
                self.R[1, 0] = K*SH
                self.R[1, 1] = CH
                self.R[2, 2] = C
                self.R[2, 3] = S/K
                self.R[3, 2] = -K*S
                self.R[3, 3] = C
            else:
                self.R = np.array([
                    [CH, SH/K, 0, 0, 0, 0],
                    [K*SH, CH, 0, 0, 0, 0],
                    [0, 0, C, S/K, 0, 0],
                    [0, 0, -K*S, C, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0, 1],
                    ])
        elif self.B > 0:
            Brho = Brho1GeV * P
            K = np.sqrt((self.B/L) / Brho)
            CH = np.cosh(K*L)
            SH = np.sinh(K*L)
            C = np.cos(K*L)
            S = np.sin(K*L)
            if hasattr(self, 'R') and not self.R==None:
                self.R -= self.R
                self.R[0, 0] = C
                self.R[0, 1] = S/K
                self.R[1, 0] = -K*S
                self.R[1, 1] = C
                self.R[2, 2] = CH
                self.R[2, 3] = SH/K
                self.R[3, 2] = K*SH
                self.R[3, 3] = CH
            else:
                self.R = np.array([
                    [C, S/K, 0, 0, 0, 0],
                    [-K*S, C, 0, 0, 0, 0],
                    [0, 0, CH, SH/K, 0, 0],
                    [0, 0, K*SH, CH, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0, 1],
                    ])

class ThinSext(BasicMag):
    """Thin sextupole class"""
    def CalcTmat(self, P=-1):
        """Calculate the T matrix for the element"""
        if P == -1:
            P = self.P
        self.T = np.zeros((6, 6, 6))
        Brho = Brho1GeV * P
        ksql = self.B / (Brho * 2)  # 2 is mysterious, but used in Lucretia!
        k2l2 = ksql * self.L
        k2l3 = k2l2 * self.L
        k2l4 = k2l3 * self.L

        self.T[0, 0, 0] = -0.5 * k2l2
        self.T[0, 0, 1] = (-1.0/3.0) * k2l3
        self.T[0, 1, 1] = (-1.0/12.0) * k2l4
        self.T[0, 2, 2] = 0.5 * k2l2
        self.T[0, 2, 3] = (1.0/3.0) * k2l3
        self.T[0, 3, 3] = (1.0/12.0) * k2l4

        self.T[1, 0, 0] = -ksql
        self.T[1, 0, 1] = -k2l2
        self.T[1, 1, 1] = -(1.0/3.0) * k2l3
        self.T[1, 2, 2] = ksql
        self.T[1, 2, 3] = k2l2
        self.T[1, 3, 3] = (1.0/3.0) * k2l3

        self.T[2, 0, 2] = k2l2
        self.T[2, 0, 3] = (1.0/3.0) * k2l3
        self.T[2, 1, 2] = (1.0/3.0) * k2l3
        self.T[2, 1, 3] = (1.0/6.0) * k2l4

        self.T[3, 0, 2] = 2.0 * ksql
        self.T[3, 0, 3] = k2l2
        self.T[3, 1, 2] = k2l2
        self.T[3, 1, 3] = (2.0/3.0) * k2l3

class Sext(ThinSext):
    """Sextupole class"""
    pass

class Oct(BasicMag):
    """Octupole class"""
    pass

class Solenoid(BasicMag):
    """Solenoid class"""
    pass

class Sbend(BasicMag):
    """Sector bend class"""
    def __init__(self, name, L=1, P=1, S=0, aper=0, apershape='ELLIPTICAL',
        is_diag=False, B=np.array([0, 0]), e_angle=0, e_curve=0, h_gap=0, 
        h_int=0):
        self.e_angle = SplitParams(e_angle)
        self.e_curve = SplitParams(e_curve)
        self.h_gap = SplitParams(h_gap)
        self.h_int = SplitParams(h_int)
        if np.array(B).size == 1:
            B = np.array([B, 0])
        BasicMag.__init__(self, name, L, P, S, aper, apershape, is_diag, B)
        self.B = B
        #try:
        #    if e_angle.shape[0] == 1:
        #        e_angle = np.array([e_angle, e_angle])
        #except AttributeError:
        #    e_angle = np.array([e_angle, e_angle])
        #try:
        #    if e_curve.shape[0] == 1:
        #        e_curve = np.array([e_curve, e_curve])
        #except AttributeError:
        #    e_curve = np.array([e_curve, e_curve])
        #try:
        #    if h_gap.shape[0] == 1:
        #        h_gap = np.array([h_gap, h_gap])
        #except AttributeError:
        #    h_gap = np.array([h_gap, h_gap])
        #try:
        #    if h_int.shape[0] == 1:
        #        h_int = np.array([h_int, h_int])
        #except AttributeError:
        #    h_int = np.array([h_int, h_int])

    def CalcRmat(self, P=-1):
        """Calculate the R matrix for the element"""
        if P == -1:
            P = self.P
        rho = Brho1GeV * P * self.L / self.B[0]
        h = 1 / rho
        n = -(self.B[1]/self.L) / (h*(self.B[0]/self.L))
        kxsq = (1-n)*h**2
        kysq = n*h**2
        kx = np.sqrt( abs(kxsq) )
        ky = np.sqrt( abs(kysq) )
        if kxsq > 0:
            Cx = np.cos(kx * self.L)
            Sx = np.sin(kx * self.L)
            Sxoverkx = Sx/kx
            signkx = 1
        elif kxsq < 0:
            Cx = np.cosh(kx * self.L)
            Sx = np.sinh(kx * self.L)
            Sxoverkx = Sx/kx
            signkx = -1
        else:
            Cx = 1
            Sx = 0
            Sxoverkx = self.L
            signkx = 0

        if kysq > 0:
            Cy = np.cos(ky * self.L)
            Sy = np.sin(ky * self.L)
            signky = 1
            Syoverky = Sy/ky
        elif kysq < 0:
            Cy = np.cosh(ky * self.L)
            Sy = np.sinh(ky * self.L)
            signky = -1
            Syoverky = Sy/ky
        else:
            Cy = 1
            Sy = 0
            Syoverky = self.L     
            signky = 0
            #gamma0 = self.gamma0
            #beta0 = self.beta0
            #alpha = self.egain / P
            #phi_perp = sqrt(abs(pi * alpha * cos(phi)) / 2.0)
            #phi_para = sqrt(abs(pi * alpha * cos(phi))) / (gamma0 * beta0)
            #C_perp = cos(phi_perp)
            #S_perp = sin(phi_perp)
            #C_para = cos(phi_para)
            #S_para = sin(phi_para)

        main_r = np.array([
            [Cx, Sxoverkx, 0, 0, 0, (h/kx**2)*(1-Cx)],
            [-signkx*kx*Sx, Cx, 0, 0, 0, (h/kx)*Sx],
            [0, 0, Cy, Syoverky, 0, 0],
            [0, 0, -signky*ky*Sy, Cy, 0, 0],
            [(h/kx)*Sx, (h/kx**2)*(1-Cx), 0, 0, 1, 
                (h**2/kx**3)*(kx*self.L - Sx)],
            [0, 0, 0, 0, 0, 1],
            ])
        phi_in = \
            self.h_int[0] * \
            self.h_gap[0] * h * \
            (1+(np.sin(self.e_angle[0])**2)) / np.cos(self.e_angle[0])
        phi_out = \
            self.h_int[1] * \
            self.h_gap[1] * h * \
            (1+(np.sin(self.e_angle[1])**2)) / np.cos(self.e_angle[1])
        edge_rin = np.eye(6)
        edge_rin[1, 0] = h * np.tan(self.e_angle[0])
        edge_rin[3, 2] = -h * np.tan(self.e_angle[0] - phi_in)
        edge_rout = np.eye(6)
        edge_rout[1, 0] = h * np.tan(self.e_angle[1])
        edge_rout[3, 2] = -h * np.tan(self.e_angle[1] - phi_out)

        self.R = np.dot(edge_rout, np.dot(main_r, edge_rin))

    def CalcTmat(self, P=1):
        """Calculate the T matrix for the element"""
        self.T = DriftTmat()
        self.Tin = DriftTmat()
        self.Tout = DriftTmat()

        rho = Brho1GeV * P * self.L / self.B[0]
        h = 1/rho
        n = -self.B[1] / (h*(self.B[0]/self.L))
        phi_in = \
            self.h_int[0] * \
            self.h_gap[0] * h * \
            (1+(np.sin(self.e_angle[0])**2)) / np.cos(self.e_angle[0])
        phi_out = \
            self.h_int[1] * \
            self.h_gap[1] * h * \
            (1+(np.sin(self.e_angle[1])**2)) / np.cos(self.e_angle[1])

        curvature = [0, 0]
        if self.e_curve[0] == 0:
            curvature[0] = np.Inf
        else:
            curvature[0] = 1.0/self.e_curve[0]

        if self.e_curve[1] == 0:
            curvature[1] = np.Inf
        else:
            curvature[1] = 1.0/self.e_curve[1]

        self.Tin[0, 0, 0] = -(h/2.0) * np.tan(self.e_angle[0])**2
        self.Tin[0, 2, 2] = (h/2.0) / np.cos(self.e_angle[0])**2
        self.Tin[1, 0, 0] = (h/(2.0*curvature[0])) / np.cos(self.e_angle[0]**3)
        self.Tin[1, 0, 1] = h * np.tan(self.e_angle[0])**2
        self.Tin[1, 0, 5] = -h * np.tan(self.e_angle[0])
        self.Tin[1, 2, 2] = h**2 * \
            (0.5 + np.tan(self.e_angle[0])**2) * np.tan(self.e_angle[0]) - \
            (h/(2*curvature[0]))/np.cos(self.e_angle[0]**3)
        self.Tin[1, 2, 3] = -h * np.tan(self.e_angle[0])**2
        self.Tin[2, 0, 2] = h * np.tan(self.e_angle[0])**2
        self.Tin[3, 0, 2] = -(h/curvature[0])/np.cos(self.e_angle[0])**3
        self.Tin[3, 0, 3] = -h * np.tan(self.e_angle[0])**2
        self.Tin[3, 1, 2] = -h / np.cos(self.e_angle[0])**2
        self.Tin[3, 2, 5] = h * np.tan(self.e_angle[0]) - \
            h * phi_in/np.cos(self.e_angle[0]-phi_in)**2
        self.Tin[0, 0, 0] = (h/2.0) * np.tan(self.e_angle[1])**2
        self.Tin[0, 2, 2] = -(h/2.0) / np.cos(self.e_angle[1])**2
        self.Tin[1, 0, 0] = (h/(2.0*curvature[1])) / \
            np.cos(self.e_angle[1]**3) - \
            h**2 * (0.5*np.tan(self.e_angle[1])**2)*np.tan(self.e_angle[1])
        self.Tin[1, 0, 1] = -h * np.tan(self.e_angle[1])**2
        self.Tin[1, 0, 5] = -h * np.tan(self.e_angle[1])
        self.Tin[1, 2, 2] = h**2 * \
            (-0.5*np.tan(self.e_angle[1])**2) * np.tan(self.e_angle[1]) - \
            (h/(2*curvature[1]))/np.cos(self.e_angle[1]**3)
        self.Tin[1, 2, 3] = h * np.tan(self.e_angle[1])**2
        self.Tin[2, 0, 2] = -h * np.tan(self.e_angle[1])**2
        self.Tin[3, 0, 2] = -(h/curvature[1])/np.cos(self.e_angle[1])**3 + \
            h**2 * np.tan(self.e_angle[1]) / np.cos(self.e_angle[1])**2
        self.Tin[3, 0, 3] = h * np.tan(self.e_angle[1])**2
        self.Tin[3, 1, 2] = h / np.cos(self.e_angle[1])**2
        self.Tin[3, 2, 5] = h * np.tan(self.e_angle[1]) - \
            h * phi_out/np.cos(self.e_angle[1]-phi_out)**2

    def TrackThruEle(self, beam_in):
        """Track a beam through the element"""
        [r_in, r_out] = RotMats(-self.tilt)
        beam_out = deepcopy(beam_in)
        if self.is_diag:
            self.DiagOut.centroid = np.mean(beam_out.x, 1)
        beam_out.x = np.dot(r_in, beam_out.x)
        doflist = range(6)
        zeromat = np.zeros((6, 6, 6))

        for partnum in range(beam_out.x.shape[1]):
            self.CalcTmat((beam_out.x[5, partnum] * self.P) + self.P)
            if not np.allclose(self.Tin, zeromat):
                for dof in doflist:
                    for Tdof1 in doflist:
                        for Tdof2 in doflist:
                            beam_out.x[dof, partnum] += \
                                self.Tin[dof, Tdof1, Tdof2] * \
                                beam_in.x[Tdof1, partnum] * \
                                beam_in.x[Tdof2, partnum]
        
        for partnum in range(beam_out.x.shape[1]):
            self.CalcRmat((beam_out.x[5, partnum] * self.P) + self.P)
            beam_out.x[:, partnum] = np.dot(self.R, beam_out.x[:, partnum])

        for partnum in range(beam_out.x.shape[1]):
            self.CalcTmat((beam_out.x[5, partnum] * self.P) + self.P)
            if not np.allclose(self.Tout, zeromat):
                for dof in doflist:
                    for Tdof1 in doflist:
                        for Tdof2 in doflist:
                            beam_out.x[dof, partnum] += \
                                self.Tout[dof, Tdof1, Tdof2] * \
                                beam_in.x[Tdof1, partnum] * \
                                beam_in.x[Tdof2, partnum]

        self.CalcRmat()
        self.CalcTmat()
        beam_out.x = np.dot(r_out, beam_out.x)

        return beam_out

# ===============================================================
# Each cavity type is a subclass of BasicCav
class AccCav(BasicCav):
    """Accelerator cavity class"""
    def CalcRmat(self, P=-1, L=-1, z=0):
        """Calculate the R matrix for the element"""
        if P == -1:
            P = self.P*1e9         
        if L == -1:
            L = self.L
        if z == 0:
            phi = self.phi
        else:
            rflambda = c_light/self.freq
            phi = self.phi + (z/rflambda)*2*np.pi 

        gamma0  = self.gamma0
        beta0   = self.beta0
        alpha = self.alpha
        phi_perp = np.sqrt(abs(np.pi * alpha * np.cos(phi)) / 2.0)
        phi_para = np.sqrt(abs(np.pi * alpha * np.cos(phi))) / (gamma0 * beta0)
        if np.cos(phi) >= 0:
            C_perp = np.cos(phi_perp)
            S_perp = np.sin(phi_perp)
            C_para = np.cos(phi_para)
            S_para = np.sin(phi_para)
        else:
            C_perp = np.cosh(phi_perp)
            S_perp = np.sinh(phi_perp)
            C_para = np.cosh(phi_para)
            S_para = np.sinh(phi_para)
    
        if not hasattr(self, 'R'):
            self.R = np.zeros((6, 6))

        self.R[0, 0] = C_perp
        self.R[0, 1] = (L/phi_perp) * S_perp
        self.R[1, 0] = -(phi_perp/L) * S_perp
        self.R[1, 1] = C_perp

        self.R[2, 2] = C_perp
        self.R[2, 3] = (L/phi_perp) * S_perp
        self.R[3, 2] = -(phi_perp/L) * S_perp
        self.R[3, 3] = C_perp

        self.R[4, 4] = C_para
        self.R[4, 5] = (1/(beta0**2 * gamma0**2)) * (L/phi_para) * S_para
        self.R[5, 4] = -(beta0**2 * gamma0**2) * (phi_para/L) * S_para
        self.R[5, 5] = C_para

    def DriftMap(self, beam_in, DistDrift, P):
        """Calculate the mapping for the drift section of the Lie map"""
        beam_out = deepcopy(beam_in)
     
        x = beam_in.x[0]
        P_x = beam_in.x[1]
        y = beam_in.x[2]
        P_y = beam_in.x[3]
        z = beam_in.x[4]
        delta = beam_in.x[5]
        
        L = self.L*DistDrift
        energy0 = np.sqrt((P*1e9)**2 + beam_in.restmass**2)
        gamma0 = energy0/beam_in.restmass
        beta0 = np.sqrt(1 - 1/(gamma0**2))
 
        P_s = np.sqrt((((1/beta0)+delta)**2) - 
            P_x**2 - P_y**2 - 
            (1/((beta0**2)*(gamma0**2))))
        
        beam_out.x[0] = x + L*(P_x/P_s)
        beam_out.x[2] = y + L*(P_y/P_s)
        beam_out.x[4] = z + L*((1/beta0)-((delta+(1/beta0))/P_s))
        return beam_out

    def KickMap(self, beam_in, DistKick, P):
        """Calculate the mapping for the kick section of the Lie map"""
        beam_out = deepcopy(beam_in)
        
        x = beam_in.x[0]
        P_x = beam_in.x[1]
        y = beam_in.x[2]
        P_y = beam_in.x[3]
        z = beam_in.x[4]
        delta = beam_in.x[5]
        
        rho = np.sqrt(x**2 + y**2)
        phi0 = self.phi
        k = self.k
        A = DistKick*(self.egain/(self.P*1e9))

        for partnum in xrange(beam_out.x.shape[1]):
            if (rho[partnum] != 0):
                beam_out.x[1] = P_x - A * (x/rho) * \
                    jn(1, (k*rho))*np.cos(phi0-(k*z))
                beam_out.x[3] = P_y - A * (y/rho) * \
                    jn(1, (k*rho))*np.cos(phi0-(k*z))
                beam_out.x[5] = delta + A*jn(0, (k*rho))*np.sin(phi0-(k*z))
            else:
                beam_out.x[5] = delta + A*jn(0, (k*rho))*np.sin(phi0-(k*z))

        P = P + A*jn(0, 0)*np.sin(phi0)
        
        return beam_out, P

    def TrackThruEle(self, beam_in):
        """Track a beam through the element"""
        intermed_beam = deepcopy(beam_in)
        try:
            DistDrift = lietrackarray[self.numdrift].dlengths
            DistKick = lietrackarray[self.numdrift].klengths    
        except IndexError:
            raise IndexError(
            'Yoshida factorisation for numdrift=%i is not defined' 
            % self.numdrift)
        except AttributeError:
            raise AttributeError(
            'Yoshida factorisation for numdrift=%i is not defined' 
            % self.numdrift)

        N = 0
        M = 0
        P = self.P
        beam_out = self.DriftMap(intermed_beam, DistDrift[N], P)
     #   beam_out = self.KickMap(beam_out, DistKick[M])
        N += 1
        while N < self.numdrift:
            beam_out, P = self.KickMap(beam_out, DistKick[M], P)
            beam_out = self.DriftMap(beam_out, DistDrift[N], P)
            N += 1
            M += 1
    
        return beam_out

class TCav(BasicCav):
    """Transverse cavity class"""
    pass

# ===============================================================
# Each diagnostic type is a subclass of BasicDiag
class BPM(BasicDiag):
    """BPM class"""
    def Processor(self, beam_in):
        """Processor class for the BPM"""
        self.DiagOut.S_pos = self.S
        if self.res[0] == 0:
            self.DiagOut.x_centroid  = np.mean(beam_in.x[0, :])
            self.DiagOut.y_centroid  = np.mean(beam_in.x[2, :])
        else:
            self.DiagOut.x_centroid = \
                np.mean(beam_in.x[0, :]) + normal(scale=self.res[0], size=1)
            self.DiagOut.y_centroid = \
                np.mean(beam_in.x[2, :]) + normal(scale=self.res[0], size=1)
        if self.res[1] == 0:
            self.DiagOut.xp_centroid = np.mean(beam_in.x[1, :])
            self.DiagOut.yp_centroid = np.mean(beam_in.x[3, :])
        else:
            self.DiagOut.xp_centroid = \
                np.mean(beam_in.x[1, :]) + normal(scale=self.res[1], size=1)
            self.DiagOut.yp_centroid = \
                np.mean(beam_in.x[3, :]) + normal(scale=self.res[1], size=1)

class Screen(BasicDiag):
    """Screen class"""
    def Processor(self, beam_in):
        """Processor class for the Screen"""
        self.DiagOut.S_pos = self.S
        self.DiagOut.beam_dist = beam_in.x[0:3:2, :]

class EmitScreen(BasicDiag):
    """EmitScreen class"""
    def Processor(self, beam_in):
        """Processor class for the EmitScreen"""
        self.DiagOut.S_pos = self.S
        self.DiagOut.beam_dist = beam_in.x

class WireScanner(BasicDiag):
    """WireScanner class"""
    pass

class OTR(Screen):
    """OTR class"""
    pass

class ICT(BasicDiag):
    """ICT class"""
    def Processor(self, beam_in):
        """Processor class for the ICT"""
        self.DiagOut.S_pos = self.S
        origparts, numparts = 0, 0
        for x in beam_in.x[0, :]:
            origparts += 1
            if not np.isnan(x):
                numparts += 1
        self.DiagOut.charge = beam_in.Q * (numparts / origparts )


# ===============================================================
# Each drift/collimator type is a subclass of BasicDrift
class Drift(BasicDrift):
    """Drift class"""
    pass
    
class Coll(BasicDrift):
    """Collimator class"""
    pass

# ===============================================================
# A mover class
class Mover:
    """Mover class"""
    def __init__(self, dof = np.array([0, 2, 5])):
        self.dof = dof
        self.act = np.zeros(len(dof))
        self.des = np.zeros(len(dof))

    def Set(self, pos):
        """Set the mover position"""
        if not len(pos) == len(self.dof):
            raise "Mover degrees of freedom do not match request"
        self.des = np.array(pos).__copy__()

    def Trim(self):
        """"Trim" to the required position"""
        self.act = self.des.__copy__()

    def GetAct(self):
        """Return the current mover position"""
        act = np.zeros(6)
        for i in range(len(self.dof)):
            act[self.dof[i]] = self.act[i]
        return act

class Wake:
    """A class to describe wakefields"""
    def __init__(self):
        self.l = None
        self.k = None

# ===============================================================
## A few useful functions
#def DriftRmat(L):
#    """Return the R matrix for a drift"""
#    R = np.array([
#        [1, L, 0, 0, 0, 0],
#        [0, 1, 0, 0, 0, 0],
#        [0, 0, 1, L, 0, 0],
#        [0, 0, 0, 1, 0, 0],
#        [0, 0, 0, 0, 1, 0],
#        [0, 0, 0, 0, 0, 1],
#        ])
#    return R
#
#def DriftTmat():
#    """Return the T matrix for a drift"""
#    T = np.zeros((6, 6, 6))
#    return T
#
#def SplitParams(param):
#    """Split a parameter into a 1x2 np.array"""
#    outparam = np.array([0.0, 0.0])
#    if np.array(param).size == 1:
#        outparam[0] = param
#        outparam[1] = param
#    elif np.array(param).size > 1:
#        outparam[0] = param[0]
#        outparam[1] = param[1]
#    return outparam
#
#def RotMats(alpha):
#    """Calculate the matrix required to rotate the beam by a
#    particular angle"""
#    C = np.cos(-alpha)
#    S = np.sin(-alpha)
#    r_inrot = np.array([
#        [C , 0 , S , 0 , 0 , 0],
#        [0 , C , 0 , S , 0 , 0],
#        [-S, 0 , C , 0 , 0 , 0],
#        [0 , -S, 0 , C , 0 , 0],
#        [0 , 0 , 0 , 0 , 1 , 0],
#        [0 , 0 , 0 , 0 , 0 , 1],
#        ])
#    C = np.cos(alpha)
#    S = np.sin(alpha)
#    r_outrot = np.array([
#        [C , 0 , S , 0 , 0 , 0],
#        [0 , C , 0 , S , 0 , 0],
#        [-S, 0 , C , 0 , 0 , 0],
#        [0 , -S, 0 , C , 0 , 0],
#        [0 , 0 , 0 , 0 , 1 , 0],
#        [0 , 0 , 0 , 0 , 0 , 1],
#        ])
#    return r_inrot, r_outrot

# ===============================================================
# Test suite
if __name__ == '__main__':
    from latticeloader import loadflatlatfile

    templine = BL.Line(loadflatlatfile('tempfile.txt')[1155:])
    mybeam = beamrep.Beam(x=np.array([0, 0, 0, 0, 0, 1.3]))
    beamout = templine.Track(mybeam)
    print beamout.x

