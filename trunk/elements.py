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
# A module defining the Element superclass, and
# individual subclasses for each element type.

from numpy import *
from matplotlib.pyplot import figure, plot
from copy import deepcopy
import accformat
from globals import *
import matplotlib.cbook as cb
import types
from numpy.random import normal
from globals import electron_mass, c_light, e_charge


# Twiss object class for addition to the Element superclass
class Twiss:
    def __init__(self):
        self.betax = None
        self.betay = None

        self.alphax = None
        self.alphay = None

        self.etax = None
        self.etay = None

        self.etapx = None
        self.etapy = None

        self.phix = None
        self.phiy = None

        self.nemitx = None
        self.nemity = None

        self.sigz = None
        self.sigP = None
        self.pz_cor = None

# 'Element' superclass, defining common qualities of all elements
# Also defines CalcRmat for a drift, which will then be overridden
# in classes where this is not true.
class Element:
    def __init__(self,name,L=1,P=1,S=0,aper=0,apershape='ELLIPTICAL',is_diag=False):
        self.name = name
        self.L = L
        self.P = P
        self.S = S
        self.aper = aper
        self.offset = array([0.0,0.0,0.0,0.0,0.0,0.0])
        self.is_diag = is_diag
        if self.is_diag:
            self.DiagOut = DiagOut()
        self.twiss = Twiss()
        self.CalcRmat()
        self.CalcTmat()

    def SetOffset(self, o) :
        if type(o) == array:
            self.offset = o
        elif type(o) == list:
            self.offset = array(o)

    def GetOffset(self) :
        return self.offset.copy()

    def SetL(self, L) :
        self.L = L

    def GetL(self) :
        return self.L 

    def SetP(self, P) :
        self.P = P

    def GetP(self) :
        return self.P

    def SetS(self,S) :
        self.S = S

    def GetS(self) :
        return S
 
    def __str__(self):
        retstr = ""
        MethodWrapperType = type(object().__hash__)
        print "Element type = " + self.__class__.__name__
        for i in dir(self):
            iattr = getattr(self,i)
            if (isinstance(iattr,types.BuiltinMethodType) or 
                    isinstance(iattr,MethodWrapperType) or
                    isinstance(iattr,types.MethodType) or
                    isinstance(iattr,types.FunctionType) or
                    isinstance(iattr,types.TypeType)):
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

    def CalcRmat(self,P=-1):
        self.R = DriftRmat(self.L)

    def CalcTmat(self,P=-1):
        self.T = DriftTmat()

    def TrackThruEle(self,beam_in):
        beam_out = deepcopy(beam_in)
        beam_out.x = dot(self.R, beam_in.x)
        return beam_out

# ===============================================================
# 'Element' has a subclass for each of the main element types
#  Magnets, cavities (not BPMs), diagnostics, and drifts & collimators
class BasicMag(Element):
    def __init__(self,name,L=1,P=1,S=0,aper=0,apershape='ELLIPTICAL',
        is_diag=False,B=0,tilt=0):
        self.B = B
        self.tilt = tilt
        Element.__init__(self, name, L, P, S, aper, apershape, is_diag)

    def SetB(self, B) :
        self.B = B
    def GetB(self) :
        return self.B
    def SetTilt(self, Tilt) :
        self.Tilt = Tilt
    def GetTilt(self) :
        return self.Tilt
    
    def TrackThruEle(self,beam_in):
        [r_in,r_out] = RotMats(-self.tilt)
        beam_out = deepcopy(beam_in)
        if self.is_diag:
            self.DiagOut.centroid = mean(beam_out.x,1)
        beam_out.x = dot(r_in,beam_out.x)
        doflist = range(6)
        zeromat = zeros((6,6,6))
        for partnum in range(beam_out.x.shape[1]):
            self.CalcRmat((beam_out.x[5,partnum] * self.P) + self.P)
            self.CalcTmat((beam_out.x[5,partnum] * self.P) + self.P)

            beam_out.x[:,partnum] = dot(self.R,beam_out.x[:,partnum])

            if not all(self.T == zeromat):
                for dof in doflist:
                    for Tdof1 in doflist:
                        for Tdof2 in doflist:
                            beam_out.x[dof,partnum] += \
                                self.T[dof,Tdof1,Tdof2] * \
                                beam_in.x[Tdof1,partnum] * \
                                beam_in.x[Tdof2,partnum]

        self.CalcRmat()
        self.CalcTmat()
        beam_out.x = dot(r_out,beam_out.x)
        return beam_out

class BasicCav(Element):
    def __init__(self,name,L=1,P=1,S=0,aper=0,apershape='ELLIPTICAL',is_diag=False,phi=0,
        freq=1e9,numdrift=2,egain=0,designQ=1e9,slices=1,loading=0,restmass=electron_mass):

        self.P        = P
        self.phi      = phi
        self.restmass = restmass
        self.freq     = freq
        self.egain    = egain
        self.numdrift = numdrift
        self.energy0  = sqrt((P*1e9)**2 + self.restmass**2)
        self.gamma0   = self.energy0/self.restmass
        self.beta0    = sqrt(1 - 1/(self.gamma0**2))
        self.alpha    = self.egain / (P*1e9)
        self.k        = self.freq*2*pi/c_light
        self.designQ  = designQ
        self.slices   = slices
        self.loading  = loading
        Element.__init__(self, name, L, P, S, aper, apershape, is_diag)
    
class BasicDiag(Element):
    def __init__(self,name,L=1,P=1,S=0,aper=0,apershape='ELLIPTICAL',
        is_diag=True,res=array([0,0])):
   
        if type(res)==int or type(res)==float or type(res)==double:
            self.res = array([res, res])
        else:
            self.res = res
        Element.__init__(self, name, L, P, S, aper, apershape, is_diag)
    
    def Processor(self, beam_in):
        pass

class DiagOut:
    pass

class BasicDrift(Element):
    pass

# ===============================================================
# Each of the magnet types is a subclass of BasicMag
class Xcor(BasicMag):
    def TrackThruEle(self,beam_in):
        beam_out = deepcopy(beam_in)
        if self.is_diag:
            self.DiagOut.centroid = mean(beam_out.x,1)
        beam_out.x = dot(DriftRmat(self.L/2.),beam_out.x)
        for partnum in range(beam_in.x.shape[1]):
            momentum = (beam_in.x[5,partnum] * self.P) + self.P
            Brho = Brho1GeV * momentum
            theta = self.B / Brho
            beam_out.x[1,partnum] += theta
        beam_out.x = dot(DriftRmat(self.L/2.),beam_out.x)
        return beam_out

class Ycor(BasicMag):
    def TrackThruEle(self,beam_in):
        beam_out = deepcopy(beam_in)
        if self.is_diag:
            self.DiagOut.centroid = mean(beam_out.x,1)
        beam_out.x = dot(DriftRmat(self.L/2.),beam_out.x)
        for partnum in range(beam_in.x.shape[1]):
            momentum = (beam_in.x[5,partnum] * self.P) + self.P
            Brho = Brho1GeV * momentum
            theta = self.B / Brho
            beam_out.x[3,partnum] += theta
        beam_out.x = dot(DriftRmat(self.L/2.),beam_out.x)
        return beam_out

class XYcor(BasicMag):
    pass

class Quad(BasicMag):
    def CalcRmat(self, P=-1):
        if P==-1:
            P=self.P
        L = self.L
        if self.B==0:
            self.R = DriftRmat(L)
        elif self.B<0:
            Brho = Brho1GeV * P
            K = sqrt(-(self.B/L) / Brho)
            CH = cosh(K*L)
            SH = sinh(K*L)
            C = cos(K*L)
            S = sin(K*L)
            if hasattr(self,'R'):
                self.R -= self.R
                self.R[0,0] = CH
                self.R[0,1] = SH/K
                self.R[1,0] = K*SH
                self.R[1,1] = CH
                self.R[2,2] = C
                self.R[2,3] = S/K
                self.R[3,2] = -K*S
                self.R[3,3] = C
            else:
                self.R = array([
                    [CH,SH/K,0,0,0,0],
                    [K*SH,CH,0,0,0,0],
                    [0,0,C,S/K,0,0],
                    [0,0,-K*S,C,0,0],
                    [0,0,0,0,1,0],
                    [0,0,0,0,0,1],
                    ])
        elif self.B>0:
            Brho = Brho1GeV * P
            K = sqrt((self.B/L) / Brho)
            CH = cosh(K*L)
            SH = sinh(K*L)
            C = cos(K*L)
            S = sin(K*L)
            if hasattr(self,'R'):
                self.R -= self.R
                self.R[0,0] = C
                self.R[0,1] = S/K
                self.R[1,0] = -K*S
                self.R[1,1] = C
                self.R[2,2] = CH
                self.R[2,3] = SH/K
                self.R[3,2] = K*SH
                self.R[3,3] = CH
            else:
                self.R = array([
                    [C,S/K,0,0,0,0],
                    [-K*S,C,0,0,0,0],
                    [0,0,CH,SH/K,0,0],
                    [0,0,K*SH,CH,0,0],
                    [0,0,0,0,1,0],
                    [0,0,0,0,0,1],
                    ])

class ThinSext(BasicMag):
    def CalcTmat(self, P=-1):
        if P==-1:
            P=self.P
        self.T = zeros((6,6,6))
        Brho = Brho1GeV * P
        ksql = self.B / (Brho * 2)  # Factor of 2 is mysterious, but used in Lucretia!
        k2l2 = ksql * self.L
        k2l3 = k2l2 * self.L
        k2l4 = k2l3 * self.L

        self.T[0,0,0] = -0.5 * k2l2
        self.T[0,0,1] = (-1.0/3.0) * k2l3
        self.T[0,1,1] = (-1.0/12.0) * k2l4
        self.T[0,2,2] = 0.5 * k2l2
        self.T[0,2,3] = (1.0/3.0) * k2l3
        self.T[0,3,3] = (1.0/12.0) * k2l4

        self.T[1,0,0] = -ksql
        self.T[1,0,1] = -k2l2
        self.T[1,1,1] = -(1.0/3.0) * k2l3
        self.T[1,2,2] = ksql
        self.T[1,2,3] = k2l2
        self.T[1,3,3] = (1.0/3.0) * k2l3

        self.T[2,0,2] = k2l2
        self.T[2,0,3] = (1.0/3.0) * k2l3
        self.T[2,1,2] = (1.0/3.0) * k2l3
        self.T[2,1,3] = (1.0/6.0) * k2l4

        self.T[3,0,2] = 2.0 * ksql
        self.T[3,0,3] = k2l2
        self.T[3,1,2] = k2l2
        self.T[3,1,3] = (2.0/3.0) * k2l3

class Sext(ThinSext):
    pass

class Oct(BasicMag):
    pass

class Solenoid(BasicMag):
    pass

class Sbend(BasicMag):
    def __init__(self,name,L=1,P=1,S=0,aper=0,apershape='ELLIPTICAL',
        is_diag=False,B=array([0,0]),e_angle=0,e_curve=0,h_gap=0,h_int=0):
        self.e_angle = SplitParams(e_angle)
        self.e_curve = SplitParams(e_curve)
        self.h_gap = S  
        self.h_int = SplitParams(h_int)
        if array(B).size==1:
            B = array([B,0])
        BasicMag.__init__(self, name, L, P, S, aper, apershape, is_diag, B)
        self.B = B
        try:
            if e_angle.shape[0] == 1:
                e_angle = array([e_angle,e_angle])
        except AttributeError:
            e_angle = array([e_angle,e_angle])
        try:
            if e_curve.shape[0] == 1:
                e_curve = array([e_curve,e_curve])
        except AttributeError:
            e_curve = array([e_curve,e_curve])
        try:
            if h_gap.shape[0] == 1:
                h_gap = array([h_gap,h_gap])
        except AttributeError:
            h_gap = array([h_gap,h_gap])
        try:
            if h_int.shape[0] == 1:
                h_int = array([h_int,h_int])
        except AttributeError:
            h_int = array([h_int,h_int])

    def CalcRmat(self, P=-1):
        if P==-1:
            P = self.P
        rho = Brho1GeV * P * self.L / self.B[0]
        h = 1 / rho
        n = -(self.B[1]/self.L) / (h*(self.B[0]/self.L))
        kxsq = (1-n)*h**2
        kysq = n*h**2
        kx = sqrt( abs(kxsq) )
        ky = sqrt( abs(kysq) )
        if kxsq>0:
            Cx = cos(kx * self.L)
            Sx = sin(kx * self.L)
            Sxoverkx = Sx/kx
            signkx = 1
        elif kxsq<0:
            Cx = cosh(kx * self.L)
            Sx = sinh(kx * self.L)
            Sxoverkx = Sx/kx
            signkx = -1
        else:
            Cx = 1
            Sx = 0
            Sxoverkx = self.L
            signkx = 0

        if kysq>0:
            Cy = cos(ky * self.L)
            Sy = sin(ky * self.L)
            signky = 1
            Syoverky = Sy/ky
        elif kysq<0:
            Cy = cosh(ky * self.L)
            Sy = sinh(ky * self.L)
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

        main_r = array([
            [Cx,Sxoverkx,0,0,0,(h/kx**2)*(1-Cx)],
            [-signkx*kx*Sx,Cx,0,0,0,(h/kx)*Sx],
            [0,0,Cy,Syoverky,0,0],
            [0,0,-signky*ky*Sy,Cy,0,0],
            [(h/kx)*Sx,(h/kx**2)*(1-Cx),0,0,1,(h**2/kx**3)*(kx*self.L - Sx)],
            [0,0,0,0,0,1],
            ])
        phi_in = \
            self.h_int[0] * \
            self.h_gap[0] * h * \
            (1+(sin(self.e_angle[0])**2)) / cos(self.e_angle[0])
        phi_out = \
            self.h_int[1] * \
            self.h_gap[1] * h * \
            (1+(sin(self.e_angle[1])**2)) / cos(self.e_angle[1])
        edge_rin = eye(6)
        edge_rin[1,0] = h * tan(self.e_angle[0])
        edge_rin[3,2] = -h * tan(self.e_angle[0] - phi_in)
        edge_rout = eye(6)
        edge_rout[1,0] = h * tan(self.e_angle[1])
        edge_rout[3,2] = -h * tan(self.e_angle[1] - phi_out)

        self.R = dot(edge_rout,dot(main_r,edge_rin))

    def CalcTmat(self, P=1):
        self.T = DriftTmat()
        self.Tin = DriftTmat()
        self.Tout = DriftTmat()

        rho = Brho1GeV * P * self.L / self.B[0]
        h = 1/rho
        n = -self.B[1] / (h*(self.B[0]/self.L))
        phi_in = \
            self.h_int[0] * \
            self.h_gap[0] * h * \
            (1+(sin(self.e_angle[0])**2)) / cos(self.e_angle[0])
        phi_out = \
            self.h_int[1] * \
            self.h_gap[1] * h * \
            (1+(sin(self.e_angle[1])**2)) / cos(self.e_angle[1])

        curvature = [0,0]
        if self.e_curve[0]==0:
            curvature[0] = Inf
        else:
            curvature[0] = 1.0/self.e_curve[0]

        if self.e_curve[1]==0:
            curvature[1] = Inf
        else:
            curvature[1] = 1.0/self.e_curve[1]

        self.Tin[0,0,0] = -(h/2.0) * tan(self.e_angle[0])**2
        self.Tin[0,2,2] = (h/2.0) / cos(self.e_angle[0])**2
        self.Tin[1,0,0] = (h/(2.0*curvature[0])) / cos(self.e_angle[0]**3)
        self.Tin[1,0,1] = h * tan(self.e_angle[0])**2
        self.Tin[1,0,5] = -h * tan(self.e_angle[0])
        self.Tin[1,2,2] = h**2 * \
            (0.5 + tan(self.e_angle[0])**2) * tan(self.e_angle[0]) - \
            (h/(2*curvature[0]))/cos(self.e_angle[0]**3)
        self.Tin[1,2,3] = -h * tan(self.e_angle[0])**2
        self.Tin[2,0,2] = h * tan(self.e_angle[0])**2
        self.Tin[3,0,2] = -(h/curvature[0])/cos(self.e_angle[0])**3
        self.Tin[3,0,3] = -h * tan(self.e_angle[0])**2
        self.Tin[3,1,2] = -h / cos(self.e_angle[0])**2
        self.Tin[3,2,5] = h * tan(self.e_angle[0]) - \
            h * phi_in/cos(self.e_angle[0]-phi_in)**2
        self.Tin[0,0,0] = (h/2.0) * tan(self.e_angle[1])**2
        self.Tin[0,2,2] = -(h/2.0) / cos(self.e_angle[1])**2
        self.Tin[1,0,0] = (h/(2.0*curvature[1])) / cos(self.e_angle[1]**3) - \
            h**2 * (0.5*tan(self.e_angle[1])**2)*tan(self.e_angle[1])
        self.Tin[1,0,1] = -h * tan(self.e_angle[1])**2
        self.Tin[1,0,5] = -h * tan(self.e_angle[1])
        self.Tin[1,2,2] = h**2 * \
            (-0.5*tan(self.e_angle[1])**2) * tan(self.e_angle[1]) - \
            (h/(2*curvature[1]))/cos(self.e_angle[1]**3)
        self.Tin[1,2,3] = h * tan(self.e_angle[1])**2
        self.Tin[2,0,2] = -h * tan(self.e_angle[1])**2
        self.Tin[3,0,2] = -(h/curvature[1])/cos(self.e_angle[1])**3 + \
            h**2 * tan(self.e_angle[1]) / cos(self.e_angle[1])**2
        self.Tin[3,0,3] = h * tan(self.e_angle[1])**2
        self.Tin[3,1,2] = h / cos(self.e_angle[1])**2
        self.Tin[3,2,5] = h * tan(self.e_angle[1]) - \
            h * phi_out/cos(self.e_angle[1]-phi_out)**2

    def TrackThruEle(self,beam_in):
        [r_in,r_out] = RotMats(-self.tilt)
        beam_out = deepcopy(beam_in)
        if self.is_diag:
            self.DiagOut.centroid = mean(beam_out.x,1)
        beam_out.x = dot(r_in,beam_out.x)
        self.alpha = self.egain / P
        doflist = range(6)
        zeromat = zeros((6,6,6))

        for partnum in range(beam_out.x.shape[1]):
            self.CalcTmat((beam_out.x[5,partnum] * self.P) + self.P)
            if not all(self.Tin == zeromat):
                for dof in doflist:
                    for Tdof1 in doflist:
                        for Tdof2 in doflist:
                            beam_out.x[dof,partnum] += \
                                self.Tin[dof,Tdof1,Tdof2] * \
                                beam_in.x[Tdof1,partnum] * \
                                beam_in.x[Tdof2,partnum]
        
        for partnum in range(beam_out.x.shape[1]):
            self.CalcRmat((beam_out.x[5,partnumrflambda] * self.P) + self.P)
            beam_out.x[:,partnum] = dot(self.R,beam_out.x[:,partnum])

        for partnum in range(beam_out.x.shape[1]):
            self.CalcTmat((beam_out.x[5,partnum] * self.P) + self.P)
            if not all(self.Tout == zeromat):
                for dof in doflist:
                    for Tdof1 in doflist:
                        for Tdof2 in doflist:
                            beam_out.x[dof,partnum] += \
                                self.Tout[dof,Tdof1,Tdof2] * \
                                beam_in.x[Tdof1,partnum] * \
                                beam_in.x[Tdof2,partnum]

        self.CalcRmat()
        self.CalcTmat()
        beam_out.x = dot(r_out,beam_out.x)

        return beam_out

# ===============================================================
# Each cavity type is a subclass of BasicCav
class AccCav(BasicCav):
    def CalcRmat(self, P=-1, L=-1, z=0):
        if P==-1: P = self.P*1e9         
        if L==-1: L = self.L
        if z==0:  phi = self.phi
        else:
            rflambda = c_light/self.freq
            phi = self.phi + (z/rflambda)*2*pi 

        energy0 = self.energy0
        gamma0  = self.gamma0
        beta0   = self.beta0
        alpha = self.alpha
        phi_perp = sqrt(abs(pi * alpha * cos(phi)) / 2.0)
        phi_para = sqrt(abs(pi * alpha * cos(phi))) / (gamma0 * beta0)
        C_perp = cos(phi_perp)
        S_perp = sin(phi_perp)
        C_para = cos(phi_para)
        S_para = sin(phi_para)
    
        if not hasattr(self,'R'):
            self.R = zeros((6,6))

        self.R[0,0] = C_perp
        self.R[0,1] = (L/phi_perp) * S_perp
        self.R[1,0] = -(phi_perp/L) * S_perp
        self.R[1,1] = C_perp

        self.R[2,2] = C_perp
        self.R[2,3] = (L/phi_perp) * S_perp
        self.R[3,2] = -(phi_perp/L) * S_perp
        self.R[3,3] = C_perp

        self.R[4,4] = C_para
        self.R[4,5] = (1/(beta0**2 * gamma0**2)) * (L/phi_para) * S_para
        self.R[5,4] = -(beta0**2 * gamma0**2) * (phi_para/L) * S_para

    def DriftMap(self,beam_in,DistDrift,P):
         beam_out = deepcopy(beam_in)
         print P
     
         x = beam_in.x[0]
         P_x = beam_in.x[1]
         y = beam_in.x[2]
         P_y = beam_in.x[3]
         z = beam_in.x[4]
         delta = beam_in.x[5]
         
         L = self.L*DistDrift
         energy0 = sqrt((P*1e9)**2 + beam_in.restmass**2)
         gamma0 = energy0/beam_in.restmass
         beta0 = sqrt(1 - 1/(gamma0**2))
 
         P_s = sqrt((((1/beta0)+delta)**2) - P_x**2 - P_y**2 - (1/((beta0**2)*(gamma0**2))))
         
         beam_out.x[0] = x + L*(P_x/P_s)
         beam_out.x[2] = y + L*(P_y/P_s)
         beam_out.x[4] = z + L*((1/beta0)-((delta+(1/beta0))/P_s))
         return beam_out

    def KickMap(self,beam_in,DistKick,P):
        from scipy.special import jn
        beam_out = deepcopy(beam_in)
        
        x = beam_in.x[0]
        P_x = beam_in.x[1]
        y = beam_in.x[2]
        P_y = beam_in.x[3]
        z = beam_in.x[4]
        delta = beam_in.x[5]
        
        rho = sqrt(x**2 + y**2)
        phi0 = self.phi
        k = self.k
        A = DistKick*(self.egain/(self.P*1e9))

        for partnum in xrange(beam_out.x.shape[1]):
            if (rho[partnum] != 0):
                beam_out.x[1] = P_x - A*(x/rho)*jn(1,(k*rho))*cos(phi0-(k*z))
                beam_out.x[3] = P_y - A*(y/rho)*jn(1,(k*rho))*cos(phi0-(k*z))
                beam_out.x[5] = delta + A*jn(0,(k*rho))*sin(phi0-(k*z))
            else:
                beam_out.x[5] = delta + A*jn(0,(k*rho))*sin(phi0-(k*z))

        P = P + A*jn(0,0)*sin(phi0)
        
        return beam_out, P

    def TrackThruEle(self,beam_in):
        from globals import lietrackarray
        intermed_beam = deepcopy(beam_in)
        try:
            DistDrift = lietrackarray[self.numdrift].dlengths
            DistKick = lietrackarray[self.numdrift].klengths    
        except IndexError:
            raise IndexError('Yoshida factorisation for numdrift=%i is not defined' % numdrift)
        except AttributeError:
            raise AttributeError('Yoshida factorisation for numdrift=%i is not defined' % numdrift)

        N = 0
        M = 0
        P = self.P
        beam_out = self.DriftMap(intermed_beam,DistDrift[N],P)
     #   beam_out = self.KickMap(beam_out,DistKick[M])
        N += 1
        while N < self.numdrift:
            beam_out,P = self.KickMap(beam_out,DistKick[M],P)
            beam_out = self.DriftMap(beam_out,DistDrift[N],P)
            N += 1
            M += 1
    
        return beam_out

class TCav(BasicCav):
    pass

# ===============================================================
# Each diagnostic type is a subclass of BasicDiag
class BPM(BasicDiag):
    def Processor(self, beam_in):
        self.DiagOut.S_pos = self.S
        if self.res[0]==0:
            self.DiagOut.x_centroid  = mean(beam_in.x[0,:])
            self.DiagOut.y_centroid  = mean(beam_in.x[2,:])
        else:
            self.DiagOut.x_centroid = \
                mean(beam_in.x[0,:]) + normal(scale=self.res[0], size=1)
            self.DiagOut.y_centroid = \
                mean(beam_in.x[2,:]) + normal(scale=self.res[0], size=1)
        if self.res[1]==0:
            self.DiagOut.xp_centroid = mean(beam_in.x[1,:])
            self.DiagOut.yp_centroid = mean(beam_in.x[3,:])
        else:
            self.DiagOut.xp_centroid = \
                mean(beam_in.x[1,:]) + normal(scale=self.res[1], size=1)
            self.DiagOut.yp_centroid = \
                mean(beam_in.x[3,:]) + normal(scale=self.res[1], size=1)

class Screen(BasicDiag):
    def Processor(self, beam_in):
        self.DiagOut.S_pos = self.S
        self.DiagOut.beam_dist = beam_in.x[0:3:2,:]

class EmitScreen(BasicDiag):
    def Processor(self, beam_in):
        self.DiagOut.S_pos = self.S
        self.DiagOut.beam_dist = beam_in.x

class WireScanner(BasicDiag):
    pass

class OTR(Screen):
    pass

class ICT(BasicDiag):
    def Processor(self, beam_in):
        self.DiagOut.S_pos = self.S
        origparts,numparts = 0,0
        for x in beam_in.x[0,:]:
            origparts += 1
            if not isnan(x):
                numparts += 1
        self.DiagOut.charge = beam_in.Q * (numparts / origparts )


# ===============================================================
# Each drift/collimator type is a subclass of BasicDrift
class Drift(BasicDrift):
    pass
    
class Coll(BasicDrift):
    pass

# ===============================================================
# A mover class
class Mover:
    def __init__(self,dof = array([0,2,5])):
        self.dof = dof
        self.act = zeros(len(dof))
        self.des = zeros(len(dof))

    def Set(self, pos):
        if not len(pos)==len(self.dof):
            raise "Mover degrees of freedom do not match request"
        self.des = array(pos).__copy__()

    def Trim(self):
        self.act = self.des.__copy__()

    def GetAct(self):
        act = zeros(6)
        for i in range(len(self.dof)):
            act[self.dof[i]] = self.act[i]
        return act

class Wake:
    def __init__(self):
        self.l = None
        self.k = None

# ===============================================================
# A few useful functions
def DriftRmat(L):
    R = array([
        [1,L,0,0,0,0],
        [0,1,0,0,0,0],
        [0,0,1,L,0,0],
        [0,0,0,1,0,0],
        [0,0,0,0,1,0],
        [0,0,0,0,0,1],
        ])
    return R

def DriftTmat():
    T = zeros((6,6,6))
    return T

def SplitParams(param):
    outparam = array([0.0,0.0])
    if array(param).size==1:
        outparam[0] = param
        outparam[1] = param
    elif array(param).size>1:
        outparam[0] = param[0]
        outparam[1] = param[1]
    return outparam

def RotMats(alpha):
    C = cos(-alpha)
    S = sin(-alpha)
    r_inrot = array([
        [C ,0 ,S ,0 ,0 ,0],
        [0 ,C ,0 ,S ,0 ,0],
        [-S,0 ,C ,0 ,0 ,0],
        [0 ,-S,0 ,C ,0 ,0],
        [0 ,0 ,0 ,0 ,1 ,0],
        [0 ,0 ,0 ,0 ,0 ,1],
        ])
    C = cos(alpha)
    S = sin(alpha)
    r_outrot = array([
        [C ,0 ,S ,0 ,0 ,0],
        [0 ,C ,0 ,S ,0 ,0],
        [-S,0 ,C ,0 ,0 ,0],
        [0 ,-S,0 ,C ,0 ,0],
        [0 ,0 ,0 ,0 ,1 ,0],
        [0 ,0 ,0 ,0 ,0 ,1],
        ])
    return r_inrot,r_outrot

# ===============================================================
# Test suite
if __name__ == '__main__':
    from latticeloader import *

    templine = Line(LoadFlatLatFile('tempfile.txt')[1155:])
    mybeam = Beam(x=array([0,0,0,0,0,1.3]))
    beamout = templine.Track(mybeam)
    print beamout.x

