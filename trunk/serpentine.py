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

"""The main class for the Serpentine code.
All analyses will use Serpentine, or a class derived from it, as their base."""

import elements as EL
import beamline
import beamrep
import latticeloader as LL
from numpy import array, zeros, std, concatenate, loadtxt
from pylab import plot, subplot, figure, xlabel, ylabel
from globals import electron_mass
from scipy import sin, sqrt

class Serpentine:
    """The primary class for this tracking code"""
    def __init__(self, line, twiss=None, beam=None, ini_mom=None, parttype='ELECTRON'):
        self.offset = zeros(6)

        if not (isinstance(line, basestring) or 
            isinstance(line, beamline.Line)):
            raise ValueError(
                "Must provide a beamline or a path to an AML file.")

        if isinstance(line, basestring):
            self.beamline, retparttype = LL.loadlatfile(line, ini_mom=ini_mom)
	    if not retparttype==parttype:
	        parttype = retparttype
        elif isinstance(line, beamline.Line):
            self.beamline = line
        self.beamline.SetSPos()

        if (isinstance(twiss, EL.Twiss) and isinstance(beam, beamrep.Beam)):
            raise ValueError("Only give twiss *or* beam, not both!")
        elif isinstance(twiss, EL.Twiss):
            self.twiss = twiss
            self.twiss_out = self.beamline.TwissProp(twiss)
            self.SingleParticle(P=self.beamline[0].P)
            self.beam_out = None
        elif isinstance(beam, beamrep.Beam):
            self.beam_in = beam
            self.beam_out = None
            self.twiss = EL.Twiss()
            self.twiss_out = EL.Twiss()
        if parttype.upper()=='PROTON':
            try:
                self.beam_in.MakeProtons()
            except AttributeError:
                pass
        elif parttype.upper()=='POSITRON':
            try:
                self.beam_in.MakePositrons()
            except AttributeError:
                pass
        elif parttype.upper()=='ELECTRON':
            pass
        else:
            raise ValueError('Unrecognised particle species: %s' % parttype)

        self.mover = EL.Mover()

    def SingleParticle(self, P=1, Q=1e-9, chargesign=-1, 
        restmass=electron_mass):
        """Create a single particle beam as an attribute of the
        Serpentine class."""
        self.beam_in = beamrep.Beam(P, Q, chargesign, restmass)

    def MakeGaussBeam(self, N=10000, Q=1e-9, pos=array([0, 0, 0, 0, 0, 1]),
        sig=array([1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 0.01]), chargesign=-1,
        restmass=electron_mass):
        """Create a multi-particle beam with a Gaussian spread in each
        degree of freedom."""
        self.beam_in = beamrep.GaussBeam(N, Q, pos, sig, chargesign, restmass)

    def TwissGaussBeam(self, N=10000, Q=1e-9, chargesign=-1, 
        restmass=electron_mass):
        """Create a multi-particle beam with a Gaussian spread in each
        degree of freedom.  The statistics """
        pos = array([0, 0, 0, 0, 0, self.beamline[0].P])
        self.beam_in = beamrep.TwissGaussBeam(self.twiss, 
            N, pos, Q, chargesign, restmass)

    def TwissProp(self):
        """Propagate the input Twiss parameters down the length of the
        lattice.  The output is stored as an attribute of self, and the Twiss
        parameters on the exit of each element are stored as attributes
        of that element."""
        self.twiss_out = self.beamline.TwissProp(self.twiss)

    def PlotTwiss(self, betax=1, betay=1):
        """Plot the Twiss parameters.  This actually calls the equivalent
        method in self.beamline."""
        self.beamline.PlotTwiss(betax=betax, betay=betay)

    def SetMomProfile(self, ini_p=None):
        """Given an initial momentum value (ini_p or the P attribute of the
        first element in the lattice), this computes and stores the expected
        design momentum at the entrance of each element.
        Currently the only element that changes the momentum are AccCav
        instances."""
        if ini_p == None:
            ini_p = self.beamline[0].P*1e9
        cum_p = ini_p
        for i in self.beamline:
            if isinstance(i, Serpentine):
                i.SetMomProfile(ini_p=cum_p)
                continue
            i.P = cum_p
            if i.__class__.__name__ == 'AccCav':
                energy = sqrt(cum_p**2 + self.beam_in.restmass**2)
                energy += i.egain * sin(i.phi)
                cum_p = sqrt(energy**2 - self.beam_in.restmass**2)

        for ele in self.beamline:
            ele.CalcRmat()

        self.TwissProp()

    def PlotMomProfile(self, formatstr='-x'):
        """Plots the momentum profile over the lattice."""
        self.beamline.PlotMomProfile(formatstr)

    def PlotEkProfile(self, formatstr='-x'):
        """Plots the kinetic energy profile over the lattice."""
        self.beamline.PlotEkProfile(restmass=self.beam_in.restmass, 
            formatstr=formatstr)

    def GetRFPhases(self):
        """Returns the phases of any RF elements in the lattice"""
        return self.beamline.GetRFPhases()
    
    def PlotRFPhases(self):
        """Plots the phases of any RF elements in the lattice"""
        self.beamline.PlotRFPhases()

    def Track(self):
        """The tracking method. Adjusts the beam to take account of
        any global mover setting, calls the Track() method of self.beamline,
        and then resets the changes made due to the mover settings."""
        offset = self.offset
        if hasattr(self, 'mover'):
            offset += self.mover.GetAct()
        self.beam_out = self.beamline.Track(self.beam_in, offset)
        if hasattr(self, 'mover'):
            offset -= self.mover.GetAct()

    def GroupElesByName(self):
        """Finds any groups of consecuctive elements that share the same
        name, and bundles these together into an internal Serpentine object."""
        group_list = list()
        for ele in self.beamline:
            if isinstance(ele, Serpentine):
                continue
            base_name = ele.name
            for newele in self.beamline[self.beamline.index(ele)+1:]:
                if ((base_name.find(newele.name)>-1) or 
                        ((newele.name).find(base_name)>-1)):
                    group_list.append(array([self.beamline.index(ele), 1 + 
                        self.beamline.index(newele)]))
                    break
        newline = beamline.Line(self.beamline[:group_list[0][0]])
        for i in xrange(len(group_list[:-1])):
            mid_line = beamline.Line(
                self.beamline[group_list[i][0]:group_list[i][1]])
            end_line = beamline.Line(
                self.beamline[group_list[i][1]:group_list[i+1][0]])
            newline.append(Serpentine(line=mid_line))
            newline = beamline.Line(newline + end_line)
        self.beamline = newline
        self.beamline.SetSPos()

    def AddMoversToGroups(self, dof=array([0, 2, 5])):
        """Adds a mover to any grouped elements (i.e. any elements
        that are Serpentine objects)."""
        for ele in self.beamline:
            if isinstance(ele, Serpentine):
                self.AddMovers(self.beamline.index(ele), dof)

    def AddMovers(self, elelist, dof=array([0, 2, 5])):
        """Add movers to a list of elements."""
        if isinstance(elelist, int) or isinstance(elelist, float):
            self.beamline[elelist].mover = EL.Mover(dof)
        else:
            for ind in elelist:
                self.beamline[ind].mover = EL.Mover(dof)

    def GetBPMReadings(self, classname='BPM'):
        """Return an array of BPM readings from the most recent tracking
        operation"""
        bpms = self.beamline.GetEleByType(classname)
        S, x, y, xp, yp = [], [], [], [], []
        for ele in bpms:
            if isinstance(ele, dict):
                line_readings = self[ele.keys()[-1]].GetBPMReadings
                print line_readings
                S  = concatenate((S, line_readings[0, :]))
                x  = concatenate((S, line_readings[1, :]))
                y  = concatenate((S, line_readings[2, :]))
                xp = concatenate((S, line_readings[3, :]))
                yp = concatenate((S, line_readings[4, :]))
                continue
            S.append(ele.DiagOut.S_pos)
            x.append(ele.DiagOut.x_centroid)
            y.append(ele.DiagOut.y_centroid)
            xp.append(ele.DiagOut.xp_centroid)
            yp.append(ele.DiagOut.yp_centroid)

        return array([S, x, y, xp, yp])
        
    def PlotBPMReadings(self, formatstr='', classname='BPM'):
        """Plot the BPM readings from the most recent tracking operation"""
        readings = self.GetBPMReadings(classname)
        subplot(211)
        plot(readings[0, :], readings[1, :], 'x'+formatstr)
        xlabel('S / m')
        ylabel('x / m')
        subplot(212)
        plot(readings[0, :], readings[2, :], 'x'+formatstr)
        xlabel('S / m')
        ylabel('y / m')

    def GetScreenReadings(self):
        """Return Screen() data from the most recent tracking operation"""
        screens = self.beamline.GetEleByType('Screen')
        beam_dist, xstd, ystd = [], [], []
        for ele in screens:
            beam_dist.append(ele.DiagOut.beam_dist)
            xstd.append(std(ele.DiagOut.beam_dist[0, :]))
            ystd.append(std(ele.DiagOut.beam_dist[1, :]))
        return beam_dist, xstd, ystd

    def PlotScreenReadings(self):
        """Plot Screen() data from the most recent tracking operation"""
        screens = self.beamline.GetEleByType('Screen')
        beam_dist = self.GetScreenReadings()[0]
        for ele in screens:
            i = screens.index(ele)
            figure(i)
            plot(beam_dist[i][0, :]*1e6, beam_dist[i][1, :], 'x')

    def AddEmitScreensToBPMs(self):
        """Adds Screen objects to any BPM objects by creating a single 
        Serpentine object containing both objects"""
        for ele in self.beamline:
            if isinstance(ele, Serpentine):
                ele.AddEmitScreensToBPMs()
            elif isinstance(ele, EL.BPM):
                newline = beamline.Line(
                    (EL.EmitScreen(name="S"+ele.name, L=0), ele))
                ind = self.beamline.index(ele)
                self.beamline[ind] = Serpentine(line=newline)
        self.beamline.SetSPos()

    def AddEleToEle(self, ele1, ele2):
        """Replace any instance of ele2 in self.beamline with a 
        Serpentine object consisting of ele1 and ele2"""
        for ele in self.beamline:
            if isinstance(ele, Serpentine):
                ele.AddEleToEle(ele1, ele2)
            elif isinstance(ele, ele2):
                newline = beamline.Line((ele1(name="S"+ele.name, L=0), ele))
                ind = self.beamline.index(ele)
                self.beamline[ind] = Serpentine(line=newline)
        self.beamline.SetSPos()

    def AddWakesToAccCavs(self, fname):
        """Extracts data from a CSV file to create a Wake() object to add
        to each AccCav element"""
        l, k = loadtxt(fname, usecols=(1, 2), skiprows=3, unpack=True)
        acccavs = self.beamline.GetEleByType('AccCav')
        for ele in acccavs:
            ele.TSRW = EL.Wake()
            ele.TSRW.l = l
            ele.TSRW.k = k
