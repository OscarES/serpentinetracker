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
import elements as BL
import beamline
import beamrep
import latticeloader as LL
from numpy import array, zeros, mean, std
from pylab import plot, subplot, figure, xlabel, ylabel
from globals import electron_mass

class Serpentine:
    def __init__(self,line, twiss=None, beam=None, P=None):
        self.offset = zeros(6)

        if not (isinstance(line, basestring) or isinstance(line, beamline.Line)):
            raise "Must provide a beamline or a path to an AML file."

        if isinstance(line, basestring):
            self.beamline,parttype = LL.LoadLatFile(line,P=P)
        elif isinstance(line, beamline.Line):
            self.beamline = line
        self.beamline.SetSPos()

        if (isinstance(twiss, BL.Twiss) and isinstance(beam, beamrep.ElectronBeam)):
            raise "Only give twiss *or* beam, not both!"
        elif isinstance(twiss, BL.Twiss):
            self.twiss = twiss
            self.twiss_out = self.beamline.TwissProp(twiss)
            self.SingleParticle(P=self.beamline[0].P)
            self.beam_out = None
        elif isinstance(beam, beamrep.ElectronBeam):
            self.beam_in = beam
            self.beam_out = None
            self.twiss = BL.Twiss()
            self.twiss_out = BL.Twiss()
        if parttype.upper()=='PROTON':
            self.beam_in.MakeProtons()
        elif parttype.upper()=='POSITRONS':
            self.beam_in.MakePositrons()
        elif parttype.upper()=='ELECTRONS':
            pass
        else:
            raise ValueError('Unrecognised particle species: %s' % parttype)

    def SingleParticle(self, P=1, Q=1e-9, chargesign=-1,restmass=electron_mass):
        self.beam_in = beamrep.Beam(P, Q, chargesign, restmass)

    def MakeGaussBeam(self,N=10000,Q=1e-9,pos=array([0,0,0,0,0,1]),
        sig=array([1e-3,1e-3,1e-3,1e-3,1e-3,0.01]),chargesign=-1,
        restmass=electron_mass):
        self.beam_in = beamrep.GaussBeam(N,Q,pos,sig,chargesign,restmass)

    def TwissGaussBeam(self, N=10000, Q=1e-9, chargesign=-1,restmass=electron_mass):
        pos = array([0,0,0,0,0,self.beamline[0].P])
        self.beam_in = beamrep.TwissGaussBeam(self.twiss,N,pos,Q,chargesign,restmass)

    def TwissProp(self):
        self.twiss_out = self.beamline.TwissProp(self.twiss)

    def PlotTwiss(self,bx=1,by=1,ax=0,ay=0,px=0,py=0):
        self.beamline.PlotTwiss(bx=bx,by=by,ax=ax,ay=ay,px=px,py=py)

    def SetMomProfile(self,ini_p=None):
        from serpentine import Serpentine
        from scipy import sin,sqrt
        if ini_p==None: ini_p = self.beamline[0].P*1e9
        cum_p = ini_p
        for i in self.beamline:
            if isinstance(i,Serpentine):
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

    def PlotMomProfile(self,formatstr='-x'):
        self.beamline.PlotMomProfile(formatstr)

    def PlotEkProfile(self,formatstr='-x'):
        self.beamline.PlotEkProfile(restmass=self.beam_in.restmass,formatstr=formatstr)
        
    def Track(self):
        self.beamline.offset = self.offset
        if hasattr(self,'Mover'):
            self.beamline.offset += self.Mover.GetAct()
        self.beam_out = self.beamline.Track(self.beam_in)
        if hasattr(self,'Mover'):
            self.beamline.offset -= self.Mover.GetAct()
        del self.beamline.offset

    def GroupElesByName(self):
        group_list = list()
        for ele in self.beamline:
            if isinstance(ele,Serpentine):
                continue
            base_name = ele.name
            for newele in self.beamline[self.beamline.index(ele)+1:]:
                if (base_name.find(newele.name)>-1) or ((newele.name).find(base_name)>-1):
                    group_list.append(array([self.beamline.index(ele),1 + self.beamline.index(newele)]))
                    break
        newline = beamline.Line(self.beamline[:group_list[0][0]])
        for i in xrange(len(group_list[:-1])):
            mid_line = beamline.Line(self.beamline[group_list[i][0]:group_list[i][1]])
            end_line = beamline.Line(self.beamline[group_list[i][1]:group_list[i+1][0]])
            newline.append(Serpentine(line=mid_line))
            newline = beamline.Line(newline + end_line)
        self.beamline = newline
        self.beamline.SetSPos()

    def AddMoversToGroups(self,dof=array([0,2,5])):
        for ele in self.beamline:
            if isinstance(ele,Serpentine):
                self.AddMovers(self.beamline.index(ele),dof)

    def AddMovers(self,elelist,dof=array([0,2,5])):
        if isinstance(elelist,int) or isinstance(elelist,float):
            self.beamline[elelist].Mover = BL.Mover(dof)
        else:
            for ind in elelist:
                self.beamline[ind].Mover = BL.Mover(dof)

    def GetBPMReadings(self):
        bpms = self.beamline.GetEleByType('BPM')
        S,x,y,xp,yp = [],[],[],[],[]
        readings = zeros((3,len(bpms)))
        for ele in bpms:
            if isinstance(ele,dict):
                line_readings = self[ele.keys()[-1]].GetBPMReadings
                print line_readings
                S  = concatenate((S,line_readings[0,:]))
                x  = concatenate((S,line_readings[1,:]))
                y  = concatenate((S,line_readings[2,:]))
                xp = concatenate((S,line_readings[3,:]))
                yp = concatenate((S,line_readings[4,:]))
                continue
            S.append(ele.DiagOut.S_pos)
            x.append(ele.DiagOut.x_centroid)
            y.append(ele.DiagOut.y_centroid)
            xp.append(ele.DiagOut.xp_centroid)
            yp.append(ele.DiagOut.yp_centroid)

        return array([S,x,y,xp,yp])
        
    def PlotBPMReadings(self,str):
        readings = self.GetBPMReadings()
        subplot(211); plot(readings[0,:],readings[1,:],'x'+str)
        xlabel('S / m')
        ylabel('x / m')
        subplot(212); plot(readings[0,:],readings[2,:],'x'+str)
        xlabel('S / m')
        ylabel('y / m')

    def GetScreenReadings(self):
        Screens = self.beamline.GetEleByType('Screen')
        beam_dist,xstd,ystd = [],[],[]
        for ele in Screens:
            beam_dist.append(ele.DiagOut.beam_dist)
            xstd.append(std(ele.DiagOut.beam_dist[0,:]))
            ystd.append(std(ele.DiagOut.beam_dist[1,:]))
        return beam_dist, xstd, ystd

    def PlotScreenReadings(self):
        Screens = self.beamline.GetEleByType('Screen')
        (beam_dist, xstd, ystd) = self.GetScreenReadings()
        for ele in Screens:
            i = Screens.index(ele)
            figure(i)
            plot(beam_dist[i][0,:]*1e6,beam_dist[i][1,:],'x')

    def AddEmitScreensToBPMs(self):
        for ele in self.beamline:
            if isinstance(ele, Serpentine):
                ele.AddEmitScreensToBPMs()
            elif isinstance(ele,BL.BPM):
                newline = beamline.Line((BL.EmitScreen(name="S"+ele.name,L=0),ele))
                ind = self.beamline.index(ele)
                self.beamline[ind] = Serpentine(line=newline)
        self.beamline.SetSPos()

    def AddEleToEle(self,ele1,ele2):
        for ele in self.beamline:
            if isinstance(ele, Serpentine):
                ele.AddEleToEle(ele1,ele2)
            elif isinstance(ele,ele2):
                newline = beamline.Line((ele1(name="S"+ele.name,L=0),ele))
                ind = self.beamline.index(ele)
                self.beamline[ind] = Serpentine(line=newline)
        self.beamline.SetSPos()

    def AddWakesToAccCavs(self,file):
        l,k = load(file, usecols=(1,2), skiprows=3, unpack=True)
        acccavs = self.beamline.GetEleByType('AccCav')
        for ele in acccavs:
            ele.TSRW = Wake()
            ele.TSRW.l = l
            ele.TSRW.k = k
