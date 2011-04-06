import pylab as pl
import scipy as sp

from serpentine import *
from elements import *

class AtfExt :
    def __init__(self) :
        print 'AtfExt:__init__'

        # set twiss parameters
        mytwiss = Twiss()
        mytwiss.betax  = 6.85338806855804
        mytwiss.alphax = 1.11230788371885
        mytwiss.etax   = 3.89188697330735e-012
        mytwiss.etaxp  = 63.1945125619190e-015
        mytwiss.betay  = 2.94129410712918
        mytwiss.alphay = -1.91105724003646
        mytwiss.etay   = 0
        mytwiss.etayp  = 0
        mytwiss.nemitx = 5.08807339588144e-006
        mytwiss.nemity = 50.8807339588144e-009
        mytwiss.sigz   = 8.00000000000000e-003
        mytwiss.sigP   = 1.03999991965541e-003
        mytwiss.pz_cor = 0

        # load beam line
        self.atfFull = Serpentine(line='newATF2lat.aml',twiss=mytwiss)
        self.atfExt  = Serpentine(line=beamline.Line(self.atfFull.beamline[947:]),twiss=mytwiss)
        
        # zero zero cors
        self.atfExt.beamline.ZeroCors()

        # Track 
        self.atfExt.Track()
        readings = self.atfExt.GetBPMReadings()
        
        figure(1)
        self.atfExt.PlotBPMReadings()
        figure(2)
        self.atfExt.PlotTwiss()

    def setMagnet(self,name, value) :
        ei = self.atfExt.beamline.FindEleByName(name)
        print ei
        e  = self.atfExt.beamline[ei[0]]
        e.B = value
        
    def plotOrbit(self) :
        self.atfExt.PlotBPMReadings()

    def plotTwiss(self) :
        self.atfExt.PlotTwiss()
        
    def run(self) :
        self.atfExt.Track()

    def jitterBeam(self) : 
        r = 1+sp.random.standard_normal()
#        self.s.beam_in.x[5,:] = (1+r/3e4)*self.nominalE
#        print r,self.s.BeamIn.x[5,:]

    
