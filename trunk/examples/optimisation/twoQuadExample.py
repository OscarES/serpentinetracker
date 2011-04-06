import pylab as pl
import scipy as sp

import serpentine as st
import latticeloader as ll
from elements import *
import beamline

class twoQuadExample :
    def __init__(self) :

        # set twiss parameters
        self.t = {}
        self.t['betax']  = 6.85338806855804
        self.t['alphax'] = 1.11230788371885
        self.t['etax']   = 3.89188697330735e-012
        self.t['etaxp']  = 63.1945125619190e-015
        self.t['betay']  = 2.94129410712918
        self.t['alphay'] = -1.91105724003646
        self.t['etay']   = 0
        self.t['etayp']  = 0
        self.t['Nemitx'] = 5.08807339588144e-006
        self.t['Nemity'] = 50.8807339588144e-009
        self.t['sigz']   = 8.00000000000000e-003
        self.t['sigP']   = 1.03999991965541e-003
        self.t['PZcor']  = 0

        # create beam line
        self.bl = beamline.Line()
        self.bl.append(Drift(name='ele1', L=0.75))
        self.bl.append(Quad(name='ele2', L=0.25, B=5))
        self.bl.append(Drift(name='ele3', L=1))
        self.bl.append(Quad(name='ele4', L=0.25, B=-5))
        self.bl.append(Drift(name='ele5',L=1))

        # create main control object
        self.s = st.Serpentine(line=self.bl,twiss=self.t);        

        # determine s locations of elements
        self.s.beamline.SetSPos()
        
        # zero zero cors
        self.s.beamline.ZeroCors()
        
        self.s.PlotTwiss()

        q1g = self.s.beamline.GetEleByName('ele2')[0].GetB
        q1s = self.s.beamline.GetEleByName('ele2')[0].SetB

        q2g = self.s.beamline.GetEleByName('ele4')[0].GetB
        q2s = self.s.beamline.GetEleByName('ele4')[0].SetB

        print q1g,q1s,q2g,q2s

