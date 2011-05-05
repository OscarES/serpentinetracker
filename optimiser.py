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
'''Classes and functions for accelerator lattice optimisation
Cass           : Optimiser, controls variables and builds objective function
Example        :
'''

# used to build objective/merit function
import string 

# serpentine
import serpentine

# numpy arrays will be useful here
from numpy import *
from matplotlib.pylab import *

# optimization library
import minuit

# global optimiser 
_go = None

class Optimiser :

    def __init__(self,s) :
        ''' Class to optimise basic lattive functions based on variable which  can be set via object function methods
        s : Serpentine or beamline (not implemented but quick) object'''
        print 'Optimiser.__init__()'
        # serpentine object (should this be copied? mmmmm)
        self.s = s
        self.variables   = []

    def AddVariable(self,name,element,className,attrName,type,value,start=0,step=1.0,max=0,min=0) :
        ''' Add variable (input/output) with constaint to Optimiser 
        elementName : name within beamline
        className   : class within element where data is stored, if list can be nested classes
        attrName    : attribute of class 
        type        : eq, lt, gt
        value       : Contraint value (equality or inequality)
        start       : Parameter start value 
        step        : Set initial step 
        max         : Maximum value 
        min         : Minimum value'''
        v = [name,element,className,attrName,type,value,start,step,max,min]
        print 'Optimizer.AddVariable',v
        self.variables.append(v)
    
        
    def TestVariables(self) :
        for v in self.variables :
            pass

    def PrintDefinition() :
        pass

    def Run(self) :
        ''' Execute minuit with Objective function '''
        print 'Optimizer.Run'
        
    def Clear() :
        ''' Clear all variables'''
        self.variables  = []

    def __repr__(self):
        pass
#        ret = '\n'.join(str(el)+" :: "+str(ele.__class__) for ele in self)        
#        return ret


def OptimiserTest() :
    import elements
    import beamline 
    import serpentine 
    import visualize

    # build simple 2 magnet system
    qf  = elements.Quad("QF",L=0.25,P=1.25,B=5)
    dr1 = elements.Drift("D1",L=0.30,P=1.25)
    qd  = elements.Quad("QD",L=0.25,P=1.25,B=-2.5)
    dr2 = elements.Drift("D2",L=1.0,P=1.25) 
    m1  = elements.BasicDiag("M1",L=0.0)

    bl = beamline.Line([qf,dr1,qd,dr2,m1])

    # make simple beam in
    # set twiss parameters
    mytwiss = elements.Twiss(betax  = 6.85338806855804,     betay  = 2.94129410712918,
                             alphax = 1.11230788371885,     alphay = -1.91105724003646,
                             etax   = 3.89188697330735e-012,etay   = 0,
                             etapx  = 63.1945125619190e-015,etapy  = 0,
                             phix   = 0,                    phiy   = 0,
                             nemitx = 5.08807339588144e-006,nemity = 50.8807339588144e-009,
                             sigz   = 8.00000000000000e-003,sigP   = 1.03999991965541e-003,
                             pz_cor = 0)

    print bl
    s = serpentine.Serpentine(bl,twiss=mytwiss)
    s.Track(); print '';

    # Visualisation check 
    visualize.matplotlib2D(s,label=True)

    # optimizer test
    o = Optimiser(s=0)
    o.AddVariable('qf1b','QF1',[''],'B','lt',50)
    o.AddVariable('qd1b','QD1',[''],'B','lt',50)
    o.AddVariable('fbetax','M1',['twiss'],'betax','eq',5)
    o.TestVariables()

    return bl



