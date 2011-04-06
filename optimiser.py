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
variable go    : global optimiser needed to pass complicated structure into objective function
Def Objective  : function created by Optimiser and passed to Minuit
Example        :
py> b = twoQuadExample()
py> o = Optimiser()
py> o.AddVariable("lt",10,b.s.beamline[0].GetB,b.s.beamline[0].SetB)
py> o.AddVariable("lt",10,b.s.beamline[2].GetB,b.s.beamline[2].SetB)
py> o.AddVariable("lt",1e-5,b.)
py> o.BuildObjective()
py> o.Run()
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
        ''' Class to optimise basic lattive functions based on variable which  can be set via object function methods'''
        print 'Optimiser.__init__()'
        # serpentine object
        self.s = s
        self.variables   = []
        self.constraints = []

    def AddVariable(self,name,type,value,getter,setter=None,start=0,step=1.0,max=0,min=0) :
        ''' Add variable (input/output) with constaint to Optimiser '''
        v = [name,type,value,getter,setter,start,step,max,min]
        print 'Optimizer.AddVariable',v
        if setter != None :
            self.variables.append(v)
        else :
            self.constraints.append(v)

    def Run(self) :
        ''' Execute minuit with Objective function '''
        print 'Optimizer.Run'
        
    def BuildObjectiveFunction(self) :
        ''' Build function from internal data '''

        # build input line
        s = "def Objective("
        i = 0
        for v in self.variables :
            s += v[0] 
            if i < len(self.variables)-1 :
                s += ","
            else :
                s += ") :\n"        
            i += 1

        # set model parameters
        i = 0
        for v in self.variables :            
            s += '\tgo.variables['+str(i)+'][3](' + v[0] +')\n'
            i += 1

        # execute track
        s += "\tgo.s.Track()\n"

        # get constraints
        i = 0
        for v in self.constraints :
            s += '\tc'+str(1)+' = go.constraints['+str(i)+'][3]()'
            i += 1

        # build objective                    

        print s 

        return 

    def Clear() :
        ''' Clear all variables'''
        pass

    def SetGlobal(self) :
        global _go 
        _go = self

class m1 :
    def __init__(self) :
        self.v1 = 2.5
        self.v2 = 1.0

    def get1(self) :
        return self.v1
    def set1(self) :
        self.v1 = v1

    def get2(self) :
        return self.v2
    def set2(self) :
        self.v2 = v2

    def getprod(self) :
        return v1*v2
    
def OptimiserTest() :
    m = m1()    
    o = Optimiser(s=0)
    o.AddVariable('v1','lt',20,m.get1,m.set1,m.get1())
    o.AddVariable('v2','lt',50,m.get2,m.set2,m.get2())
    o.AddVariable('pr','eq',5,m.getprod)    
    o.BuildObjectiveFunction()
    



