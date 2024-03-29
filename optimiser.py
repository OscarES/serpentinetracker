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

class FncClass(object) :
    def __init__(self, narg, varnames) :
        self.co_argcount = narg
        self.co_varnames = tuple(varnames)

class Optimiser(object) :

    def __init__(self,s, debug=False) :
        ''' Class to optimise basic lattive functions based on variable which  can be set via object function methods
        s : Serpentine or beamline (not implemented but quick) object'''
        print 'Optimiser.__init__()'
        # debug variable for verbose output
        self.debug = debug

        # serpentine object (should this be copied? mmmmm)
        self.s = s
        self.variables   = {}
        self.constraints = {}

    def AddVariable(self,name,element,attrName,step=1.0,min=0,max=0) :
        ''' Add variable (input/output) with range 
        elementName : name within beamline
        className   : class within element where data is stored, if list can be nested classes
        attrName    : attribute of class 
        type        : eq, lt, gt
        value       : Contraint value (equality or inequality)
        start       : Parameter start value 
        step        : Set initial step 
        max         : Maximum value 
        min         : Minimum value'''
        self.variables[name] = {'data':{'ele':element,'attr':attrName,'step':step,'min':min,'max':max}}
    
    def AddConstraint(self,name,element,attrName,type,value, weight) :
        ''' Add constraint (output)'''
#        v = {'name':name, data:{'ele':element,'attr':attrName,'type':type,'value':value}}
        self.constraints[name] = {'data':{'ele':element,'attr':attrName,'type':type,'value':value,'weight':weight}}

    def TestVariables(self) :
        ''' Test we can extract and set variables in the variables list'''
        e = None 
        for v in self.variables :
            print 'Optimiser.TestVariables  >',v,self.variables[v],self.GetValue(v)

    def TestConstraints(self) :
        ''' Test if we can extract the output '''
        for v in self.constraints :
            print 'Optimiser.TestConstraints>',v,self.constraints[v],self.GetValue(v)

    def GetValue(self,v) :                   
        if type(v) == str :
            try :
                v = self.variables[v]['data']
            except :
                v = self.constraints[v]['data']
                
        e = self.s.beamline.GetEleByName(v['ele'])[0]
        for a in v['attr'] :
            e = e.__getattribute__(a)        
        return e

    def SetValue(self, v, value) :
        if type(v) == str :
            try :
                v = self.variables[v]['data']
            except KeyError :
                print 'Optimiser : no such variables' 
                return 

        e = self.s.beamline.GetEleByName(v['ele'])[0]
        for a in v['attr'][:-2] :
            e = e.__getattribute__(a)        
        e.__setattr__(v['attr'][-1],value)

    def EvalConstraint(self, v) :
        v = self.constraints[v]['data'] 

        val = self.GetValue(v)
        if v['type'] == 'eq' : 
            ret = v['weight']*(val - v['value'])**2
        elif v['type'] == 'lt' :
            if val > v['value'] :
                ret = v['weight']*(val - v['value'])**2
            else :
                ret = 0.0
        elif v['type'] == 'gt' :
            if val < v['value'] :
                ret = v['weight']*(val - v['value'])**2
            else :
                ret = 0.0             
        
        return ret

    def PrintDefinition() :
        pass

    def Run(self) :
        ''' Execute minuit with Objective function '''
        print 'Optimizer.Run>'
        
        # check for a system to optimise
        if len(self.variables) == 0 or len(self.constraints) == 0 :
            print 'Optimiser.Run> Insufficent variables/constraints'

        # make object callable to minuit
        self.func_code = FncClass(len(self.variables),self.variables.keys())
               
        # create minuit object
        try :
            self.m = minuit.Minuit(self)
        except MinuitError : 
            pass

        # set starting values and limits 
        for v in self.variables : 
            self.m.values[v] = self.GetValue(v)
            self.m.limits[v] = (self.variables[v]['data']['min'],self.variables[v]['data']['max'])
        self.m.printMode = 0
        self.m.migrad()
        print self.m.values        
        for v in self.constraints :
            print v,self.GetValue(v)
        
    def Clear(self) :
        ''' Clear all variables and constraints'''
        print 'Optimizer.Clear>'
        self.variables  = []
        self.constaints = []

    def __call__(self,*args) :
        ''' Objective function''' 
        if self.debug :
            print 'Optimiser.Objective>',args

        # set variables
        i = 0
        for v in self.variables :
            self.SetValue(v,args[i])
            i += 1

        # Twiss progagation
        for ele in self.s.beamline:
            ele.CalcRmat()
        self.s.TwissProp();
            
        # Track in serpentine object
        # self.s.Track(); print '\n'

        # loop constraints 
        sum = 0 
        for v in self.constraints :
            if self.debug : 
                print 'Optimiser.Objective> val=', self.GetValue(v)
            contrib = self.EvalConstraint(v)
            sum += contrib
        if self.debug :
            print 'Optimiser.Objective> sum=',sum,'\n'
        return sum

    def __repr__(self):
        return ''


def OptimiserTest() :
    ''' Simple two magnet focus optimisation :
    To run 
    python>import visualize 
    python>[s,o] = optimiser.OptimiserTest()
    python>o.Run()'''


    import elements
    import beamline 
    import serpentine 
    import visualize

    # build simple 2 magnet system
    qf  = elements.Quad("QF",L=0.25,P=1.25,B=5)
    dr1 = elements.Drift("D1",L=0.30,P=1.25)
    qd  = elements.Quad("QD",L=0.25,P=1.25,B=-5)
    dr2 = elements.Drift("D2",L=1.0,P=1.25) 
    m1  = elements.BasicDiag("M1",L=0.0)

    bl = beamline.Line([qf,dr1,qd,dr2,m1])

    # make simple beam in
    # set twiss parameters
    mytwiss = elements.Twiss(betax  = 6.85338806855804,     betay  = 2.94129410712918,
                             alphax = 1.11230788371885,     alphay = -1.91105724003646,
                             etax   = 3.89188697330735e-012,etay   = 0,
                             etaxp  = 63.1945125619190e-015,etayp  = 0,
                             phix   = 0,                    phiy   = 0,
                             nemitx = 5.08807339588144e-006,nemity = 50.8807339588144e-009,
                             sigz   = 8.00000000000000e-003,sigP   = 1.03999991965541e-003,
                             pz_cor = 0)

    print bl
    s = serpentine.Serpentine(bl,twiss=mytwiss)
    # Visualisation check 
    vis = visualize.Visualize()
    figure(1)
    vis.PlotTwiss(s)

    # optimizer test
    o = Optimiser(s)
    o.AddVariable('qf1b','QF',['B'],0.1,-50,50)
    o.AddVariable('qd1b','QD',['B'],0.1,-50,50)
    o.AddConstraint('fbetax','M1',['twiss','betax'],'lt',1.5,1.0)
    o.AddConstraint('fbetay','M1',['twiss','betay'],'eq',0.3,20.0)
    o.TestVariables()
    o.TestConstraints()
    o.Run()

    # Visualisation check
    vis.PlotTwiss(o.s)
    
    return [s,o]



