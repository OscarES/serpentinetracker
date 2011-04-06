#!/usr/bin/python

from serpentine import *
import numpy as np
from scipy import optimize
from pylab import figure, plot, title, show, ylim, xlabel, ylabel

# A function to scan a corrector, and record BPM readings
class CorSweep:
    def __init__(self,serpobj,vals,corind,bpmind):
        self.serpobj = serpobj
        self.sweepvals = vals
        self.corind = corind
        self.quadbpm = self.serpobj.beamline[bpmind[0]]
        self.witbpm = self.serpobj.beamline[bpmind[1]]

        # Get the initial corr value
        self.orig_curr = self.serpobj.beamline[self.corind].GetB()

        # Get a 'setter' method for the corr current
        self.setter = self.serpobj.beamline[self.corind].SetB
        # And a 'track' method
        self.track = self.serpobj.Track

    def do_sweep(self):
        # Loop around the values, setting the corr,
        # tracking, and recording the BPM value
        witbpmpos,quadbpmpos = [],[]
        for cor_curr in self.sweepvals:
            self.setter(cor_curr)
            self.track()
            quadbpmpos.append(self.quadbpm.DiagOut.x_centroid)
            witbpmpos.append(self.witbpm.DiagOut.x_centroid)
        # Reset the corrector
        self.setter(self.orig_curr)
        self.bpmpos = array([[quadbpmpos], [witbpmpos]]).squeeze()

class BBA:
    def __init__(self,corsweep,quadind,shuntval):
        self.corsweep = corsweep
        self.quadind = quadind
        self.shuntval = shuntval
        self.guess = [0.,1.]

        self.setter = self.corsweep.serpobj.beamline[self.quadind].SetB
        self.getter = self.corsweep.serpobj.beamline[self.quadind].GetB
        self.origval = self.getter()

    def shunt(self):
        """Shunt the quad"""
        self.setter(self.origval * self.shuntval)
    
    def unshunt(self):
        """Restore the quad to its nominal value"""
        self.setter(self.origval)

    def run(self):
        """Run the BBA:
           1: Corrector sweep
           2: Quad shunt
           3: Corrector sweep
           4: Restore the quad"""
        self.corsweep.do_sweep()
        self.ini_sweep = self.corsweep.bpmpos
        self.shunt()
        self.corsweep.do_sweep()
        self.fin_sweep = self.corsweep.bpmpos
        self.unshunt()

# Perform a linear fit
def fitline(x, y, guess, fitfunc):
    errfunc = lambda guess, x0, y0: (y0 - fitfunc(guess, x0))**2
    p1, success = optimize.leastsq(errfunc, guess[:], args=(x,y))
    return p1

# Set variable 'mytwiss' from a file
execfile('ATF2extTwiss.py')

# Create the Serpentine object for entire ATF2
fullATF2 = Serpentine(line='newATF2lat.aml',twiss=mytwiss)

# Isolate the extraction from the ring using the element name,
# and create a new Serpentine object
ext_start = fullATF2.beamline.FindEleByName('KEX1A')
ATF2ext = Serpentine(line=beamline.Line(fullATF2.beamline[ext_start[0]:]),twiss=mytwiss)

# Plot the Twiss parameters
figure(10);ATF2ext.PlotTwiss()

# Zero all correctors to nominal, track a single particle, and plot results
ATF2ext.beamline.ZeroCors()
ATF2ext.Track()
figure(1); ATF2ext.PlotBPMReadings('b')

# Add an offset in 'x' to a quad, re-track, and re-plot
ATF2ext.beamline[194].SetOffset([0.1e-3,0,0,0,0,0])
ATF2ext.Track()
figure(1); ATF2ext.PlotBPMReadings('r')

# Generate new CorSweep and BBA objects
cor_currs = np.arange(-2e-3,2e-3,step=0.25e-3)
corsweep = CorSweep(ATF2ext, cor_currs, 190, [196, 202])
atfbba = BBA(corsweep,194,0.5)

# Do the BBA
atfbba.run()

# Define a function to fit to (a straight line!)
fitfunc = lambda guess, x0: guess[0] + guess[1] * x0

# Perform the fit to the two corrector sweeps (using the function
# defined earlier).
p1 = fitline(atfbba.ini_sweep[0,:], atfbba.ini_sweep[1,:], array([0,0]), fitfunc)
p2 = fitline(atfbba.fin_sweep[0,:], atfbba.fin_sweep[1,:], array([0,0]), fitfunc)

# Calculate the crossing point of the two straight lines
fit_offs = (p1[0] - p2[0]) / (p2[1] - p1[1])

# Calculate new data based on the fit
newX = np.linspace(atfbba.ini_sweep[0,0],atfbba.ini_sweep[0,-1],num=1000)
newY1 = fitfunc(p1,newX)
newY2 = fitfunc(p2,newX)

# Plot the data, the fit, and the crossing point
figure(2)
plot(atfbba.ini_sweep[0,:]*1e3,atfbba.ini_sweep[1,:]*1e3,'bx')
plot(atfbba.fin_sweep[0,:]*1e3,atfbba.fin_sweep[1,:]*1e3,'rx')
plot(newX*1e3,newY1*1e3,'b')
plot(newX*1e3,newY2*1e3,'r')
plot([fit_offs*1e3,fit_offs*1e3],ylim())
title('Measured offset = %0.3f mm' % (fit_offs*1e3))
xlabel('Quad BPM / mm')
ylabel('Witness BPM / mm')
show()

