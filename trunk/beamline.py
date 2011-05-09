#!/usr/bin/python
#
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
"""Define the physics functions and classes 
(e.g. tracking, Rmat calcs, etc.)"""

# numpy arrays will be useful here
#from numpy import *
import numpy as np
from matplotlib.pylab import plot, subplot, xlabel, ylabel, legend
#from elements import *
#from scipy import weave
from utilities import RotMats
import beamadjust
import copy
from itertools import repeat
import re

# ===============================================================
# Lists of beamline components are almost, but not quite, the right
# way to define the lattice.  Here we define a class to inherit from
# list, but with the multiplication operator redefined to do what
# we want
# The physics tools will be added on as methods of this new class
class Line(list):
    """class Line:  A class designed to hold the list of elements that define
    the beamline lattice.  This class inherits from Python's built-in 'list' 
    class."""

    def __mul__(self, fact):
        """Allows multiplication of a small lattice subset by an integer in 
        order to easily define a repeated section"""
        new_line = Line()
        copyfunc = lambda x: new_line.extend(copy.deepcopy(x))
        for rep in repeat(copyfunc, rep):
            rep(self)
        #for rep in range(fact):
        #    new_line.extend(copy.deepcopy(self))
        return new_line

    def __repr__(self):
        def namecatch(inst):
            try: return str(inst.name)
            except AttributeError: return "No name attr"

        ret = '\n'.join(namecatch(ele)+" :: "+str(ele.__class__) for ele in self)
        return ret

    def FindEleByName(self, name):
        """Returns the indices at which elements named 'name' can be found in 
        self."""
        borkedlist = self._FindEleByName(name)
        fixedlist = list()
        for i in borkedlist:
            fixedlist.append(fixborkedlist(i))
        return fixedlist

    def _FindEleByName(self, name):
        """A method to help find elements by their name.
        This method is only to be called by FindEleByName()."""
        p = re.compile("^" + name +"$")
        indlist = list()
        for i in range(len(self)):
            if self[i].__class__.__name__ == 'Serpentine':
                intern_list = self[i].beamline.FindEleByName(name)
                try:
                    for int_i in intern_list:
                        indlist.append([i, int_i])
                except TypeError:
                    indlist.append(intern_list)
            elif p.match(self[i].name):
                indlist.append(i)
        if indlist:
            return indlist
        else:
            raise ValueError(name + ": Not found.")

    def FindEleByType(self, classname):
        """Returns the indices at which elements of class 'classtype' can be 
        found in self."""
        borkedlist = self.FindEleByTypeinit(classname)
        fixedlist = list()
        for i in borkedlist:
            fixedlist.append(fixborkedlist(i))
        return fixedlist

    def FindEleByTypeinit(self, classname):
        """A method to help find elements by their type.
        This method is only to be called by FindEleByType()."""
        p = re.compile("^" + classname + "$")
        indlist = list()
        for i in range(len(self)):
            if self[i].__class__.__name__ == 'Serpentine':
                try: intern_list = self[i].beamline.FindEleByTypeinit(classname)
                except ValueError: pass
                try: [indlist.append([i, int_i]) for int_i in intern_list]
                except UnboundLocalError: pass
            elif p.match(self[i].__class__.__name__):
                indlist.append(i)
        if indlist:
            return indlist
        else:
            raise ValueError(classname + ": Not found.")

    def GetEleByType(self, classname):
        """Returns a list of elements of class 'classtype' from self.  This 
        returns the elements themselves, not their indices."""
        elems = list() 
        indlist = self.FindEleByType(classname)
        for i in indlist:
            if not isinstance(i, list):
                elems.append(self[i])
            else:
                obj = self
                for depth in xrange(len(i)-1):
                    obj = obj[i[depth]].beamline
                obj = obj[i[-1]]
                elems.append(obj)
        return elems
            
    def FindEleByObj(self, obj):
        """Returns the index at which the object 'obj' can be found in self."""
        for i in range(len(self)):
            if self[i].__class__.__name__ == 'Serpentine':
                intern_list = self[i].beamline.FindEleByObj(obj)
                eledict = dict()
                eledict[i] = intern_list
                return eledict
            if obj == self[i] :
                return i
        return -1

    def GetEleByName(self, name):
        """Returns a list of elements named 'name' from self.  This returns 
        the elements themselves, not their indices."""
        elems = list() 
        indlist = self.FindEleByName(name)
        for i in indlist:
            if not isinstance(i, list):
                elems.append(self[i])
            else:
                obj = self
                for depth in xrange(len(i)-1):
                    obj = obj[i[depth]].beamline
                obj = obj[i[-1]]
                elems.append(obj)
        return elems

    def RmatAtoB(self, first, last):
        """Returns the 6D R-matrix between the entrance of self[first], and 
        the exit of self[last]."""
        rmat = np.eye(6)
        for i in self[first:last+1]:
            rmat = np.dot(i.R, rmat)
        return rmat
    
    def Track(self, beam, offset=np.array([0, 0, 0, 0, 0, 0])):
        """The primary tracking method for the Line class.
        
        It loops around each element of self, calculates the offsets due to
        movers, offsets, etc., recalculates the energy variable of the beam
        being tracked to delta_E/E, and then calls the 'TrackThruEle' method
        of the element in question.
        
        Once tracking is complete for that element, the offset and the beam's
        energy variable are reset to their original values.
        
        The loop then continues on to the next element.
        
        The 'beam' input should be an object of class 'ElectronBeam' (or one 
        which inherits from that class).
        
        Track returns the beam that results from tracking through the lattice.
        """
        prog = ProgressBar(0, len(self), 77)
        beam_out = copy.deepcopy(beam)
        # Beam momentum is defined as absolute, but R matrices expect delta_P/P
        for ele in self:
            if ele.__class__.__name__ == 'Serpentine':
                ele.beam_in = beam_out
                ele.Track()
                beam_out = ele.beam_out
                continue
            if sum(offset**2)>0:
                beam_out.x = self._AdjustBeamByLineOffset(ele, beam_out, offset)
            try:
                beam_out.x = beamadjust.AdjustBeamWithMover(ele, beam_out)
            except AttributeError: pass
            if sum(ele.offset**2)>0:
                beam_out.x = beamadjust.AdjustBeamByOffset(ele, beam_out)
            try:
                ele.Processor(beam_out)
            except AttributeError: pass
            beam_out.x[5, :] = (beam_out.x[5, :] - ele.P) / ele.P
            beam_out = ele.TrackThruEle(beam_out)
            beam_out.x[5, :] = (beam_out.x[5, :] * ele.P) + ele.P
            if sum(ele.offset**2):
                beam_out.x = beamadjust.ReAdjustBeamByOffset(ele, beam_out)
            if hasattr(ele, 'Mover'):
                beam_out.x = beamadjust.ReAdjustBeamWithMover(ele, beam_out)
            if sum(offset**2)>0:
                beam_out.x = self._ReAdjustBeamByLineOffset(
                    ele, beam_out, offset
                    )
            prog.updateAmount(self.index(ele))
            print prog,  "\r",
        return beam_out

    def _AdjustBeamByLineOffset(self, ele, beam_out, offset):
        """Correct the beam position by the offset specified for the entire
        beamline before the call to Track()"""
        r_in = RotMats(-offset[5])[0]
        line_length = self[-1].S - self[0].S
        dist_along_line = ele.S - self[0].S
        dist_normed = dist_along_line - (line_length/2) # dist from line centre
        delta_x = (dist_normed * offset[1]) + offset[0]
        delta_y = (dist_normed * offset[3]) + offset[2]
        delta_xp = offset[1]
        delta_yp = offset[3]
        beam_out.x[0, :] -= delta_x
        beam_out.x[1, :] -= delta_xp
        beam_out.x[2, :] -= delta_y
        beam_out.x[3, :] -= delta_yp
        beam_out.x = np.dot(r_in, beam_out.x)
        return beam_out.x

    def _ReAdjustBeamByLineOffset(self, ele, beam_out, offset):
        """Reset the beam position by the offset specified for the entire
        beamline after the call to Track()"""
        r_out = RotMats(-offset[5])[1]
        line_length = self[-1].S - self[0].S
        dist_along_line = ele.S - self[0].S
        dist_normed = dist_along_line - (line_length/2) # dist from line centre
        delta_x = (dist_normed * offset[1]) + offset[0]
        delta_y = (dist_normed * offset[3]) + offset[2]
        delta_xp = offset[1]
        delta_yp = offset[3]
        beam_out.x[0, :] += delta_x
        beam_out.x[1, :] += delta_xp
        beam_out.x[2, :] += delta_y
        beam_out.x[3, :] += delta_yp
        beam_out.x = np.dot(r_out, beam_out.x)
        return beam_out.x

    def SetSPos(self, ini_s=0):
        """Sets the longitudinal position of each element based on an initial 
        value that defines the location of the upstream end of the first 
        element (ini_s), and the length of each subsequent element."""
        cum_s = ini_s
        for i in self:
            if i.__class__.__name__ == 'Serpentine':
                i.beamline.SetSPos(ini_s=cum_s)
                continue
            i.S = cum_s
            cum_s += i.L

    def TwissProp(self, ini_twiss):
        """Propagates an initial twiss object ('ini_twiss') through the 
        lattice.

        For each element, the twiss calculated at its downstream end are 
        stored as an attribute of that element.  The twiss output at the 
        end of the lattice are returned from this function."""
        sum_phix, sum_phiy = 0, 0
        final_twiss = copy.deepcopy(ini_twiss)
        finalgammax = (1+ini_twiss.alphax**2) / ini_twiss.betax
        finalgammay = (1+ini_twiss.alphay**2) / ini_twiss.betay
        
        for ele in self:
            ele.twiss = copy.deepcopy(final_twiss)

            if ele.__class__.__name__ == 'Serpentine':
                ele.TwissProp()
                continue

            det_x = np.linalg.det(ele.R[0:2, 0:2])
            det_y = np.linalg.det(ele.R[2:4, 2:4])

            deltaphix = np.arctan(ele.R[0, 1] / \
                (final_twiss.betax*ele.R[0, 0] - 
                    final_twiss.alphax*ele.R[0, 1]))
            deltaphiy = np.arctan(ele.R[2, 3] / \
                (final_twiss.betay*ele.R[2, 2] - 
                    final_twiss.alphay*ele.R[2, 3]))
            sum_phix += deltaphix
            sum_phiy += deltaphiy

            betax  = final_twiss.betax
            alphax = final_twiss.alphax
            gammax = finalgammax
            betay  = final_twiss.betay
            alphay = final_twiss.alphay
            gammay = finalgammay

            final_twiss.betax = (
                (ele.R[0, 0]**2 * betax)  +  
                (-2*ele.R[0, 0]*ele.R[0, 1] * alphax)  +  
                (ele.R[0, 1]**2 * gammax)
                )  /  det_x
            final_twiss.alphax = (
                (-ele.R[0, 0]*ele.R[1, 0] * betax)  +  
                ((ele.R[0, 0]*ele.R[1, 1] + ele.R[0, 1]*ele.R[1, 0]) * 
                    alphax)  +  
                (-ele.R[0, 1]*ele.R[1, 1] * gammax)
                )  /  det_x
            finalgammax = (1 + final_twiss.alphax**2) / final_twiss.betax

            final_twiss.betay = (
                (ele.R[2, 2]**2 * betay)  +  
                (-2*ele.R[2, 2]*ele.R[2, 3] * alphay)  +  
                (ele.R[2, 3]**2 * gammay)
                )  /  det_y
            final_twiss.alphay = (
                (-ele.R[2, 2]*ele.R[3, 2] * betay)  +  
                ((ele.R[2, 2]*ele.R[3, 3] + ele.R[2, 3]*ele.R[3, 2]) * 
                    alphay)  +  
                (-ele.R[2, 3]*ele.R[3, 3] * gammay)
                )  /  det_y
            finalgammay = (1 + final_twiss.alphay**2) / final_twiss.betay

            final_twiss.phix = sum_phix
            final_twiss.phiy = sum_phiy

        return final_twiss

    def ZeroCors(self):
        """Sets the field of all correctors in the lattice to zero.
        This is useful for reverting to the default lattice after a
        steering operation has been performed."""
        import elements
        for ele in self:
            if (type(ele) == elements.Xcor or 
                type(ele) == elements.Ycor or
                type(ele) == elements.XYcor):
                ele.B = 0

    def SingleRmat(self, i):
        """Returns the already calculated R-matrix for beamline[i].
        i.e. it returns beamline[i].R."""
        return self[i].R

    def PlotRparam(self, param1=1, param2=1):
        """Plots the value of the R matrix element R[param1,param2] vs S 
        
        Note that param1 and param2 use 'Matlab-style' indexing, rather 
        than 'Python-style'. i.e. they can be any integer between 1 and 
        6 inclusive."""
        spos = np.zeros(len(self))
        rparam = np.ones(len(self))
        for ele in self:
            spos[self.index(ele)] = ele.S
            rparam[self.index(ele)] = ele.R[param1-1, param2-1]
        plot(spos, rparam, '-x')

    def GetMomProfile(self):
        """Returns the momentum profile of the reference particle"""
        spos = [ele.S for ele in self]
        mom  = [ele.P for ele in self]
        return (spos, mom)
    
    def PlotMomProfile(self, formatstr='-x'):
        """Plots the momentum profile of the reference particle"""
        spos, mom = self.GetMomProfile()
        plot(spos, mom, formatstr)

    def GetEkProfile(self, restmass):
        """Returns the kinetic energy profile of the reference particle"""
        spos    = [ele.S for ele in self]
        kenergy = [np.sqrt(ele.P**2+restmass**2)-restmass for ele in self]
        return (spos, kenergy)

    def PlotEkProfile(self, restmass, formatstr='-x'):
        """Plots the kinetic energy profile of the reference particle"""
        spos, kenergy = self.GetEkProfile(restmass)
        plot(spos, kenergy, formatstr)

    def GetRFPhases(self):
        """Returns the RF phases of the AccCav objects in beamline."""
        acccavs = self.GetEleByType('AccCav')
        return [ele.phi for ele in acccavs]
        
    def PlotRFPhases(self):
        """Plots the RF phases of the AccCav objects in beamline."""
        plot(self.GetRFPhases(), 'x')

    def XRmat(self, ind=0):
        """Print the 2x2 block of the R matrix corresponding to the 
        horizontal transverse space. 'ind' is the element for which the 
        value is printed."""
        print self[ind].name + " x matrix:"
        print self[ind].R[0:2, 0:2]

    def YRmat(self, ind=0):
        """Print the 2x2 block of the R matrix corresponding to the 
        vertical transverse space. 'ind' is the element for which the 
        value is printed."""
        print self[ind].name + " y matrix:"
        print self[ind].R[2:4, 2:4]

    def LongRmat(self, ind=0):
        """Print the 2x2 block of the R matrix corresponding to the 
        longitudinal space. 'ind' is the element for which the value is 
        printed."""
        print self[ind].name + " longitudinal matrix:"
        print self[ind].R[4:6, 4:6]

    def GetTwiss(self):
        """Returns a dictionary object containing the Twiss paramters
        calculated for the beamline."""
        twiss_dict = {}
        twiss_dict['S'] = []
        twiss_dict['betax'] = []
        twiss_dict['betay'] = []
        twiss_dict['alphax'] = []
        twiss_dict['alphay'] = []
        twiss_dict['phix'] = []
        twiss_dict['phiy'] = []
        twiss_dict['etax'] = []
        twiss_dict['etay'] = []
        twiss_dict['etapx'] = []
        twiss_dict['etapy'] = []
        for ele in self:
            if ele.__class__.__name__ == 'Serpentine':
                subtwiss_dict = ele.beamline.GetTwiss()
                twiss_dict['S'].extend(subtwiss_dict['S'])
                twiss_dict['betax'].extend(subtwiss_dict['betax'])
                twiss_dict['betay'].extend(subtwiss_dict['betay'])
                twiss_dict['alphax'].extend(subtwiss_dict['alphax'])
                twiss_dict['alphay'].extend(subtwiss_dict['alphay'])
                twiss_dict['phix'].extend(subtwiss_dict['phix'])
                twiss_dict['phiy'].extend(subtwiss_dict['phiy'])
                twiss_dict['etax'].extend(subtwiss_dict['etax'])
                twiss_dict['etay'].extend(subtwiss_dict['etay'])
                twiss_dict['etapx'].extend(subtwiss_dict['etapx'])
                twiss_dict['etapy'].extend(subtwiss_dict['etapy'])
            else:
                twiss_dict['S'].append(ele.S)
                twiss_dict['betax'].append(ele.twiss.betax)
                twiss_dict['betay'].append(ele.twiss.betay)
                twiss_dict['alphax'].append(ele.twiss.alphax)
                twiss_dict['alphay'].append(ele.twiss.alphay)
                twiss_dict['phix'].append(ele.twiss.phix)
                twiss_dict['phiy'].append(ele.twiss.phiy)
                twiss_dict['etax'].append(ele.twiss.etax)
                twiss_dict['etay'].append(ele.twiss.etay)
                twiss_dict['etapx'].append(ele.twiss.etapx)
                twiss_dict['etapy'].append(ele.twiss.etapy)
        return twiss_dict

    def PlotTwiss(self, betax=1, betay=1):
        """PlotTwiss(self, betax=1, betay=1)
        Plot the twiss parameters.
        if betax: plot Beta_x versus S
        if betay: plot Beta_y versus S"""
        twiss_dict = self.GetTwiss()

        numplots = (betax or betay)
        for i in range(numplots):
            xstr = ''
            if (betax or betay):
                if betax:
                    subplot(numplots, 1, i+1)
                    plot(twiss_dict['S'], twiss_dict['betax'], '-bx')
                    xstr = 'Beta_x / m  '
                    betax = 0
                if betay:
                    subplot(numplots, 1, i+1)
                    plot(twiss_dict['S'], twiss_dict['betay'], '-rx')
                    xstr = xstr + '&  Beta_y / m'
                    betay = 0
                xlabel('S / m')
                ylabel(xstr)
                legend(('Beta_x', 'Beta_y'), 0)
                continue

class ProgressBar:
    """A class to display a progress bar when tracking through a beamline."""
    def __init__(self, minvalue = 0, maxvalue = 10, totalwidth=12):
        self.progbar = "[]"   # This holds the progress bar string
        self.min = minvalue
        self.max = maxvalue
        self.span = maxvalue - minvalue
        self.width = totalwidth
        self.amount = 0       # When amount == max, we are 100% done 
        self.progbar = ""
        self.percentdone = 0
        self.updateAmount(0)  # Build progress bar string

    def updateAmount(self, new_amount = 0):
        """Calculate the percentage compled, and update the progbar string."""
        if new_amount < self.min:
            new_amount = self.min
        if new_amount > self.max:
            new_amount = self.max
        self.amount = new_amount
        self.percentDone()
        self.makestr()

    def percentDone(self):
        """Figure out the new percent done, round to an integer"""
        difffrommin = float(self.amount - self.min)
        percentdone = (difffrommin / float(self.span)) * 100.0
        self.percentdone = int(round(percentdone))

    def makestr(self):
        """Figure out how many hash bars the percentage should be"""
        allfull = self.width - 2
        numhashes = (self.percentdone / 100.0) * allfull
        numhashes = int(round(numhashes))

        # build a progress bar with hashes and spaces
        self.progbar = "[" + '#'*numhashes + ' '*(allfull-numhashes) + "]"

        # figure out where to put the percentage, roughly centered
        percentplace = (len(self.progbar) / 2) - len(str(self.percentdone)) 
        percentstring = str(self.percentdone) + "%"

        # slice the percentage into the bar
        self.progbar = self.progbar[0:percentplace] + percentstring + \
            self.progbar[percentplace+len(percentstring):]

    def __str__(self):
        return str(self.progbar)

def fixborkedlist(borkedlist):
    """A method to repair the broken lists returned by the find methods.
    This function should not be called by users."""
    buildlist = list()
    if isinstance(borkedlist, int):
        return borkedlist
    for i in borkedlist:
        if isinstance(i, int):
            buildlist.append(i)
        else:
            newlist = fixborkedlist(i)
            for newi in newlist:
                buildlist.append(newi)
    return buildlist
    
# A test suite
if __name__ == '__main__':
    from elements import Drift, Quad
    import beamrep
    import matplotlib.pylab as plt

    Shortline = Line()
    Shortline.append(Drift(name='ele1', L=0.75))
    Shortline.append(Quad(name='ele2', L=0.25, B=5))
    Shortline.append(Drift(name='ele3', L=1))
    Shortline.append(Quad(name='ele4', L=0.25, B=-5))
    beamline = Shortline * 5

#     print "="*20
#     print "    SingleRmat"
#     print "="*20
#     for i in range(0, len(beamline)):
#         print "%s: " % i,
#         print SingleRmat(beamline, i)

#     print
#     print "="*20
#     print "    RmatAtoB"
#     print "="*20
#     for i in range(0, len(beamline)):
#         print "%s: " % i,
#         print RmatAtoB(beamline, 0, i)

    print
    print "="*20
    print "    TwissProp"
    print "="*20
    i_twiss = {}
    i_twiss['betax'] = 1
    i_twiss['alphax'] = 0
    i_twiss['betay'] = 2
    i_twiss['alphay'] = 0

    f_twiss = beamline.TwissProp(i_twiss)
    plt.figure(1)
    beamline.PlotTwiss(f_twiss, ax=1, ay=1, px=1, py=1)

    print "Assigning beam..."
    beamin = beamrep.GaussBeam(N=1e4)
    print "Starting tracking..."
    # profile.run('beamout = elements.Tracking(beamin)')
    # profile.run('beamout = beamline.Track(beamin)')
    beamout = beamline.Track(beamin)
    print "Done.  Now printing figures."
    plt.figure(2)
    plt.subplot(121)
    plt.plot(beamin.x[0, :], beamin.x[1, :], 'bx')
    plt.subplot(121)
    plt.plot(beamout.x[0, :], beamout.x[1, :], 'r.')
    plt.subplot(122)
    plt.plot(beamin.x[2, :], beamin.x[3, :], 'bx')
    plt.subplot(122)
    plt.plot(beamout.x[2, :], beamout.x[3, :], 'r.')

    plt.show()
