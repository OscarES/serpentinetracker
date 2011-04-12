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
# Define the physics functions and classes (e.g. tracking, Rmat calcs, etc.)

# numpy arrays will be useful here
from numpy import *
from matplotlib.pylab import *
from elements import *
from scipy import weave

# ===============================================================
# Lists of beamline components are almost, but not quite, the right
# way to define the lattice.  Here we define a class to inherit from
# list, but with the multiplication operator redefined to do what
# we want
# The physics tools will be added on as methods of this new class
class Line(list):
    """class Line:  A class designed to hold the list of elements that define
    the beamline lattice.  This class inherits from Python's built-in 'list' class."""

    def __mul__(self,fact):
        """Allows multiplication of a small lattice subset by an integer in order to easily
        define a repeated section"""
        new_Line = Line()
        for f in range(fact): new_Line.extend(deepcopy(self))
        return new_Line

    def __repr__(self):
        return '\n'.join(str(ele.name)+" :: "+str(ele.__class__) for ele in self)

    def FindEleByName(self,name):
        """Returns the indices at which elements named 'name' can be found in self."""
        borkedlist = self._FindEleByName(name)
        fixedlist = list()
        for i in borkedlist:
            fixedlist.append(FixBorkedList(i))
        return fixedlist

    def _FindEleByName(self,name):
        from serpentine import Serpentine
        indlist = list()
        for i in range(len(self)):
            if isinstance(self[i],Serpentine):
                intern_list = self[i].beamline.FindEleByName(name)
                try:
                    for int_i in intern_list:
                        indlist.append([i,int_i])
                except TypeError:
                    indlist.append(inter_list)
            elif self[i].name == name:
                indlist.append(i)
        return indlist

    def FindEleByType(self,classname):
        """Returns the indices at which elements of class 'classtype' can be found in self."""
        borkedlist = self._FindEleByType(classname)
        fixedlist = list()
        for i in borkedlist:
            fixedlist.append(FixBorkedList(i))
        return fixedlist

    def _FindEleByType(self,classname):
        from serpentine import Serpentine
        indlist = list()
        for i in range(len(self)):
            if isinstance(self[i],Serpentine):
                intern_list = self[i].beamline._FindEleByType(classname)
                for int_i in intern_list:
                    indlist.append([i,int_i])
            elif self[i].__class__.__name__ == classname:     
                indlist.append(i)
        return indlist

    def GetEleByType(self,classname):
        """Returns a list of elements of class 'classtype' from self.  This returns the elements
        themselves, not their indices."""
        elems = list() 
        indlist = self.FindEleByType(classname)
        for i in indlist:
            if not isinstance(i,list):
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
        from serpentine import Serpentine
        for i in range(len(self)):
            if isinstance(self[i],Serpentine):
                intern_list = self[i].beamline.FindEleByObj(obj)
                eledict = dict()
                eledict[i] = intern_list
                return eledict
            if obj == self[i] :
                return i
        return -1

    def GetEleByName(self,name):
        """Returns a list of elements named 'name' from self.  This returns the elements
        themselves, not their indices."""
        elems = list() 
        indlist = self.FindEleByName(name)
        for i in indlist:
            if not isinstance(i,list):
                elems.append(self[i])
            else:
                obj = self
                for depth in xrange(len(i)-1):
                    obj = obj[i[depth]].beamline
                obj = obj[i[-1]]
                elems.append(obj)
        return elems

    def RmatAtoB(self, first, last):
        """Returns the 6D R-matrix between the entrance of self[first], and the exit of self[last]."""
        R = eye(6)
        for i in self[first:last+1]:
            R = dot(i.R, R)
        return R
    
    def Track(self,Beam):
        """The primary tracking method for the Line class.
        
        It loops around each element of self, calculates the offsets due to
        movers, offsets, etc., recalculates the energy variable of the beam
        being tracked to delta_E/E, and then calls the 'TrackThruEle' method
        of the element in question.
        
        Once tracking is complete for that element, the offset and the beam's
        energy variable are reset to their original values.
        
        The loop then continues on to the next element.
        
        The 'Beam' input should be an object of class 'ElectronBeam' (or one which
        inherits from that class).
        
        Track returns the beam that results from tracking through the lattice."""
        from serpentine import Serpentine
        prog = progressBar(0,len(self),77)
        beam_out = deepcopy(Beam)
        # Beam momentum is defined as absolute, whereas the R matrices expect delta_P/P
        for ele in self:
            if isinstance(ele, Serpentine):
                ele.beam_in = beam_out
                ele.Track()
                beam_out = ele.beam_out
                continue
            if sum(self.offset**2)>0:
                beam_out.x = self._AdjustBeamByLineOffset(ele,beam_out)
            if hasattr(ele,'Mover'):
                beam_out.x = self._AdjustBeamWithMover(ele,beam_out)
            if sum(ele.offset**2)>0:
                beam_out.x = self._AdjustBeamByOffset(ele,beam_out)
            if ele.is_diag:
                ele.Processor(beam_out)
            beam_out.x[5,:] = (beam_out.x[5,:] - ele.P) / ele.P
            beam_out = ele.TrackThruEle(beam_out)
            beam_out.x[5,:] = (beam_out.x[5,:] * ele.P) + ele.P
            if sum(ele.offset**2):
                beam_out.x = self._ReAdjustBeamByOffset(ele,beam_out)
            if hasattr(ele,'Mover'):
                beam_out.x = self._ReAdjustBeamWithMover(ele,beam_out)
            if sum(self.offset**2)>0:
                beam_out.x = self._ReAdjustBeamByLineOffset(ele,beam_out)
            prog.updateAmount(self.index(ele))
            print prog, "\r",
        return beam_out

    def _AdjustBeamByLineOffset(self,ele,beam_out):
        [r_in,r_out] = RotMats(-self.offset[5])
        line_length = self[-1].S - self[0].S
        dist_along_line = ele.S - self[0].S
        dist_normed = dist_along_line - (line_length/2) # dist from centre of the line
        delta_x = (dist_normed * self.offset[1]) + self.offset[0]
        delta_y = (dist_normed * self.offset[3]) + self.offset[2]
        delta_xp = self.offset[1]
        delta_yp = self.offset[3]
        beam_out.x[0,:] -= delta_x
        beam_out.x[1,:] -= delta_xp
        beam_out.x[2,:] -= delta_y
        beam_out.x[3,:] -= delta_yp
        beam_out.x = dot(r_in,beam_out.x)
        return beam_out.x

    def _ReAdjustBeamByLineOffset(self,ele,beam_out):
        [r_in,r_out] = RotMats(-self.offset[5])
        line_length = self[-1].S - self[0].S
        dist_along_line = ele.S - self[0].S
        dist_normed = dist_along_line - (line_length/2) # dist from centre of the line
        delta_x = (dist_normed * self.offset[1]) + self.offset[0]
        delta_y = (dist_normed * self.offset[3]) + self.offset[2]
        delta_xp = self.offset[1]
        delta_yp = self.offset[3]
        beam_out.x[0,:] += delta_x
        beam_out.x[1,:] += delta_xp
        beam_out.x[2,:] += delta_y
        beam_out.x[3,:] += delta_yp
        beam_out.x = dot(r_out,beam_out.x)
        return beam_out.x

    def _AdjustBeamByOffset(self,ele,beam_out):
        [r_in,r_out] = RotMats(-ele.offset[5])
        beam_out.x[0,:] = beam_out.x[0,:] - ele.offset[0]
        beam_out.x[1,:] = beam_out.x[1,:] - ele.offset[1]
        beam_out.x[2,:] = beam_out.x[2,:] - ele.offset[2]
        beam_out.x[3,:] = beam_out.x[3,:] - ele.offset[3]
        beam_out.x = dot(r_in,beam_out.x)
        return beam_out.x
    
    def _ReAdjustBeamByOffset(self,ele,beam_out):
        [r_in,r_out] = RotMats(-ele.offset[5])
        beam_out.x[0,:] = beam_out.x[0,:] + ele.offset[0]
        beam_out.x[1,:] = beam_out.x[1,:] + ele.offset[1]
        beam_out.x[2,:] = beam_out.x[2,:] + ele.offset[2]
        beam_out.x[3,:] = beam_out.x[3,:] + ele.offset[3]
        beam_out.x = dot(r_out,beam_out.x)
        return beam_out.x
    
    def _AdjustBeamWithMover(self,ele,beam_out):
        moverpos = ele.Mover.GetAct()
        [r_in,r_out] = RotMats(-moverpos[5])
        beam_out.x[0,:] = beam_out.x[0,:] - moverpos[0]
        beam_out.x[1,:] = beam_out.x[1,:] - moverpos[1]
        beam_out.x[2,:] = beam_out.x[2,:] - moverpos[2]
        beam_out.x[3,:] = beam_out.x[3,:] - moverpos[3]
        beam_out.x = dot(r_in,beam_out.x)
        return beam_out.x
    
    def _ReAdjustBeamWithMover(self,ele,beam_out):
        moverpos = ele.Mover.GetAct()
        [r_in,r_out] = RotMats(-moverpos[5])
        beam_out.x[0,:] = beam_out.x[0,:] + moverpos[0]
        beam_out.x[1,:] = beam_out.x[1,:] + moverpos[1]
        beam_out.x[2,:] = beam_out.x[2,:] + moverpos[2]
        beam_out.x[3,:] = beam_out.x[3,:] + moverpos[3]
        beam_out.x = dot(r_out,beam_out.x)
        return beam_out.x
    
    def SetSPos(self,ini_s=0):
        """Sets the longitudinal position of each element based on an initial value
        that defines the location of the upstream end of the first element (ini_s),
        and the length of each subsequent element."""
        from serpentine import Serpentine
        cum_s = ini_s
        for i in self:
            if isinstance(i,Serpentine):
                i.beamline.SetSPos(ini_s=cum_s)
                continue
            i.S = cum_s
            cum_s += i.L

    def TwissProp(self, ini_twiss):
        """Propagates an initial twiss object ('ini_twiss') through the lattice.

        For each element, the twiss calculated at its downstream end are stored as an
        attribute of that element.  The twiss output at the end of the lattice are
        returned from this function."""
        from serpentine import Serpentine
        from elements import Twiss
        import copy

        sum_phix,sum_phiy = 0,0
        final_twiss = copy.deepcopy(ini_twiss)
        finalgammax = (1+ini_twiss.alphax**2) / ini_twiss.betax
        finalgammay = (1+ini_twiss.alphay**2) / ini_twiss.betay
        
        for ele in self:
            ele.twiss = copy.deepcopy(final_twiss)

            if isinstance(ele, Serpentine):
                ele.TwissProp()
                continue

            det_x = det(ele.R[0:2,0:2])
            det_y = det(ele.R[2:4,2:4])

            deltaphix = arctan(ele.R[0,1] / \
                (final_twiss.betax*ele.R[0,0] - final_twiss.alphax*ele.R[0,1]))
            deltaphiy = arctan(ele.R[2,3] / \
                (final_twiss.betay*ele.R[2,2] - final_twiss.alphay*ele.R[2,3]))
            sum_phix += deltaphix
            sum_phiy += deltaphiy

            betax,alphax,gammax = final_twiss.betax,final_twiss.alphax,finalgammax
            betay,alphay,gammay = final_twiss.betay,final_twiss.alphay,finalgammay

            final_twiss.betax = (
                (ele.R[0,0]**2 * betax)  +  
                (-2*ele.R[0,0]*ele.R[0,1] * alphax)  +  
                (ele.R[0,1]**2 * gammax)
                )  /  det_x
            final_twiss.alphax = (
                (-ele.R[0,0]*ele.R[1,0] * betax)  +  
                ((ele.R[0,0]*ele.R[1,1] + ele.R[0,1]*ele.R[1,0]) * alphax)  +  
                (-ele.R[0,1]*ele.R[1,1] * gammax)
                )  /  det_x
            finalgammax = (1 + final_twiss.alphax**2) / final_twiss.betax

            final_twiss.betay = (
                (ele.R[2,2]**2 * betay)  +  
                (-2*ele.R[2,2]*ele.R[2,3] * alphay)  +  
                (ele.R[2,3]**2 * gammay)
                )  /  det_y
            final_twiss.alphay = (
                (-ele.R[2,2]*ele.R[3,2] * betay)  +  
                ((ele.R[2,2]*ele.R[3,3] + ele.R[2,3]*ele.R[3,2]) * alphay)  +  
                (-ele.R[2,3]*ele.R[3,3] * gammay)
                )  /  det_y
            finalgammay = (1 + final_twiss.alphay**2) / final_twiss.betay

            final_twiss.phix = sum_phix
            final_twiss.phiy = sum_phiy

        return final_twiss

    def ZeroCors(self):
        """Sets the field of all correctors in the lattice to zero.

        This is useful for reverting to the default lattice after a
        steering operation has been performed."""
        for ele in self:
            if (str(ele.__class__)=='elements.Xcor' or 
                str(ele.__class__)=='elements.Ycor' or
                str(ele.__class__)=='elements.XYcor'):
                ele.B = 0

    def SingleRmat(self, i):
        """Returns the already calculated R-matrix for beamline[i].
        i.e. it returns beamline[i].R."""
        return self[i].R

    def PlotRparam(self,param1=1,param2=1):
        """Plots the value of the R matrix element R[param1,param2] vs S position
        
        Note that param1 and param2 use 'Matlab-style' indexing, rather than 'Python-style'.
        i.e. they can be any integer between 1 and 6 inclusive."""
        S = zeros(len(self))
        Rparam = ones(len(self))
        for ele in self:
            S[self.index(ele)] = ele.S
            Rparam[self.index(ele)] = ele.R[param1-1,param2-1]
        plot(S,Rparam,'-x')

    def GetMomProfile(self):
        S = [ele.S for ele in self]
        P = [ele.P for ele in self]
        return (S,P)
    
    def PlotMomProfile(self,formatstr='-x'):
        S,P = self.GetMomProfile()
        plot(S,P,formatstr)

    def GetEkProfile(self,restmass):
        S  = [ele.S for ele in self]
        Ek = [sqrt(ele.P**2+restmass**2)-restmass for ele in self]
        return (S,Ek)

    def PlotEkProfile(self,restmass,formatstr='-x'):
        S,Ek = self.GetEkProfile(restmass)
        plot(S,Ek,formatstr)

    def XRmat(self,ind=0):
        """Print the 2x2 block of the R matrix corresponding to the horizontal transverse space.
        'ind' is the element for which the value is printed."""
        print self[ind].name + " x matrix:"
        print self[ind].R[0:2,0:2]

    def YRmat(self,ind=0):
        """Print the 2x2 block of the R matrix corresponding to the vertical transverse space.
        'ind' is the element for which the value is printed."""
        print self[ind].name + " y matrix:"
        print self[ind].R[2:4,2:4]

    def LongRmat(self,ind=0):
        """Print the 2x2 block of the R matrix corresponding to the longitudinal space.
        'ind' is the element for which the value is printed."""
        print self[ind].name + " longitudinal matrix:"
        print self[ind].R[4:6,4:6]

    def GetTwiss(self):
        from serpentine import Serpentine
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
            if isinstance(ele,Serpentine):
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

    def PlotTwiss(self,bx=1,by=1,ax=0,ay=0,px=0,py=0):
        twiss_dict = self.GetTwiss()

        numplots = (bx or by) + (ax or ay) + (px or py)
        for i in range(numplots):
            xstr = ''
            if (bx or by):
                if bx:
                    subplot(numplots,1,i+1); plot(twiss_dict['S'],twiss_dict['betax'],'-bx')
                    xstr = 'Beta_x / m  '
                    bx = 0
                if by:
                    subplot(numplots,1,i+1); plot(twiss_dict['S'],twiss_dict['betay'],'-rx')
                    xstr = xstr + '&  Beta_y / m'
                    by = 0
                xlabel('S / m')
                ylabel(xstr)
                legend(('Beta_x','Beta_y'),0)
                continue
            if (ax or ay):
                if ax:
                    subplot(numplots,1,i+1); plot(twiss_dict['S'],twiss_dict['alphax'],'-bx')
                    ax = 0
                if ay:
                    subplot(numplots,1,i+1); plot(twiss_dict['S'],twiss_dict['alphay'],'-rx')
                    ay = 0
                xlabel('S / m')
                continue
            if (px or py):
                if px:
                    subplot(numplots,1,i+1); plot(twiss_dict['S'],twiss_dict['phix'],'-bx')
                    px = 0
                if py:
                    subplot(numplots,1,i+1); plot(twiss_dict['S'],twiss_dict['phiy'],'-rx')
                    py = 0
                xlabel('S / m')
                continue

class progressBar:
    def __init__(self, minValue = 0, maxValue = 10, totalWidth=12):
        self.progBar = "[]"   # This holds the progress bar string
        self.min = minValue
        self.max = maxValue
        self.span = maxValue - minValue
        self.width = totalWidth
        self.amount = 0       # When amount == max, we are 100% done 
        self.updateAmount(0)  # Build progress bar string

    def updateAmount(self, newAmount = 0):
        if newAmount < self.min: newAmount = self.min
        if newAmount > self.max: newAmount = self.max
        self.amount = newAmount

        # Figure out the new percent done, round to an integer
        diffFromMin = float(self.amount - self.min)
        percentDone = (diffFromMin / float(self.span)) * 100.0
        percentDone = round(percentDone)
        percentDone = int(percentDone)

        # Figure out how many hash bars the percentage should be
        allFull = self.width - 2
        numHashes = (percentDone / 100.0) * allFull
        numHashes = int(round(numHashes))

        # build a progress bar with hashes and spaces
        self.progBar = "[" + '#'*numHashes + ' '*(allFull-numHashes) + "]"

        # figure out where to put the percentage, roughly centered
        percentPlace = (len(self.progBar) / 2) - len(str(percentDone)) 
        percentString = str(percentDone) + "%"

        # slice the percentage into the bar
        self.progBar = self.progBar[0:percentPlace] + percentString + \
            self.progBar[percentPlace+len(percentString):]

    def __str__(self):
        return str(self.progBar)

def FixBorkedList(borkedlist):
    buildlist = list()
    if isinstance(borkedlist,int):
        return borkedlist
    for i in borkedlist:
        if isinstance(i,int):
            buildlist.append(i)
        else:
            newlist = FixBorkedList(i)
            for newi in newlist:
                buildlist.append(newi)
    return buildlist
    
# A test suite
if __name__ == '__main__':
    from elements import *
    from beamrep import *
    from matplotlib.pylab import *
    import profile

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
    ini_twiss = {}
    ini_twiss['betax'] = 1
    ini_twiss['alphax'] = 0
    ini_twiss['betay'] = 2
    ini_twiss['alphay'] = 0

    final_twiss = beamline.TwissProp(ini_twiss)
    figure(1);PlotTwiss(final_twiss,ax=1,ay=1,px=1,py=1)

    print "Assigning beam..."
    beamin = GaussBeam(N=1e4)
    print "Starting tracking..."
    # profile.run('beamout = elements.Tracking(beamin)')
    profile.run('beamout = beamline.Track(beamin)')
    print "Done.  Now printing figures."
    figure(2)
    subplot(121);plot(beamin.x[0,:],beamin.x[1,:],'bx')
    subplot(121);plot(beamout.x[0,:],beamout.x[1,:],'r.')
    subplot(122);plot(beamin.x[2,:],beamin.x[3,:],'bx')
    subplot(122);plot(beamout.x[2,:],beamout.x[3,:],'r.')

    show()
