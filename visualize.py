#    Copyright 2009,  Stephen Molloy
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

from matplotlib.path import Path
from matplotlib.patches import PathPatch
'''Module to control of the visualisation of beamline, beamrep and ...'''

from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate



class Visualize(object) :
    def __init__(self) :
        '''Class to visualize numerical data as function of beamline distance "S" '''
        self.xlim = None
        self.plotNames = []
        self.axes = []
        self.acb  = None

    # from serpentine
    def PlotBPMReadings(self,so, formatstr='', classname='BPM'):
        """Plot the BPM readings from the most recent tracking operation"""
        readings = so.GetBPMReadings(classname)
        plt.plot(readings[0, :], readings[1, :], '-rx')
        plt.plot(readings[0, :], readings[2, :], '-xb')
        plt.ylabel('x/y / m')

        self.PostPlot("PlotBPMReadings")
        
    # From beamline 
    def PlotRparam(self, so, param1=1, param2=1):
        """Plots the value of the R matrix element R[param1,param2] vs S 
        
        Note that param1 and param2 use 'Matlab-style' indexing, rather 
        than 'Python-style'. i.e. they can be any integer between 1 and 
        6 inclusive."""
        spos = np.zeros(len(so.beamline))
        rparam = np.ones(len(so.beamline))
        for ele in so.beamline:
            spos[so.beamline.index(ele)] = ele.S
            rparam[so.beamline.index(ele)] = ele.R[param1-1, param2-1]
        plt.plot(spos, rparam, '-x')
        plt.ylabel('R_'+str(param1)+str(param2))

        self.PostPlot("PlotRParam")

    # From beamline 
    def PlotMomProfile(self, so, formatstr='-x'):
        """Plots the momentum profile of the reference particle"""
        spos, mom = so.beamline.GetMomProfile()
        plt.plot(spos, mom, formatstr)
        plt.ylabel('P / GeV/c')

        self.PostPlot("PlotMomProfile")

    # From beamline (get data from serpentine and beam line... )
    def PlotEkProfile(self, so, formatstr='-x'):
        """Plots the kinetic energy profile of the reference particle"""
        spos, kenergy = so.beamline.GetEkProfile(so.beam_in.restmass)
        plt.plot(spos, kenergy, formatstr)

        self.PostPlot("PlotEkProfile")

    def PlotRFPhases(self, so):
        """Plots the RF phases of the AccCav objects in beamline."""
        so.plot(so.beamline.GetRFPhases(), 'x')

        self.PostPlot("PlotRFPhases")

    # From beamline 
    def PlotTwiss(self,so, betax=True, betay=True, spline=False) :
        """PlotTwiss(self, betax=True, betay=True, spline=False)
        Plot the twiss parameters.
        if betax: plot Beta_x versus S
        if betay: plot Beta_y versus S"""
        twiss_dict = so.beamline.GetTwiss()
    
        if betax :
            xarr = np.array(twiss_dict['S'])
            yarr = np.array(twiss_dict['betax'])
            if spline == False:
                plt.plot(xarr, yarr, '-bx')
            else :
                pass
            f= interpolate.InterpolatedUnivariateSpline(xarr,yarr,k=3)
            xarrp = np.linspace(xarr.min(),xarr.max(),1000)
        #                plt.plot(xarrp,f(xarrp),'-kx')
            xstr = 'Beta_x / m  '
            plt.plot(xarr,yarr,'-bx',label='Beta_x')
        if betay :
            xarr = np.array(twiss_dict['S'])
            yarr = np.array(twiss_dict['betay'])
            if spline == False :
                plt.plot(xarr, yarr, '-rx')
            else :
                pass
            f= interpolate.InterpolatedUnivariateSpline(xarr,yarr,k=3)
            xarrp = np.linspace(xarr.min(),xarr.max(),1000)
        #                plt.plot(xarrp,f(xarrp),'-kx')
            xstr = xstr + '&  Beta_y / m'
            plt.plot(xarr,yarr,'-rx',label='Beta_y')
             
        plt.ylabel('beta_{x,y}')
        plt.legend(loc=0)

        self.PostPlot("PlotTwiss")

    def Matplotlib2D(self,so, projection='sx', options = '', labelmag = False, labeldiag = False) :    
        '''Draw matplotlib representation of beamline.  
        so         : serpentine object (could be beamline)
        projection : 'sx','sy (no implemented yet'
        options    : undefined as yet
        label      : mark each element with its name
        return     : none'''
        
    #####################################
    # Draw beam line 
    #####################################
        bl_verticies = []
        bl_codes     = [] 

        eheight = 0.025

    # first point
        bl_codes.append(Path.MOVETO)
        bl_verticies.append((so.beamline[0].S,0))

        for e in so.beamline[1:] :
            bl_codes.append(Path.LINETO)
            bl_verticies.append((e.S,0))

        # last point 
        bl_codes.append(Path.CLOSEPOLY)
        bl_verticies.append((so.beamline[-1].S,0))

        # make path patch
        bl_verticies = np.array(bl_verticies,float)
        bl_path      = Path(bl_verticies,bl_codes)
        bl_pathpatch = PathPatch(bl_path, facecolor='None', edgecolor = 'green')

        # plot and update
        axe = plt.gca()
        axe.add_patch(bl_pathpatch)
        axe.dataLim.update_from_data_xy(bl_verticies)
        axe.autoscale_view()

        # set ranges 
        xmin = bl_verticies[:,0].min()
        xmax = bl_verticies[:,0].max()
        ymin = bl_verticies[:,1].min()
        ymax = bl_verticies[:,1].max()    
        xdiff = xmax-xmin
        axe.set_xlim(xmin-0.05*xdiff,xmax+0.05*xdiff)
        axe.set_ylim(ymin-eheight*4.5,ymax+eheight*4.5)

        #####################################
        # Draw beam elements
        #####################################    
        for e in so.beamline :
            # swtich on each element type
            textsloc = e.S+e.L/2.0
            if   e.__class__.__name__ == "Quad" :
                if e.B > 0 :
                    rect = Rectangle((e.S,0),e.L,eheight)
                    if labelmag : 
                        plt.text(textsloc,1.75*eheight,e.name,size=12, rotation=-90, ha="center",va="center", clip_on=True)
                else :
                    rect = Rectangle((e.S,0),e.L,-eheight)
                    if labelmag : 
                        plt.text(textsloc,-1.75*eheight,e.name,size=12, rotation=-90, ha="center",va="center", clip_on=True)
                axe.add_patch(rect)            
            elif e.__class__.__name__ == "Sext" :
                if e.B > 0 :
                    rect = Rectangle((e.S,0),e.L,eheight,facecolor='green')
                    if labelmag : 
                        plt.text(textsloc,3.5*eheight,e.name,size=12, rotation=-90, ha="center",va="center", clip_on=True)
                else :
                    rect = Rectangle((e.S,0),e.L,-eheight,facecolor='green')
                    if labelmag : 
                        plt.text(textsloc,-3.5*eheight,e.name,size=12, rotation=-90, ha="center",va="center", clip_on=True)
                axe.add_patch(rect)            
            elif e.__class__.__name__ == "Sbend" :
                rect = Rectangle((e.S,-eheight/2.0),e.L,eheight,facecolor='red')
                axe.add_patch(rect)
                if labelmag : 
                    plt.text(textsloc,3.5*eheight,e.name,size=12, rotation=-90, ha="center",va="center", clip_on=True)
            elif e.__class__.__name__ == "BasicDiag" :
                rect = Rectangle((e.S,-eheight/2.0),e.L,eheight,fill=False,ls='dashed')
                axe.add_patch(rect)            
                if labeldiag : 
                    plt.text(textsloc,-3.5*eheight,e.name,size=12, rotation=-90, ha="center",va="center", clip_on=True)
            elif e.__class__.__name__ == "BPM" :
                rect = Rectangle((e.S,-eheight/2.0),e.L,eheight,fill=False,ls='dashed')
                axe.add_patch(rect)            
                if labeldiag : 
                    plt.text(textsloc,-3.5*eheight,e.name,size=12, rotation=-90, ha="center",va="center", clip_on=True)
            elif e.__class__.__name__ == "Screen" :
                rect = Rectangle((e.S,-eheight/2.0),e.L,eheight,fill=False,ls='dashed')
                axe.add_patch(rect)            
                if labeldiag : 
                    plt.text(textsloc,-3.5*eheight,e.name,size=12, rotation=-90, ha="center",va="center", clip_on=True)
            elif e.__class__.__name__ == "EmitScreen" :
                rect = Rectangle((e.S,-eheight/2.0),e.L,eheight,fill=False,ls='dashed')
                axe.add_patch(rect)            
                if labeldiag : 
                    plt.text(textsloc,-3.5*eheight,e.name,size=12, rotation=-90, ha="center",va="center", clip_on=True)
            elif e.__class__.__name__ == "OTR" :
                rect = Rectangle((e.S,-eheight/2.0),e.L,eheight,fill=False,ls='dashed')
                axe.add_patch(rect)            
                if labeldiag : 
                    plt.text(textsloc,-3.5*eheight,e.name,size=12, rotation=-90, ha="center",va="center", clip_on=True)
            elif e.__class__.__name__ == "ICT" :
                rect = Rectangle((e.S,-eheight/2.0),e.L,eheight,fill=False,ls='dashed')
                axe.add_patch(rect)            
                if labeldiag : 
                    plt.text(textsloc,-3.5*eheight,e.name,size=12, rotation=-90, ha="center",va="center")
            elif e.__class__.__name__ == "Xcor" :
                pass 
            elif e.__class__.__name__ == "Ycor" : 
                pass
                
        # set axis labels etc
        axe.yaxis.set_ticklabels("")

        self.PostPlot("Matplotlib2D")

    def XAxesLabel(self) :
        plt.xlabel('S / m')
    
    def Update(self, cb=False) :
        """Update all axes limits from stored values""" 
        for a in self.axes :
            if a != self.acb :
                a.set_xlim(self.xlim)
            
        plt.show()

        # loop over all figures 
#        fnums = plt.get_fignums()
#        for f in fnums : 
#            f = plt.figure(f)
            # loop over all subplots
#            for a in f.get_axes() :
#                if a != self.acb :
#                    a.set_xlim(self.xlim)        
#                elif cb == False : 
                    # remove callback make the change and then reinstall callback.
#                    pass 
            
    def PostPlot(self, plotName = '') :
        # keep list of plots
        self.plotNames.append(plotName)

        # if no limits set some    
        if self.xlim == None : 
            self.UpdateLimits()        
        # apply consistent limits
        self.SetLimits()
        # keep axes for redrawing later
        self.AddAxes()
        # update all plots 
        self.Update()        

    def UpdateLimits(self) : 
        """Get current plot limits and store locally"""
        a = plt.gca()
        self.xlim = a.get_xlim()
        print self.xlim

    def SetLimits(self) :
        """Set the visualisation limits from the current axis"""
        a = plt.gca()
        a.set_xlim(self.xlim)

    def AddAxes(self) :
        self.axes.append(plt.gca())
    
    def ObserveAxes(self) :
        """Function to install axes change callback"""        
        self.acb = plt.gca()
        self.acb.callbacks.connect('xlim_changed',self.CallbackUpdate)
        
    def CallbackUpdate(self,ra) :
        self.xlim  = self.acb.get_xlim()
        self.Update(True)

def VisualizeTestRecursion() :
    import visualize
    import elements
    import serpentine
    import beamline 

    print 'visualize.VisualizeTestRecursion'
    
    # set twiss parameters
    mytwiss = elements.Twiss()
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
    
    qf  = elements.Quad("QF",L=0.25,P=1.25,B=5)
    dr1 = elements.Drift("D1",L=0.50,P=1.25)
    qd  = elements.Quad("QD",L=0.25,P=1.25,B=-5)
    dr2 = elements.Drift("D2",L=0.5,P=1.25) 
    m1  = elements.BasicDiag("M1",L=0.0)    

    fodo = beamline.Line([qf,dr1,qd,dr2,m1])
    fodo_sim = serpentine.Serpentine(line=fodo,twiss=mytwiss)
    
    vis = visualize.Visualize()
    plt.figure(1)
    plt.subplot(3,1,1)
    vis.Matplotlib2D(fodo_sim,labelmag=False, labeldiag=False)    
    vis.ObserveAxes()
    plt.subplot(3,1,2)
    vis.PlotTwiss(fodo_sim)    

    return fodo

def VisualizeTestATF2() :
    import visualize
    import elements
    import serpentine
    import beamline

    print 'visualize.VisualizeTestATF2()'
    
    # set twiss parameters
    mytwiss = elements.Twiss()
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
    atfFull = serpentine.Serpentine(line='./examples/atf/newATF2lat.aml',twiss=mytwiss)
    atfExt  = serpentine.Serpentine(line=beamline.Line(atfFull.beamline[947:]),twiss=mytwiss)
    
    # zero zero cors
    atfExt.beamline.ZeroCors()

    # Track 
    atfExt.Track()
    readings = atfExt.GetBPMReadings()
        
    vis = visualize.Visualize()
    plt.figure(1)
    plt.subplot(3,1,1)
    vis.Matplotlib2D(atfExt,labelmag=False, labeldiag=False)
    vis.ObserveAxes()
    plt.subplot(3,1,2)
    vis.PlotTwiss(atfExt)
    plt.subplot(3,1,3)
    vis.PlotBPMReadings(atfExt,'b')
    vis.XAxesLabel()

    plt.figure(2)
    plt.subplot(3,1,1)
    vis.PlotRparam(atfExt,1,1)
    plt.subplot(3,1,2)
    vis.PlotRparam(atfExt,2,2)
    plt.subplot(3,1,3)
    vis.PlotMomProfile(atfExt)
    vis.XAxesLabel()

    return vis

    
    
