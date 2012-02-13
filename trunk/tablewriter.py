import numpy as np
import serpentine
from elements import Twiss

class TableWriter:
    """Class useful for repeat tracking and writing the results to a data file
    required input:
    *sp* - Serpentine object
    
    keyword arguments
    *pramlist* - list of parameters to record every pulse
    *eltlist* - type of element for which to record the parameters in pramlist
    *pulses* - number of pulses"""

    xdict = {'x':0,'xp':1,'y':2,'yp':3,'z':4,'P':5}
    tlist = Twiss().__dict__.keys()

    def __init__(self,sp,pramlist=[],eltlist=['BPM'],pulses=10):

        self.serp = sp
        self.params = pramlist
        self.beam = self.serp.beam_in
        self.np = pulses
        sortlist = []
        for e in eltlist:
            sortlist.extend(self.serp.beamline.FindEleByType(e))
        self.eles = sortlist[:]
        for (i,s) in enumerate(sortlist):
            self.eles[i] = self.serp.beamline[s]

    def fill(self,callback=lambda:None):
        """Class method: perform tracking and collect data"""

        l = '# N_train\tN_bunch\tindex\ttype'
        for p in self.params: l += '\t'+p
        l += '\n'
        self.lines = [l]
        for n in range(self.np):
            self.serp.RefreshBeam()
            self.serp.Track()
            for b in range(len(self.serp.beam_in.x[0])):
                for e in self.eles:
                    l = str(n)+'\t'+str(b)+'\t'+str(self.serp.beamline.index(e))+'\t'+e.__class__.__name__+'\t'
                    for p in self.params:
                        if p in self.xdict: 
                            l += '\t'+str(e.x[self.xdict[p]][b])
                        elif p in self.tlist:
                            l += '\t'+str(e.twiss.__getattribute__(p))
                        else: 
                            l += '\t'+str(e.__getattribute__(p)) 
                    l += '\n'
                    self.lines.append(l)
            callback()
            
    def writeData(self,filename=''):
        """Write data to file
        *filename* - output filename including path"""
        self.filename = filename
        of = open(self.filename, 'w')
        of.writelines(self.lines)
        of.close()

