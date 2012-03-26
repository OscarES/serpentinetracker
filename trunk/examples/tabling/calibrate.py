import tablewriter
import beamline
from serpentine import Serpentine
import numpy as np
from elements import BPM
from elements import Drift
from elements import REF

class MoverCalibration:
    """Calibration class used to write a table of calibration data
    
    arguments:
    *beamline* - Serpentine Beamline instance
    *bpm* - name of BPM to be calibrated
    
    keyword arguments:
    *dir* - direction of calibration: 'x' or 'y'
    *range* - half range of mover scan
    *pulses* - number of pulses per step in the calibration
    *steps* - number of steps in the mover scan"""

    def __init__(self,beamline,bpm,dir='x',range=250.e-6,pulses=10,steps=5):

        self.bpm = beamline.GetEleByName(bpm)[0]
        try: self.quad = beamline.GetEleByName(bpm[1:])[0]
        except AttributeError: self.quad = beamline.GetEleByType('Drift')[0]
        self.range = range
        self.pulses = pulses
        self.pcount = 0
        self.scount = 0
        self.steps = steps
        self.dir = 'xy'.find(dir)
        self.mpos = np.array([-1*self.range]
                             +(self.steps-1)*[2*self.range/float(self.steps-1)]
                             +[-1*self.range])

    def calibrate(self):

        table = tablewriter.TableWriter(ATF2ext,eltlist=['REF','BPM'],
                                        pramlist=['name','S','x','xp','y','yp','sigz','betax','offset_x','offset_y','DiagOut.x_centroid','DiagOut.y_centroid'],
                                        pulses = self.steps*self.pulses)
        self.makeMove()
        table.fill(self.makeMove)
        self.makeMove()
        table.writeData(self.bpm.name+'cal.dat')

    def makeMove(self):

        self.pcount += 1
        if self.pcount > self.scount*self.pulses:
            self.bpm.offset[2*self.dir] += self.mpos[self.scount]
            self.quad.offset[2*self.dir] += self.mpos[self.scount]
            self.scount += 1            

# Create the Serpentine object for entire ATF2
fullATF2 = Serpentine(line='ATF2.aml',twiss='atftwiss.txt')

# Isolate the extraction from the ring using the element name,
# and create a new Serpentine object
ext_start = fullATF2.beamline.FindEleByName('KEX1A')
ATF2ext = Serpentine(line=beamline.Line(fullATF2.beamline[ext_start[0]:]),twiss='atftwiss.txt')

# Replace all ICTs with REFs
for n in ATF2ext.beamline.GetEleByType('ICT'):
    index = ATF2ext.beamline.index(n)
    if 'REF' not in n.name: continue
    ATF2ext.beamline[index] = REF(name=n.name,L=n.L,P=n.P)

# Make drift smaller, Add IPBPMs, Add new drifts
idrift = ATF2ext.beamline.FindEleByName('L114C')[0]
ATF2ext.beamline[idrift].L = 1.444665
ipT1 = BPM('IPT1', 0)
ipD1 = Drift('IPD1', 0.076000)
ipT2 = BPM('IPT2', 0)
ipD2 = Drift('IPD2', 0.164000)
ipT3 = BPM('IPT3', 0)
ipD3 = Drift('IPD3', 0.076000)
ipT4 = BPM('IPT4', 0)
ipD4 = Drift('IPD4', 1.444665)

ATF2ext.beamline.insert(idrift+1, ipD4)
ATF2ext.beamline.insert(idrift+1, ipT4)
ATF2ext.beamline.insert(idrift+1, ipD3)
ATF2ext.beamline.insert(idrift+1, ipT3)
ATF2ext.beamline.insert(idrift+1, ipD2)
ATF2ext.beamline.insert(idrift+1, ipT2)
ATF2ext.beamline.insert(idrift+1, ipD1)
ATF2ext.beamline.insert(idrift+1, ipT1)

ATF2ext.beamline.SetSPos()
ATF2ext.TwissProp()

mc = MoverCalibration(ATF2ext.beamline,'MQM16FF','y')
mc.calibrate()
