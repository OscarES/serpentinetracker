import tablewriter
import beamline
from serpentine import Serpentine
from elements import BPM
from elements import Drift

# Create the Serpentine object for entire ATF2
fullATF2 = Serpentine(line='ATF2.aml',twiss='atftwiss.txt')

# Isolate the extraction from the ring using the element name,
# and create a new Serpentine object
ext_start = fullATF2.beamline.FindEleByName('KEX1A')
ATF2ext = Serpentine(line=beamline.Line(fullATF2.beamline[ext_start[0]:]),twiss='atftwiss.txt')

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

table = tablewriter.TableWriter(ATF2ext,eltlist=['BPM'],pramlist=['name','S','x','xp','y','yp','sigz','betax'])
table.fill()
table.writeData('atf2res.dat')
