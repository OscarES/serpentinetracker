#!/opt/local/bin/python
import beamline
import beamrep
import serpentine as st
from elements import *

def printres(beamline, inds):
    try:
        for i in inds:
            if type(i)==int:
                print beamline[i].name
            else:
                printres(beamline[i[0]].beamline, i[1])
    except TypeError:
        print beamline[inds].name

bl = beamline.Line()
bl.append(Drift(name='drift0', L=1))
ind = 1
bpmgrd0 = beamline.Line()
bpmgrd1 = beamline.Line()
for d in range(3):
    bpmgrd0.append(BPM(name='bpm'+str(ind),L=0,res=1e-7))
    bpmgrd0.append(Drift(name='drift'+str(ind),L=1))
    ind += 1

for d in range(3):
    bpmgrd1.append(Drift(name='drift'+str(ind),L=1))
    bpmgrd1.append(BPM(name='bpm'+str(ind),L=0,res=1e-7))
    ind += 1

bl.append(st.Serpentine(line=bpmgrd0))
bl.append(Quad(name='clicquad',L=0.5,B=50))
bl.append(st.Serpentine(line=bpmgrd1))

bm = beamrep.Beam(P=11)
s = st.Serpentine(line=bl,beam=bm)

qind = s.beamline.FindEleByType('Drift')
printres(s.beamline, qind)

