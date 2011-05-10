#!/opt/local/bin/python
import beamline
import beamrep
import serpentine as st
from elements import *

def printres(beamline, inds):
    def printele(beamline, i):
        if type(i)==int:
            print beamline[i].name
        elif type(i[1])==int:
            print beamline[i[0]].beamline[i[1]].name
        else:
            printele(beamline[i[0]].beamline, i[1])
    for i in inds:
        printele(beamline, i)

def printbeamline(beamline, level=0):
    print "Level %i:" % level
    for ele in beamline:
        if type(ele)==st.Serpentine:
            level = printbeamline(ele.beamline, level=level-1)
        else:
            print "    " + ele.name
    print "Level %i:" % (level+1)
    return level+1

bl = beamline.Line()
bl.append(Drift(name='drift0', L=1))
ind = 1
bpmgrd0 = beamline.Line()
bpmgrd1 = beamline.Line()
bpmgrd2 = beamline.Line()
bpmgrd3 = beamline.Line()
for d in range(3):
    bpmgrd0.append(BPM(name='bpm'+str(ind),L=0,res=1e-7))
    bpmgrd0.append(Drift(name='drift'+str(ind),L=1))
    ind += 1

for d in range(3):
    bpmgrd1.append(Drift(name='drift'+str(ind),L=1))
    bpmgrd1.append(BPM(name='bpm'+str(ind),L=0,res=1e-7))
    ind += 1

for d in range(3):
    bpmgrd2.append(BPM(name='bpm'+str(ind),L=0,res=1e-7))
    bpmgrd2.append(Drift(name='drift'+str(ind),L=1))
    ind += 1

for d in range(3):
    bpmgrd3.append(BPM(name='bpm'+str(ind),L=0,res=1e-7))
    bpmgrd3.append(Drift(name='drift'+str(ind),L=1))
    ind += 1

bpmgrd2.append(st.Serpentine(line=bpmgrd3))
bpmgrd2.append(BPM(name='bpmtmp', L=0., res=1e-7))
bpmgrd1.append(st.Serpentine(line=bpmgrd2))

bl.append(st.Serpentine(line=bpmgrd0))
bl.append(Quad(name='clicquad',L=0.5,B=50))
bl.append(st.Serpentine(line=bpmgrd1))

bm = beamrep.Beam(P=11)
s = st.Serpentine(line=bl,beam=bm)

printbeamline(s.beamline)
print " "

qind = s.beamline.FindEleByName('drift11')
print qind
printres(s.beamline, qind)

qind = s.beamline.FindEleByName('bpmtmp')
print qind
printres(s.beamline, qind)

print s.beamline.GetEleByName('drift7')
print s.beamline.GetEleByName('bpmtmp')

