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
from beamline import *
from elements import *
from beamrep import *
import xml.dom.minidom as minidom
import accformat
from globals import Brho1GeV
from numpy import radians

def LoadLatFile(name=None, P=None):
    if name==None:
        raise "Must specify a filename."
    try:
        filestr = accformat.latticereader(name)
    except IOError:
        if name.split('.')[-1]=='dat':
            filestr = TraceWinStr(name).printstr()
        else: raise
    if not P==None:
        fileroot = minidom.parseString(filestr)
        machineroot = fileroot.getElementsByTagName('machine')[-1]
        latticeroot = machineroot.getElementsByTagName('lattice')[-1]
        pcnode = minidom.Document().createElement('pc')
        pcnode.setAttribute('design',str(P))
        latticeroot.appendChild(pcnode)
        filestr = fileroot.toprettyxml()

    beamline,parttype = LoadLat(filestr)
    return (beamline,parttype)

def LoadFlatLatFile(name=None):
    if name==None:
        raise "Must specify a filename."
    filehand = open(name)
    filestr = filehand.read()
    filehand.close()
    beamline = LoadLat(filestr)
    return beamline

def LoadLat(filestr):
    fileroot = minidom.parseString(filestr)
    positionroot,beamroot,latticeroot = None,None,None

    # Get the machine node.  It is an error if this doesn't exist
    try:
        machineroot = fileroot.getElementsByTagName('machine')[-1]
    except IndexError:
        raise "No machine node found in this file."

    # Get the tracking_lattice.  It is an error if this doesn't exist
    try:
        trackingroot = machineroot.getElementsByTagName('tracking_lattice')[-1]
    except IndexError:
        raise "No tracking_lattice node found in this file."

    # Get the beam node.  Resort to defaults if non-existent
    try:
        beamroot = machineroot.getElementsByTagName('beam')[-1]
        parttype = beamroot.getElementsByTagName('particle')[-1].getAttribute('type')
    except IndexError:
        parttype = 'ELECTRON'
        print "No beam node.  Using defaults."

    # Get the lattice node.  Resort to defaults if non-existent
    try:
        latticeroot = machineroot.getElementsByTagName('lattice')[-1]
        P = float(GetDesignFromEle('pc',latticeroot))/1e9
    except IndexError:
        print "No lattice node.  Using P = 1 GeV."
        P = 1.0

    # Now loop around the childNodes and extract each element
    beamline = Line()
    for i1 in trackingroot.childNodes:
        if i1.nodeType == i1.ELEMENT_NODE:
            beamline = MakeElement(i1,P,beamline,parttype)

    return (beamline,parttype)

def MakeElement(node,P,beamline,parttype):
    numeles = len(beamline)
    Brho = Brho1GeV * P
    elename = node.getAttribute('name')
    for i1 in node.childNodes:
        if not i1.nodeType==i1.ELEMENT_NODE:
            i1.parentNode.removeChild(i1)

    try: L = float(GetDesignFromEle('length',node))
    except IndexError: L = 0

    # Now loop over the elements and try to assign them
    for i1 in node.childNodes:
        if i1.localName=='quadrupole':
            B = ExtractB(i1,P)
            beamline.append(Quad(name=elename,L=L,P=P,B=B*L))
        elif i1.localName=='bend':
            B = ExtractSbend(i1,P,L)
            beamline.append(Sbend(name=elename,L=L,P=P,B=B*L))
        elif i1.localName=='kicker':
            xnode,ynode = None,None
            try:
                xnode = float(GetDesignFromEle('x_kick',i1))
                xnode *= Brho
                xnode *= L
            except IndexError:
                pass
            try:
                ynode = float(GetDesignFromEle('y_kick',i1))
                ynode *= Brho
                ynode *= L
            except IndexError:
                pass
            if not xnode==None and not ynode==None:
                beamline.append(XYcor(name=elename,L=L,P=P,B=([xnode,ynode])))
            elif not xnode==None:
                beamline.append(Xcor(name=elename,L=L,P=P,B=xnode))
            elif not ynode==None:
                beamline.append(Ycor(name=elename,L=L,P=P,B=ynode))
        elif i1.localName=='multipole': pass
        elif i1.localName=='octupole': pass
        elif i1.localName=='linac_cavity':
            gradient = float(GetDesignFromEle('gradient',i1))
            phase    = radians(float(GetDesignFromEle('phase0',i1)) + 90)
            freq     = float(GetDesignFromEle('rf_freq',i1))
            if parttype.upper()   =='PROTON'  : restmass = proton_mass
            elif parttype.upper() =='ELECTRON': restmass = electron_mass
            beamline.append(AccCav(name=elename,L=L,P=P,freq=freq,phi=phase,numdrift=2,egain=gradient*L,restmass=restmass))
        elif i1.localName=='sextupole':
            B = ExtractB(i1,P)
            B *= L
            beamline.append(Sext(name=elename,L=L,P=P,B=B))
        elif i1.localName=='solenoid': pass
        elif i1.localName=='wiggler': pass
        elif i1.localName=='instrument':
            if i1.getAttribute('type')=='INST' or i1.getAttribute('type')=='MONI':
                beamline.append(BPM(name=elename,L=L,P=P))
            elif i1.getAttribute('type')=='IMON':
                beamline.append(ICT(name=elename,L=L,P=P))
            elif i1.getAttribute('type')=='PROF':
                beamline.append(Screen(name=elename,L=L,P=P))
            elif i1.getAttribute('type')=='WIRE':
                beamline.append(WireScanner(name=elename,L=L,P=P))
    
    # If we haven't assigned anything yet, assign a drift
    if len(beamline)==numeles: beamline.append(Drift(name=elename,L=L,P=P))

    return beamline

def ExtractB(node,P):
    from globals import Brho1GeV
    k,ku = None,None
    try:
        k = GetDesignFromEle('k',node)
    except IndexError:
        try:
            ku = GetDesignFromEle('k_u',node)
        except IndexError:
            pass
    Brho = Brho1GeV * P
    if not k==None:
        B = float(k) * Brho
    elif not ku==None:
        B = float(ku)
    return B

def ExtractSbend(node,P,L):
    from globals import Brho1GeV
    B = array([0.0,0.0])
    g,gu = None,None
    g_node = node.getElementsByTagName('g')
    try:
        g = g_node[-1].getAttribute('design')
    except IndexError:
        gu_node = node.getElementsByTagName('g_u')
        try:
            gu = gu_node[-1].getAttribute('design')
        except IndexError:
            pass
    Brho = Brho1GeV * P
    if not g==None:
        B[0] = float(g) * Brho
    elif not gu==None:
        B[0] = float(gu)
    mult_node = node.getElementsByTagName('scaled_multipole')
    try:
        k_coef_node = mult_node[-1].getElementsByTagName('k_coef')
        k_coef = k_coef_node[-1].getAttribute('design')
        B[1] = float(k_coef) * B[0]
    except IndexError:
        pass
    return B

def GetDesignFromEle(elename,node):
    return node.getElementsByTagName(elename)[-1].getAttribute('design')

class TraceWinStr:
    def __init__(self,name):
        doc = minidom.Document()

        f = open(name)
        contents = f.readlines()
        data = [i.replace('\n','').split() for i in contents if not (i[0]==';' or i[0]=='\n')]

        machine = doc.createElement("machine")
        machine.appendChild(doc.createElement("tracking_lattice"))
        machine.appendChild(doc.createElement("master_list"))
        machine.appendChild(doc.createElement("beam"))
        machine.appendChild(doc.createElement("lattice"))
        self.machine = machine

        self.makebeam()

        self.freq = None

        for ele in data:
            if ele[0].upper()   == 'DRIFT':  self.makedrift(ele)
            elif ele[0].upper() == 'QUAD':   self.makequad(ele)
            elif ele[0].upper() == 'NCELLS': self.makecav(ele)
            elif ele[0].upper() == 'FREQ':   self.freq = float(ele[1]) * 1e6 # MHz-->Hz
            elif ele[0].upper() == 'END':    break
            elif ele[0].upper() == 'SET_ADV': pass
            elif ele[0].upper() == 'LATTICE': pass
            else: self.makedrift(ele)

    def makebeam(self):
        doc = minidom.Document()
        beamnode = self.machine.getElementsByTagName('beam')[-1]
        particlenode = doc.createElement("particle")
        particlenode.setAttribute('type','PROTON')
        beamnode.appendChild(particlenode)

    def makedrift(self,ele):
        tlatticenode = self.machine.getElementsByTagName("tracking_lattice")[-1]
        elenode = minidom.Document().createElement("element")
        elenode.setAttribute('name','drift')

        lengthnode = minidom.Document().createElement("length")
        lengthnode.setAttribute('design',str(float(ele[1])*1e-3))
        lengthnode.setAttribute('actual',str(float(ele[1])*1e-3))
        elenode.appendChild(lengthnode)

        tlatticenode.appendChild(elenode)

    def makequad(self,ele):
        self.makedrift(ele)
        tlatticenode = self.machine.getElementsByTagName("tracking_lattice")[-1]
        elenode = tlatticenode.childNodes[-1]
        elenode.setAttribute('name','quad')

        quadnode = minidom.Document().createElement("quadrupole")
        knode    = minidom.Document().createElement("k_u")
        knode.setAttribute('design',ele[1])
        knode.setAttribute('actual',ele[1])

        quadnode.appendChild(knode)
        elenode.appendChild(quadnode)
        tlatticenode.appendChild(elenode)
        
    def makecav(self,ele):
        if not ele[1]=='1': raise ValueError("Only pi-mode cavities are supported")
        doc = minidom.Document()

        Bgeo = float(ele[3])
        halfbetalambda = 0.5 * Bgeo * (c_light/self.freq)
        length = (halfbetalambda * int(ele[2])) + float(ele[10])/1e3 + float(ele[11])/1e3
        egain = float(ele[4])
        phi = float(ele[5])
        aper = float(ele[6])
        self.makedrift(ele)
        tlatticenode = self.machine.getElementsByTagName("tracking_lattice")[-1]
        elenode = tlatticenode.childNodes[-1]
        elenode.setAttribute('name','TM010cav')

        lengthnode = elenode.getElementsByTagName("length")[-1]
        lengthnode.setAttribute('design',str(length))
        lengthnode.setAttribute('actual',str(length))

        gradnode = doc.createElement("gradient")
        gradnode.setAttribute('design',str(egain/length))
        gradnode.setAttribute('actual',str(egain/length))

        freqnode = doc.createElement("rf_freq")
        freqnode.setAttribute('design',str(self.freq))
        freqnode.setAttribute('actual',str(self.freq))

        phasenode = doc.createElement("phase0")
        phasenode.setAttribute('design',str(phi))
        phasenode.setAttribute('actual',str(phi))

        betageonode = doc.createElement("beta_geo")
        betageonode.setAttribute('design',str(Bgeo))
        betageonode.setAttribute('actual',str(Bgeo))

        betapartnode = doc.createElement("beta_particle")
        betapartnode.setAttribute('design',ele[12])
        betapartnode.setAttribute('actual',ele[12])

        cavnode = doc.createElement("linac_cavity")
        cavnode.appendChild(freqnode)
        cavnode.appendChild(gradnode)
        cavnode.appendChild(phasenode)
        cavnode.appendChild(betageonode)
        cavnode.appendChild(betapartnode)

        elenode.appendChild(cavnode)

    def printstr(self):
        return self.machine.toprettyxml()

if __name__ == '__main__':
    from serpentine import Serpentine
    import matplotlib.pyplot as plt
    from scipy import sqrt
    execfile('ESSTwiss.py')
    in_energy = 50e6
    inP = sqrt((in_energy+proton_mass)**2 - proton_mass**2)
    blah = Serpentine(line="C.CQ.DCR.5.Step.dat",twiss=mytwiss,P=inP)
    blah.SetMomProfile()
    plt.figure()
    blah.PlotMomProfile(formatstr='-rx')
    blah.PlotEkProfile(formatstr='-bo')
    plt.figure()
    blah.PlotTwiss()
    plt.show()

