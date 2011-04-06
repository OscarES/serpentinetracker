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

def LoadLatFile(name=None):
    if name==None:
        raise "Must specify a filename."
    filestr = accformat.latticereader(name)
    beamline = LoadLat(filestr)
    return beamline

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
        positionroot = beamroot.getElementsByTagName('position')[-1]
    except IndexError:
        print "No beam node.  Using defaults."

    # Get the lattice node.  Resort to defaults if non-existent
    try:
        latticeroot = machineroot.getElementsByTagName('lattice')[-1]
        P = float(GetDesignFromEle('pc',latticeroot))/1e9
        print P
    except IndexError:
        print "No lattice node.  Using P = 1 GeV."
        P = 1

    # Now loop around the childNodes and extract each element
    beamline = Line()
    for i1 in trackingroot.childNodes:
        if i1.nodeType == i1.ELEMENT_NODE:
            beamline = MakeElement(i1,P,beamline)

    return beamline

def MakeElement(node,P,beamline):
    from globals import Brho1GeV
    numeles = len(beamline)
    Brho = Brho1GeV * P
    elename = node.getAttribute('name')
    for i1 in node.childNodes:
        if not i1.nodeType==i1.ELEMENT_NODE:
            i1.parentNode.removeChild(i1)

    try:
        L = float(GetDesignFromEle('length',node))
    except IndexError:
        L = 0

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
        elif i1.localName=='multipole':
            pass
        elif i1.localName=='octupole':
            pass
        elif i1.localName=='rf_cavity':
            pass
        elif i1.localName=='sextupole':
            B = ExtractB(i1,P)
            B *= L
            beamline.append(Sext(name=elename,L=L,P=P,B=B))
        elif i1.localName=='solenoid':
            pass
        elif i1.localName=='wiggler':
            pass
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
    if len(beamline)==numeles:
        beamline.append(Drift(name=elename,L=L,P=P))

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

if __name__ == '__main__':
    beamline = LoadFlatLatFile('tempfile.txt')
    print beamline
