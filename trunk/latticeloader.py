#!/usr/bin/python
"""A module to allow importing/exporting of various lattice formats
into the base format of Serpentine."""

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
from beamline import Line
import elements
import xml.dom.minidom as minidom
import accformat
from globals import Brho1GeV, electron_mass, proton_mass, c_light
from numpy import radians, array

def loadlatfile(name=None, ini_mom=None):
    """def loadlatfile(name=None, ini_mom=None)
    Load a lattice file in any format understood by Serpentine.
    name is the filename of the lattice
    ini_mom is the initial momentum of the lattice, unless it is specified
        in the lattice file."""
    if name == None:
        raise ValueError("Must specify a filename.")
    try:
        filestr = accformat.latticereader(name)
    except IOError:
        if name.split('.')[-1] == 'dat':
            filestr = TraceWinStr(name).printstr()
        else: raise
    if not ini_mom == None:
        fileroot = minidom.parseString(filestr)
        machineroot = fileroot.getElementsByTagName('machine')[-1]
        latticeroot = machineroot.getElementsByTagName('lattice')[-1]
        pcnode = minidom.Document().createElement('pc')
        pcnode.setAttribute('design', str(ini_mom))
        latticeroot.appendChild(pcnode)
        filestr = fileroot.toprettyxml()

    beamline, parttype = loadlat(filestr)
    return (beamline, parttype)

def loadflatlatfile(name=None):
    """def loadflatlatfile(name=None)
    Loads a "flat file".  i.e. one that has already been expanded by UAP.
    name is the filename of the flat lattice."""
    if name == None:
        raise ValueError("Must specify a filename.")
    filehand = open(name)
    filestr = filehand.read()
    filehand.close()
    beamline = loadlat(filestr)
    return beamline

def loadlat(filestr):
    """Takes a lattice as an expanded AML string, and returns
    a beamline and a string specifying the particle type."""
    fileroot = minidom.parseString(filestr)
    beamroot, latticeroot = None, None

    # Get the machine node.  It is an error if this doesn't exist
    try:
        machineroot = fileroot.getElementsByTagName('machine')[-1]
    except IndexError:
        raise ValueError("No machine node found in this file.")

    # Get the tracking_lattice.  It is an error if this doesn't exist
    try:
        trackingroot = machineroot.getElementsByTagName('tracking_lattice')[-1]
    except IndexError:
        raise ValueError("No tracking_lattice node found in this file.")

    # Get the beam node.  Resort to defaults if non-existent
    try:
        beamroot = machineroot.getElementsByTagName('beam')[-1]
        partroot = beamroot.getElementsByTagName('particle')[-1]
        parttype = partroot.getAttribute('type')
    except IndexError:
        parttype = 'ELECTRON'
        print "No beam node.  Using defaults."

    # Get the lattice node.  Resort to defaults if non-existent
    try:
        latticeroot = machineroot.getElementsByTagName('lattice')[-1]
        mom = float(get_design('pc', latticeroot))/1e9
    except IndexError:
        print "No lattice node.  Using P = 1 GeV."
        mom = 1.0

    # Now loop around the childNodes and extract each element
    beamline = Line()
    for child in trackingroot.childNodes:
        if child.nodeType == child.ELEMENT_NODE:
            beamline = make_element(child, mom, beamline, parttype)

    return (beamline, parttype)

def make_element(node, mom, beamline, parttype):
    """make_element(node, mom, beamline, parttype)
    Takes a single XML node, a momentum (mom), a beamline, and a string
    indicating the particle type.
    The data in the node is then extracted to create a new element.  This
    is then appended to beamline.  The new beamline is returned."""
    numeles = len(beamline)
    rigidity = Brho1GeV * mom
    elename = node.getAttribute('name')
    for child in node.childNodes:
        if not child.nodeType == child.ELEMENT_NODE:
            child.parentNode.removeChild(child)

    try:
        length = float(get_design('length', node))
    except IndexError:
        length = 0

    # Now loop over the elements and try to assign them
    for child in node.childNodes:
        if child.localName == 'quadrupole':
            bfield = extract_bfield(child, mom)
            beamline.append(
                elements.Quad(name=elename, L=length, P=mom, B=bfield*length)
                )
        elif child.localName == 'bend':
            bfield = extract_sbend(child, mom)
            beamline.append(
                elements.Sbend(name=elename, L=length, P=mom, B=bfield*length)
                )
        elif child.localName == 'kicker':
            xnode, ynode = None, None
            try:
                xnode = float(get_design('x_kick', child))
                xnode *= rigidity
                xnode *= length
            except IndexError:
                pass
            try:
                ynode = float(get_design('y_kick', child))
                ynode *= rigidity
                ynode *= length
            except IndexError:
                pass
            if not xnode == None and not ynode == None:
                beamline.append(
                    elements.XYcor(
                        name=elename, L=length, P=mom, B=([xnode, ynode])
                        )
                    )
            elif not xnode == None:
                beamline.append(
                    elements.Xcor(name=elename, L=length, P=mom, B=xnode)
                    )
            elif not ynode == None:
                beamline.append(
                    elements.Ycor(name=elename, L=length, P=mom, B=ynode)
                    )
        elif child.localName == 'multipole':
            pass
        elif child.localName == 'octupole':
            pass
        elif child.localName == 'linac_cavity':
            if parttype.upper()   == 'PROTON'  :
                restmass = proton_mass
            elif parttype.upper() == 'ELECTRON':
                restmass = electron_mass
            beamline.append(
                elements.AccCav(
                    name     = elename,
                    L        = length,
                    P        = mom,
                    freq     = float(get_design('rf_freq', child)),
                    phi      = radians(float(get_design('phase0', child)) + 90),
                    numdrift = 2,
                    egain    = float(get_design('gradient', child)) * length,
                    restmass = restmass
                    )
                )
        elif child.localName == 'sextupole':
            bfield = extract_bfield(child, mom)
            bfield *= length
            beamline.append(
                elements.Sext(name=elename, L=length, P=mom, B=bfield)
                )
        elif child.localName == 'solenoid':
            pass
        elif child.localName == 'wiggler':
            pass
        elif child.localName == 'instrument':
            if child.getAttribute('type') == 'INST':
                beamline.append(elements.BPM(name=elename, L=length, P=mom))
            elif child.getAttribute('type') == 'MONI':
                beamline.append(elements.BPM(name=elename, L=length, P=mom))
            elif child.getAttribute('type') == 'IMON':
                beamline.append(elements.ICT(name=elename, L=length, P=mom))
            elif child.getAttribute('type') == 'PROF':
                beamline.append(elements.Screen(name=elename, L=length, P=mom))
            elif child.getAttribute('type') == 'WIRE':
                beamline.append(
                    elements.WireScanner(name=elename, L=length, P=mom)
                    )
    
    # If we haven't assigned anything yet, assign a drift
    if len(beamline) == numeles:
        beamline.append(elements.Drift(name=elename, L=length, P=mom))

    return beamline

def extract_bfield(node, mom):
    """Extracts the magnetic field from node, and performs the appropriate
    calculations to convert it to a Serpentine B field.
    This value is returned."""
    kval, kuval = None, None
    try:
        kval = get_design('k', node)
    except IndexError:
        try:
            kuval = get_design('k_u', node)
        except IndexError:
            pass
    rigidity = Brho1GeV * mom
    if not kval == None:
        bfield = float(kval) * rigidity
    elif not kuval == None:
        bfield = float(kuval)
    return bfield

def extract_sbend(node, mom):
    """Extracts the field strength of a sector bend."""
    bfield = array([0.0, 0.0])
    gval, guval = None, None
    g_node = node.getElementsByTagName('g')
    try:
        gval = g_node[-1].getAttribute('design')
    except IndexError:
        gu_node = node.getElementsByTagName('g_u')
        try:
            guval = gu_node[-1].getAttribute('design')
        except IndexError:
            pass
    rigidity = Brho1GeV * mom
    if not gval == None:
        bfield[0] = float(gval) * rigidity
    elif not guval == None:
        bfield[0] = float(guval)
    mult_node = node.getElementsByTagName('scaled_multipole')
    try:
        k_coef_node = mult_node[-1].getElementsByTagName('k_coef')
        k_coef = k_coef_node[-1].getAttribute('design')
        bfield[1] = float(k_coef) * bfield[0]
    except IndexError:
        pass
    return bfield

def get_design(elename, node):
    """Returns the design attribute from an element."""
    return node.getElementsByTagName(elename)[-1].getAttribute('design')

class TraceWinStr:
    """A class to allow extraction of a lattice from a TraceWine file."""
    def __init__(self, name):
        doc = minidom.Document()

        fobj = open(name)
        contents = fobj.readlines()
        data = [
            i.replace('\n', '').split() 
            for i in contents 
            if not (i[0] == ';' or i[0] == '\n')]

        machine = doc.createElement("machine")
        machine.appendChild(doc.createElement("tracking_lattice"))
        machine.appendChild(doc.createElement("master_list"))
        machine.appendChild(doc.createElement("beam"))
        machine.appendChild(doc.createElement("lattice"))
        self.machine = machine

        self.makebeam()

        self.freq = None

        for ele in data:
            if ele[0].upper()   == 'DRIFT':
                self.makedrift(ele)
            elif ele[0].upper() == 'QUAD':
                self.makequad(ele)
            elif ele[0].upper() == 'NCELLS':
                self.makecav(ele)
            elif ele[0].upper() == 'FREQ':
                self.freq = float(ele[1]) * 1e6 # MHz-->Hz
            elif ele[0].upper() == 'END':
                break
            elif ele[0].upper() == 'SET_ADV':
                pass
            elif ele[0].upper() == 'LATTICE':
                pass
            else:
                self.makedrift(ele)

    def makebeam(self):
        """Since this is for TraceWin files, assume protons.
        Make a beam node specifying this."""
        doc = minidom.Document()
        beamnode = self.machine.getElementsByTagName('beam')[-1]
        particlenode = doc.createElement("particle")
        particlenode.setAttribute('type', 'PROTON')
        beamnode.appendChild(particlenode)

    def makedrift(self, ele):
        """Makes a drift node in AML format"""
        tlatticenode = self.machine.getElementsByTagName("tracking_lattice")[-1]
        elenode = minidom.Document().createElement("element")
        elenode.setAttribute('name', 'drift')

        lengthnode = minidom.Document().createElement("length")
        lengthnode.setAttribute('design', str(float(ele[1])*1e-3))
        lengthnode.setAttribute('actual', str(float(ele[1])*1e-3))
        elenode.appendChild(lengthnode)

        tlatticenode.appendChild(elenode)

    def makequad(self, ele):
        """Makes a quad node in AML format"""
        self.makedrift(ele)
        tlatticenode = self.machine.getElementsByTagName("tracking_lattice")[-1]
        elenode = tlatticenode.childNodes[-1]
        elenode.setAttribute('name', 'quad')

        quadnode = minidom.Document().createElement("quadrupole")
        knode    = minidom.Document().createElement("k_u")
        knode.setAttribute('design', ele[1])
        knode.setAttribute('actual', ele[1])

        quadnode.appendChild(knode)
        elenode.appendChild(quadnode)
        tlatticenode.appendChild(elenode)
        
    def makecav(self, ele):
        """Makes a cavity node in AML format"""
        if not ele[1] == '1':
            raise ValueError("Only pi-mode cavities are supported")
        doc = minidom.Document()

        length = (
            (0.5 * float(ele[3]) * (c_light/self.freq) * int(ele[2])) + 
            float(ele[10])/1e3 + 
            float(ele[11])/1e3
            )
        egain = float(ele[4])
        self.makedrift(ele)
        tlatticenode = self.machine.getElementsByTagName("tracking_lattice")[-1]
        elenode = tlatticenode.childNodes[-1]
        elenode.setAttribute('name', 'TM010cav')

        lengthnode = elenode.getElementsByTagName("length")[-1]
        lengthnode.setAttribute('design', str(length))
        lengthnode.setAttribute('actual', str(length))

        apernode = doc.createElement("aperture")
        apernode.setAttribute("at", "BOTH")
        apernode.setAttribute("shape", "CIRCLE")
        apernode.appendChild(createdesignnode("xy_limit", ele[6]))

        cavnode = doc.createElement("linac_cavity")
        cavnode.appendChild(createdesignnode("rf_freq", self.freq))
        cavnode.appendChild(createdesignnode("gradient", egain/length))
        cavnode.appendChild(createdesignnode("phase0", ele[5]))
        cavnode.appendChild(createdesignnode("beta_geo", ele[3]))
        cavnode.appendChild(createdesignnode("beta_particle", ele[12]))

        elenode.appendChild(cavnode)

    def printstr(self):
        """Returns the machine node of AML (basically the whole thing)
        as a pretty-printed XML string."""
        return self.machine.toprettyxml()

def createdesignnode(name, val):
    """Returns an AML node with the specified name, and the design|actual
    attributes set to val."""
    doc = minidom.Document()
    node = doc.createElement(name)
    node.setAttribute("design", str(val))
    node.setAttribute("actual", str(val))
    return node

if __name__ == '__main__':
    from serpentine import Serpentine
    from elements import Twiss
    import matplotlib.pyplot as plt
    from scipy import sqrt, pi
    mytwiss = Twiss()
    mytwiss.betax  = 14.4755754890
    mytwiss.alphax = -4.3155742015
    mytwiss.etax   = 0
    mytwiss.etaxp  = 0
    
    mytwiss.betay  = 5.1520133054
    mytwiss.alphay = -0.4908428984
    mytwiss.etay   = 0
    mytwiss.etayp  = 0
    
    Ek = 49.5e6 # MeV
    m0 = proton_mass
    Et = Ek + m0
    gamma0 = Et/m0
    betagamma = sqrt(gamma0**2 - 1)
    
    mytwiss.nemitx = pi * 0.23e-006# * betagamma
    mytwiss.nemity = pi * 0.23e-006# * betagamma
    
    betaz  = 6.6758214228
    alphaz = 0.0163358142
    nemitz = pi * 0.31e-006# * betagamma
    
    mytwiss.sigz = sqrt((nemitz/betagamma) * betaz)
    mytwiss.sigP = sqrt((nemitz/betagamma) / betaz)
    mytwiss.pz_cor = 0
    
    lambda352 = c_light/352.21e6
    onedegat352 = lambda352/360

    IEnergy = 50e6
    inP = sqrt((IEnergy+proton_mass)**2 - proton_mass**2)
    blah = Serpentine(line="C.CQ.DCR.5.Step.dat", twiss=mytwiss, P=inP)
    blah.SetMomProfile()
    plt.figure()
    blah.PlotMomProfile(formatstr='-rx')
    blah.PlotEkProfile(formatstr='-bo')
    plt.figure()
    blah.PlotTwiss()
    plt.show()

