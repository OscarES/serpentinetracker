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
from serpentine import *
execfile('ATF2extTwiss.py')
fullATF2 = Serpentine(line='newATF2lat.aml',twiss=mytwiss)
ATF2ext = Serpentine(line=beamline.Line(fullATF2.beamline[947:]),twiss=mytwiss)
ATF2ext.TwissGaussBeam(N=1000)
ATF2ext.beamline.ZeroCors()
ATF2ext.Track()
#ATF2ext.beamline.FindEleByType('BPM')

#newline    = beamline.Line(ATF2ext.beamline[:26])
#small_line = beamline.Line(ATF2ext.beamline[26:29])
#end_line   = beamline.Line(ATF2ext.beamline[29:])
#
#small_serp = Serpentine(line=small_line)
#newline.append(small_serp)
#newline = beamline.Line(newline+end_line)
#
#test_serp = Serpentine(line = newline,twiss=mytwiss)
#bpminds = test_serp.beamline.FindEleByType('BPM')
#test_serp.TwissGaussBeam(N=1000)
#test_serp.beamline.ZeroCors()
#test_serp.Track()
#readings = test_serp.GetBPMReadings()
#
