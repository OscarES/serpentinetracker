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
# Test each element, in order to compare with Lucretia

from beamrep import *
from beamline import *
from elements import *
from matplotlib.pylab import *

mytwiss = Twiss()

mytwiss.betax = 17.2390
mytwiss.alphax = -3.2950
mytwiss.etax = 0
mytwiss.etapx = 0
mytwiss.phix = 0

mytwiss.betay = 17.1690
mytwiss.alphay = -3.2780
mytwiss.etay = 0
mytwiss.etapy = 0
mytwiss.phiy = 0

mytwiss.nemitx = 1e-6
mytwiss.nemity = 1e-6

mytwiss.sigz = 1e-1
mytwiss.sigP = 1.28e-6
mytwiss.pz_cor = 0

# mybeam = TwissGaussBeam(mytwiss, N=200, pos=array([50e-3,0,0,0,0,1.3]), Q=8e-10)
mybeam = TwissGaussBeam(mytwiss, N=5000, pos=array([0,0,0,0,0,1.3]), Q=8e-10)
# figure(1);mybeam.PlotXSpace()
# figure(2);mybeam.PlotYSpace()
# figure(3);mybeam.PlotLongSpace()
# show()

# mybeam = Beam(P=1.3)
# mybeam.x = array([
#     [0,  -1e-4,1e-4,0,   0,   0,   0,   0,   0,   0,   0   ],
#     [0,   0,   0,  -1e-4,1e-4,0,   0,   0,   0,   0,   0   ],
#     [0,   0,   0,   0,   0,  -1e-4,1e-4,0,   0,   0,   0   ],
#     [0,   0,   0,   0,   0,   0,   0,  -1e-4,1e-4, 0,  0   ],
#     [0,   0,   0,   0,   0,   0,   0,   0,   0,  -1e-4,1e-4],
#     [1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3, 1.3 ]
#     ])
# #mybeam.x = array([
#    [0,  -1e-4,1e-4],
#    [0,   0,   0  ],
#    [0,   0,   0  ],
#    [0,   0,   0  ],
#    [0,   0,   0  ],
#    [1.3, 1.3, 1.3]
#    ])
# mybeam.Q = ones(mybeam.x.shape[1]) * mybeam.Q/mybeam.x.shape[1]
# 
# beamline = Line([
#     Sbend(Name='sbend',P=1.3,B=0.2,L=1),
#     Sbend(Name='sbend_eangle46.12',P=1.3,B=0.2,L=1,EAngle=46.12e-3*ones(2)),
#     Sbend(Name='sbend_ECurve0.1',P=1.3,B=0.2,L=1,ECurve=0.1*ones(2)),
#     Quad(Name='quad_normal',P=1.3,B=1),
#     Quad(Name='quad_skew',P=1.3,B=1,Tilt=pi/4),
#     Xcor(Name='xcor',P=1.3,B=0.002,L=0.1),
#     Ycor(Name='ycor',P=1.3,B=0.002,L=0.1),
#     ThinSext(Name='thinsext_normal',P=1.3,B=1000000,L=0.1),
#     ThinSext(Name='thinsext_skew',P=1.3,B=1000000,L=0.1,Tilt=pi/8),
#     ])
# 
beamline = Line([
    AccCav(name='acc_cav',L=1.5,P=1.3,S=0,egain=20e6,designQ=1e9,freq=1.3e9,phi=pi/2)
    ])
# beamline[0].TSRW = Wake()
# l,k = load('srwf_lcls_tran.dat', usecols=(1,2), skiprows=3, unpack=True)
# beamline[0].TSRW.l = l
# beamline[0].TSRW.k = k

counter = 0
for ele in beamline:
    counter += 1
    newline = Line([ele])
    newline.offset = zeros(6)
    beamout = newline.Track(mybeam)
    
    figure(counter)
    subplot(121)
    plot(mybeam.x[0,:]*1e3,mybeam.x[1,:]*1e3,'rx')
    plot(beamout.x[0,:]*1e3,beamout.x[1,:]*1e3,'bo')
    title('x Phase space')
    xlabel('x / mm')
    ylabel("x' / m.rad")

    subplot(122)
    plot(mybeam.x[2,:]*1e3,mybeam.x[3,:]*1e3,'rx')
    plot(beamout.x[2,:]*1e3,beamout.x[3,:]*1e3,'bo')
    title('y Phase space')
    xlabel('y / mm')
    ylabel("y' / m.rad")

    figure(counter+100)
    plot(mybeam.x[4,:],mybeam.x[5,:],'rx')
    plot(beamout.x[4,:],beamout.x[5,:],'bo')

    show()
  # savefig(ele.Name + '.pdf')
    # close('all')

