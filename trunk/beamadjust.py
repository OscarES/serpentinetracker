#!/usr/bin/python
#
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
from utilities import RotMats
import numpy as np

def AdjustBeamByOffset(ele, beam_out):
    [r_in, r_out] = RotMats(-ele.offset[5])
    beam_out.x[0, :] = beam_out.x[0, :] - ele.offset[0]
    beam_out.x[1, :] = beam_out.x[1, :] - ele.offset[1]
    beam_out.x[2, :] = beam_out.x[2, :] - ele.offset[2]
    beam_out.x[3, :] = beam_out.x[3, :] - ele.offset[3]
    beam_out.x = np.dot(r_in, beam_out.x)
    return beam_out.x

def ReAdjustBeamByOffset(ele, beam_out):
    [r_in, r_out] = RotMats(-ele.offset[5])
    beam_out.x[0, :] = beam_out.x[0, :] + ele.offset[0]
    beam_out.x[1, :] = beam_out.x[1, :] + ele.offset[1]
    beam_out.x[2, :] = beam_out.x[2, :] + ele.offset[2]
    beam_out.x[3, :] = beam_out.x[3, :] + ele.offset[3]
    beam_out.x = np.dot(r_out, beam_out.x)
    return beam_out.x

def AdjustBeamWithMover(ele, beam_out):
    moverpos = ele.Mover.GetAct()
    [r_in, r_out] = RotMats(-moverpos[5])
    beam_out.x[0, :] = beam_out.x[0, :] - moverpos[0]
    beam_out.x[1, :] = beam_out.x[1, :] - moverpos[1]
    beam_out.x[2, :] = beam_out.x[2, :] - moverpos[2]
    beam_out.x[3, :] = beam_out.x[3, :] - moverpos[3]
    beam_out.x = np.dot(r_in, beam_out.x)
    return beam_out.x

def ReAdjustBeamWithMover(ele, beam_out):
    moverpos = ele.Mover.GetAct()
    [r_in, r_out] = RotMats(-moverpos[5])
    beam_out.x[0, :] = beam_out.x[0, :] + moverpos[0]
    beam_out.x[1, :] = beam_out.x[1, :] + moverpos[1]
    beam_out.x[2, :] = beam_out.x[2, :] + moverpos[2]
    beam_out.x[3, :] = beam_out.x[3, :] + moverpos[3]
    beam_out.x = np.dot(r_out, beam_out.x)
    return beam_out.x

