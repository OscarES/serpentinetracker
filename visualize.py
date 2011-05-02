#    Copyright 2009,  Stephen Molloy
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

from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import numpy as np

def matplotlib2D(so, projection='sx', options = '') :    
    #####################################
    # Draw beam line 
    #####################################
    bl_verticies = []
    bl_codes     = [] 

    # first point
    bl_codes.append(Path.MOVETO)
    bl_verticies.append((so.beamline[0].S,0))
    
    for e in so.beamline[1:-2] :
        bl_codes.append(Path.LINETO)
        bl_verticies.append((e.S,0))
    
    # last point 
    bl_codes.append(Path.CLOSEPOLY)
    bl_verticies.append((so.beamline[-1].S,0))

    # make path patch
    bl_verticies = np.array(bl_verticies,float)
    bl_path      = Path(bl_verticies,bl_codes)
    bl_pathpatch = PathPatch(bl_path, facecolor='None', edgecolor = 'green')

    # plot and update
    axe = plt.gca()
    axe.add_patch(bl_pathpatch)
    axe.dataLim.update_from_data_xy(bl_verticies)
    axe.autoscale_view()
    plt.show()

    #####################################
    # Draw beam elements
    #####################################    
    for e in so.beamline :
        # swtich on each element type
        if   e.__class__.__name__ == "Quad" :
            if e.B > 0 :
                rect = Rectangle((e.S,0),e.L,0.05)
            else :
                rect = Rectangle((e.S,0),e.L,-0.05)
            axe.add_patch(rect)
        elif e.__class__.__name__ == "Sext" :
            if e.B > 0 :
                rect = Rectangle((e.S,0),e.L,0.05,facecolor='green')
            else :
                rect = Rectangle((e.S,0),e.L,-0.05,facecolor='green')
            axe.add_patch(rect)            
        elif e.__class__.__name__ == "Sbend" :
            rect = Rectangle((e.S,-0.025),e.L,0.05,facecolor='red')
            axe.add_patch(rect)

    plt.show()

    
