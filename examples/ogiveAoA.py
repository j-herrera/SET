import sys
sys.path.append('../')

import numpy as np

from SET.shock import shock
from SET.io import read_ply, write_ply


filename = '../geometries/ogive.ply'
geo = read_ply(filename)

flow = {'mach': 9, 'gamma': 1.4, 'vector':np.array([np.cos(10*np.pi/180),np.sin(10*np.pi/180),0])}

sh = shock('ogiveAoA', geo, flow)

for i in range(len(sh.shocks)):
    write_ply(sh.name + "_" + str(i) + ".ply", sh.shocks[i])
