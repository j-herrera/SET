import sys
sys.path.append('../')

import numpy as np

from SET.shock import shock
from SET.io import read_ply, write_ply


filename = '../geometries/blunted_elliptic_ogive.ply'
geo = read_ply(filename)

flow = {'mach': 9, 'gamma': 1.4, 'vector':np.array([1,0,0])}

sh = shock('blunted_elliptic_ogive', geo, flow, 1)

for i in range(len(sh.shocks)):
    write_ply(sh.name + "_" + str(i) + ".ply", sh.shocks[i])
