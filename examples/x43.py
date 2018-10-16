import sys
sys.path.append('../')

import numpy as np

from SET.shock import shock
from SET.io import read_ply, write_ply


filename = '../geometries/x43.ply'
geo = read_ply(filename)

flow = {'mach': 10, 'gamma': 1.4, 'vector':np.array([1,0,0])}

sh = shock('x43', geo, flow, model=2, multi_shock=1, cut_angle=5, eps=0.04, minpts=10)

for i in range(len(sh.shocks)):
    write_ply(sh.name + "_" + str(i) + ".ply", sh.shocks[i])
