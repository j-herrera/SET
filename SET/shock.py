import numpy as np

import SET.wrapper as wp

class shock(object):
    def __init__(self, name, geo, flow, model=0, cut_angle=90, multi_shock=0, eps=0.05, minpts=10):
        # creates cgal shock
        inputs = wp.optMap()
        inputs['gamma'] = flow['gamma']
        inputs['Mach'] = flow['mach']
        inputs['model'] = model
        inputs['cut_angle'] = cut_angle
        inputs['multi_shock'] = multi_shock
        inputs['eps'] = eps
        inputs['minpts'] = minpts

        cgal_geo = self.vf2cgal(geo)

        sh = cgal_geo.shock(inputs, (-flow['vector']).tolist())
        self._cpp_struct = sh

        # extracts faces and vertices
        nclusters = sh.return_nclusters()
        self._shocks = []

        for i in range(nclusters):
            fcs = np.array(sh.return_faces(i))[:, np.newaxis]
            vcs = np.array(sh.return_vertices(i))[:, np.newaxis]
            self._shocks.append(np.array([fcs.astype(int).reshape((int(fcs.size / 3), 3)), vcs.reshape((int(vcs.size / 3), 3))]))

        # defines class properties
        self._name = name + '_shock'
        self._flow = flow

    # --- defines getter methods ---
    @property
    def shocks(self):
        return self._shocks

    @property
    def name(self):
        return self._name

    @property
    def flow(self):
        return self._flow

    @property
    def cpp_struct(self):
        return self._cpp_struct

    def vf2cgal(self, geo):
        # creates a list of vertices and faces for CGAL
        cgal_vertices = geo['vertices'].reshape((1, geo['vertices'].size))[0].tolist()
        cgal_faces = geo['faces'].reshape((1, geo['faces'].size))[0].tolist()

        return wp.CGAL_geometry(cgal_vertices, cgal_faces)
