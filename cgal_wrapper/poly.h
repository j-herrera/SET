// Copyright (c) 2015 Javier Herrera All Rights Reserved.

#ifndef POLY_H
#define POLY_H

#include <boost/python.hpp>

#include "./types.h"
#include "./PolyhedralSurf.h"
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include "./shock.h"

#include <vector>
#include <array>

template <class HDS>
class Build_from_list : public CGAL::Modifier_base<HDS> {
typedef std::vector<std::array<double, 3> > Points_3;
typedef std::array<int,3> Facet;
typedef std::vector<Facet> Surface;
typedef typename HDS::Vertex::Point Point_3;

const Points_3 *pts;
const Surface *fcs;
std::map<int, Facet_handle> *map;

public:
Build_from_list(const Points_3 *plst, const Surface *flst, std::map<int, Facet_handle> *_map) : pts(plst), fcs(flst), map(_map) {
}

void operator() ( HDS& hds) {
  CGAL::Polyhedron_incremental_builder_3<HDS> B(hds);
  B.begin_surface( (*pts).size(), (*fcs).size());
  typedef typename Points_3::size_type size_type;

  for(size_type i=0; i < (*pts).size(); i++) {
    B.add_vertex(
      Point_3((*pts)[i][0], (*pts)[i][1], (*pts)[i][2])
      );
  }

  for(size_type i=0; i < (*fcs).size(); i++) {
    Facet_handle tfh;

    tfh = B.begin_facet();
    B.add_vertex_to_facet( (*fcs)[i][0]);
    B.add_vertex_to_facet( (*fcs)[i][1]);
    B.add_vertex_to_facet( (*fcs)[i][2]);
    B.end_facet();

    map->emplace(i, tfh);
  }
  if(B.error())
  {
    std::cerr << "An error occured while creating a Polyhedron" << std::endl;
    B.rollback();
  }

  B.end_surface();
}
};

struct shock_struct
{
  std::vector<std::vector<double> > vertices;
  std::vector<std::vector<double> > faces;
  unsigned int nclusters;
  std::vector<Polyhedron> structure;

  std::vector<double> return_vertices(unsigned int i) {
    return this->vertices[i];
  }

  std::vector<double> return_faces(unsigned int i) {
    return this->faces[i];
  }

  unsigned int return_nclusters() {
    return this->nclusters;
  }

  std::vector<Polyhedron> return_structure() {
    return this->structure;
  }
};

typedef Polyhedron::HalfedgeDS HalfedgeDS;
class CGAL_geometry {
public:
CGAL_geometry(const boost::python::list & vertex, const boost::python::list & faces){
  std::vector<std::array<double, 3> > v;
  std::vector<std::array<int, 3> > f;

  unsigned int nr = boost::python::len(vertex) / 3;
  for (unsigned int i = 0; i < nr; ++i) {
    std::array<double, 3> pts { {0., 0., 0.} };
    for (unsigned int j = 0; j < 3; ++j) {
      pts[j] = boost::python::extract<double>(vertex[3 * i +j]);
    }
    v.push_back(pts);
  }

  nr = boost::python::len(faces) / 3;
  for (unsigned int i = 0; i < nr; ++i) {
    std::array<int, 3> fcs { {0, 0, 0} };
    for (unsigned int j = 0; j < 3; ++j) {
      fcs[j] = boost::python::extract<int>(faces[3 * i +j]);
    }
    f.push_back(fcs);
  }

  Build_from_list<HalfedgeDS> builder(&v, &f, &m);
  G.delegate( builder);

}

shock_struct shock(std::map<std::string, double> & opt, const boost::python::list & ul){
  double gamma =  opt["gamma"];
  double Mach = opt["Mach"];
  double fmd = opt["model"];
  double cut_angle = opt["cut_angle"];
  double msh = opt["multi_shock"];
  double eps = opt["eps"];
  int minpts = opt["minpts"];

  std::array<double, 3> fl;
  for (unsigned int j = 0; j < 3; ++j) {
    fl[j] = boost::python::extract<double>(ul[j]);
  }

  Vector u(fl[0], fl[1], fl[2]);

  Polyhedron Ss = G;

  int model = 0;
  if (fmd > 0.5 and fmd < 1.4) {
    model = 1;
  } else if (fmd > 1.5){
    model = 2;
  }

  bool multi_shock = false;
  if (msh > 0.5) {
    multi_shock = true;
  }

  std::vector<Polyhedron> S = compute_shock(Ss, gamma, Mach, u, model, cut_angle, multi_shock, eps, minpts);

  shock_struct tout;
  tout.nclusters = S.size();
  tout.structure = S;

  for(unsigned int j=0; j< S.size(); j++) {
    std::vector<double> out;
    std::map<Vertex_handle, int> mv;
    int vn = 0;
    for ( Vertex_iterator i = S[j].vertices_begin(); i != S[j].vertices_end(); ++i ) {
      mv.emplace(i, vn);
      out.push_back(i->point().x());
      out.push_back(i->point().y());
      out.push_back(i->point().z());
      vn++;
    }

    std::vector<double> out2;
    for ( Facet_iterator i = S[j].facets_begin(); i != S[j].facets_end(); ++i ) {
      Halfedge_handle h = i->halfedge();
      out2.push_back(mv[h->vertex()]);
      out2.push_back(mv[h->next()->vertex()]);
      out2.push_back(mv[h->prev()->vertex()]);
    }

    tout.vertices.push_back(out);
    tout.faces.push_back(out2);
  }

  return tout;
}

Polyhedron G;
std::map<int, Facet_handle> m;
};

#endif /* POLY_H */
