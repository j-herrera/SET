// Copyright (c) 2015 Javier Herrera All Rights Reserved.

#ifndef PATH_H
#define PATH_H

#include <math.h>
#include "kdtree.h"

#include <iostream>
#include <fstream>

#include <utility>
#include <vector>
#include <list>
#include <map>

#include "./types.h"
#include "./PolyhedralSurf.h"
#include "./interp.h"
#include "./ode.h"

struct output
{
  std::map<Facet_handle, std::pair<double, Point> > maps;
  std::vector<std::vector<Point> > streams;
};

output path(Polyhedron & G, Vector & Free_Stream, double max_angle) {
  // Parameters
  double kt = 0.1; // k mean min triangle edge / max_vel dt

  // Compute per VERTEX velocities and max velocity
  // Vector Free_Stream(-1.0, 0.0, 0.0);
  std::map< Point, Vector > vel;
  double max_vel = 0;

  for ( Vertex_iterator v = G.vertices_begin(); v != G.vertices_end(); ++v ) {
    HV_circulator j = v->vertex_begin();

    Vector velT(0.0, 0.0, 0.0);

    int i = 0;
    do {
      if (j->is_border()) {
        continue;
      }
      Halfedge_handle h = j->facet()->halfedge();

      // Whyyyy?
      if (CGAL::collinear(h->vertex()->point(),
                          h->next()->vertex()->point(),
                          h->prev()->vertex()->point())) {
        continue;
      }

      Vector norm = -CGAL::unit_normal(h->vertex()->point(),
                                       h->next()->vertex()->point(),
                                       h->prev()->vertex()->point());
      Vector vel_xyz = Free_Stream - (Free_Stream * norm) * norm;

      velT = velT - vel_xyz;
      i++;
    } while ( ++j != v->vertex_begin());

    vel.emplace(v->point(), (1.0/i) * velT);
    if (sqrt(((1.0/i) * velT).squared_length()) > max_vel) {
      max_vel = sqrt(((1.0/i) * velT).squared_length());
    }
  }

  // Declare output
  output result;

  // Compute centroids.  Filter back-face
  std::vector<Point> cent;
  std::vector<Facet_handle> fh;

  for ( Facet_iterator i = G.facets_begin(); i != G.facets_end(); ++i ) {
    Halfedge_handle h = i->halfedge();
    if (CGAL::collinear(h->vertex()->point(),
                        h->next()->vertex()->point(),
                        h->prev()->vertex()->point())) {
      result.maps.emplace(h->facet(), std::make_pair(-1, h->vertex()->point()));
      continue;
    }

    Vector norm = -CGAL::unit_normal(h->vertex()->point(),
                                     h->next()->vertex()->point(),
                                     h->prev()->vertex()->point());
    double dir = Free_Stream * norm;

    Point c = CGAL::centroid(h->vertex()->point(),
                             h->next()->vertex()->point(),
                             h->prev()->vertex()->point());

    if (dir < -0.8) { //-1e-4
      result.maps.emplace(h->facet(), std::make_pair(-1, c));
      continue;
    }

    cent.push_back(c);
    fh.push_back(h->facet());
  }

  // Compute mean min triangle edge
  int nedg = 0;
  double tlen = 0.0;
  for ( Facet_iterator i = G.facets_begin(); i != G.facets_end(); ++i ) {
    Halfedge_handle h = i->halfedge();
    std::vector<double> edgs;
    edgs.push_back((h->vertex()->point() - h->opposite()->vertex()->point()).squared_length());
    edgs.push_back((h->next()->vertex()->point() - h->next()->opposite()->vertex()->point()).squared_length());
    edgs.push_back((h->prev()->vertex()->point() - h->prev()->opposite()->vertex()->point()).squared_length());

    tlen += sqrt(*std::min_element(edgs.begin(), edgs.end()));
    nedg++;
  }

  tlen /= nedg;

  double dt = kt * tlen / max_vel;

  // Compute AABB tree
  Tree tree(faces(G).first, faces(G).second, G);
  tree.accelerate_distance_queries();

  // Compute nanoflann KD-tree
  PointCloud cloud;
  cloud.pts.resize(cent.size());

  for (unsigned int i = 0; i < cent.size(); i++) {
    cloud.pts[i] = cent[i];
  }

  typedef KDTreeSingleIndexAdaptor<L2_Simple_Adaptor<double, PointCloud>,PointCloud,3> kd_tree_t;

  kd_tree_t treeSI(3 /*dim*/, cloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
  treeSI.buildIndex();

  nanoflann::SearchParams params;
  //params.sorted = false;

  // Integrate streamlines
  using namespace boost::numeric::odeint;
  state_type x(3);

  stream strm(&tree, &vel);

  // Fancy lambda comparator (Flow vector base) // -Wfloat-equal
  Vector iv = Free_Stream;
  Vector jv = CGAL::cross_product( iv, Vector(std::rand(), std::rand(), std::rand()));
  jv = jv / sqrt(jv.squared_length());
  Vector kv = CGAL::cross_product( iv, jv);

  auto comp = [&](const Point &a, const Point &b) -> bool {
                if ( (a-Point(0.0, 0.0, 0.0)) * iv != (b-Point(0.0, 0.0, 0.0)) * iv) {
                  return (a-Point(0.0, 0.0, 0.0)) * iv > (b-Point(0.0, 0.0, 0.0)) * iv;
                }

                if ( (a-Point(0.0, 0.0, 0.0)) * jv != (b-Point(0.0, 0.0, 0.0)) * jv) {
                  return (a-Point(0.0, 0.0, 0.0)) * jv > (b-Point(0.0, 0.0, 0.0)) * jv;
                }

                return (a-Point(0.0, 0.0, 0.0)) * kv > (b-Point(0.0, 0.0, 0.0)) * kv;
              };

  std::map<Point, Facet_handle, decltype(comp)> candidates(comp);       //Kernel::Less_xyz_3
  std::map<Facet_handle, Point> centroids;
  for (size_t i = 0; i < cent.size(); i++) {
    candidates.emplace(cent[i], fh[i]);
    centroids.emplace(fh[i], cent[i]);
  }

  // Stagnation point close filtering
  std::vector<Point> stgp;

  while (!candidates.empty()) {
    x[0] = candidates.begin()->first.x();
    x[1] = candidates.begin()->first.y();
    x[2] = candidates.begin()->first.z();

    std::vector<Point> stline;
    stline.push_back(Point(x[0], x[1], x[2]));

    runge_kutta4< state_type > stepper; //runge_kutta4, euler ...

    Point old(x[0], x[1], x[2]);

    double clen = 0.0;
    Vector s2;
    for ( double t=0.0; t < 10000.0 * dt; t+= dt ) {
      stepper.do_step(strm, x, t, dt);

      Point ptt;
      try {
        ptt = tree.closest_point(Point(x[0], x[1], x[2]));
      } catch (...) {
        ptt = Point(x[0], x[1], x[2]);
      }

      x[0] = ptt.x();
      x[1] = ptt.y();
      x[2] = ptt.z();
      clen += sqrt((ptt - old).squared_length());

      for (std::vector<Point>::iterator it = stgp.begin(); it !=  stgp.end(); it++) {
        if ((*it - ptt).squared_length() <= 1.e-6) {
          stline.push_back(*it);
          old = *it;
          break;
        }
      }

      if ((old - ptt).squared_length() <= 1.e-12) {
        stgp.push_back(ptt);
        break;
      }

      Vector s1 = ptt - old;
      double angle;
      if (t < dt) {
        angle = 0;
      } else {
        angle = acos( (s1 * s2) / (sqrt(s1.squared_length()) * sqrt(s2.squared_length())));
      }
      if ( angle > max_angle * M_PI / 180.0 ) {
        break;
      }
      s2 = s1;

      stline.push_back(ptt);
      old = ptt;
    }

    result.streams.push_back(stline);

    std::vector<double> stlen;
    double len = 0;
    stlen.push_back(len);
    for (unsigned int i = 1; i < stline.size(); i++) {
      len += sqrt((stline[i]-stline[i-1]).squared_length());
      stlen.push_back(len);
      // std::cout << stline[i] << std::endl;
    }

    Point_and_primitive_id pp;
    Polyhedron::Face_handle f;
    try{
      pp = tree.closest_point_and_primitive(stline[stline.size()-1]);
      f = pp.second;
    } catch (...) {
      // break;
    }

    try {
      result.maps.emplace(candidates.at(centroids.at(f)), std::make_pair(0, old));
      candidates.erase(centroids.at(f));
    } catch (...) {

    }

    std::list<Point> LVstp;
    double rad = (stline[1]-stline[0]).squared_length();
    Point qp = stline[0];
    std::vector<double> query_pt = {qp.x(), qp.y(), qp.z()};
    std::vector<std::pair<size_t,double> >   ret_matches;

    size_t nMatches = treeSI.radiusSearch(&query_pt[0], rad, ret_matches, params);

    for (size_t i=0; i<nMatches; i++) {
      LVstp.push_back(cloud.pts[ret_matches[i].first]);
    }

    for (std::list<Point>::const_iterator it=LVstp.begin(); it != LVstp.end(); it++) {
      try {
        result.maps.emplace(candidates.at(*it), std::make_pair(len, old));
        candidates.erase((*it));
      } catch (...) {

      }
    }

    for (unsigned int j = 1; j < stline.size()-1; j++) {
      LVstp.clear();
      rad = (stline[j]-stline[j+1]).squared_length()/2;
      qp = stline[j] + (stline[j+1]-stline[j])/2;
      query_pt = {qp.x(), qp.y(), qp.z()};

      ret_matches.clear();
      nMatches = treeSI.radiusSearch(&query_pt[0], rad, ret_matches, params);

      for (size_t i=0; i<nMatches; i++) {
        LVstp.push_back(cloud.pts[ret_matches[i].first]);
      }

      for (std::list<Point>::const_iterator it=LVstp.begin(); it != LVstp.end(); it++) {
        try {
          result.maps.emplace(candidates.at(*it), std::make_pair(len - (stlen[j]+stlen[j+1])/2, old));
          candidates.erase((*it));
        } catch (...) {

        }
      }
    }

  }

  return result;
}

#endif /* PATH_H */
