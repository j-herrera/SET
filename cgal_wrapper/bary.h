// Copyright (c) 2015 Javier Herrera All Rights Reserved.

#ifndef BARY_H
#define BARY_H

#include "./types.h"

Vector to_bary(Facet_handle& fh, Point& p, std::map< Point, Vector >& velocity){

  Point a = fh->halfedge()->vertex()->point(); //u
  Point b = fh->halfedge()->next()->vertex()->point(); //v
  Point c = fh->halfedge()->prev()->vertex()->point(); //w

  Vector v0 = b - a, v1 = c - a, v2 = p - a;

  Vector cvc = CGAL::cross_product(v2, v0);
  Vector cvb = CGAL::cross_product(v1, v2);
  Vector ar = CGAL::cross_product(v1, v0);

  double v = sqrt(cvb.squared_length() / ar.squared_length());
  double w = sqrt(cvc.squared_length() / ar.squared_length());
  double u = 1.0 - v - w;

  return u * velocity[a] + v * velocity[b] + w * velocity[c];
}

std::pair<std::array<double, 2>, std::vector<Vector> > to_bary(Facet_handle& fh, Point& p, std::map< Point, std::pair<std::array<double, 2>, std::vector<Vector> > >& curvature){

  Point a = fh->halfedge()->vertex()->point(); //u
  Point b = fh->halfedge()->next()->vertex()->point(); //v
  Point c = fh->halfedge()->prev()->vertex()->point(); //w

  Vector v0 = b - a, v1 = c - a, v2 = p - a;

  Vector cvc = CGAL::cross_product(v2, v0);
  Vector cvb = CGAL::cross_product(v1, v2);
  Vector ar = CGAL::cross_product(v1, v0);

  double v = sqrt(cvb.squared_length() / ar.squared_length());
  double w = sqrt(cvc.squared_length() / ar.squared_length());
  double u = 1.0 - v - w;

  std::pair<std::array<double, 2>, std::vector<Vector> > ka = curvature[a];
  std::pair<std::array<double, 2>, std::vector<Vector> > kb = curvature[b];
  std::pair<std::array<double, 2>, std::vector<Vector> > kc = curvature[c];

  double k1 = u * ka.first[0] + v * kb.first[0] + w * kc.first[0];
  double k2 = u * ka.first[1] + v * kb.first[1] + w * kc.first[1];

  std::array<double, 2> tmpk = {k1, k2};

  std::vector<Vector> tmpv;
  tmpv.push_back(u * ka.second[0] + v * kb.second[0] + w * kc.second[0]);
  tmpv.push_back(u * ka.second[1] + v * kb.second[1] + w * kc.second[1]);
  tmpv.push_back(u * ka.second[2] + v * kb.second[2] + w * kc.second[2]);

  std::pair<std::array<double, 2>, std::vector<Vector> > tmpp;
  tmpp =  std::make_pair(tmpk, tmpv);
  return tmpp;
}

#endif /* BARY_H */
