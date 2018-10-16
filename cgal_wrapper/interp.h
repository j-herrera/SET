// Copyright (c) 2015 Javier Herrera All Rights Reserved.

#ifndef INTERP_H
#define INTERP_H

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include "./types.h"
#include "./bary.h"
#include "./curvature.h"

typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;

typedef std::vector< std::pair< Point, Kernel::FT  > > Point_coordinate_vector;

Vector get_interp(const Point& p, Tree * tree, std::map< Point, Vector >& velocity){
  try {
    Point_and_primitive_id pp = (*tree).closest_point_and_primitive(p);
    Point closest_point = pp.first;
    Polyhedron::Face_handle f = pp.second;
    return to_bary(f, closest_point, velocity);
  } catch (...) {
    std::cout << "No closest point" << std::endl;
    return Vector(0.0, 0.0, 0.0);
  }
}

std::pair<std::array<double, 2>, std::vector<Vector> > get_interp(const Point& p, Tree * tree, std::map< Point, std::pair<std::array<double, 2>, std::vector<Vector> > >& curvature){
  try {
    Point_and_primitive_id pp = (*tree).closest_point_and_primitive(p);
    Point closest_point = pp.first;
    Polyhedron::Face_handle f = pp.second;
    return to_bary(f, closest_point, curvature);
  } catch (...) {
    std::cout << "No closest point" << std::endl;

    std::array<double, 2> tmpk = {0, 0};
    std::vector<Vector> tmpv;
    tmpv.push_back(Vector(0.0, 0.0, 0.0));
    tmpv.push_back(Vector(0.0, 0.0, 0.0));
    tmpv.push_back(Vector(0.0, 0.0, 0.0));

    std::pair<std::array<double, 2>, std::vector<Vector> > tmpp;
    tmpp =  std::make_pair(tmpk, tmpv);
    return tmpp;
  }
}

#endif /* INTERP_H */
