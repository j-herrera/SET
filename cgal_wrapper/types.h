// Copyright (c) 2015 Javier Herrera All Rights Reserved.

#ifndef TYPES_H
#define TYPES_H

#define CGAL_EIGEN3_ENABLED

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
typedef typename CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef CGAL::Vector_3<Kernel> Vector;
typedef CGAL::Vector_2<Kernel> Vector2;
typedef Kernel::Point_3 Point;
typedef Kernel::Point_2 Point2;
typedef Kernel::Triangle_3 Triangle;
typedef Kernel::Segment_3 Segment;
typedef Kernel::Plane_3 Plane;
typedef Kernel::Ray_3 Ray;
typedef Kernel::Line_3 Line;
typedef Kernel::FT FT;

#include <CGAL/Aff_transformation_3.h>
typedef CGAL::Aff_transformation_3<Kernel> Transf;

#include <CGAL/Delaunay_triangulation_3.h>
typedef CGAL::Delaunay_triangulation_3<Kernel>           Triangulation;

typedef Triangulation::Cell_handle Cell_handle;
typedef Triangulation::Vertex_handle TVertex_handle;
typedef Triangulation::Locate_type Locate_type;
typedef Triangulation::Point TPoint;

// #include <CGAL/Polyhedron_3.h>
//
// typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
// typedef Polyhedron::Facet_iterator Facet_iterator;
// typedef Polyhedron::Facet_handle Facet_handle;
// typedef Polyhedron::Vertex_iterator Vertex_iterator;
// typedef Polyhedron::Halfedge_handle Halfedge_handle;
// typedef Polyhedron::Halfedge_iterator Halfedge_iterator;
// typedef Polyhedron::Vertex_handle Vertex_handle;
// typedef Polyhedron::Halfedge_around_vertex_circulator HV_circulator;
// typedef Polyhedron::Halfedge_around_facet_circulator HF_circulator;
//
// #include <boost/timer/timer.hpp>

#endif /* TYPES_H */
