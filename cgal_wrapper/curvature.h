#ifndef CURVATURE_H
#define CURVATURE_H

#include <CGAL/Monge_via_jet_fitting.h>

#include <fstream>
#include <cassert>

#include <CGAL/property_map.h>


using namespace std;

#include "./PolyhedralSurf.h"
#include "./PolyhedralSurf_operations.h"
#include "./PolyhedralSurf_rings.h"

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;

// #include "./vtkexport.h"

//Kernel of the Polyhedron
typedef Kernel::Point_3 DPoint;
typedef Kernel::Vector_3 DVector;

//HDS
typedef Polyhedron::Vertex_handle Vertex_handle;
typedef Polyhedron::Vertex Vertex;
typedef Polyhedron::Halfedge_handle Halfedge_handle;
typedef Polyhedron::Halfedge Halfedge;
typedef Polyhedron::Vertex_iterator Vertex_iterator;
typedef Polyhedron::Facet_handle Facet_handle;
typedef Polyhedron::Facet Facet;

struct Hedge_cmp {
  bool operator() (Halfedge_handle a,  Halfedge_handle b) const {
    return &*a < &*b;
  }
};

struct Facet_cmp {
  bool operator() (Facet_handle a, Facet_handle b) const {
    return &*a < &*b;
  }
};

//Vertex property map, with std::map
typedef std::map<Vertex*, int> Vertex2int_map_type;
typedef boost::associative_property_map< Vertex2int_map_type > Vertex_PM_type;
typedef T_Polyhedron_rings<Polyhedron, Vertex_PM_type > Poly_rings;

//Hedge property map, with enriched Halfedge with its length
// typedef HEdge_PM<Polyhedron> Hedge_PM_type;
// typedef T_Polyhedron_hedge_ops<Polyhedron, Hedge_PM_type> Poly_hedge_ops;
//Hedge property map, with std::map
typedef std::map<Halfedge_handle, double, Hedge_cmp> Hedge2double_map_type;
typedef boost::associative_property_map<Hedge2double_map_type> Hedge_PM_type;
typedef T_Polyhedron_hedge_ops<Polyhedron, Hedge_PM_type> Poly_hedge_ops;

// //Facet property map with enriched Facet with its normal
// typedef Facet_PM<Polyhedron> Facet_PM_type;
// typedef T_Polyhedron_facet_ops<Polyhedron, Facet_PM_type> Poly_facet_ops;
//Facet property map, with std::map
typedef std::map<Facet_handle, Vector, Facet_cmp> Facet2normal_map_type;
typedef boost::associative_property_map<Facet2normal_map_type> Facet_PM_type;
typedef T_Polyhedron_facet_ops<Polyhedron, Facet_PM_type> Poly_facet_ops;

typedef CGAL::Monge_via_jet_fitting<Kernel> My_Monge_via_jet_fitting;
typedef My_Monge_via_jet_fitting::Monge_form My_Monge_form;


// default parameter values and global variables
unsigned int d_fitting = 2;
unsigned int d_monge = 2;
unsigned int nb_rings = 0; //seek min # of rings to get the required #pts
unsigned int nb_points_to_use = 0; //
bool verbose = false;
unsigned int min_nb_points = (d_fitting + 1) * (d_fitting + 2) / 2;


//gather points around the vertex v using rings on the
//Polyhedron. the collection of points resorts to 3 alternatives:
// 1. the exact number of points to be used
// 2. the exact number of rings to be used
// 3. nothing is specified
void gather_fitting_points(Vertex* v,
                           std::vector<DPoint> &in_points,
                           Vertex_PM_type& vpm)
{
  //container to collect vertices of v on the Polyhedron
  std::vector<Vertex*> gathered;
  //initialize
  in_points.clear();

  //OPTION -p nb_points_to_use, with nb_points_to_use != 0. Collect
  //enough rings and discard some points of the last collected ring to
  //get the exact "nb_points_to_use"
  if ( nb_points_to_use != 0 ) {
    Poly_rings::collect_enough_rings(v, nb_points_to_use, gathered, vpm);
    if ( gathered.size() > nb_points_to_use ) gathered.resize(nb_points_to_use);
  }
  else { // nb_points_to_use=0, this is the default and the option -p is not considered;
    // then option -a nb_rings is checked. If nb_rings=0, collect
    // enough rings to get the min_nb_points required for the fitting
    // else collect the nb_rings required
    if ( nb_rings == 0 )
      Poly_rings::collect_enough_rings(v, min_nb_points, gathered, vpm);
    else Poly_rings::collect_i_rings(v, nb_rings, gathered, vpm);
  }

  //store the gathered points
  std::vector<Vertex*>::iterator
  itb = gathered.begin(), ite = gathered.end();
  CGAL_For_all(itb,ite) in_points.push_back((*itb)->point());
}


std::map< Point, std::pair<std::array<double, 2>, std::vector<Vector>> > curvatures(Polyhedron & P) {
  //create property maps
  //-----------------------------
  //Vertex, using a std::map
  Vertex2int_map_type vertex2props;
  Vertex_PM_type vpm(vertex2props);

  //Hedge, with enriched hedge
  //HEdgePM_type hepm = get_hepm(boost::edge_weight_t(), P);
  //Hedge, using a std::map
  Hedge2double_map_type hedge2props;
  Hedge_PM_type hepm(hedge2props);

  //Facet PM, with enriched Facet
  //FacetPM_type fpm = get_fpm(boost::vertex_attribute_t(), P);
  //Facet PM, with std::map
  Facet2normal_map_type facet2props;
  Facet_PM_type fpm(facet2props);

  //initialize Polyhedral data : length of edges, normal of facets
  Poly_hedge_ops::compute_edges_length(P, hepm);
  Poly_facet_ops::compute_facets_normals(P, fpm);

  //MAIN LOOP: perform calculation for each vertex
  //----------------------------------------------
  std::vector<DPoint> in_points;  //container for data points
  Vertex_iterator vitb, vite;

  //initialize the tag of all vertices to -1
  vitb = P.vertices_begin(); vite = P.vertices_end();
  CGAL_For_all(vitb,vite) put(vpm, &(*vitb), -1);

  vitb = P.vertices_begin(); vite = P.vertices_end();
  std::map< Point, std::pair<std::array<double, 2>, std::vector<Vector>> > k;
  for (; vitb != vite; vitb++) {
    //initialize
    Vertex* v = &(*vitb);
    in_points.clear();
    My_Monge_form monge_form;

    //gather points around the vertex using rings
    gather_fitting_points(v, in_points, vpm);

    //skip if the nb of points is to small
    if ( in_points.size() < min_nb_points )
    {std::cerr << "not enough pts for fitting this vertex" << in_points.size() << std::endl;
     continue; }

    // perform the fitting
    My_Monge_via_jet_fitting monge_fit;
    monge_form = monge_fit(in_points.begin(), in_points.end(),
                           d_fitting, d_monge);
    //switch min-max ppal curv/dir wrt the mesh orientation
    const DVector normal_mesh = Poly_facet_ops::compute_vertex_average_unit_normal(v, fpm);
    monge_form.comply_wrt_given_normal(normal_mesh);

    std::array<double, 2> tmpk = {monge_form.principal_curvatures(0), monge_form.principal_curvatures(1)};
    std::vector<Vector> tmpv;
    tmpv.push_back(monge_form.maximal_principal_direction());
    tmpv.push_back(monge_form.minimal_principal_direction());
    tmpv.push_back(monge_form.normal_direction());

    std::pair<std::array<double, 2>, std::vector<Vector>> tmpp;
    tmpp =  std::make_pair(tmpk, tmpv);

    k.emplace(vitb->point(), tmpp);
  } //all vertices processed
    // toVTK(P, curv, "curv.vtk", "min_curvature");
  return k;
}

std::map< Point, std::pair<std::array<double, 2>, std::vector<Vector>> > curvatures(Polyhedron & P, std::list<Point> &st_points) {
  //create property maps
  //-----------------------------
  //Vertex, using a std::map
  Vertex2int_map_type vertex2props;
  Vertex_PM_type vpm(vertex2props);

  //Hedge, with enriched hedge
  //HEdgePM_type hepm = get_hepm(boost::edge_weight_t(), P);
  //Hedge, using a std::map
  Hedge2double_map_type hedge2props;
  Hedge_PM_type hepm(hedge2props);

  //Facet PM, with enriched Facet
  //FacetPM_type fpm = get_fpm(boost::vertex_attribute_t(), P);
  //Facet PM, with std::map
  Facet2normal_map_type facet2props;
  Facet_PM_type fpm(facet2props);

  //initialize Polyhedral data : length of edges, normal of facets
  Poly_hedge_ops::compute_edges_length(P, hepm);
  Poly_facet_ops::compute_facets_normals(P, fpm);

  //MAIN LOOP: perform calculation for each vertex
  //----------------------------------------------
  std::vector<DPoint> in_points;  //container for data points
  Vertex_iterator vitb, vite;

  //initialize the tag of all vertices to -1
  vitb = P.vertices_begin(); vite = P.vertices_end();
  CGAL_For_all(vitb,vite) put(vpm, &(*vitb), -1);

  // Compute AABB tree
  Tree tree(faces(P).first, faces(P).second, P);
  tree.accelerate_distance_queries();

  std::map< Point, std::pair<std::array<double, 2>, std::vector<Vector>> > k;
  while (!st_points.empty()) {
    //initialize
    Point tmpP = st_points.back();
    st_points.pop_back();

    Point_and_primitive_id pp;
    Polyhedron::Face_handle f;
    try {
      pp = tree.closest_point_and_primitive(tmpP);
      f = pp.second;
    } catch (...) {
      continue;
    }

    Vertex* v = &(*(f->halfedge()->vertex()));
    in_points.clear();
    My_Monge_form monge_form;

    //gather points around the vertex using rings
    gather_fitting_points(v, in_points, vpm);

    //skip if the nb of points is to small
    if ( in_points.size() < min_nb_points )
    {std::cerr << "not enough pts for fitting this vertex" << in_points.size() << std::endl;
     continue; }

    // perform the fitting
    My_Monge_via_jet_fitting monge_fit;
    monge_form = monge_fit(in_points.begin(), in_points.end(),
                           d_fitting, d_monge);
    //switch min-max ppal curv/dir wrt the mesh orientation
    const DVector normal_mesh = Poly_facet_ops::compute_vertex_average_unit_normal(v, fpm);
    monge_form.comply_wrt_given_normal(normal_mesh);

    std::array<double, 2> tmpk = {monge_form.principal_curvatures(0), monge_form.principal_curvatures(1)};
    std::vector<Vector> tmpv;

    Vector dmt;
    dmt = monge_form.maximal_principal_direction();
    dmt = dmt / sqrt(dmt.squared_length());
    tmpv.push_back(dmt);

    dmt = monge_form.minimal_principal_direction();
    dmt = dmt / sqrt(dmt.squared_length());
    tmpv.push_back(dmt);

    dmt = monge_form.normal_direction();
    dmt = dmt / sqrt(dmt.squared_length());
    tmpv.push_back(dmt);

    std::pair<std::array<double, 2>, std::vector<Vector>> tmpp;
    tmpp =  std::make_pair(tmpk, tmpv);

    k.emplace(tmpP, tmpp);
  } //all vertices processed
    // toVTK(P, curv, "curv.vtk", "min_curvature");
  return k;
}

#endif /* CURVATURE_H */
