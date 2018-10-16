// Copyright (c) 2015 Javier Herrera All Rights Reserved.

#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include "poly.h"

BOOST_PYTHON_MODULE(wrapper)
{
  boost::python::class_<std::vector<double> >("std::vector<double>")
  .def(boost::python::vector_indexing_suite<std::vector<double> >() )
  ;

  boost::python::class_<std::map<std::string, double> >("optMap")
  .def(boost::python::map_indexing_suite<std::map<std::string, double> >() )
  ;

  boost::python::class_<shock_struct>("shock_struct")
  .def("return_vertices", &shock_struct::return_vertices)
  .def("return_faces", &shock_struct::return_faces)
  .def("return_nclusters", &shock_struct::return_nclusters)
  .def("return_structure", &shock_struct::return_structure)
  ;

  boost::python::class_<CGAL_geometry>("CGAL_geometry", boost::python::init<const boost::python::list &, const boost::python::list & >())
  .def("shock", &CGAL_geometry::shock)
  ;
}
