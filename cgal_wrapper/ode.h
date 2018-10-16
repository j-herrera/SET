// Copyright (c) 2015 Javier Herrera All Rights Reserved.

#ifndef ODE_H
#define ODE_H

#include <boost/numeric/odeint.hpp>
#include "./types.h"
#include "./interp.h"

typedef std::vector< double > state_type;

class stream {

  Tree * tree;
  std::map< Point, Vector > * velocity;

public:
  stream( Tree * _tree, std::map< Point, Vector > * _velocity ) : tree(_tree), velocity(_velocity) {
  }

  void operator() ( const state_type &x, state_type &dxdt, const double /* t */ )
  {
    // boost::timer::auto_cpu_timer t;
    Point proy = Point(x[0], x[1], x[2]);
    Vector vel = get_interp(proy, tree, (*velocity));
    dxdt[0] = vel.x();
    dxdt[1] = vel.y();
    dxdt[2] = vel.z();
  }
};

#endif /* ODE_H */
