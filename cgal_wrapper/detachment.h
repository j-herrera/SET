// Copyright (c) 2015 Javier Herrera All Rights Reserved.

#include "wedge.h"
#include <math.h>
#include <algorithm>
#include <boost/math/tools/roots.hpp>
#include <boost/numeric/odeint.hpp>

using namespace boost::numeric::odeint;
using namespace boost::math::tools;

template <class T>
struct hyp_functor
{
  hyp_functor(T const& _e, T const& _a, T const& _b) : e(_e), a(_a), b(_b){
  }
  T operator() (T const& t)
  {

    double r = a * (pow(e, 2) - 1.0) / (1.0 + e * cos(t));
    double dr = a * e * (pow(e, 2) - 1.0) * sin(t) / pow(e * cos(t) + 1.0, 2);
    double dydx = (dr * sin(t) + r * cos(t)) / (dr * cos(t) - r * sin(t));
    double ang = b - atan(dydx);

    return ang;
  }

private:
  T e;
  T a;
  T b;
};

Point compute_point(double e, double a, Point f, double b) {
  double guess = 1;
  double factor = 1.5;

  const boost::uintmax_t maxit = 40;
  boost::uintmax_t it = maxit;
  bool is_rising = false;
  int digits = std::numeric_limits<double>::digits;
  int get_digits = digits - 3;
  eps_tolerance<double> tol(get_digits);
  std::pair<double, double> r = bracket_and_solve_root(hyp_functor<double>(e, a, b), guess, factor, is_rising, tol, it);
  double theta = r.first + (r.second - r.first)/2;

  double rt = a * (pow(e, 2) - 1.0) / (1.0 + e * cos(theta));
  Point pt(-rt * cos(theta) + f.x(),  -rt * sin(theta) + f.y(), 0);
  return pt;
}

double k_factor(double mach, double gamma) {
  double d = 0.386 * exp(4.67 / pow(mach,2));
  double rc = 1.386 * exp(1.8 / pow(mach - 1, 0.75));

  double bth = asin(1.0/mach);
  double e = 1.0/cos(bth);
  double a = rc / (pow(e, 2) - 1.0);
  double r0 = a * (pow(e, 2) - 1.0) / (1.0 + e);

  Point f(-1.0 - d + r0, 0.0, 0.0);

  std::pair<double, double> maxA = sonicWedgeAngle(mach, gamma);

  double pangle = M_PI/2 - maxA.first;

  Point p(-cos(pangle), sin(pangle), 0);
  Point s = compute_point(e, a, f, maxA.second);

  double tmu = M_PI/2.0;

  Vector v = s - p;
  v = v/sqrt(v.squared_length());

  double rtmu = atan2(v.y(), v.x()) - maxA.first;
  double kons = rtmu / tmu;

  return kons;
}
