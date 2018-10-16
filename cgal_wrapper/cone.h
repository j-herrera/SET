// Copyright (c) 2015 Javier Herrera All Rights Reserved.

#ifndef CONE_H
#define CONE_H

#include <math.h>
#include <algorithm>
#include <boost/math/tools/roots.hpp>
#include <boost/numeric/odeint.hpp>

using namespace boost::numeric::odeint;
using namespace boost::math::tools;


typedef std::vector< double > state_type;

// Taylor Maccoll equation
class TaylorMaccoll {

  double gamma;

public:
  TaylorMaccoll( double _gamma ) : gamma(_gamma) {
  }

  void operator() ( const state_type &x, state_type &dxdt, const double t)
  {
    double A = (gamma-1)/2 * (1 - pow(x[0],2) - pow(x[1],2));

    dxdt[0] = x[1];
    dxdt[1] = (x[1] * x[0] * x[1] - A * (2 * x[0] + x[1]/tan(t)) ) /
              (A - x[1] * x[1]);
  }
};

// Stop condition
struct cond {
  bool operator() (const state_type &xtt) {
    return xtt[1] >= 0.0;
  }
};

// Newton Rhapson iterator
template <class T>
struct mbt_functor
{ // Mach-beta-theta oblique shock and first derivative
  mbt_functor(T const& _c, T const& _Ma, T const& _gamma, TaylorMaccoll const& _TM) : c(_c), Ma(_Ma), gamma(_gamma), TM(_TM){
  }
  T operator() (T const& x)
  {
    double d = atan(2 / tan(x) * (pow(Ma,2) * pow(sin(x),2) - 1) / (pow(Ma,2) * (gamma + cos(2*x)) +2) );
    double Ma2 = 1/sin(x - d) * sqrt((1+(gamma-1)/2 * pow(Ma,2) * pow(sin(x),2)) / (gamma * pow(Ma,2) * pow(sin(x),2) - (gamma-1)/2 ));
    double V = 1/sqrt(2/((gamma-1)*pow(Ma2,2))+1);
    double Vr = V * cos(x-d);
    double Vt = -(V * sin(x-d));

    state_type xt(2);
    xt[0] = Vr;
    xt[1] = Vt;


    auto stepper = make_dense_output(1.0e-4, 1.0e-4, runge_kutta_dopri5<state_type>());

    const double dt = -0.001;
    auto ode_range = make_adaptive_range(boost::ref(stepper), TM, xt, x, 0.0, dt);

    cond stop;
    auto found_iter = std::find_if(ode_range.first, ode_range.second, stop);

    if(found_iter == ode_range.second)
    {
      T fx = 0.0 + dt;
      return fx;
    }

    double t0 = stepper.previous_time();
    double t1 = stepper.current_time();
    double t_m;
    state_type x_m(2);

    while (t0 - t1 > 1e-4) {
      t_m = 0.5 * (t0 + t1);
      stepper.calc_state(t_m, x_m);
      if (stop(x_m))
        t1 = t_m;
      else
        t0 = t_m;
    }
    t_m = 0.5 * (t0 + t1);
    T fx = t_m - c;

    return fx;
  }

private:
  T c;
  T Ma;
  T gamma;
  TaylorMaccoll TM;
};

template <class T>
T mbt(T Ma, T c, T gm, TaylorMaccoll TM)
{
  // T guess = M_PI/4;
  // T factor = 2;
  const boost::uintmax_t maxit = 20;
  boost::uintmax_t it = maxit;
  // bool is_rising = true;
  int digits = std::numeric_limits<T>::digits;
  int get_digits = digits - 3;
  eps_tolerance<T> tol(get_digits);

  // std::pair<T, T> r = bracket_and_solve_root(mbt_functor<T>(c, Ma, gm, TM), guess, factor, is_rising, tol, it);
  try {
    std::pair<T, T> r = bisect(mbt_functor<T>(c, Ma, gm, TM), c, 70*M_PI/180, tol, it);
    return r.first + (r.second - r.first)/2;
  } catch(...){
    //std::cerr << "Bisection for shock angle failed" << std::endl;
    return -1;
  }
}

// Actual shock functions
double obliqueShockAngle(double Mach, double c, double g) {
  TaylorMaccoll TM(g);
  return mbt(Mach, c, g, TM);
}

double afterMach(double Ma, double x, double c, double gamma){
  if(x <= asin(1.0 / Ma)) {
    return Ma;
  }

  TaylorMaccoll TM(gamma);

  double d = atan(2 / tan(x) * (pow(Ma,2) * pow(sin(x),2) - 1) / (pow(Ma,2) * (gamma + cos(2*x)) +2) );
  double Ma2 = 1/sin(x - d) * sqrt((1+(gamma-1)/2 * pow(Ma,2) * pow(sin(x),2)) / (gamma * pow(Ma,2) * pow(sin(x),2) - (gamma-1)/2 ));
  double V = 1/sqrt(2/((gamma-1)*pow(Ma2,2))+1);
  double Vr = V * cos(x-d);
  double Vt = -(V * sin(x-d));

  V = sqrt(pow(Vr,2) + pow(Vt,2));

  state_type xt(2);
  xt[0] = Vr;
  xt[1] = Vt;

  runge_kutta4< state_type > stepper;
  integrate_const( stepper, TM, xt, x, c, -0.001);

  V = sqrt(pow(xt[0],2) + pow(xt[1],2));
  return sqrt(2 * pow(V,2) / (1 - pow(V,2)) / (gamma-1));
}

double ppo(double g, double m){
  return pow(1+(g-1)/2 * pow(m,2), -g/(g-1));
}

double afterPressure(double Ma, double x,  double c, double gamma){
  double m1n = Ma * sin(x);
  double p2p1 = 1 + 2*gamma/(gamma+1) * (pow(m1n,2)-1);

  double mc = afterMach( Ma, x,  c, gamma);

  double m2n = sqrt(((gamma-1) * pow(m1n,2) + 2) / (2 * gamma * pow(m1n,2) - (gamma-1)));

  double d = atan(2 / tan(x) * (pow(Ma,2) * pow(sin(x),2) - 1) / (pow(Ma,2) * (gamma + cos(2*x)) +2) );

  double m22 = m2n/sin(x-d);

  return ppo(gamma, mc)/ppo(gamma, m22) * p2p1;
}

template <class T>
struct mbt2_functor
{
  mbt2_functor(T const& _Ma, T const& _gamma) : M(_Ma), g(_gamma){
  }
  T operator() (T const& x)
  {

    T b = obliqueShockAngle(M, x, g);
    if (b == -1) {
      return -1;
    }
    T fx = afterMach(M, b, x, g)  - 1;

    return fx;
  }

private:
  T M;
  T g;
};

std::pair<double, double> sonicAngle(double M, double g) {
  double guess = 0.1;
  double factor = 2;

  const boost::uintmax_t maxit = 20;
  boost::uintmax_t it = maxit;
  bool is_rising = false;
  int digits = std::numeric_limits<double>::digits;
  int get_digits = digits - 3;
  eps_tolerance<double> tol(get_digits);
  std::pair<double, double> r = bracket_and_solve_root(mbt2_functor<double>(M, g), guess, factor, is_rising, tol, it);
  double max_angle = r.first + (r.second - r.first)/2;
  return std::make_pair (max_angle, obliqueShockAngle(M, max_angle, g));
}

#endif /* CONE_H */
