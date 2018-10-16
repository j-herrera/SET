// "Copyright (c) 2015 Javier Herrera All Rights Reserved."

#ifndef WEDGE_H
#define WEDGE_H

#include <boost/math/tools/roots.hpp>
using boost::math::tools::newton_raphson_iterate;
#include <math.h>

using namespace boost::math::tools;

template <class T>
struct mbt_wedge_functor
{ // Mach-beta-theta oblique shock and first derivative
  mbt_wedge_functor(T const& Mach, T const& theta, T const& gamma) : Ma(Mach), th(theta), gm(gamma){
  }
  std::pair<T, T> operator() (T const& x)
  {
    T fx =  (2/tan(x) * (pow(Ma, 2) * pow(sin(x), 2)-1)) /
           (pow(Ma, 2) * (gm + cos(2*x))+2)-tan(th);  // Diff(estimate - value)

    T dx = (4 * pow(Ma, 2) * sin(2*x) / tan(x) * (pow(Ma, 2)* pow(sin(x), 2)-1))/
           (pow(pow(Ma, 2) * (cos(2*x)+gm)+2, 2))
           +(4 * pow(Ma, 2) * pow(cos(x), 2) - 2*pow(1/sin(x), 2)*
             (pow(Ma, 2) * pow(sin(x), 2)-1))/(pow(Ma, 2)*(cos(2*x)+gm) +2);    // 1st derivative
    return std::make_pair(fx, dx);             // 'return' both fx and dx.
  }

private:
  T Ma;
  T th;
  T gm;
};

template <class T>
T mbt_wedge(T Ma, T th, T gm)
{ // return root of x using 1st derivative and Newton_Raphson.
  T guess = th;
  T min = 0;
  T max = M_PI / 2;
  int digits = std::numeric_limits<T>::digits;
  boost::uintmax_t maxit = 10;

  T result = newton_raphson_iterate(mbt_wedge_functor<T>(Ma, th, gm), guess, min, max, digits, maxit);
  return result;
}


double obliqueShockWedgeAngle(double infMach, double theta, double gamma) {
  return mbt_wedge(infMach, theta, gamma);
}

double afterMachWedge(double M, double b, double t, double gamma){
  return 1.0 / sin(b - t) * sqrt((1.0 + (gamma - 1.0) / 2.0 * pow(M, 2) * pow(sin(b), 2)) / (gamma * pow(M, 2) * pow(sin(b), 2) - (gamma - 1.0) / 2.0));
}

double tbm(double M, double b, double g) {
  return atan(1.0/(tan(b) * (((g + 1.0) * pow(M, 2)) / (2.0 * (pow(M, 2) * pow(sin(b), 2) - 1.0)) - 1.0)));
}

template <class T>
struct mbt2_wedge_functor
{
  mbt2_wedge_functor(T const& _Ma, T const& _gamma) : M(_Ma), g(_gamma){
  }
  T operator() (T const& x)
  {

    T t = tbm(M, x, g);
    T fx = afterMachWedge(M, x, t, g)  - 1;
    // std::cout << "fx: " << fx << "  x: " << x*180/M_PI << std::endl;

    return fx;
  }

private:
  T M;
  T g;
};

std::pair<double, double> sonicWedgeAngle(double M, double g) {
  double guess = (M_PI/2 - asin(1 / M)) / 2 + asin(1 / M);
  double factor = 2;

  const boost::uintmax_t maxit = 20;
  boost::uintmax_t it = maxit;
  bool is_rising = false;
  int digits = std::numeric_limits<double>::digits;
  int get_digits = digits - 3;
  eps_tolerance<double> tol(get_digits);

  std::pair<double, double> r = bracket_and_solve_root(mbt2_wedge_functor<double>(M, g), guess, factor, is_rising, tol, it);

  double max_shock = r.first + (r.second - r.first)/2;

  return std::make_pair (tbm(M,max_shock, g), max_shock);
}

#endif /* WEDGE_H */
