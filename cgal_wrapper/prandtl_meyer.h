// "Copyright (c) 2016 Javier Herrera All Rights Reserved."

#ifndef PM_H
#define PM_H

#include <boost/math/tools/roots.hpp>
using boost::math::tools::newton_raphson_iterate;
#include <math.h>

double pm(double M, double gamma){
  return sqrt((gamma + 1) / (gamma - 1)) * atan(sqrt((gamma - 1) /(gamma + 1) * (pow(M, 2) - 1))) - atan(sqrt(pow(M, 2) - 1));
}

template <class T>
struct pm_functor
{
  pm_functor(T const& Mach, T const& theta, T const& gamma) : Ma(Mach), th(theta), gm(gamma){
  }
  std::pair<T, T> operator() (T const& x)
  {
    T fx = pm(x, gm) - pm(Ma, gm) - th;  // Diff(estimate - value)

    T dx = ((gm-1)*sqrt((gm+1)/(gm-1)) * Ma) / ((gm+1)*sqrt((gm-1)*(pow(Ma,2)-1)/(gm+1)) * ((gm-1)*(pow(Ma,2)-1)/(gm+1)+1)) - 1/(Ma*sqrt(pow(Ma,2)-1)); // 1st derivative
    return std::make_pair(fx, dx);             // 'return' both fx and dx.
  }

private:
  T Ma;
  T th;
  T gm;
};

template <class T>
T pm_M2(T Ma, T th, T gm)
{ // return root of x using 1st derivative and Newton_Raphson.
  T guess = 2*Ma;
  T min = 1;
  T max = 1e3;
  int digits = std::numeric_limits<T>::digits;
  boost::uintmax_t maxit = 10;

  T M2 = newton_raphson_iterate(pm_functor<T>(Ma, th, gm), guess, min, max, digits, maxit);

  return M2;
}


std::pair<double, double> pm_cond(double Mach, double theta, double p1, double gamma) {
  double M2 = pm_M2(Mach, theta, gamma);
  double p2 = p1 *  pow((1+(gamma-1)/2 * pow(Mach,2))/(1+(gamma-1)/2 * pow(M2,2)), gamma/(gamma-1));

  return std::make_pair(M2,p2);
}

#endif /* PM_H */
