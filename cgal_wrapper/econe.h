// Copyright (c) 2015 Javier Herrera All Rights Reserved.

#include "cone.h"
#include "wedge.h"

double ellipticShockAngle(double mach, double e, double angle, double gmm, double cas){
  if (e == 0){
    return obliqueShockAngle(mach, angle, gmm);
  } else if(e == 1){
    return obliqueShockWedgeAngle(mach, angle, gmm);
  } else{
    double k = sqrt((e+1)/(1-e));

    double b = 0;
    if ( cas == 1 ) {
      b = tan(angle) * k;
    } else if (cas == 0) {
      b = tan(angle);
    }

    double thm = atan(b * sqrt(1-e));
    double dl = thm + (pow(e,2)/32.0)*(3.0-2.0*pow(sin(thm),2))*sin(2*thm);

    double beta = obliqueShockAngle(mach, dl, gmm);

    double sg = beta / dl;

    double ep = (e/4.0)*(1+(pow(e,2)/32.0)*(15.0-20.0*pow(sin(thm),2)+8.0*pow(sin(thm),4)))*sin(2*thm);
    double g = 6.0 * pow(sg,3) / (3.0/(sqrt(pow(sg,2)-1.0) * cos(1.0/sg)) + 6.0/(gmm+1.0) * (pow(sg,6)+pow(sg,2)) + 3.0*pow(sg,4) - pow(sg,2) - 5.0);

    double theta = 0;
    if ( cas == 1 ) {
      theta = beta - ep * g;
    } else if (cas == 0) {
      theta = beta + ep * g;
    }

    if(theta < angle) {
      return 1.1*angle;
    }
    return theta;
  }
}

double ellipticAfterMach(double mach, double e, double angle, double gmm, double cas) {
  if (e == 1){
    double beta = obliqueShockWedgeAngle(mach, angle, gmm);
    return afterMachWedge(mach, beta, angle, gmm);
  } else {
    double beta = obliqueShockAngle(mach, angle, gmm);
    return afterMach(mach, beta, angle, gmm);
  }
}

std::pair<double, double> maxAngle(int model, double M, double g) {
  if (model == 1){
    return sonicWedgeAngle( M,  g);
  } else {
    return sonicAngle( M,  g);
  }
}
