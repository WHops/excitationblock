#include <stdio.h>
#include <math.h>
#include "auto_f2c.h"
#include "tm.h"

int func (integer ndim, const double *u, const integer *icp, const double *par, integer ijac, double *f, double *dfdu, double *dfdp) {
  f[0]= (350.0*h*pow(-100.0*v - 3.5, 3)*(-v + 0.055)/(pow(-1 + 0.0301973834223185*exp(-100.0*v), 3)*pow(4*exp(-55.5555555555556*v - 10.0L/3.0L) + (-100.0*v - 3.5)/(-1 + 0.0301973834223185*exp(-100.0*v)), 3)) - 90.0*pow(n, 4)*(-v - 0.09) + 1.0*v + 0.065 - I)/Cm;
  f[1]= -625.0*n*exp(-12.5*v - 11.0L/20.0L) + 1000.0*(-n + 1)*(-50.0*v - 1.7)/(-1 + 0.0333732699603261*exp(-100.0*v));
  f[2]= -5000.0*h/(1 + 0.0608100626252179*exp(-100.0*v)) + 1000.0*(-0.35*h + 0.35)*exp(-50.0*v - 29.0L/10.0L);

  return 0;
}

int stpnt (integer ndim, double t, double *u, double *par) {
  par[0]= 0.0;
  par[1]= 0.01;
  par[10]= 0;
  u[0]= 0;
  u[1]= 0;
  u[2]= 0;
  return 0;
}

int icnd (integer ndim, const double *par, const integer *icp, integer nint, const double *u, const double *uold, const double *udot, const double *upold, integer ijac, double *fi, double *dint) {

  return 0;
}

int bcnd (integer ndim, const double *par, const integer *icp, integer nbc, const double *u0, const double *u1, integer ijac, double *fb, double *dbc) {
  fb[0]= v_left - v_right;
  fb[1]= n_left - n_right;
  fb[2]= h_left - h_right;
  fb[3]= (350.0*h_left*pow(-100.0*v_left - 3.5, 3)*(-v_left + 0.055)/(pow(-1 + 0.0301973834223185*exp(-100.0*v_left), 3)*pow(4*exp(-55.5555555555556*v_left - 10.0L/3.0L) + (-100.0*v_left - 3.5)/(-1 + 0.0301973834223185*exp(-100.0*v_left)), 3)) - 90.0*pow(n_left, 4)*(-v_left - 0.09) + 1.0*v_left + 0.065 - I)/Cm;
  return 0;
}

int fopt (integer ndim, const double *u, const integer *icp, const double *par, integer ijac, double *fs, double *dfdu, double *dfdp) {
  return 0;
}

int pvls (integer ndim, const double *u, double *par) {
  return 0;
}
