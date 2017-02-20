#include <stdio.h>
#include <math.h>
#include "auto_f2c.h"
#include "nak1.h"

int func (integer ndim, const double *u, const integer *icp, const double *par, integer ijac, double *f, double *dfdu, double *dfdp) {
  f[0]= (gL*(-v - 65) + 35*h*(-v + 55)*pow(-0.1*v - 3.5, 3)/(pow(-1 + 0.0301973834223185*exp(-0.1*v), 3)*pow(4*exp(-1.0L/18.0L*v - 10.0L/3.0L) + (-0.1*v - 3.5)/(-1 + 0.0301973834223185*exp(-0.1*v)), 3)) + 9*pow(n, 4)*(-v - 90) + I)/Cm;
  f[1]= -5*h/(1 + 0.0608100626252179*exp(-0.1*v)) + 0.35*(-h + 1)*exp(-1.0L/20.0L*v - 29.0L/10.0L);
  f[2]= -0.625*n*exp(-1.0L/80.0L*v - 11.0L/20.0L) + 5*(-n + 1)*(-0.01*v - 0.34)/(-1 + 0.0333732699603261*exp(-0.1*v));

  return 0;
}

int stpnt (integer ndim, double t, double *u, double *par) {
  par[0]= 0;
  par[1]= 0.1;
  par[2]= 1;
  par[3]= 1;
  par[4]= 0.1;
  par[5]= 1;
  par[6]= 1;
  par[7]= 1;
  par[8]= 1;
  par[9]= 0.1;
  par[10]= 0;
  par[11]= 1;
  par[12]= 1;
  par[13]= 1;
  par[14]= 1;
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
  fb[1]= h_left - h_right;
  fb[2]= n_left - n_right;
  fb[3]= (gL*(-v_left - 65) + 35*h_left*(-v_left + 55)*pow(-0.1*v_left - 3.5, 3)/(pow(-1 + 0.0301973834223185*exp(-0.1*v_left), 3)*pow(4*exp(-1.0L/18.0L*v_left - 10.0L/3.0L) + (-0.1*v_left - 3.5)/(-1 + 0.0301973834223185*exp(-0.1*v_left)), 3)) + 9*pow(n_left, 4)*(-v_left - 90) + I)/Cm;
  return 0;
}

int fopt (integer ndim, const double *u, const integer *icp, const double *par, integer ijac, double *fs, double *dfdu, double *dfdp) {
  return 0;
}

int pvls (integer ndim, const double *u, double *par) {
  return 0;
}
