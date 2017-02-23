#include <stdio.h>
#include <math.h>
#include "auto_f2c.h"
#include "wbtemp.h"

int func (integer ndim, const double *u, const integer *icp, const double *par, integer ijac, double *f, double *dfdu, double *dfdp) {
  f[0]= (pow(1.2, (1.0L/10.0L)*dT)*pow(2.0, 0.3*dT)*gNa*h*pow(-0.1*v - 3.5, 3)*((55.0L/293.0L)*dT - v + 55)/(pow(-1 + 0.0301973834223185*exp(-0.1*v), 3)*pow(0.14269597338901*pow(2.0, 0.1*dT)*exp(-1.0L/18.0L*v) + pow(2.0, 0.1*dT)*(-0.1*v - 3.5)/(-1 + 0.0301973834223185*exp(-0.1*v)), 3)) + pow(1.2, (1.0L/10.0L)*dT)*gK*pow(n, 4)*(-90.0L/293.0L*dT - v - 90) + 0.1*pow(1.2, (1.0L/10.0L)*dT)*(-65.0L/293.0L*dT - v - 65) + I)/Cm;
  f[1]= -0.360593631487804*pow(2.0, 0.1*dT)*n*exp(-0.0125*v) + 5*pow(2.0, 0.1*dT)*(-n + 1)*(-0.01*v - 0.34)/(-1.0 + 0.0333732699603261*exp(-0.1*v));
  f[2]= -5*pow(2.0, 0.1*dT)*h/(1 + 0.0608100626252179*exp(-0.1*v)) + 0.0192581270197425*pow(2.0, 0.1*dT)*(-h + 1)*exp(-0.05*v);

  return 0;
}

int stpnt (integer ndim, double t, double *u, double *par) {
  par[0]= 0.0;
  par[1]= 0.0;
  par[2]= 350.0;
  par[3]= 90.0;
  par[4]= 0.01;
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
  fb[3]= (pow(1.2, (1.0L/10.0L)*dT)*pow(2.0, 0.3*dT)*gNa*h_left*pow(-0.1*v_left - 3.5, 3)*((55.0L/293.0L)*dT - v_left + 55)/(pow(-1 + 0.0301973834223185*exp(-0.1*v_left), 3)*pow(0.14269597338901*pow(2.0, 0.1*dT)*exp(-1.0L/18.0L*v_left) + pow(2.0, 0.1*dT)*(-0.1*v_left - 3.5)/(-1 + 0.0301973834223185*exp(-0.1*v_left)), 3)) + pow(1.2, (1.0L/10.0L)*dT)*gK*pow(n_left, 4)*(-90.0L/293.0L*dT - v_left - 90) + 0.1*pow(1.2, (1.0L/10.0L)*dT)*(-65.0L/293.0L*dT - v_left - 65) + I)/Cm;
  return 0;
}

int fopt (integer ndim, const double *u, const integer *icp, const double *par, integer ijac, double *fs, double *dfdu, double *dfdp) {
  return 0;
}

int pvls (integer ndim, const double *u, double *par) {
  return 0;
}
