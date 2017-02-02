#include <stdio.h>
#include <math.h>
#include "auto_f2c.h"
#include "tm.h"

int func (integer ndim, const double *u, const integer *icp, const double *par, integer ijac, double *f, double *dfdu, double *dfdp) {
  f[0]= (-350.0*h*pow(-100.0*v - 3.5, 3)*(v - 0.055)/(pow(-1 + 0.0301973834223185*exp(-100.0*v), 3)*pow(4*exp(-55.5555555555556*v - 10.0L/3.0L) + (-100.0*v - 3.5)/(-1 + 0.0301973834223185*exp(-100.0*v)), 3)) - 1.0*v - 10.0*(v - 0.055)*(pow(Nap, 3) + Nap1 + Nap2 + Nap3) - 90.0*(v + 0.09)*(Kdr1 + Kdr2 + Kdr3 + Kdr4 + pow(n, 4)) - 0.065 + I)/Cm;
  f[1]= 1000.0*(-5.0*n + (-50.0*v - 1.7)/((-1 + 0.0333732699603261*exp(-100.0*v))*(0.125*exp(-12.5*v - 11.0L/20.0L) + (-10.0*v - 0.34)/(-1 + 0.0333732699603261*exp(-100.0*v)))))*(0.125*exp(-12.5*v - 11.0L/20.0L) + (-10.0*v - 0.34)/(-1 + 0.0333732699603261*exp(-100.0*v)));
  f[2]= 1000.0*(-5.0*h + 0.35*exp(-50.0*v - 29.0L/10.0L)/(0.07*exp(-50.0*v - 29.0L/10.0L) + 1.0/(1 + 0.0608100626252179*exp(-100.0*v))))*(0.07*exp(-50.0*v - 29.0L/10.0L) + 1.0/(1 + 0.0608100626252179*exp(-100.0*v)));
  f[3]= -3750.0*Nap3;
  f[4]= -2500.0*Nap2;
  f[5]= -1250.0*Nap1;
  f[6]= -1250.0*Nap + 1250.0/(1 + 0.263597138115727*exp(-66.6666666666667*v));
  f[7]= -4000.0*Kdr4*(0.125*exp(-12.5*v - 11.0L/20.0L) + (-10.0*v - 0.34)/(-1 + 0.0333732699603261*exp(-100.0*v)));
  f[8]= -3000.0*Kdr3*(0.125*exp(-12.5*v - 11.0L/20.0L) + (-10.0*v - 0.34)/(-1 + 0.0333732699603261*exp(-100.0*v)));
  f[9]= -2000.0*Kdr2*(0.125*exp(-12.5*v - 11.0L/20.0L) + (-10.0*v - 0.34)/(-1 + 0.0333732699603261*exp(-100.0*v)));
  f[10]= -1000.0*Kdr1*(0.125*exp(-12.5*v - 11.0L/20.0L) + (-10.0*v - 0.34)/(-1 + 0.0333732699603261*exp(-100.0*v)));

  return 0;
}

int stpnt (integer ndim, double t, double *u, double *par) {
  par[0]= 0.0;
  par[1]= 0.01;
  par[10]= 0;
  u[0]= 0;
  u[1]= 0;
  u[2]= 0;
  u[3]= 0;
  u[4]= 0;
  u[5]= 0;
  u[6]= 0;
  u[7]= 0;
  u[8]= 0;
  u[9]= 0;
  u[10]= 0;
  return 0;
}

int icnd (integer ndim, const double *par, const integer *icp, integer nint, const double *u, const double *uold, const double *udot, const double *upold, integer ijac, double *fi, double *dint) {

  return 0;
}

int bcnd (integer ndim, const double *par, const integer *icp, integer nbc, const double *u0, const double *u1, integer ijac, double *fb, double *dbc) {
  fb[0]= v_left - v_right;
  fb[1]= n_left - n_right;
  fb[2]= h_left - h_right;
  fb[3]= Nap3_left - Nap3_right;
  fb[4]= Nap2_left - Nap2_right;
  fb[5]= Nap1_left - Nap1_right;
  fb[6]= Nap_left - Nap_right;
  fb[7]= Kdr4_left - Kdr4_right;
  fb[8]= Kdr3_left - Kdr3_right;
  fb[9]= Kdr2_left - Kdr2_right;
  fb[10]= Kdr1_left - Kdr1_right;
  fb[11]= (-350.0*h_left*pow(-100.0*v_left - 3.5, 3)*(v_left - 0.055)/(pow(-1 + 0.0301973834223185*exp(-100.0*v_left), 3)*pow(4*exp(-55.5555555555556*v_left - 10.0L/3.0L) + (-100.0*v_left - 3.5)/(-1 + 0.0301973834223185*exp(-100.0*v_left)), 3)) - 1.0*v_left - (10.0*v_left - 0.55)*(Nap1_left + Nap2_left + Nap3_left + pow(Nap_left, 3)) - (90.0*v_left + 8.1)*(Kdr1_left + Kdr2_left + Kdr3_left + Kdr4_left + pow(n_left, 4)) - 0.065 + I)/Cm;
  return 0;
}

int fopt (integer ndim, const double *u, const integer *icp, const double *par, integer ijac, double *fs, double *dfdu, double *dfdp) {
  return 0;
}

int pvls (integer ndim, const double *u, double *par) {
  return 0;
}
