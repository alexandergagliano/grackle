#ifndef RK45_H
#define RK45_H
#include <stdlib.h>
#include <stdio.h>
//#include <parameters.h>
//#include <reactions.h>
#include <chemistry.h>

int check_solution(double *Y, double *Ylst, double *err, double *tol, int ispecies, int water_rates, int NR);

int rk45_integrate(double *r, double *Y, double dt, double *dtrcmd, int ispecies, double UV, int water_rates);
#endif
