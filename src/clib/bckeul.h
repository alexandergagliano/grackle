#ifndef BCKEUL_H
#define BCKEUL_H
#include <stdlib.h>
#include <stdio.h>

int backward_euler_step(double *r, double *Y, double dt, double *dtrcmd, int newton, double *Ylst, int ispecies, double UV, int water_rates);
int backward_euler_NR(double *r, double *Y, double dt, double *dtrcmd, int ispecies, double UV, int water_rates);
#endif
