#ifndef FWDEUL_H
#define FWDEUL_H
#include <stdlib.h>
#include <stdio.h>

int forward_euler_integrate(double *r,double *Y,double dt,int ispecies,double UV,int water_rates);

#endif
