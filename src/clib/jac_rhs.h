#ifndef JAC_RHS
#define JAC_RHS
#include <stdlib.h>
#include <stdio.h>

void clip_jacobian(double **J, int ispecies, int water_rates);
void get_f_vector(double *rate, double *F, double *Y, int ispecies, int water_rates, int UV);
void get_jacobian(double rate[], double **J, double Y[], int ispecies, int water_rates, int UV);

#endif
