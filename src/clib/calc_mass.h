#ifndef CALC_MASS_H
#define CALC_MASS_H
#include <stdlib.h>
#include <stdio.h>
#include <chemistry.h>


double calculate_mass(double *Y, int ispecies, int water_rates);
double calculate_metl_mass(double *Y, int ispecies, int water_rates);
#endif
