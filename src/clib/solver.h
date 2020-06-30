#ifndef SOLVER_H
#define SOLVER_H
#include <stdlib.h>
#include <stdio.h>

#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"

int do_network_step(double *Y, double *r, double dt, double *dtrcmd, int ispecies, double UV, int water_rates) ;
int integrate_network(int water_rate, double *Y, double T0, double T1, double n, double metl, double UV, int CRX, double dtwant, int *nstp, code_units *my_units, int ispecies, int H2_shield, double crsHI, double k24, int water_only, chemistry_data_storage *my_rates);

#endif
