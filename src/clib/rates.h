#include <stdlib.h>
#include <stdio.h>

#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"

double calc_H4(double alpha[], double nH, double temp);

int get_rates(int water_rate, double *rate, double temp, double n, double t_dust, double Z, double UV, int CRX, code_units *my_units, int ispecies, double Y[], int H2_shield, double crsHI, double k24, int water_only, chemistry_data_storage *my_rates);

