#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <jac_rhs.h>
#include <chemistry.h>
#include <gsl/gsl_math.h>

// clip magnitude of jacobian elements to 1e-99 to avoid underflow error
void clip_jacobian(double **J, int ispecies, int water_rates){
  int i, j;
  double jabs, sign;

  for (i = 0; i < nSpecies; i++) {
    for (j = 0; j < nSpecies; j++) {
      //sign = (double) GSL_SIGN(J[i][j]);
      jabs = fabs(J[i][j]);
      if (jabs < 1.e-50){
        J[i][j] = 0.;
      }
    }
  }
}

