#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fwdeul.h>
#include <const.h>
#include <rates.h>
#include <chemistry.h>
#include <stdbool.h>
#include <err.h>

extern reaction_t *my_reactions;

// An implementation of the forward euler method. Because it's so simple, it's very fast. This also
// means it's wildly inaccurate and should probably only be used for a first approximation in
// a more stable or accurate scheme.
int forward_euler_integrate(double *r,double *Y,double dt,int ispecies,double UV,int water_rates){
  double *rhs;
  rhs = (double*)malloc(nSpecies * sizeof(double));

  // Initialize our right-hand side.
  memset(rhs,0,nSpecies*sizeof(double));
  cal_rhs(rhs, Y, my_reactions, r);


  // Apply the Euler step.
  for(int i = 0; i < nSpecies; i++){
    Y[i] = Y[i] + dt*rhs[i];
  }

  free(rhs);
  return ALLOK;
}
