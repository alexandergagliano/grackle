#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <solver.h>
#include <rk45.h>
#include <bckeul.h>
#include <err.h>
#include <const.h>
#include <chemistry.h>
#include <calc_mass.h>
#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"
#include "rates.h"

int do_network_step(double *Y, double *r, double dt, double *dtrcmd, int ispecies, double UV, int water_rates) {

  // backward-Euler solver for Waternet
  int ierr = backward_euler_NR(r, Y, dt, dtrcmd, ispecies, UV, water_rates);
  //int ierr = backward_euler_step(r, Y, dt, dtrcmd, 0, Y, ispecies, UV, water_rates);
  //int ierr = rk45_integrate(r, Y, dt, dtrcmd, ispecies, UV, water_rates);
  
  return ierr;
}

int integrate_network(int water_rates, double *Y, double T0, double T1, double n, double metl, double UV, int UV_molec, int CRX, double dtwant, int *nstp, code_units *my_units, int ispecies, int H2_shield, double crsHI, double k24, int water_only, chemistry_data_storage *my_rates)
  {

  double *Ylst = (double *) malloc(nSpecies * sizeof(double)); // last iter abundances
  int s;
  memcpy(Ylst, Y, sizeof(double) * nSpecies);
  double dttry, dtdone, dtrcmd;
  int max_substp;
  max_substp = 5.e3; // The maximum number of iterations - if we fail after this, we're in trouble.

  int j, ierr;
  double *Y_old = (double *) malloc(nSpecies * sizeof(double));

  // initialize the array of reactions
  double *r = (double *) malloc(nReactions * sizeof(double));

  // We start by trying to jump the entire time interval, 
  // instead of assuming subcycling from the get go
  dttry =  dtwant ;
  dtdone = 0.; 
  dtrcmd = 1.e99; 

  // clip abundances on input to tiny (non-zero) number
  for (int i = 1; i < nSpecies; i++) {
    Y[i] = fmax(Y[i], TINY);
  }
  memcpy(Y_old, Y, sizeof(double) * nSpecies);

    // the rates are updated based on the density and temp
    // of the cell
    // note: now we don't change the rates within the substeps!
    int rerr = get_rates(water_rates, r, T0, n, T0, metl, UV, UV_molec, CRX, my_units, ispecies, Y, H2_shield, crsHI, k24, water_only, my_rates);
    if (rerr != 0){
       free(r);
       free(Ylst);
       free(Y_old);
       return rerr;
    }

  // substep loop
  for (j = 1; j < max_substp; j++) {

    // perform one step integration
    ierr = do_network_step(Y, r, dttry, &dtrcmd, ispecies, UV, water_rates);

    // Floor the species that make up less than 1.e-33% of the 
    // total number density
    double noise_tol = 1.e-35;
    double sumY = 0;
    for (int i = 0; i < nSpecies; i++) {
      sumY += Y[i];
    }
    
    for (int i = 0; i < nSpecies; i++) {
      // If species are below the noise, then just floor them
      if ((Y[i]/sumY) < noise_tol) {
        Y[i] = 0.; //TINY;
      }
    }
    
    //check success of step
    if (ierr == 0){
      dtdone += dttry;

      if (dtdone >= dtwant ) {
        //cycle complete
        ierr = ALLOK;
        break;
      }
      memcpy(Y_old, Y, sizeof(double) * nSpecies);
    }
    else {
      //failure
      memcpy(Y, Y_old, sizeof(double) * nSpecies);
    } 

    // set dttry to dtrcmd, unless the remaining time is 
    // even smaller
    dttry = fmin(dtrcmd, dtwant - dtdone);

    // check for time step underflow
    if (dttry < SMALL) {
        dttry = SMALL;
        ierr = SMLDT;
        break;
    }

    /* Is this needed? 
    if ( j == max_substp - 1) {
      printf("Max Waternet steps reached!\n");
      exit(99);
      ierr = MXSTP;
      break;
    }*/
  } 

  //free the malloc'd arrays
  free(r);
  free(Ylst);
  free(Y_old);
  return ierr;
}
