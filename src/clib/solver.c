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

  // backward-Euler
  int ierr = backward_euler_NR(r, Y, dt, dtrcmd, ispecies, UV, water_rates);
  //int ierr = backward_euler_step(r, Y, dt, dtrcmd, 0, Y, ispecies, UV, water_rates);
  //int ierr = rk45_integrate(r, Y, dt, dtrcmd, ispecies, UV, water_rates);
  
  return ierr;
}

int integrate_network(int water_rates, double *Y, double T0, double T1, double n, double metl, double UV, int CRX, double dtwant, int *nstp, code_units *my_units, int ispecies, int H2_shield, double crsHI, double k24, int water_only, chemistry_data_storage *my_rates)
  {

  double *Ylst = (double *) malloc(nSpecies * sizeof(double)); // last iter abundances
  int s;
  memcpy(Ylst, Y, sizeof(double) * nSpecies);
  double dttry, dtdone, dtrcmd;
  int max_substp;
  if(water_only){
    max_substp = 5000; // max. number of substeps
  }
  else{
    max_substp = 500;
  }
  int j, ierr;
  double *Y_old = (double *) malloc(nSpecies * sizeof(double));

  // init
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

  // substep loop
  for (j = 1; j < max_substp; j++) {

    // the rates are updated based on the density and temp
    // of the cell
    int rerr = get_rates(water_rates, r, T0, n, T0, metl, UV, CRX, my_units, ispecies, Y, H2_shield, crsHI, k24, water_only, my_rates);
    if (rerr != 0){return rerr;}

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
    
    /*
    if (j%200 == 0) {
      // print status
      printf("substp %6d  dtwant = %9.2e  dtdone = %9.2e dttry = %9.2e  dtrcmd = %9.2e  err = %5s\n",j,dtwant,dtdone,dttry,dtrcmd,errmsg[ierr]);
    }
    */

    //check success of step
    if (ierr == 0){
      //printf("dt = %g\n",dttry);
      // success
      dtdone += dttry;

      // compute greatest change in chemical over time step
      /*
      double tmp, drel, dydt;
      int idx = 0;
      drel = 1.e-99;
      for (int i = 0; i < nSpecies; i++) {
        //printf("cvg x, all abunds: %3d: Ylst = %15.6e  Y = %15.6e\n", i, Ylst[i], Y[i]);
//	printf("Y[%i] = %.2e\n", i, Y[i]);
        //tmp = fabs(Y[i] - Y_old[i])/Y_old[i];
        dydt = (Y[i] - Y_old[i])/dttry;
        printf("dY[%d]/dt = %g\n",i,dydt);
        if ((Y[i] / sumY > noise_tol) && (tmp > drel)) {
          idx = i;
          drel = tmp;
        }

       }
      printf("\n");
      printf("Greatest Change: Y_old[%2d] = %.5e  Y[%2d] = %.5e  DY[%d] = %.5e\n",idx,Y_old[idx],idx,Y[idx],idx,drel);
      */

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
  //      printf("time step underflow :/ \n"); 
        dttry = SMALL;
        ierr = SMLDT;
        break;
    }

    //maximum number of iterations reached
    /*
    if ( j == max_substp - 1) {
      // SJ: NO!!! DO NOT DO THIS!!!
      printf("max sub time steps reached\n");
      exit(99);
      ierr = MXSTP;
      break;
    }
    */
  } 

  //free the malloc'd arrays
  free(r);
  free(Ylst);
  free(Y_old);
  return ierr;
}
