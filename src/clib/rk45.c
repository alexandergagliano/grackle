#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <err.h>
#include <bckeul.h>
#include <const.h>
#include <rates.h>
#include <linalg.h>
#include <chemistry.h>
#include "calc_mass.h"
#include "calc_charge.h"

extern reaction_t *my_reactions;

// rk4/5 cash-karp coefficients
double b21 = 1. / 5., 
       b31 = 3. / 40., 
       b32 = 9. / 40., 
       b41 = 3. / 10., 
       b42 = -9. / 10., 
       b43 = 6. / 5., 
       b51 = -11. / 54., 
       b52 = 5. * 0.5, 
       b53 = -70. / 27., 
       b54 = 35. / 27., 
       b61 = 1631. / 55296., 
       b62 = 175 / 512., 
       b63 = 575. / 13824., 
       b64 = 44275. / 110592., 
       b65 = 253. / 4096., 
       c1  = 37. / 378.,
       c1s = 2825. / 27648., 
       c2  = 0.,
       c2s = 0., 
       c3  = 250. / 621.,
       c3s = 18575. / 48384., 
       c4  = 125. / 594.,
       c4s = 13525. / 55296., 
       c5  = 0.,
       c5s = 277. / 14336., 
       c6  = 512. / 1771.,
       c6s = 1. / 4.;

double tiny = 1.e-35;
double maxfac = 2., minfac = 0.5;

int check_solution(double *Y, double *Ylst, double *err, double *tol, int ispecies, int water_rates, int NR) {

  double minY = Y[0];
  double sumM0 = 0;
  double sumM1 = 0;
  double sumY = 0;
  double sumYlst = 0;
  int ineg = 0;
  for (int i = 1; i < nSpecies; i++)
  {
     sumY = sumY + fabs(Y[i]);
     sumYlst = sumYlst + fabs(Ylst[i]);
     if (Y[i] < minY) {
       minY = Y[i];
       ineg = i;
     }
  }
  double pos_tol = 1.e-3;
  double neg_tol = 1.e-15 * ispecies * pos_tol; // 1e-15
  double Ynegmax = (minY < 0.) ? fabs(minY) : 0.;
  int ierr;

  ierr = ALLOK;

  // check for negative abundances
  //if (Ynegmax >= neg_tol && (!NR)) {
  if (Ynegmax/sumY >= neg_tol) {
    *err = Ynegmax/sumY;
    *tol = neg_tol;
    ierr = NEG_X;

    //if (NR) {
    //  printf("neg x: %3d Ylst = %15.6e  Y = %15.6e\n", ineg, Ylst[ineg], Y[ineg]);
    //}
    //for (int i = 0; i < nSpecies; i++) {
    //  printf("neg x, all abunds: %3d: Ylst = %15.6e  Y = %15.6e\n", i, Ylst[i], Y[i]);
    //}
    
    goto error;
  }

  // check for NaNs
  for (int i = 0; i < nSpecies; i++) {
     if (Y[i] != Y[i] || isnan(Y[i])) {
       *err = 1.e10;
       *tol = TINY;
       // abundance is NaN, reject time step and try with smaller one
       ierr = NAN_X;
       goto error;
     } 
  }

  // check for mass conservation
  sumM0 = calculate_mass(Ylst, ispecies, water_rates);
  sumM1 = calculate_mass(Y, ispecies, water_rates);
  if (fabs(sumM0 - sumM1)/fabs(sumM0) > pos_tol)
  {
  //  printf("Mass Conservation violation! pChange(x) = %.5e \n", sumFrac);
    *err = fabs(sumM0 - sumM1)/fabs(sumM0);
    *tol = pos_tol;
    ierr = MASSC;
    goto error;
  }

  // Check for charge conservation
  double sumQ0 = calculate_charge(Ylst, ispecies, water_rates);
  double sumQ1 = calculate_charge(Y, ispecies, water_rates);
  if (fabs(sumQ0 - sumQ1)/fabs(sumQ0) > pos_tol){
    *err = fabs(sumQ0 - sumQ1)/fabs(sumQ0);
    *tol = pos_tol;
    ierr = CHRGC;
    goto error;
  }

   // check for convergence of abundances
   double noise_tol, eps = tiny, my_eps = tiny;
   int i, i_grd;
   my_eps = 0.0;
   i_grd = 0;
   for (i = 0; i < nSpecies; i++) {
    noise_tol = 1.e-35;
    //noise_tol = 1.e-18; // 1.e-15;

    if (((Ylst[i]/sumYlst) > noise_tol) || ((Y[i]/sumY) > noise_tol)) {
      /*my_eps = fabs((Y[i] - Ylst[i]) / fmax(Ylst[i], SMALL));*/
      my_eps = fabs((Y[i]/sumY - Ylst[i]/sumYlst) / fmax(Ylst[i]/sumYlst, SMALL));
      //printf("my_eps = %.2e\n", my_eps);
      if (my_eps > eps) {
        eps = my_eps;
        i_grd = i;
      }
    }

    *err = eps; *tol = pos_tol;
    ierr = (*err < *tol) ? 0 : E_TOL;

    /*
    if (ierr == E_TOL){
    	printf("convergence: Y5[%d] = %.6e  Y4[%d] = %.6e\n", i_grd, Ylst[i_grd], i_grd, Y[i_grd]);
    	printf("convergence: eps = %.6e\n", eps);
    }*/

    if (ierr == E_TOL){goto error;}
  }

  error: 
     return ierr;
}


int rk45_integrate(double *r, double *Y, double dt, double* dtrcmd, int ispecies, double UV, int water_rates) {

  // RHS and initialize to zero
  double rmax = 2; double rmin = 0.2;
  //static double *dydt; // RHS (dY/dt)
  double *dydt;

  //static first = 1;
  //static double *k1, *k2, *k3, *k4, *k5, *k6, *Y0, *Y4, *Y5;
  double *k1, *k2, *k3, *k4, *k5, *k6, *Y0, *Y4, *Y5;

  // k1, k2, k3, k4, k5 and k6 or the cash-carp method - intermediate solutions
 // if(first){
    dydt = (double*) malloc(nSpecies * sizeof(double));
    k1 = (double*) malloc(nSpecies*sizeof(double));
    k2 = (double*) malloc(nSpecies*sizeof(double));
    k3 = (double*) malloc(nSpecies*sizeof(double));
    k4 = (double*) malloc(nSpecies*sizeof(double));
    k5 = (double*) malloc(nSpecies*sizeof(double));
    k6 = (double*) malloc(nSpecies*sizeof(double));
    Y0 = (double*) malloc(nSpecies*sizeof(double));
    Y4 = (double*) malloc(nSpecies*sizeof(double));
    Y5 = (double*) malloc(nSpecies*sizeof(double));
 // }
  memset(dydt, 0., nSpecies*sizeof(double*));

  // entry abundances Y0, 4th order solution Y4 and 5th order solution Y5
  memcpy(Y0, Y, sizeof(double)*nSpecies);

  // get RHS (dY/dt) at entry abundances Y0 and calculate k1 from it, and new
  // abundances from that
  //calculate_dydt(r, dydt, Y0, p);
  //get_f_vector(r, dydt, Y0, ispecies, water_rates, (int) UV);

  cal_rhs(dydt, Y0, my_reactions, r);
  for (int i = 0; i < nSpecies; i++) {
    //printf("i = %i, dydt[i] = %.2e\n", i, dydt[i]);
    k1[i] = dt * dydt[i];
    Y[i] = Y0[i] + (b21 * k1[i]);
  }

  // use k1 to get new abundances; get k2 from those
  //calculate_dydt(r, dydt, Y, p);
  //get_f_vector(r, dydt, Y, ispecies, water_rates, (int) UV);
  cal_rhs(dydt, Y, my_reactions, r);
  for (int i = 0; i < nSpecies; i++) {
    k2[i] = dt * dydt[i];
    Y[i] = Y0[i] + (b31 * k1[i] + b32 * k2[i]);
  }

  // use k2 to get k3
  //calculate_dydt(r, dydt, Y, p);
  //get_f_vector(r, dydt, Y, ispecies, water_rates, (int) UV);
  cal_rhs(dydt, Y, my_reactions, r);
  for (int i = 0; i < nSpecies; i++) {
    k3[i] = dt * dydt[i];
    Y[i] = Y0[i] + (b41*k1[i] + b42*k2[i] + b43*k3[i]);
  }

  // use k3 to get k4
  //calculate_dydt(r, dydt, Y, p);
  //get_f_vector(r, dydt, Y, ispecies, water_rates, (int) UV);
  cal_rhs(dydt, Y, my_reactions, r);
  for (int i = 0; i < nSpecies; i++) {
    k4[i] = dt * dydt[i];
    Y[i] = Y0[i] + (b51*k1[i] + b52*k2[i] + b53*k3[i] + b54*k4[i]);
  }

  // use k4 to get k5
  //calculate_dydt(r, dydt, Y, p);
  //get_f_vector(r, dydt, Y, ispecies, water_rates, (int) UV);
  cal_rhs(dydt, Y, my_reactions, r);
  for (int i = 0; i < nSpecies; i++) {
    k5[i] = dt * dydt[i];
    Y[i] = Y0[i] + (b61*k1[i] + b62*k2[i] + b63*k3[i] + b64*k4[i] + b65*k5[i]);
  }

  // get k6
  //calculate_dydt(r, dydt, Y, p);
  //get_f_vector(r, dydt, Y, ispecies, water_rates, (int) UV);
  cal_rhs(dydt, Y, my_reactions, r);
  for (int i = 0; i < nSpecies; i++) {
    k6[i] = dt * dydt[i];
  }

  double err, tol;
  int ierr;
  // calculate 4th and 5th order solutions and set result Y to 5th order solution
  for (int i = 0; i < nSpecies; i++) {
    Y5[i] = Y0[i] + c1 *k1[i] + c2 *k2[i] + c3 *k3[i] + c4 *k4[i] + c5 *k5[i] + c6 *k6[i];
    Y4[i] = Y0[i] + c1s*k1[i] + c2s*k2[i] + c3s*k3[i] + c4s*k4[i] + c5s*k5[i] + c6s*k6[i];
    Y[i] = Y5[i];

//   printf("i = %i, Previous Y[i]: %.2e; Final Y[i] = %.2e\n; Change in Y[i] is %.2e\n", i, Y0[i], Y[i], (c1 *k1[i] + c2 *k2[i] + c3 *k3[i] + c4 *k4[i] + c5 *k5[i] + c6 *k6[i]));
//   printf("Y[i]/Y[HI] = %.2e\n", Y[i]/Y[H]);
  }

  ierr = check_solution(Y4, Y5, &err, &tol, ispecies, water_rates, 0); 

  // recommend time step
  double fac;

//take the maximum value in Y5
  if (ierr == ALLOK){
    fac = abs( tol / fmax(err,SMALL) );
    fac = pow(fac,0.2);
  }
  else{
    fac = abs( tol / fmax(err,SMALL) );
    fac = 0.5*pow(fac,0.2);
  }
  *dtrcmd = dt * fmax(fmin(fac,maxfac),minfac);
  if (ierr == E_TOL){
	  printf("dttry = %.2e   dtrcmd = %.2e\n", dt, *dtrcmd);
  }
  /*
  if (ierr != ALLOK){
  printf("RK45: dttry = %.2e  msg = %5s  err = %.2e  tol = %.2e  dtrcmd = %.2e\n", dt, errmsg[ierr], err, tol, *dtrcmd);
  }
  */

  return ierr;
}
