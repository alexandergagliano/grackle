#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <err.h>
#include <bckeul.h>
#include <const.h>
#include <chemistry.h>
#include <stdbool.h>
#include "jac_rhs_builder.h"
#include "linalg.h"
#include "fwdeul.h"
#include "rk45.h"

extern reaction_t *my_reactions;

// compute a backward-euler time step. if newton is true (non-zero), this is a
// step of a newton-raphson iteration, in which case the abundances for the
// previous iteration should also be provided in the variable Ylst
int backward_euler_step(double *r, double *Y, double dt, double *dtrcmd, int newton, double *Ylst, int ispecies, double UV, int water_rates) {

  //const double rmax = 2, rmin = 0.2;
  //double err = TINY;
  //double tol = TINY;

  //remove this stuff for now
  //static int first = 1;
  //static double *F, **J;
  double *F, **J;
  //static gsl_matrix *A;
  gsl_matrix *A;
  //static gsl_vector *b;
  gsl_vector *b;
  //if(first){
    F = (double *) malloc(nSpecies * sizeof(double));
    J = (double **) malloc(nSpecies * sizeof(double *));
    for (int i = 0; i < nSpecies; i++){
      J[i] = (double *) malloc(nSpecies * sizeof(double));
    }
    b = gsl_vector_alloc(nSpecies);
    A = gsl_matrix_alloc(nSpecies,nSpecies);
  //  first = 0;
  //}

  memset(F, 0., nSpecies*sizeof(double));
  for(int i = 0; i < nSpecies; i++){
    memset(J[i],0,nSpecies*sizeof(double));
  }


  double dtinv = 1. / dt;
  int i, j, ierr;

  // We use Ylst for NR and Y from last time step for BE.
  if(newton){
    // For NR, do a quick floor on all the abundances; we shouldn't be guessing
    // negative numbers.
    for(int i = 0; i < nSpecies; i++){
      Ylst[i] = fmax(0.0,Ylst[i]);
    }
    cal_jacobian(J, Ylst, my_reactions, r);
    cal_rhs(F, Ylst, my_reactions, r);
  }
  else{
    cal_jacobian(J, Y, my_reactions, r);
    cal_rhs(F, Y, my_reactions, r);
  }

  // clip abs. values to 1e-99
  clip_jacobian(J, ispecies, water_rates);

  // get RHS (dY/dt): use Y from last time step if BE or last iteration (Ylst) if NR
  (newton) ? cal_rhs(F, Ylst, my_reactions, r) : cal_rhs(F, Y, my_reactions, r);


/*
  for (int k = 0; k < nSpecies; k++){
     printf("dydt[%i] = %.2e\n", k, F[k]);
  }
  */
  
  // set up matrix A and vector b for the Ax = b linear solve
  if (newton) {
    // NR: copy -dt*J to A and Ylst - Y_previous - dt * F to b
    for (i = 0; i < nSpecies; i++) {
      gsl_vector_set(b, i, Ylst[i] - Y[i] - dt * F[i]);
      for (j = 0; j < nSpecies; j++) {
        // if N-R multiply J by dt
        gsl_matrix_set(A, i, j, -dt *J[i][j]);
      }
    }
  }
  else {
    // BE: copy -J to A and F to b
    for (i = 0; i < nSpecies; i++) {
      gsl_vector_set(b, i, F[i]);
      for (j = 0; j < nSpecies; j++) {
        // if N-R multiply J by dt
        gsl_matrix_set(A, i, j, -J[i][j]);
      }
    }
  }

  /*printf("max/min A: %.5e, %.5e\n", gsl_matrix_max(A), gsl_matrix_min(A));*/

  // BE: A = I/dt - J; NR: A = I - dt*J
  if (newton) {
    for (i = 0; i < nSpecies; i++) {
      gsl_matrix_set(A, i, i, 1. + gsl_matrix_get(A, i, i));
    }
  }
  else {
    for (i = 0; i < nSpecies; i++) {
      gsl_matrix_set(A, i, i, dtinv + gsl_matrix_get(A, i, i));
    }
  }


  // solve system Ax = b for x (returned as b)
  ierr = solve_Ax_equals_b(A, b);
  if (ierr != 0) {
    goto error;
  }

  // update the abundances
  if (newton) {
    for (i = 0; i < nSpecies; i++) {
      Ylst[i] -= gsl_vector_get(b, i);
    }
  }
  else {
    for (i = 0; i < nSpecies; i++) {
      Y[i] += gsl_vector_get(b, i);
    }
  }
  ierr = ALLOK;
  goto error;
  //removing this stuff for now;
  error: 

    gsl_vector_free(b);
    gsl_matrix_free(A);
    free(F);
    for (int i = 0; i < nSpecies; i++){
      free(J[i]);
    }
    free(J);
    return ierr;
   // return ALLOK;

//error:
 // return ierr;
}

// perform Newton-Raphson iterations using backward Euler to converge to solution
int backward_euler_NR(double *r, double *Y, double dt, double *dtrcmd, int ispecies, double UV, int water_rates) {

  int iter, ierr, i;
  const int maxiter = 7;
  //don't make it static! we're in a loop
  //static int first_call = 1;
  double *Yold, *Ycur, *Ylst;
  //if(first_call){
    Yold = (double *) malloc(nSpecies * sizeof(double)); // abundances at prev. time step
    Ycur = (double *) malloc(nSpecies * sizeof(double)); // abundances at current iteration
    Ylst = (double *) malloc(nSpecies * sizeof(double)); // abundances at last iteration
  //  first_call = 0;
  //}
  const double rmax = 2, rmin = 0.2;
  double err = TINY, tol;
  

  // initialize arrays
  memcpy(Yold, Y   , sizeof(double)*nSpecies);
  memcpy(Ycur, Yold, sizeof(double)*nSpecies);
  forward_euler_integrate(r,Ycur,dt,ispecies,UV,water_rates);
  memcpy(Ylst, Ycur, sizeof(double)*nSpecies);

  // perform newton iterations
  for (iter = 0; iter < maxiter; iter++) {
    //printf("Starting iteration number %i\n", iter + 1);

    // perform one backward euler step
    ierr = backward_euler_step(r, Yold, dt, dtrcmd, 1, Ycur, ispecies, UV, water_rates);
    if (ierr != 0) break;

    ierr = check_solution(Ycur, Ylst, &err, &tol, ispecies, water_rates, 0);

    if (ierr == ALLOK) { break; }
    if (ierr == NAN_X) { break; }

    if (iter == maxiter) { break; }

    // copy Ycur (current abundances) to Ylst (abundances at last iteration)
    memcpy(Ylst, Ycur, sizeof(double)*nSpecies);
  }

  //ierr = check_solution(Ycur, Ylst, &err, &tol, ispecies, water_rates, 1);

  memcpy(Y, Ylst, sizeof(double)*nSpecies);

  // recommended time step to try next
  *dtrcmd = fmax(fmin(pow(( tol / err ), HALF), rmax), rmin) * dt;
  //printf("dttry = %.2e  dtrcmd = %.2e  err = %.2e  tol = %.2e  ierr %d  err = %s\n", dt, *dtrcmd, err, tol, ierr, errmsg[ierr]);
  free(Yold);
  free(Ycur);
  free(Ylst);
  return ierr;
}
