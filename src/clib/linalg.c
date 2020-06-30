#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

int solve_Ax_equals_b(gsl_matrix *A, gsl_vector *b) {
  int s, ierr;
  static gsl_vector *x;
  static gsl_permutation *p;
  static int first = 1;
  if(first){
    p = gsl_permutation_alloc (b->size);
    x = gsl_vector_alloc(b->size);
    first = 0;
  }

  // do the LU decomposition
  ierr = gsl_linalg_LU_decomp(A, p, &s);
  if (ierr != 0) {
    printf("error in LU decomposition");
    goto error;
  }

  // solve the system
  ierr = gsl_linalg_LU_solve(A, p, b, x);
  if (ierr != 0) {
    printf("error in LU solve");
    goto error;
  }

  // write result back to vector b
  gsl_vector_memcpy(b, x);

  //clean up
  //gsl_permutation_free (p);
  //gsl_vector_free      (x);

  return 0;

error:
  //gsl_permutation_free (p);
  //gsl_vector_free      (x);
  return ierr;
}


