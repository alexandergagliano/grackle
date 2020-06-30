#ifndef LINALG
#define LINALG
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

int solve_Ax_equals_b(gsl_matrix *A, gsl_vector *b);
#endif
