#ifndef JAC_RHS_BUILDER_H
#define JAC_RHS_BUILDER_H
#include <math.h>
#include "grackle_macros.h"
// This header file defines the structures and functions necessary for chemical reactions so that
// we can build the right-hand side and Jacobian automatically.

// Define min and max macros
/*#define max(a,b) \
  ({ __typeof__ (a) _a = (a); \
     __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
#define min(a,b) \
  ({ __typeof__ (a) _a = (a); \
     __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })*/

// reaction_t defines a type that can be used to store chemical reactions. It assumes two
// different reactants and up to three products. Unused reactant and product indices can be
// expressed by setting the abundance coefficient to zero and setting the species to < 0.
typedef struct {
  int rate; // The index for the corresponding reaction rate.
  // Note that the abundance coefficients are typically integers. We force this so that
  // we can get some speedup in our power calculations.
  int nra; // The coefficient attached to the first reactant.
  int nrb; // The coefficient attached to the second reactant.
  int ra; // The index for the first reactant species.
  int rb; // The index for the second reactant species.

  int npa; // The abundance coefficient for the first product species.
  int npb; // The abundance coefficient for the second produt species.
  int npc; // The abundance coefficient for the third product species.
  int npd; // The abundance coefficient for the fourth product species.
  int sa; // The index for the first product species.
  int sb; // The index for the second product species.
  int sc; // The index for the third product species.
  int sd; // The index for the fourth product species.
} reaction_t;

// A simple function to build a new reaction. Basically it's just a setter, but it makes life easier.
// build_reaction {{{
void build_reaction(reaction_t* restrict reaction, int rate, int nra, int nrb, int ra, int rb,
                           int npa, int npb, int npc, int npd, 
                           int sa, int sb, int sc, int sd);
// }}}

// A function to calculate the righthand side from a specific set of reactions.
void cal_rhs(double* restrict F,double* restrict Y, reaction_t* restrict reactions,
             double* restrict rates);

// A function to calculate the jacobian from a specific set of reactions.
void cal_jacobian(double** restrict J, double* restrict Y, reaction_t* restrict reactions, 
                  double* restrict rates);

#endif
