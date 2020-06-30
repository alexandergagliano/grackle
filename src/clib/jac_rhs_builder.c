#include "jac_rhs_builder.h"
#include <math.h>
#include <chemistry.h>
#include <stdlib.h>

void build_reaction(reaction_t* restrict reaction, int rate, int nra, int nrb, int ra, int rb,
                           int npa, int npb, int npc, int npd,
                           int sa, int sb, int sc, int sd){
  // Make sure everything conforms to the correct bounds so we don't screw up the Jacobian.
  // Rate
  reaction->rate = max(rate,0);
  
  // Reactants
  reaction->nra = (ra >= 0) ? fmax(nra,0.0) : 0.0;
  reaction->nrb = (rb >= 0) ? fmax(nrb,0.0) : 0.0;
  if(reaction->nra == 0){
    reaction->ra = -1;
  }
  else{
    reaction->ra = ra;
  }
  if(reaction->nrb == 0){
    reaction->rb = -1;
  }
  else{
    reaction->rb = rb;
  }

  // Products
  reaction->npa = (sa >= 0) ? fmax(npa, 0) : 0;
  reaction->npb = (sb >= 0) ? fmax(npb, 0) : 0;
  reaction->npc = (sc >= 0) ? fmax(npc, 0) : 0;
  reaction->npd = (sd >= 0) ? fmax(npd, 0) : 0;
  if(reaction->npa == 0){
    reaction->sa = -1;
  }
  else{
    reaction->sa = sa;
  }
  if(reaction->npb == 0){
    reaction->sb = -1;
  }
  else{
    reaction->sb = sb;
  }
  if(reaction->npc == 0){
    reaction->sc = -1;
  }
  else{
    reaction->sc = sc;
  }
  if(reaction->npd == 0){
    reaction->sd = -1;
  }
  else{
    reaction->sd = sd;
  }
}

// A function to calculate the righthand side from a specific set of reactions.
// Arguments: double *F - The righthand side array of size nSpecies
//            double *Y - The species abundance array of size nSpecies
//            reaction_t *reactions - The array of reactions of size nReactions
//            double *rates - The array of rate coefficients of size nReactions
void cal_rhs(double* restrict F, double* restrict Y, reaction_t* restrict reactions, 
             double* restrict rates){
  // Loop through all the species.
  for(int i = 0; i < nSpecies; i++){
    // Make sure the righthand side is initialized to zero.
    F[i] = 0.0;
  }
  // Loop through all reactions, looking for the reactants and products.
  for(int i = 0; i < nReactions; i++){
    // Grab the indices of interest for our reaction.
    int ra = reactions[i].ra;
    int rb = reactions[i].rb;
    int sa = reactions[i].sa;
    int sb = reactions[i].sb;
    int sc = reactions[i].sc;
    int sd = reactions[i].sd;

    double r;
    // Get the reaction rate.
    if (reactions[i].npa == 0 && reactions[i].npb  == 0 
        && reactions[i].npc == 0 && reactions[i].npd == 0){
        continue;
    } else{
	r = rates[reactions[i].rate];
    }

    /*double ya = (ra >= 0) ? pow(Y[ra],reactions[i].nra) : 1;
    double yb = (rb >= 0) ? pow(Y[rb],reactions[i].nrb) : 1;*/
    double ya = 1.0;
    double yb = 1.0;
    // NOTE: hypothetically, any negative species will be forced to have
    // nra or nrb = 0, so these checks shouldn't be needed. Uncomment them
    // if there are problems with seg faults and bad memory accesses.
    //if(ra >= 0){
      for(int p = 0; p < reactions[i].nra; p++){
        ya *= Y[ra];
      }
    //}
    //if(rb >= 0){
      for(int p = 0; p < reactions[i].nrb; p++){
        yb *= Y[rb];
      }
    //}

    // Calculate the rhs terms for each product and reactant
    // Reactants
    if(ra > -1 && ra < nSpecies){
      F[ra] -= reactions[i].nra*ya*yb*r;
    }
    if(rb > -1 && rb < nSpecies){
      F[rb] -= reactions[i].nrb*ya*yb*r;
    }
    // Products
    if (sa > -1){
      F[sa] += reactions[i].npa*ya*yb*r;
    }
    if (sb > -1){
      F[sb] += reactions[i].npb*ya*yb*r;
    }
    if (sc > -1){
      F[sc] += reactions[i].npc*ya*yb*r;
    }
    if (sd > -1){
      F[sd] += reactions[i].npd*ya*yb*r;
    }
  }
}


// A function to calculate the jacobian from a specific set of reactions.
// Arguments: double **J - The Jacobian matrix of size nSpecies x nSpecies
//            double *Y - The species abundance array of size nSpecies
//            reaction_t *reactions - The array of reactions of size nReactions
//            double *rates - The array of rate coefficients of size nReactions
void cal_jacobian(double** restrict J, double* restrict Y, reaction_t* restrict reactions, 
                  double* restrict rates){
  // Initialize the Jacobian.
  for(int i = 0; i < nSpecies; i++){
    for(int j = 0; j < nSpecies; j++){
      J[i][j] = 0.0;
    }
  }

  // Loop over all the reactions and compute their Jacobian terms.
  for(int i = 0; i < nReactions; i++){
    // Start collecting some frequently used information.
    int ra = reactions[i].ra;
    int rb = reactions[i].rb;
    double nra = reactions[i].nra;
    double nrb = reactions[i].nrb;
    int sa = reactions[i].sa;
    int sb = reactions[i].sb;
    int sc = reactions[i].sc;
    int sd = reactions[i].sd;
    //double ya, yda, yb, ydb;
    double r;
    if (reactions[i].npa == 0 && reactions[i].npb  == 0 
        && reactions[i].npc == 0 && reactions[i].npd == 0){
      continue;
    } else{
        r = rates[reactions[i].rate];
    }
    
    // We perform these checks a lot, so let's save them.
    int ura = ra > -1 && ra < nSpecies;
    int urb = rb > -1 && rb < nSpecies;

    double ya = 1.0;
    double yda = nra;
    double yb = 1.0;
    double ydb = nrb;
    /*if(ura){
      ya = pow(Y[ra],nra);
      yda = nra*pow(Y[ra],nra-1.0);
    }
    else{
      ya = 1.0;
      yda = 0.0;
    }
    if(urb){
      yb = pow(Y[rb],nrb);
      ydb = nrb*pow(Y[rb],nrb-1.0);
    }
    else{
      yb = 1.0;
      ydb = 0.0;
    }*/
    for(int p = 0; p < nra; p++){
      ya *= Y[ra];
    }
    for(int p = 0; p < nra-1; p++){
      yda *= Y[ra];
    }
    for(int p = 0; p < nrb; p++){
      yb *= Y[rb];
    }
    for(int p = 0; p < nrb-1; p++){
      ydb *= Y[rb];
    }

    // Calculate Jacobian terms for the reactants.
    /*if(ura){
      J[ra][ra] -= nra*yda*yb*r;
      if(urb){
        J[ra][rb] -= nra*ya*ydb*r;
      }
    }
    if(urb){
      J[rb][rb] -= nrb*ya*ydb*r;
      if(ura){
        J[rb][ra] -= nrb*yda*yb*r;
      }
    }
    
    // Calculate Jacobian terms for the products.
    if(sa > -1){
      if(ura){
        J[sa][ra] += reactions[i].npa*yda*yb*r;
      }
      if(urb){
        J[sa][rb] += reactions[i].npa*ya*ydb*r;
      }
    }
    if(sb > -1){
      if(ura){
        J[sb][ra] += reactions[i].npb*yda*yb*r;
      }
      if(urb){
        J[sb][rb] += reactions[i].npb*ya*ydb*r;
      }
    }
    if(sc > -1){
      if(ura){
        J[sc][ra] += reactions[i].npc*yda*yb*r;
      }
      if(urb){
        J[sc][rb] += reactions[i].npc*ya*ydb*r;
      }
    }
    if(sd > -1){
      if(ura){
        J[sd][ra] += reactions[i].npd*yda*yb*r;
      }
      if(urb){
        J[sd][rb] += reactions[i].npd*ya*ydb*r;
      }
    }*/

    if(ura){
      J[ra][ra] -= nra*yda*yb*r;
      if(urb){
        J[ra][rb] -= nra*ya*ydb*r;
      }
      if(sa > -1){
        J[sa][ra] += reactions[i].npa*yda*yb*r;
      }
      if(sb > -1){
        J[sb][ra] += reactions[i].npb*yda*yb*r;
      }
      if(sc > -1){
        J[sc][ra] += reactions[i].npc*yda*yb*r;
      }
      if(sd > -1){
        J[sd][ra] += reactions[i].npd*yda*yb*r;
      }
    }

    if(urb){
      J[rb][rb] -= nrb*ya*ydb*r;
      if(ura){
        J[rb][ra] -= nrb*yda*yb*r;
      }
      if(sa > -1){
        J[sa][rb] += reactions[i].npa*ya*ydb*r;
      }
      if(sb > -1){
        J[sb][rb] += reactions[i].npb*ya*ydb*r;
      }
      if(sc > -1){
        J[sc][rb] += reactions[i].npc*ya*ydb*r;
      }
      if(sd > -1){
        J[sd][rb] += reactions[i].npd*ya*ydb*r;
      }
    }
  }
}
