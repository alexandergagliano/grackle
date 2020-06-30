#include <stdlib.h>
#include <stdio.h>
#include <chemistry.h>
#include "calc_charge.h"

double calculate_charge(double *Y, int ispecies, int water_rates){
  double charges = Y[Hplus] - Y[el] + Y[Oplus] + Y[OHplus] + Y[H2Oplus] +
                   Y[H3Oplus] + Y[O2plus] + Y[Cplus] + Y[COplus];
  if(ispecies > 1){
    charges -= Y[Hmin];
    if(ispecies > 2){
      charges += Y[Dplus];
    }
  }
  if(water_rates == 3){
    charges += Y[CHplus] + Y[CH2plus] + Y[Heplus] + Y[H3plus] + Y[HCOplus] +
               Y[H2plus] + Y[HeHplus] + Y[CH3plus] + Y[CH4plus] + Y[CH5plus] + 
	       Y[O2Hplus];
  }
  return charges;
}
