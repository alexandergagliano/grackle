#include <stdlib.h>
#include <stdio.h>
#include <chemistry.h>

double calculate_mass(double *Y, int ispecies, int water_rates)
{
     double weight[nSpecies];

     weight[H] = 1.008;
     weight[el] = 0.000000;
     weight[Hplus] = 1.008;
     weight[O] = 15.999;
     weight[OH] = weight[O] + weight[H];
     weight[OHplus] = weight[OH];
     weight[H2O] = 2.0*weight[H] + weight[O];
     weight[O2] = 2.0*weight[O];
     weight[Oplus] = weight[O];
     weight[H2Oplus] = weight[H2O];
     weight[H3Oplus] = 3.0*weight[H] + weight[O];
     weight[O2plus] = 2.0*weight[O];
     weight[Cplus] = 12.011;
     weight[C] = weight[Cplus];
     weight[CH] = weight[C] + weight[H];
     weight[CH2] = weight[C] + 2.0*weight[H];
     weight[CH3] = weight[C] + 3.0*weight[H];
     weight[CH4] = weight[C] + 4.0*weight[H];
     weight[CO] = weight[C] + weight[O];
     weight[COplus] = weight[CO];
     weight[CO2] = weight[CO] + weight[O];

     if (ispecies > 1)
     {
       weight[Hmin] = weight[H];
       weight[H2m] = 2.0*weight[H]; 
       if (ispecies > 2)
       {
         weight[D] = 2.0*weight[H];
         weight[Dplus] = weight[D];
         weight[HD] = 3.0*weight[H];
       }
     }
  
     if (water_rates ==3)
     {
         weight[CHplus] = weight[CH];
         weight[CH2plus] = weight[CH2];
         weight[He] = 4.003;
         weight[Heplus] = weight[He];
         weight[HeHplus] = weight[He] + weight[H];
         weight[H2plus] = 2.0*weight[H];
         weight[H3plus] = 3.0*weight[H];
         weight[HCOplus] = weight[CH] + weight[O];
         weight[CH3plus] = weight[CH3];
         if (ispecies > 1){weight[H2plus] = weight[H2m];}
	 weight[CH4plus] = weight[CH3] + weight[H];
         weight[CH5plus] = weight[CH3] + 2.0*weight[H];
       	 weight[O2Hplus] = weight[O2] + weight[H];
     }

     double sum = 0.0; int i;
     for (i = 0; i < nSpecies; i++){
           sum += Y[i]*weight[i];
     }
     return sum;
}

double calculate_metl_mass(double *Y, int ispecies, int water_rates)
{
     double weight[nSpecies];

     weight[H] = 0.00;
     weight[el] = 0.000000;
     weight[Hplus] = 0.00;
     weight[O] = 16.0;
     weight[OH] = 16.0 + 1.00;
     weight[OHplus] = weight[OH];
     weight[H2O] = 2.0*1.00 + 16.0;
     weight[O2] = 32.0;
     weight[Oplus] = 16.0;
     weight[H2Oplus] = weight[H2O];
     weight[H3Oplus] = 3.0 + 16.0;
     weight[O2plus] = 32.0;
     weight[Cplus] = 12.0;
     weight[C] = 12.0;
     weight[CH] = 13.0;
     weight[CH2] = 14.0;
     weight[CH3] = 15.0;
     weight[CH4] = 16.0;
     weight[CO] = 28.0;
     weight[COplus] = 28.0;
     weight[CO2] = 44.0;

     if (ispecies > 1)
     {
       weight[Hmin] = 0.00;
       weight[H2m] = 0.00;
       if (ispecies > 2)
       {
         weight[D] = 0.00;
         weight[Dplus] = 0.00;
         weight[HD] = 0.00;
       }
     }

     if (water_rates ==3)
     {
         weight[CHplus] = weight[CH];
         weight[CH2plus] = weight[CH2];
         weight[He] = 0.00;
         weight[Heplus] = 0.00;
         weight[H2plus] = 0.00;
         weight[H3plus] = 0.00;
         weight[HeHplus] = 0.00;
         weight[HCOplus] = weight[CH] + weight[O];
         weight[CH3plus] = weight[CH3];
         if (ispecies > 1){weight[H2plus] = 0.00;}
         weight[CH4plus] = weight[CH3] + weight[H];
         weight[CH5plus] = weight[CH3] + 2.0*weight[H];
         weight[O2Hplus] = weight[O2] + weight[H];
     }

     double sum = 0.0; int i;
     for (i = 0; i < nSpecies; i++){
           sum += ((Y[i]*weight[i]));
     }
     return sum/(6.022e23);
}

