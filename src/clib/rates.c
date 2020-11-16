#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "phys_constants.h"
#include <gsl/gsl_math.h>
#include <chemistry.h>
#include <rates.h>
#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"


/* ---------------------------------------------------
 *
 *                   Water Network
 *               Last Modified: 05/08/19
 *                Author: Alex Gagliano
 *
 * ---------------------------------------------------
 *
 * The following lists the chemical reactions 
 * implemented in waternet, including those originally
 * defined in Omukai et al. (2005) and Bialy & 
 * Sternberg (2015 and 2019). Some of these are already 
 * defined in Grackle; where this is the case, we 
 * floor the reaction rates in Grackle and use the 
 * rates in this file. The table is constructed as 
 * follows:
 *
 * ------- 01/02/03:    a    + b   ->   c   +  d
 *
 * Where the reaction occurs between a and b
 * and produces c and d.
 *
 * 01 gives the listing of the reaction in 
 * Omukai et al.(2005),
 * 02 gives the reaction number in Grackle (2016), 
 * and 03 gives the reaction number in Bialy &
 * Sternberg (2019). 
 * If an entry is missing, it means this
 * reaction is not used or given in the paper. 
 *
 * Note: reaction 8 is H2 formation on dust grains!
 *
 * -------------------------------------------------
 *
 *  Species: H, H+, H2, H-, D, D+, HD, electrons, O, 
 *          OH, H2O, O2, O+, OH+, H2O+, H3O+, O2+, 
 *          C+, C, CH, CH2, CH3, CH4, CO, CO+, CO2
 *          
 * --------------Hydrogen Reactions------------------
 *
 * ------- H1/k02/R290:    HII   + e   -> HI    + phot
 * ------- H2/k07/R287:    HI    + e   -> HM   + phot
 * ------- H3/k08/R6:      HM    + HI  -> H2I   + e
 * ------- H4/k13/R294:    H2I   + HI  -> HI    + 2HI
 * ------- H5/k22/:       HI    + 2HI -> H2I   + HI
 * ------- H6/k21/:       2HI   + H2I -> H2I   + H2I
 * ------- H7//R10:        H2I   + H2I -> 2HI   + H2I
 * ------- H8/h2dust/R295: 2HI   + gr  -> H2I
 *
 * --------------Deuterium Reactions------------------
 *
 * ------- D1/k50/:  DI    + HII -> DII   + HI
 * ------- D2/k51/:  DII   + HI  -> DI    + HII
 * ------- D3/k54/:  DI    + H2I -> HI    + HDI
 * ------- D4/k52/:  DII   + H2I -> HII   + HDI
 * ------- D5/k55/:  HDI   + HI  -> H2I   + DI
 * ------- D6/k53/:  HDI   + HII -> H2I   + DII
 *
 * --------------Metal Reactions------------------
 *
 * ------- Z1//R57:    OII   + HI   -> HII   + OI
 * ------- Z2//R42:    HII   + OI   -> OII   + HI
 * ------- Z3//R170:   OII   + H2I  -> OHII  + HI
 * ------- Z4//R172:   OHII  + H2I  -> H2OII + HI
 * ------- Z5//R167:   H2OII + H2I  -> H3OII + HI
 * ------- Z6//RR88:   H2OII + e    -> OI    + H2I
 * ------- Z7//R90:    H2OII + e    -> OHI   + HI
 * ------- Z8//R96:    H3OII + e    -> OHI   + 2HI
 * ------- Z9//R93:    H3OII + e    -> H2OI  + HI
 * ------- Z10//R284:  OI    + HI   -> OHI   + phot
 * ------- Z11//R259:  OI    + H2I  -> OHI   + HI
 * ------- Z12//R260:  H2I   + OHI  -> H2OI  + HI
 * ------- Z13//R272:  OHI   + OHI  -> H2OI  + OI
 * ------- Z14//R43:   HII   + OHI  -> OHII  + HI 
 * ------- Z15//R40:   HII   + H2OI -> H2OII + HI
 * ------- Z16//R268:  HI    + OHI  -> H2I   + OI
 * ------- Z17//R266:  HI    + H2OI -> OHI   + H2I
 * ------- Z18//R286:  OI    + OI   -> O2I   + phot
 * ------- Z19//R271:  OI    + OHI  -> O2I   + HI
 * ------- Z20//R267:  HI    + O2I  -> OHI   + OI
 * ------- Z21//R41:   HII   + O2I  -> O2II  + HI
 * ------- Z22//R99:   O2II  + e    -> OI    + OI
 * ------- Z23//R252:  OI    + CHI  -> COI   + HI
 * ------- Z24//R234:  CI    + OHI  -> COI   + HI
 * ------- Z25//R233:  CI    + O2I  -> COI   + OI
 * ------- Z26//R104:  CII   + O2I  -> OII   + CO
 * ------- Z27//:      OHI   + COI  -> CO2I  + HI
 * ------- Z28//R288:  CII   + e    -> CI    + phot
 * ------- Z29//R105:  CII   + OHI  -> COII  + HI
 * ------- Z30//R54:   COII  + HI   -> HII   + COI
 * ------- Z31//R283:  CI    + HI   -> CHI   + phot
 * ------- Z32//R254:  CI    + H2I  -> CHI   + HI
 * ------- Z33//R264:  HI    + CHI  -> CI    + H2I
 * ------- Z34//R257:  H2I   + CHI  -> CH2I + HI
 * ------- Z35//R261:  HI    + CH2I -> CHI  + H2I
 * ------- Z36//R255:  H2I   + CH2I -> CH3I + HI
 * ------- Z37//R262:  HI    + CH3I -> CH2I + H2I
 * ------- Z38//R256:  H2I   + CH3I -> CH4I + HI
 * ------- Z39//R263:  HI    + CH4I -> H2I  + CH3I 
 * ------- Z40//R279:  H2I   + CI   -> CH2I + phot
 *
 * If water_rates = 3, we additionally consider
 * 9 species and 66 reactions that were not in
 * Omukai (2005). These are given as 
 *
 * ------- 03/04:    a    + b   ->   c   +  d
 *
 * Where 03 is the reaction in Bialy & Sternberg
 * (2019) and 04 is the reaction number in waternet,
 * based on an extension of Omukai's original 
 * numbering scheme.
 *
 * -------------------------------------------------
 * Species: CH+, CH2+, He, He+, H3+, HCO+, CH3+, H2+, 
 *          HeH+ 
 *
 * --------------Hydrogen Reactions------------------
 *
 * ------- R156/H9:   H2II  + H2I -> H3II  + HI
 * ------- R92/H10:   H3II  + e   -> HI    + 2HI
 * ------- R91/H11:   H3II  + e   -> H2I   + HI
 * ------- R56/H12:   HeII  + HI  -> HeI   + HII 
 * ------- R51/H13:   HeII  + H2I -> HeI   + H2II
 * ------- R168/H14:  HeII  + H2I -> He    + HII + HI
 * ------- R158/H15:  H2II  + HeI -> HeHII + HI
 * ------- R169/H16:  HeHII + H2I -> H3II  + HeI
 * ------- R194/H17:  HeHII + HI  -> H2II  + HeI
 * ------- R277/H18:  HeI   + HII -> HeHII + p
 *
 * --------------Metal Reactions------------------
 * 
 * ------- R239/Z43:  CH2I  + OI  -> COI   + H2I
 * ------- R240/Z44:  CH2I  + OI  -> COI   + 2HI
 * ------- R203/Z45:  HeII  + COI -> HeI   + CII + OI
 * ------- R103/Z46:  O2I   + CII -> COII  + OI
 * ------- R187/Z48:  H3II  + OI  -> OHII  + H2I
 * ------- R95/Z49:   H3OII + e   -> OHI   + H2I
 * ------- R94/Z50:   H3OII + e   -> OI    + H2I
 * ------- R89/Z51:   H2OII + e   -> OI    + 2HI
 * ------- R101/Z52:  OHII  + e   -> OI    + HI
 * ------- R183/Z53:  H3II  + COI -> HCOII + H2I
 * ------- R96/Z54:   H3OII + e   -> OHI   + 2HI
 * ------- R234/Z55:  OHI   + CI  -> COI   + HI
 * ------- R105/Z56:  OHI   + CII -> COII  + HI
 * ------- R43/Z57:   OHI   + HII -> OHII  + HI
 * ------- R102/Z58:  H2OI  + CII -> HCOII + HI
 * ------- R176/Z59:  HCOII + H2OI-> H3OII + COI
 * ------- R184/Z60:  H2OI  + H3II-> H3OII + H2I
 * ------- R178/Z61:  H3II  + CI  -> CHII  + H2I 
 * ------- R163/Z62:  CHII  + H2I -> CH2II + HI
 * ------- R278/Z63:  CII   + H2I -> CH2II + ph
 * ------- R164/Z64:  CH2II + H2I -> CH3II + HI
 * ------- R76/Z65:   CH3II + e   -> CH2I  + HI
 * ------- R77/Z66:   CH3II + e   -> CHI   + H2I
 * ------- R78/Z67:   CH3II + e   -> CHI   + 2HI
 * ------- R72/Z68:   CHII  + e   -> CI    + HI
 * ------- R75/Z69:   CH2II + e   -> CHI   + HI
 * ------- R74/Z70:   CH2II + e   -> CI    + 2HI
 * ------- R73/Z71:   CH2II + e   -> CI    + H2I
 * ------- R189/Z72:  CHII  + e   -> CII   + H2I
 * ------- R264/Z73:  CHI   + HI  -> CI    + H2I
 * ------- R166/Z74:  COII  + H2I -> HCOII + HI
 * ------- R97/Z75:   HCOII + e   -> COI   + HI
 * ------- R1/Z76:    OI    + CHI -> HCOII + e
 *  TODO: Add Z77 - Z281 Reactions!! 
 *
 * ---------------UV Reactions---------------------
 *
 * ------- PH16/UV1:   OHI   + p  -> OI    + HI
 * ------- PH11/UV2:   H2OI  + p  -> OI    + H2I
 * ------- PH10/UV3:   H2OI  + p  -> OHI   + HI
 * ------- PH14/UV4:   O2I   + p  -> OI    + OI
 * ------- PH6/UV5:    CHI   + p  -> CI    + HI 
 * ------- PH8/UV6:    COI   + p  -> CI    + OI
 * ------- PH20/UV7:   CI    + p  -> CII   + e
 * -------     /UV8:   OI    + p  -> OII   + e 
 * ------- PH23/UV9:   HM    + p  -> HI    + e
 * ------- PH19/UV10:  CHI   + p  -> CHII  + e 
 *  TODO: Add UV11 - UV40 Reactions!!
 *
 * ---------------CRX-Reactions--------------------
 *
 * ------- CRP1/CRX1:    CI  + crp   -> CII   + e
 * ------- CRP2/CRX2:    COI + crp   -> CI    + OI
 * ------- CRP3/CRX3:    CHII+ crp   -> CI    + HII
 * ------- CRP4/CRX4:    CHI + crp   -> CHI   + HI
 * ------- CRP5/CRX5:    H2OI+ crp   -> OHI   + HI
 * ------- CRP6/CRX6:    OHI + crp   -> OI    + HI
 * ------- CRP7/CRX7:    CH4I+ crp   -> CH2I  + H2I
 * ------- CRP8/CRX8:    O2I + crp   -> O2II  + e
 * ------- CRP9/CRX9:    O2I + crp   -> OI    + OI
 * ------- CRP10/CRX13:  CH2I+ crp   -> CH2II + e
 * ------- CRP11/CRX14:  CH2I+ crp   -> CHI   + HI
 * ------- CRP13/CRX15:  CH3I+ crp   -> CH2I  + HI
 * ------- CRP14/CRX16:  CH3I+ crp   -> CHI   + H2I
 * ------- CRX4/CRX17:   HeI + crx   -> HeII  + e 
 * ------- CRX1/CRX18:   HI  + crx   -> HII   + e
 * ------- CRX2/CRX19:   H2I + crx   -> H2II  + e
 * ------- CRX3/CRX20:   H2I + crx   -> HII   + e + HI
 * ------- CRX12/CRX21:  CH3I + crp  -> CH3II + e
 *
 * ---------------------------------------------------
 */

double calc_H4(double alpha[], double nH, double temp){
      double tinv = 1.0/temp;

      double lg_h1 = alpha[0] + alpha[1]*log10(temp) + alpha[2]*pow(log10(temp), 2.0) + alpha[3]*pow(log10(temp), 3.0) + alpha[4]*log10(1 + alpha[5] * tinv);
      double lg_h2 = alpha[6] * tinv;
      double lg_l1 = alpha[7] + alpha[8]*log10(temp) + alpha[9]*pow(log10(temp), 2.0) + alpha[10]*log10(1 + alpha[11] * tinv);
      double lg_l2 = alpha[12] * tinv;
      double ln_c1 = alpha[13] + alpha[14]*log10(temp) + alpha[15]*pow(log10(temp), 2.0) + alpha[16] * tinv;
      double n_c1 = pow(10.0, ln_c1);
      double p = alpha[18] + alpha[19]*exp(-temp/1850.0) + alpha[20] * exp(-temp/440.0);
      double ln_c2 = alpha[17] + ln_c1;
      double n_c2 = pow(10.0, ln_c2);

      double lg = lg_h1 - (lg_h1 - lg_l1)/(1 + pow((nH/n_c1), p)) + lg_h2 - (lg_h2 - lg_l2)/(1 + pow((nH/n_c2), p));

      double gamma = pow(10.0, lg);
      return gamma;
}

int get_rates(int water_rates, double *rate, double temp, double n, double t_dust, double Z, double UV, int UVbackground_molec_redshift_on, int CRX, code_units *my_units, int ispecies, double Y[], int H2_shield, double crsHI, double k24, int water_only, chemistry_data_storage *my_rates) 
{
  double uxyz = my_units->length_units; 
  double utim = my_units->time_units;
  double urho = my_units->density_units;
  double uaye = my_units->a_units;
  double aye = my_units->a_value;
  double tbase1 = utim;
  double xbase1 = uxyz/(aye*uaye);  // uxyz is [x]*a     
  double dbase1 = urho * pow((aye*uaye), 3); // urho is [dens]/a^3
  double kunit   = (pow(uaye, 3) * mh) / (dbase1 * tbase1);

//================================================================================================//

  // DEBUGGING TEMPERATURE FOR NOW
  //if (water_only){temp = 100.0;}

  double Tev = temp * 8.621738e-5;
  //double XRTev = XRT * 8.621738e-5;
  double logTev = log(Tev);
  double tinv = 1.0 / temp;
  double G = 6.67e-8; // cm^3/s^2/g
  double Z_solar = 0.01295;
  double dom = urho*pow(aye, 3)/mh;
  /* Note: The forula for Jeans length is sqrt((15*kboltz*T)/(4*pi*G*mu*rho)) */
  /*       obviously mu is the mean mass per particle and rho is the density*/
  // n gets passed in not as a number density, but as a mass density in code units
  // We need to convert to physical units and then to number density  

  n = n*urho/mh;
  double tov300 = temp / 300.;
  double T2 = temp / 100.;

  double mu = 2*mh;
  double rho_H2 = n*mu; //note that we're assuming nearly everything here is molecular hydrogen (n is actually total number density)
  double jeans = sqrt(15.0*kboltz* temp/(4.0*pi*G*mu*rho_H2)); //jeans length
  double Nh2 = n*jeans; // estimating H2 column density
  //double xx = Nh2/(5.0e14); // dimensionless
  //double b5 = sqrt( 2.0 * kboltz*temp/mh)/100000.0; //doppler parameter, from cm/s to km/s
  // Self-shielding of molecular Hydrogen
  //double shield = 0.965/pow(1.0 + xx/b5, 1.1) + 0.035/pow(1.0+xx, 0.5) * exp(-pow(1180.0,-1)*pow(1.0 + xx, 0.5)); // eqn (26) in the textbook given below

   //these coefficients are taken from Bialy & Sternberg 2015
   
   // if shielding is high, assume we use the LW-blocked rates instead of the thin rates
   // If column density of H2 is greater than 1.e22 cm^-2, then we'll have the LW-blocked rates

   rate[H1] = 0.0;
   rate[Z1] = 0.0;
   rate[Z2] = 0.0;
   rate[Z7] = 0.0;
   rate[Z8] = 0.0;
   rate[Z9] = 0.0;
   rate[Z10] = 0.0;
   rate[Z13] = 0.0;
   rate[Z14] = 0.0;
   rate[Z15] = 0.0;
   rate[Z18] = 0.0;
   rate[Z19] = 0.0;
   rate[Z20] = 0.0;
   rate[Z21] = 0.0;
   rate[Z22] = 0.0;
   rate[Z23] = 0.0;
   rate[Z24] = 0.0;
   rate[Z25] = 0.0;
   rate[Z26] = 0.0;
   if (water_rates != 3){
   	rate[Z27] = 0.0;
   }
   rate[Z28] = 0.0;
   rate[Z29] = 0.0;
   rate[Z30] = 0.0;
   rate[Z31] = 0.0;

   if (ispecies > 1){
      rate[H2] = 0.0;
      rate[H3] = 0.0;
      rate[H4] = 0.0;
      rate[H5] = 0.0;
      rate[H6] = 0.0;
      rate[H7] = 0.0;
      rate[H8] = 0.0;
      rate[Z3] = 0.0;
      rate[Z4] = 0.0;
      rate[Z5] = 0.0;
      rate[Z6] = 0.0;
      rate[Z11] = 0.0;
      rate[Z12] = 0.0;
      rate[Z16] = 0.0;
      rate[Z17] = 0.0;
      rate[Z32] = 0.0;
      rate[Z33] = 0.0;
      rate[Z34] = 0.0;
      rate[Z35] = 0.0;
      rate[Z36] = 0.0;
      rate[Z37] = 0.0;
      rate[Z38] = 0.0;
      rate[Z39] = 0.0;
      rate[Z40] = 0.0;

      if (ispecies > 2){
         rate[D1]  = 0.0;
         rate[D2]  = 0.0;
         rate[D3]  = 0.0;
         rate[D4]  = 0.0;
         rate[D5]  = 0.0;
         rate[D6]  = 0.0;
      }
   }

   if ((int) UV) {
      rate[UV1] = 0.0;
      rate[UV2] = 0.0;
      rate[UV3] = 0.0;
      rate[UV4] = 0.0;
      rate[UV5] = 0.0;
      rate[UV6] = 0.0;
      rate[UV7] = 0.0;
      rate[UV8] = 0.0;
      rate[UV9] = 0.0;
      rate[UV10] = 0.0;
      rate[UV11] = 0.0;
      rate[UV12] = 0.0;
      rate[UV13] = 0.0;
      rate[UV14] = 0.0;
      rate[UV15] = 0.0;
      rate[UV16] = 0.0;
      rate[UV17] = 0.0;
      rate[UV18] = 0.0;
      rate[UV19] = 0.0;
      rate[UV20] = 0.0;
      rate[UV21] = 0.0;
      rate[UV22] = 0.0;
      rate[UV23] = 0.0;
      rate[UV24] = 0.0;
      rate[UV25] = 0.0;
      rate[UV26] = 0.0;
      rate[UV27] = 0.0;
      rate[UV28] = 0.0;
      rate[UV29] = 0.0;
      rate[UV30] = 0.0;
      rate[UV31] = 0.0;
      rate[UV32] = 0.0;
      rate[UV33] = 0.0;
      rate[UV34] = 0.0;
      rate[UV35] = 0.0;
      rate[UV36] = 0.0;
      rate[UV37] = 0.0;
      rate[UV38] = 0.0;
      rate[UV39] = 0.0;
      rate[UV40] = 0.0;
   }

   if (water_rates == 3){
      rate[H9] = 0.0;
      rate[H10] = 0.0;
      rate[H11] = 0.0;
      rate[H12] = 0.0;
      rate[H13] = 0.0;
      rate[H14] = 0.0;
      rate[H15] = 0.0;
      rate[H16] = 0.0;
      rate[H17] = 0.0;
      rate[H18] = 0.0;
      rate[H19] = 0.0;
      rate[H20] = 0.0;
      rate[H21] = 0.0;
      rate[H22] = 0.0;
      rate[H23] = 0.0;
      rate[H24] = 0.0;
      rate[H25] = 0.0;

      rate[Z43] = 0.0;
      rate[Z44] = 0.0;
      rate[Z45] = 0.0;
      rate[Z46] = 0.0;
      rate[Z48] = 0.0;
      rate[Z49] = 0.0;
      rate[Z50] = 0.0;
      rate[Z51] = 0.0;
      rate[Z52] = 0.0;
      rate[Z53] = 0.0;
      rate[Z54] = 0.0;
      rate[Z55] = 0.0;
      rate[Z56] = 0.0;
      rate[Z57] = 0.0;
      rate[Z58] = 0.0;
      rate[Z59] = 0.0;
      rate[Z60] = 0.0;
      rate[Z61] = 0.0;
      rate[Z62] = 0.0;
      rate[Z63] = 0.0;
      rate[Z64] = 0.0;
      rate[Z65] = 0.0;
      rate[Z66] = 0.0;
      rate[Z67] = 0.0;
      rate[Z68] = 0.0;
      rate[Z69] = 0.0;
      rate[Z70] = 0.0;
      rate[Z71] = 0.0;
      rate[Z72] = 0.0;
      rate[Z73] = 0.0;
      rate[Z74] = 0.0;
      rate[Z75] = 0.0;
      rate[Z76] = 0.0;
      rate[Z77] = 0.0;
      rate[Z78] = 0.0;
      rate[Z79] = 0.0;
      rate[Z80] = 0.0;
      rate[Z81] = 0.0;
      rate[Z82] = 0.0;
      rate[Z83] = 0.0;
      rate[Z84] = 0.0;
      rate[Z85] = 0.0;
      rate[Z86] = 0.0;
      rate[Z87] = 0.0;
      rate[Z88] = 0.0;
      rate[Z89] = 0.0;
      rate[Z90] = 0.0;
      rate[Z91] = 0.0;
      rate[Z92] = 0.0;
      rate[Z93] = 0.0;
      rate[Z94] = 0.0;
      rate[Z95] = 0.0;
      rate[Z96] = 0.0;
      rate[Z97] = 0.0;
      rate[Z98] = 0.0;
      rate[Z99] = 0.0;
      rate[Z100] = 0.0;
      rate[Z101] = 0.0;
      rate[Z102] = 0.0;
      rate[Z103] = 0.0;
      rate[Z104] = 0.0;
      rate[Z105] = 0.0;
      rate[Z106] = 0.0;
      rate[Z107] = 0.0;
      rate[Z108] = 0.0;
      rate[Z109] = 0.0;
      rate[Z110] = 0.0;
      rate[Z111] = 0.0;
      rate[Z112] = 0.0;
      rate[Z113] = 0.0;
      rate[Z114] = 0.0;
      rate[Z115] = 0.0;
      rate[Z116] = 0.0;
      rate[Z117] = 0.0;
      rate[Z118] = 0.0;
      rate[Z119] = 0.0;
      rate[Z120] = 0.0;
      rate[Z121] = 0.0;
      rate[Z122] = 0.0;
      rate[Z123] = 0.0;
      rate[Z124] = 0.0;
      rate[Z125] = 0.0;
      rate[Z126] = 0.0;
      rate[Z127] = 0.0;
      rate[Z128] = 0.0;
      rate[Z129] = 0.0;
      rate[Z130] = 0.0;
      rate[Z131] = 0.0;
      rate[Z132] = 0.0;
      rate[Z133] = 0.0;
      rate[Z134] = 0.0;
      rate[Z135] = 0.0;
      rate[Z136] = 0.0;
      rate[Z137] = 0.0;
      rate[Z138] = 0.0;
      rate[Z139] = 0.0;
      rate[Z140] = 0.0;
      rate[Z141] = 0.0;
      rate[Z142] = 0.0;
      rate[Z143] = 0.0;
      rate[Z144] = 0.0;
      rate[Z145] = 0.0;
      rate[Z146] = 0.0;
      rate[Z147] = 0.0;
      rate[Z148] = 0.0;
      rate[Z149] = 0.0;
      rate[Z150] = 0.0;
      rate[Z151] = 0.0;
      rate[Z152] = 0.0;
      rate[Z153] = 0.0;
      rate[Z154] = 0.0;
      rate[Z155] = 0.0;
      rate[Z156] = 0.0;
      rate[Z157] = 0.0;
      rate[Z158] = 0.0;
      rate[Z159] = 0.0;
      rate[Z160] = 0.0;
      rate[Z161] = 0.0;
      rate[Z162] = 0.0;
      rate[Z163] = 0.0;
      rate[Z164] = 0.0;
      rate[Z165] = 0.0;
      rate[Z166] = 0.0;
      rate[Z167] = 0.0;
      rate[Z168] = 0.0;
      rate[Z169] = 0.0;
      rate[Z170] = 0.0;
      rate[Z171] = 0.0;
      rate[Z172] = 0.0;
      rate[Z173] = 0.0;
      rate[Z174] = 0.0;
      rate[Z175] = 0.0;
      rate[Z176] = 0.0;
      rate[Z177] = 0.0;
      rate[Z178] = 0.0;
      rate[Z179] = 0.0;
      rate[Z180] = 0.0;
      rate[Z181] = 0.0;
      rate[Z182] = 0.0;
      rate[Z183] = 0.0;
      rate[Z184] = 0.0;
      rate[Z185] = 0.0;
      rate[Z186] = 0.0;
      rate[Z187] = 0.0;
      rate[Z188] = 0.0;
      rate[Z189] = 0.0;
      rate[Z190] = 0.0;
      rate[Z191] = 0.0;
      rate[Z192] = 0.0;
      rate[Z193] = 0.0;
      rate[Z194] = 0.0;
      rate[Z195] = 0.0;
      rate[Z196] = 0.0;
      rate[Z197] = 0.0;
      rate[Z198] = 0.0;
      rate[Z199] = 0.0;
      rate[Z200] = 0.0;
      rate[Z201] = 0.0;
      rate[Z202] = 0.0;
      rate[Z203] = 0.0;
      rate[Z204] = 0.0;
      rate[Z205] = 0.0;
      rate[Z206] = 0.0;
      rate[Z207] = 0.0;
      rate[Z208] = 0.0;
      rate[Z209] = 0.0;
      rate[Z210] = 0.0;
      rate[Z211] = 0.0;
      rate[Z212] = 0.0;
      rate[Z213] = 0.0;
      rate[Z214] = 0.0;
      rate[Z215] = 0.0;
      rate[Z216] = 0.0;
      rate[Z217] = 0.0;
      rate[Z218] = 0.0;
      rate[Z219] = 0.0;
      rate[Z220] = 0.0;
      rate[Z221] = 0.0;
      rate[Z222] = 0.0;
      rate[Z223] = 0.0;
      rate[Z224] = 0.0;
      rate[Z225] = 0.0;
      rate[Z226] = 0.0;
      rate[Z227] = 0.0;
      rate[Z228] = 0.0;
      rate[Z229] = 0.0;
      rate[Z230] = 0.0;
      rate[Z231] = 0.0;
      rate[Z232] = 0.0;
      rate[Z233] = 0.0;
      rate[Z234] = 0.0;
      rate[Z235] = 0.0;
      rate[Z236] = 0.0;
      rate[Z237] = 0.0;
      rate[Z238] = 0.0;
      rate[Z239] = 0.0;
      rate[Z240] = 0.0;
      rate[Z241] = 0.0;
      rate[Z242] = 0.0;
      rate[Z243] = 0.0;
      rate[Z244] = 0.0;
      rate[Z245] = 0.0;
      rate[Z246] = 0.0;
      rate[Z247] = 0.0;
      rate[Z248] = 0.0;
      rate[Z249] = 0.0;
      rate[Z250] = 0.0;
      rate[Z251] = 0.0;
      rate[Z252] = 0.0;
      rate[Z253] = 0.0;
      rate[Z254] = 0.0;
      rate[Z255] = 0.0;
      rate[Z256] = 0.0;
      rate[Z257] = 0.0;
      rate[Z258] = 0.0;
      rate[Z259] = 0.0;
      rate[Z260] = 0.0;
      rate[Z261] = 0.0;
      rate[Z262] = 0.0;
      rate[Z263] = 0.0;
      rate[Z264] = 0.0;
      rate[Z265] = 0.0;
      rate[Z266] = 0.0;
      rate[Z267] = 0.0;
      rate[Z268] = 0.0;
      rate[Z269] = 0.0;
      rate[Z270] = 0.0;
      rate[Z271] = 0.0;
      rate[Z272] = 0.0;
      rate[Z273] = 0.0;
      rate[Z274] = 0.0;
      rate[Z275] = 0.0;
      rate[Z276] = 0.0;
      rate[Z277] = 0.0;
      rate[Z278] = 0.0;
      rate[Z279] = 0.0;
      rate[Z280] = 0.0;
      rate[Z281] = 0.0;

      if (CRX > 0){
      rate[CRX1] = 0.0; 
      rate[CRX2] = 0.0;
      rate[CRX3] = 0.0;
      rate[CRX4] = 0.0;
      rate[CRX5] = 0.0;
      rate[CRX6] = 0.0;
      rate[CRX7] = 0.0;
      rate[CRX8] = 0.0;
      rate[CRX9] = 0.0;
      rate[CRX13] = 0.0;
      rate[CRX14] = 0.0;
      rate[CRX15] = 0.0;
      rate[CRX16] = 0.0;
      rate[CRX17] = 0.0;
      rate[CRX18] = 0.0;
      rate[CRX19] = 0.0;
      rate[CRX20] = 0.0;
      rate[CRX21] = 0.0;
      }
   }

  /* Hydrogen formation reactions */ 

  /* NOTE: H1-H6 rates are already calculated in calc_rates_g.F
           with updated values, so we provide the original 
           Omukai (2005) rates here only for the purposes of the 
           One-zone freefall test.*/

  if (water_rates == 1 && water_only)
  {
        //Deuterium reactions
    rate[D1] = 3.7e-10 * pow( temp, 0.28) * exp( -43.0 * tinv);
    rate[D2] = 3.7e-10 * pow( temp, 0.28);
    rate[D3] = 9.0e-11 * exp( -3876.0 * tinv);
    rate[D4] = 2.1e-9;
    rate[D5] = 3.2e-11 * exp( -3624.0 * tinv);
    rate[D6] = 1.0e-9 * exp( -464.0 * tinv); 
  }

   double Zdash = Z/Z_solar;
   double IUV;
   double IUV_shield;
   if ((int) UV) {
      IUV = 1.0;
      IUV_shield = 0.0;
      //if (water_only){
      //   IUV = 1.0;
      //} else {
      // double IUV_0 = 1.0;
      //   IUV = IUV_0 * exp( - 38.0 * Zdash);
      //}

      if (!UVbackground_molec_redshift_on){
      //if (water_only){
         rate[UV1] = 2.527e-10*IUV;
         rate[UV2] = 3.221e-11*IUV;
         rate[UV3] = 4.803e-10*IUV;
         rate[UV4] = 5.126e-10*IUV;
         rate[UV5] = 4.748e-10*IUV;
         rate[UV6] = 2.04e-10*IUV_shield; //off onless we get to 1.e5K, then we explode
         rate[UV7] = 3.15e-10*IUV_shield;
         rate[UV8] = 1.e-12*IUV_shield; //hokey addition for now (same rate as UV7)
         rate[UV9] = 2.07e-09*IUV;
         rate[UV10] = 1.389e-10*IUV;
         rate[UV11] = 3.36e-10*IUV;
         rate[UV12] = 4.989e-11*IUV;
         rate[UV13] = 2.343e-10*IUV;
         rate[UV14] = 9.492e-10*IUV;
         rate[UV15] = 1.104e-10*IUV;
         rate[UV16] = 7.363e-12*IUV;
         rate[UV17] = 1.087e-10*IUV;
         rate[UV18] = 4.31e-10*IUV;
         rate[UV19] = 5.39e-12*IUV_shield;
         rate[UV20] = 2.631e-11*IUV;
         rate[UV21] = 9.215e-16*IUV;
         rate[UV22] = 6.82e-12*IUV_shield;
         rate[UV23] = 3.14e-11*IUV_shield;
         rate[UV24] = 7.60e-11*IUV_shield;
         rate[UV25] = 5.00e-10*IUV;
         rate[UV26] = 3.00e-10*IUV;
         rate[UV34] = 2.80e-10*IUV;
         rate[UV37] = 5.60e-11*IUV;
         rate[UV38] = 5.60e-11*IUV;
      } else{
	   //Use the UV rates from the table! And interpolate where necessary
	   double Redshift = 1.0 / (my_units->a_value * my_units->a_units) - 1;
	   if (Redshift <= 0){ printf("Error! To run the water network with UV, start a cosmology run.\n"); }
	     // find interpolation index
	    if (Redshift <= my_rates->UVbackground_table.zmax_molec){
            double *zvec = my_rates->UVbackground_table.z_molec;
            double slope;
            int index=0;
            while (Redshift > zvec[index])
              index++;
            if(index == 0) index=1;
            if(index == my_rates->UVbackground_table.Nz_molec) index--;

                         /* *** UV1 *** */
            slope = (my_rates->UVbackground_table.UV1[index] -
           my_rates->UVbackground_table.UV1[index-1]) / (zvec[index] - zvec[index-1]);
           rate[UV1] = (Redshift - zvec[index-1]) * slope +
                my_rates->UVbackground_table.UV1[index-1];

	                /* *** UV2 *** */
            slope = (my_rates->UVbackground_table.UV2[index] -
           my_rates->UVbackground_table.UV2[index-1]) / (zvec[index] - zvec[index-1]);
           rate[UV2] = (Redshift - zvec[index-1]) * slope +
                my_rates->UVbackground_table.UV2[index-1];

	                /* *** UV3 *** */
            slope = (my_rates->UVbackground_table.UV3[index] -
           my_rates->UVbackground_table.UV3[index-1]) / (zvec[index] - zvec[index-1]);
           rate[UV3] = (Redshift - zvec[index-1]) * slope +
                my_rates->UVbackground_table.UV3[index-1];

	                /* *** UV4 *** */
            slope = (my_rates->UVbackground_table.UV4[index] -
           my_rates->UVbackground_table.UV4[index-1]) / (zvec[index] - zvec[index-1]);
           rate[UV4] = (Redshift - zvec[index-1]) * slope +
                my_rates->UVbackground_table.UV4[index-1];

	                /* *** UV5 *** */
            slope = (my_rates->UVbackground_table.UV5[index] -
           my_rates->UVbackground_table.UV5[index-1]) / (zvec[index] - zvec[index-1]);
           rate[UV5] = (Redshift - zvec[index-1]) * slope +
                my_rates->UVbackground_table.UV5[index-1];

	                /* *** UV6 *** */
            slope = (my_rates->UVbackground_table.UV6[index] -
           my_rates->UVbackground_table.UV6[index-1]) / (zvec[index] - zvec[index-1]);
           rate[UV6] = (Redshift - zvec[index-1]) * slope +
                my_rates->UVbackground_table.UV6[index-1];

	                /* *** UV7 *** */
            slope = (my_rates->UVbackground_table.UV7[index] -
           my_rates->UVbackground_table.UV7[index-1]) / (zvec[index] - zvec[index-1]);
           rate[UV7] = (Redshift - zvec[index-1]) * slope +
                my_rates->UVbackground_table.UV7[index-1];

	                /* *** UV8 *** */
            slope = (my_rates->UVbackground_table.UV8[index] -
           my_rates->UVbackground_table.UV8[index-1]) / (zvec[index] - zvec[index-1]);
           rate[UV8] = (Redshift - zvec[index-1]) * slope +
                my_rates->UVbackground_table.UV8[index-1];

	                /* *** UV9 *** */
	   /* Because UV9 is already in Grackle
            slope = (my_rates->UVbackground_table.UV9[index] -
           my_rates->UVbackground_table.UV9[index-1]) / (zvec[index] - zvec[index-1]);
           rate[UV9] = (Redshift - zvec[index-1]) * slope +
                my_rates->UVbackground_table.UV9[index-1];
		*/

	                /* *** UV10 *** */
            slope = (my_rates->UVbackground_table.UV10[index] -
           my_rates->UVbackground_table.UV10[index-1]) / (zvec[index] - zvec[index-1]);
           rate[UV10] = (Redshift - zvec[index-1]) * slope +
                my_rates->UVbackground_table.UV10[index-1];

	                /* *** UV11 *** */
            slope = (my_rates->UVbackground_table.UV11[index] -
           my_rates->UVbackground_table.UV11[index-1]) / (zvec[index] - zvec[index-1]);
           rate[UV11] = (Redshift - zvec[index-1]) * slope +
                my_rates->UVbackground_table.UV11[index-1];

	                /* *** UV12 *** */
            slope = (my_rates->UVbackground_table.UV12[index] -
           my_rates->UVbackground_table.UV12[index-1]) / (zvec[index] - zvec[index-1]);
           rate[UV12] = (Redshift - zvec[index-1]) * slope +
                my_rates->UVbackground_table.UV12[index-1];

	                /* *** UV13 *** */
            slope = (my_rates->UVbackground_table.UV13[index] -
           my_rates->UVbackground_table.UV13[index-1]) / (zvec[index] - zvec[index-1]);
           rate[UV13] = (Redshift - zvec[index-1]) * slope +
                my_rates->UVbackground_table.UV13[index-1];

	                /* *** UV14 *** */
            slope = (my_rates->UVbackground_table.UV14[index] -
           my_rates->UVbackground_table.UV14[index-1]) / (zvec[index] - zvec[index-1]);
           rate[UV14] = (Redshift - zvec[index-1]) * slope +
                my_rates->UVbackground_table.UV14[index-1];

	                /* *** UV15 *** */
            slope = (my_rates->UVbackground_table.UV15[index] -
           my_rates->UVbackground_table.UV15[index-1]) / (zvec[index] - zvec[index-1]);
           rate[UV15] = (Redshift - zvec[index-1]) * slope +
                my_rates->UVbackground_table.UV15[index-1];

	                /* *** UV16 *** */
            slope = (my_rates->UVbackground_table.UV16[index] -
           my_rates->UVbackground_table.UV16[index-1]) / (zvec[index] - zvec[index-1]);
           rate[UV16] = (Redshift - zvec[index-1]) * slope +
                my_rates->UVbackground_table.UV16[index-1];

	                /* *** UV17 *** */
            slope = (my_rates->UVbackground_table.UV17[index] -
           my_rates->UVbackground_table.UV17[index-1]) / (zvec[index] - zvec[index-1]);
           rate[UV17] = (Redshift - zvec[index-1]) * slope +
                my_rates->UVbackground_table.UV17[index-1];

	                /* *** UV18 *** */
	   /* commenting out because UV9, UV18 are definitely already in grackle
            slope = (my_rates->UVbackground_table.UV18[index] -
           my_rates->UVbackground_table.UV18[index-1]) / (zvec[index] - zvec[index-1]);
           rate[UV18] = (Redshift - zvec[index-1]) * slope +
                my_rates->UVbackground_table.UV18[index-1];
		*/

	                /* *** UV19 *** */
            slope = (my_rates->UVbackground_table.UV19[index] -
           my_rates->UVbackground_table.UV19[index-1]) / (zvec[index] - zvec[index-1]);
           rate[UV19] = (Redshift - zvec[index-1]) * slope +
                my_rates->UVbackground_table.UV19[index-1];

	                /* *** UV20 *** */
            slope = (my_rates->UVbackground_table.UV20[index] -
           my_rates->UVbackground_table.UV20[index-1]) / (zvec[index] - zvec[index-1]);
           rate[UV20] = (Redshift - zvec[index-1]) * slope +
                my_rates->UVbackground_table.UV20[index-1];

	                /* *** UV21 *** */
            slope = (my_rates->UVbackground_table.UV21[index] -
           my_rates->UVbackground_table.UV21[index-1]) / (zvec[index] - zvec[index-1]);
           rate[UV21] = (Redshift - zvec[index-1]) * slope +
                my_rates->UVbackground_table.UV21[index-1];

	                /* *** UV22 *** */
            slope = (my_rates->UVbackground_table.UV22[index] -
           my_rates->UVbackground_table.UV22[index-1]) / (zvec[index] - zvec[index-1]);
           rate[UV22] = (Redshift - zvec[index-1]) * slope +
                my_rates->UVbackground_table.UV22[index-1];

	                /* *** UV23 *** */
            slope = (my_rates->UVbackground_table.UV23[index] -
           my_rates->UVbackground_table.UV23[index-1]) / (zvec[index] - zvec[index-1]);
           rate[UV23] = (Redshift - zvec[index-1]) * slope +
                my_rates->UVbackground_table.UV23[index-1];

	                /* *** UV24 *** */
            slope = (my_rates->UVbackground_table.UV24[index] -
           my_rates->UVbackground_table.UV24[index-1]) / (zvec[index] - zvec[index-1]);
           rate[UV24] = (Redshift - zvec[index-1]) * slope +
                my_rates->UVbackground_table.UV24[index-1];

	                /* *** UV25 *** */
            slope = (my_rates->UVbackground_table.UV25[index] -
           my_rates->UVbackground_table.UV25[index-1]) / (zvec[index] - zvec[index-1]);
           rate[UV25] = (Redshift - zvec[index-1]) * slope +
                my_rates->UVbackground_table.UV25[index-1];

	                /* *** UV26 *** */
            slope = (my_rates->UVbackground_table.UV26[index] -
           my_rates->UVbackground_table.UV26[index-1]) / (zvec[index] - zvec[index-1]);
           rate[UV26] = (Redshift - zvec[index-1]) * slope +
                my_rates->UVbackground_table.UV26[index-1];

	                /* *** UV34 *** */
            slope = (my_rates->UVbackground_table.UV34[index] -
           my_rates->UVbackground_table.UV34[index-1]) / (zvec[index] - zvec[index-1]);
           rate[UV34] = (Redshift - zvec[index-1]) * slope +
                my_rates->UVbackground_table.UV34[index-1];

	                /* *** UV37 *** */
            slope = (my_rates->UVbackground_table.UV37[index] -
           my_rates->UVbackground_table.UV37[index-1]) / (zvec[index] - zvec[index-1]);
           rate[UV37] = (Redshift - zvec[index-1]) * slope +
                my_rates->UVbackground_table.UV37[index-1];

	                /* *** UV38 *** */
            slope = (my_rates->UVbackground_table.UV38[index] -
           my_rates->UVbackground_table.UV38[index-1]) / (zvec[index] - zvec[index-1]);
           rate[UV38] = (Redshift - zvec[index-1]) * slope +
                my_rates->UVbackground_table.UV38[index-1];

	   //estimating H2 shielding
//	   if ((H2_shield > 0) & (Nh2 > 1.e22)){
           if (H2_shield > 0){
              rate[UV6] = 0.0;
              rate[UV7] = 0.0;
              rate[UV8] = 0.0;
              rate[UV19] = 0.0;
              rate[UV22] = 0.0;
              rate[UV23] = 0.0;
              rate[UV24] = 0.0;
           }

	 }
      }

   }

   if (water_rates == 3)
   {
     /* Hydrogen formation reactions */
     /* NOTE: H1-H6 rates are already calculated in calc_rates_g.F
           with updated values, so we provide the original 
           Omukai (2005) rates here only for the purposes of the 
           One-zone freefall test. */

       double xH2,zeta, zeta_16; 
       double zeta_set [4];
       zeta_16 = 0.0;
       if (CRX > 0){
	  // make sure our CRX flag is feasible
	  if (CRX < 5){
      	     xH2 = (Y[H2m] + Y[H2plus])/n;
	     double zeta_set[] = {1.e-20, 1.e-16, 1.e-15, 5.e-15};
             zeta = zeta_set[(CRX-1)]; //s^-1
	  } else {
	     printf("Error: CRX flag not understood.\n");
	     return -1;
	  }
	  zeta_16 = zeta/1.e-16; //s^-1
       }
       
       /* Not implemented yet - added physics to be considerred
       // normalized dust-to-gas ratio (Eqn. 5 
       // of Bialy & Sternberg, 2019)
       double Zd_dash;
       double Z0 = 0.2;
       if (Zdash >= Z0) {
          Zd_dash = Zdash;  
       } else { 
          Zd_dash = Z0 * pow((Zdash/Z0),3.0);
       } 
       */
       

       /****** BIALY RATES **********/
      //implementing H4, Formula A1 from Martin et al. (1996)
      //NOTE: we take the full dissociation rate as the 
      // sum of contributions from collisional-induced 
      // excitations and dissociative tunneling rates

      double nH = Y[H] + Y[Hmin];
      double alpha[21];

      // Starting off with collisional-induced coefficients
      alpha[0]  =  -1.784239e2;
      alpha[1]  =  -6.842243e1;
      alpha[2]  =   4.320243e1;
      alpha[3]  =  -4.633167e0;
      alpha[4]  =   6.970086e1;
      alpha[5]  =   4.087038e4;
      alpha[6]  =  -2.370570e4;
      alpha[7]  =   1.288953e2;
      alpha[8]  =  -5.391334e1;
      alpha[9]  =  5.315517e0;
      alpha[10] = -1.973427e1;
      alpha[11] =  1.678095e4;
      alpha[12] = -2.578611e4;
      alpha[13] =  1.482123e1;
      alpha[14] = -4.890915e0;
      alpha[15] =  4.749030e-1;
      alpha[16] = -1.338283e2;
      alpha[17] = -1.164408e0;
      alpha[18] =  8.227443e-1;
      alpha[19] =  5.864073e-1;
      alpha[20] = -2.056313e0;

      double gamma_CIE = calc_H4(alpha, nH, temp);

      // Moving on to the dissociative tunneling coefficients 
      alpha[0]  =  -1.427664e2;
      alpha[1]  =  4.270741e1;
      alpha[2]  =  -2.027365e0;
      alpha[3]  =  -2.582097e-1;
      alpha[4]  =  2.136094e1;
      alpha[5]  =  2.753531e4;
      alpha[6]  =  -2.146779e4;
      alpha[7]  =  6.034928e1;
      alpha[8]  =  -2.743096e1;
      alpha[9]  =  2.676150e0;
      alpha[10] =  -1.128215e1;
      alpha[11] =  1.425455e4;
      alpha[12] =  -2.312520e4;
      alpha[13] =  9.305564e0;
      alpha[14] =  -2.464009e0;
      alpha[15] =  1.985955e-1;
      alpha[16] =  7.430600e2;
      alpha[17] =  -1.174242e0;
      alpha[18] =  7.502286e-1;
      alpha[19] =  2.358848e-1;
      alpha[20] =  2.937507e0;

      double gamma_dt = calc_H4(alpha, nH, temp);

//      if (water_only){
      // We combine these for the full collisional rate
//         rate[H4] = gamma_CIE + gamma_dt;
 //     }

      if(water_only){
        rate[H1] = 3.50e-12 * pow(tov300, -0.75);
        rate[H2] = 3.37e-16 * pow(tov300, 0.64) * exp(-9.2 * tinv);
        rate[H3] = 4.32e-09 * pow(tov300, -0.39) * exp(-39.4 * tinv);
        rate[UV9] = 2.07e-09*IUV;
        rate[UV18] = 4.31e-10*IUV; //these two are already in Grackle! So turn them on only if running tests
         rate[H4] = gamma_CIE + gamma_dt;
         //Here we assume that Zd_dash = Zdash, 
         //not the broken pow law used in Bialy & Sternberg (2019)
         rate[H8] = 3.00e-17 * pow(T2, 0.5) * Zdash;
         rate[H8] *= n / Y[H];
      } 
      //rates H5 and H6 not present in Bialy's expanded 2019 network... 
      rate[H7] = 1.0e-8 * exp(-84100.0 * tinv);

      rate[H9] = 2.08e-09;
      rate[H10] = 4.36e-08*pow(tov300, -0.52);
      rate[H11] = 2.34e-08*pow(tov300, -0.52);
      rate[H12] = 1.20e-15*pow(tov300, 0.25);
      rate[H13] = 7.20e-15;
      rate[H14] = 3.70e-14* exp(-3.50e+01 * tinv);
      rate[H15] = 1.30e-10;
      rate[H16] = 1.50e-09;
      rate[H17] = 9.10e-10;
      rate[H18] = 5.26e-20*pow(tov300, -0.51);
      rate[H19] = 7.51e-08*pow(tov300, -0.50);
      rate[H20] = 1.15e-18*pow(tov300, 1.49) * exp(-2.28e+02 * tinv);
      rate[H21] = 6.40e-10;
      rate[H22] = 8.41e-09*pow(tov300, 0.11 )* exp(-1.02e+05 * tinv);
      rate[H23] = 7.51e-08*pow(tov300, -0.50);
      rate[H24] = 7.51e-08*pow(tov300, -0.50);
      rate[H25] = 7.51e-08*pow(tov300, -0.50);

      //Metallicity Reactions
      rate[Z1]  = 5.66e-10 * pow(tov300, 0.36) * exp(8.6 * tinv);
      rate[Z2]  = 6.86e-10 * pow(tov300, 0.26) * exp(-224.0 * tinv);
      rate[Z3]  = 1.70e-9;
      rate[Z4]  = 1.01e-9;
      rate[Z5]  = 6.40e-10;
      rate[Z6]  = 3.90e-8   * pow(tov300, -0.50);
      rate[Z7]  = 8.60e-8   * pow(tov300, -0.50);
      rate[Z8]  = 3.05e-7   * pow(tov300, -0.50);
      rate[Z9]  = 7.09e-8   * pow(tov300, -0.50);
      rate[Z10] = 9.90e-19 * pow(tov300, -0.38);
      rate[Z11] = 3.14e-13 * pow(tov300, 2.70)   * exp (-3150.0 * tinv);
      rate[Z12] = 2.05e-12 * pow(tov300, 1.52)   * exp (-1740.0 * tinv);
      rate[Z13] = 1.65e-12 * pow(tov300, 1.14)   * exp (-50.0   * tinv);
      rate[Z14] = 2.10e-9 * pow(tov300, -0.50);
      rate[Z15] = 6.90e-9 * pow(tov300, -0.50);
      rate[Z16] = 6.99e-14 * pow(tov300, 2.80)   * exp (-1950.0  * tinv);
      rate[Z17] = 1.59e-11 * pow(tov300, 1.20)   * exp (-9610.0 * tinv);
      rate[Z18] = 4.90e-20 * pow(tov300, 1.58);
      rate[Z19] = 3.69e-11 * pow(tov300, -0.27)   * exp (-12.9 * tinv);
      rate[Z20] = 2.61e-10 * exp(-8160.0 * tinv);
      rate[Z21] = 2.00e-9;
      rate[Z22] = 1.95e-7  * pow (tov300, -0.70);
      rate[Z23] = 6.02e-11 * pow (tov300, 0.10) * exp(4.50 * tinv);
      rate[Z24] = 1.10e-10;
      rate[Z25] = 5.56e-11 * pow (tov300, 0.41) * exp(26.9 * tinv);
      rate[Z26] = 4.54e-10;
      rate[Z28] = 2.36e-12 * pow(tov300, -0.29) * exp(17.6 * tinv);
      rate[Z29] = 7.70e-10 * pow(tov300, -0.50);
      rate[Z30] = 7.50e-10;
      rate[Z31] = 1.00e-17;
      rate[Z32] = 6.64e-10 * exp(-11700.0 * tinv);
      rate[Z33] = 1.31e-10 * exp(-80.0 * tinv);
      rate[Z34] = 5.46e-10 * exp(-1940.0 * tinv);
      rate[Z35] = 2.20e-10;
      rate[Z36] = 5.18e-11 * pow(tov300, 0.17) * exp(-6400.0 * tinv);
      rate[Z37] = 1.00e-10 * exp(-7600.0 * tinv);
      rate[Z38] = 6.86e-14 * pow(tov300, 2.74) * exp(-4740.0 * tinv);
      rate[Z39] = 5.94e-13 * pow(tov300, 3.00) * exp(-4040.0 * tinv);
      rate[Z40] = 1.00e-17;
      rate[Z43] = 8.00e-11;
      rate[Z44] = 1.33e-10;
      rate[Z45] = 1.60e-09;
      rate[Z46] = 3.42e-10; 
      rate[Z48] = 7.98e-10 * pow(tov300, -0.16) * exp(-1.40 * tinv);
      rate[Z49] = 5.37e-08 * pow(tov300, -0.50);
      rate[Z50] = 5.60e-09 * pow(tov300, -0.50);
      rate[Z51] = 3.05e-07 * pow(tov300, -0.50);
      rate[Z52] = 3.75e-08 * pow(tov300, -0.50);
      rate[Z53] = 1.36e-09 * pow(tov300, -0.14)* exp(3.40e+00 * tinv);

      rate[Z58] = 9.00e-10*pow(tov300, -0.50);
      rate[Z59] = 2.50e-09*pow(tov300, -0.50);
      rate[Z60] = 5.90e-09*pow(tov300, -0.50);
      rate[Z61] = 2.00e-09;
      rate[Z62] = 1.20e-09;
      rate[Z63] = 2.00e-16*pow(tov300, -1.30)* exp(-2.30e+01 * tinv);
      rate[Z64] = 1.60e-09;
      rate[Z65] = 7.75e-08*pow(tov300, -0.50);
      rate[Z66] = 1.95e-07*pow(tov300, -0.50);
      rate[Z67] = 2.00e-07*pow(tov300, -0.40);
      rate[Z68] = 1.50e-07*pow(tov300, -0.42);
      rate[Z69] = 1.60e-07*pow(tov300, -0.60);
      rate[Z70] = 4.03e-07*pow(tov300, -0.60);
      rate[Z71] = 7.68e-08*pow(tov300, -0.60);
      rate[Z72] = 9.06e-10*pow(tov300, -0.37)* exp(-2.91e+01 * tinv);

      rate[Z74] = 7.50e-10;
      rate[Z75] = 2.40e-07*pow(tov300, -0.69);
      rate[Z76] = 1.09e-11*pow(tov300, -2.19)* exp(-1.65e+02 * tinv);

      rate[Z77] = 6.44e-10;
      rate[Z78] = 4.90e-12*pow(tov300, 0.50)* exp(-4.58e03 * tinv);
      rate[Z79] = 2.16e-09;
      rate[Z80] = 5.00e-10;
      rate[Z81] = 1.05e-09;
      rate[Z82] = 1.10e-10 *pow(tov300, 0.50) * exp(-7.77e04 * tinv);
      rate[Z83] = 7.51e-08 * pow(tov300, -0.50);
      rate[Z84] = 1.00e-09;
      rate[Z85] = 1.00e-09;
      rate[Z86] = 1.00e-09;
      rate[Z87] = 1.00e-10;
      rate[Z88] = 1.00e-09;
      rate[Z89] = 1.00e-10;
      rate[Z90] = 6.00e-09* exp(-4.02e+04 * tinv);
      rate[Z91] = 5.80e-09* exp(-5.29e+04 * tinv);
      rate[Z92] = 6.00e-09*exp(-5.23e+04 * tinv); 
      rate[Z93] = 6.00e-09* exp(-5.09e+04 * tinv); 
      rate[Z94] = 6.00e-09* exp(-4.02e+04 * tinv);
      rate[Z95] = 5.80e-09* exp(-5.29e+04 * tinv); 
      rate[Z96] = 6.00e-09* exp(-5.23e+04 * tinv); 
      rate[Z97] = 6.00e-09* exp(-5.09e+04 * tinv); 
      rate[Z98] = 5.20e-10;
      rate[Z99] = 3.80e-10*pow(tov300, -0.50 );
      rate[Z100] = 1.10e-10;
      rate[Z101] = 5.20e-11;
      rate[Z102] = 4.30e-10;
      rate[Z103] = 4.70e-10;
      rate[Z104] = 9.70e-10;
      rate[Z105] = 4.30e-10;
      rate[Z106] = 4.80e-10;
      rate[Z107] = 3.90e-10;
      rate[Z108] = 7.93e-10;
      rate[Z109] = 3.20e-10*pow(tov300, -0.50);
      rate[Z110] = 3.40e-10*pow(tov300, -0.50);
      rate[Z111] = 3.50e-10*pow(tov300, -0.50);
      rate[Z112] = 3.10e-10*pow(tov300, -0.50);
      rate[Z113] = 3.50e-10*pow(tov300, -0.50);
      rate[Z114] = 1.20e-10;
      rate[Z115] = 1.40e-09;
      rate[Z116] = 3.40e-09;
      rate[Z117] = 1.50e-09;
      rate[Z118] = 1.90e-09*pow(tov300, -0.50);
      rate[Z119] = 1.00e-09;
      rate[Z120] = 1.40e-09;
      rate[Z121] = 7.10e-10*pow(tov300, -0.50);
      rate[Z122] = 3.90e-09*pow(tov300, -0.50);
      rate[Z123] = 8.00e-10;
      rate[Z124] = 7.60e-10*pow(tov300, -0.50);
      rate[Z125] = 4.60e-10;
      rate[Z126] = 1.72e-09*pow(tov300, -0.50);
      rate[Z127] = 6.30e-15*pow(tov300, 0.75);
      rate[Z128] = 5.10e-11;
      rate[Z129] = 5.00e-10*pow(tov300, -0.50);
      rate[Z130] = 6.05e-11*pow(tov300, -0.50);
      rate[Z131] = 3.30e-11;
      rate[Z132] = 8.90e-10;
      rate[Z133] = 3.20e-09*pow(tov300, -0.50);
      rate[Z134] = 1.90e-11;
      rate[Z135] = 3.60e-10*pow(tov300, -0.50);
      rate[Z136] = 1.40e-10;
      rate[Z137] = 1.59e-09*pow(tov300, -0.50);
      rate[Z138] = 5.90e-10;
      rate[Z139] = 3.10e-10*pow(tov300, -0.50);
      rate[Z140] = 1.75e-07*pow(tov300, -0.50);
      rate[Z141] = 1.75e-07*pow(tov300, -0.50);
      rate[Z142] = 4.76e-08*pow(tov300, -0.52);
      rate[Z143] = 1.40e-08*pow(tov300, -0.52);
      rate[Z144] = 1.96e-07*pow(tov300, -0.52);
      rate[Z145] = 1.40e-08*pow(tov300, -0.52);
      rate[Z146] = 8.40e-09*pow(tov300, -0.52);
      rate[Z147] = 2.00e-07*pow(tov300, -0.48);
      rate[Z148] = 1.60e-08*pow(tov300, -0.43);
      rate[Z149] = 1.00e-08*pow(tov300, -0.60);
      rate[Z150] = 3.00e-07*pow(tov300, -0.50);
      rate[Z151] = 1.20e-09;
      rate[Z152] = 1.10e-09;
      rate[Z153] = 1.00e-11;
      rate[Z154] = 1.10e-09;
      rate[Z155] = 5.20e-11;
      rate[Z156] = 1.00e-09;
      rate[Z157] = 1.20e-09;
      rate[Z158] = 5.80e-10*pow(tov300, -0.50);
      rate[Z159] = 2.90e-09*pow(tov300, -0.50);
      rate[Z160] = 1.00e-11;
      rate[Z161] = 9.70e-10;
      rate[Z162] = 3.50e-10;
      rate[Z163] = 7.50e-10*pow(tov300, -0.50);
      rate[Z164] = 9.10e-10;
      rate[Z165] = 7.50e-10;
      rate[Z166] = 9.60e-10;
      rate[Z167] = 4.30e-10;
      rate[Z168] = 4.70e-10;
      rate[Z169] = 9.40e-10;
      rate[Z170] = 8.60e-10;
      rate[Z171] = 8.50e-10;
      rate[Z172] = 4.80e-10;
      rate[Z173] = 4.00e-10;
      rate[Z174] = 1.50e-09;
      rate[Z175] = 1.40e-09;
      rate[Z176] = 2.60e-09*pow(tov300, -0.50);
      rate[Z177] = 4.55e-10;
      rate[Z178] = 1.40e-09;
      rate[Z179] = 1.95e-10;
      rate[Z180] = 1.31e-09;
      rate[Z181] = 1.00e-09;
      rate[Z182] = 3.70e-09*pow(tov300, -0.50);
      rate[Z183] = 6.90e-10*pow(tov300, -0.50);
      rate[Z184] = 3.20e-10*pow(tov300, -0.50);
      rate[Z185] = 3.40e-10*pow(tov300, -0.50);
      rate[Z186] = 6.80e-10*pow(tov300, -0.50);
      rate[Z187] = 6.30e-10*pow(tov300, -0.50);
      rate[Z188] = 3.50e-10*pow(tov300, -0.50);
      rate[Z189] = 3.10e-10*pow(tov300, -0.50);
      rate[Z190] = 6.20e-10*pow(tov300, -0.50);
      rate[Z191] = 3.50e-10*pow(tov300, -0.50);
      rate[Z192] = 8.40e-10;
      rate[Z193] = 1.40e-09;
      rate[Z194] = 2.30e-09;
      rate[Z195] = 2.40e-09;
      rate[Z196] = 1.00e-09;
      rate[Z197] = 2.30e-09;
      rate[Z198] = 1.14e-10;
      rate[Z199] = 7.10e-10*pow(tov300, -0.50);
      rate[Z200] = 3.40e-09*pow(tov300, -0.50);
      rate[Z201] = 1.90e-09;
      rate[Z202] = 1.50e-09;
      rate[Z203] = 7.60e-10*pow(tov300, -0.50);
      rate[Z204] = 1.00e-10* exp(-4.64e+03 * tinv);
      rate[Z205] = 4.89e-11*pow(tov300, -0.14)* exp(3.61e+01 * tinv);
      rate[Z206] = 6.40e-10;
      rate[Z207] = 2.10e-09*pow(tov300, -0.50);
      rate[Z208] = 8.84e-10*pow(tov300, -0.50);
      rate[Z209] = 8.20e-10*pow(tov300, -0.50);
      rate[Z210] = 1.70e-09;
      rate[Z211] = 2.10e-09;
      rate[Z212] = 2.40e-09;
      rate[Z213] = 1.20e-09*pow(tov300, -0.50);
      rate[Z214] = 9.30e-10* exp(-1.00e+02 * tinv);
      rate[Z215] = 3.42e-10*pow(tov300, -0.16)* exp(-1.40e+00 * tinv);
      rate[Z216] = 1.30e-09*pow(tov300, -0.50);
      rate[Z217] = 1.00e-09* exp(-7.08e+03 * tinv);
      rate[Z218] = 7.00e-10* exp(-1.06e+04 * tinv);
      rate[Z219] = 1.00e-11;
      rate[Z220] = 1.50e-10;
      rate[Z221] = 7.50e-10;
      rate[Z222] = 7.50e-10;
      rate[Z223] = 1.80e-09;
      rate[Z224] = 2.40e-10;
      rate[Z225] = 9.50e-10;
      rate[Z226] = 8.50e-11;
      rate[Z227] = 4.80e-10;
      rate[Z228] = 1.10e-09*pow(tov300, -0.50);
      rate[Z229] = 2.86e-10*pow(tov300, -0.50);
      rate[Z230] = 2.04e-10*pow(tov300, -0.50);
      rate[Z231] = 1.10e-09;
      rate[Z232] = 1.10e-09*pow(tov300, -0.50);
      rate[Z233] = 1.10e-10;
      rate[Z234] = 3.60e-10*pow(tov300, -0.50);
      rate[Z235] = 1.00e-09;
      rate[Z236] = 2.20e-10;
      rate[Z237] = 4.00e-11;
      rate[Z238] = 6.20e-10;
      rate[Z239] = 7.10e-10;
      rate[Z240] = 1.30e-09*pow(tov300, -0.50);
      rate[Z241] = 7.00e-10*pow(tov300, -0.50);
      rate[Z242] = 7.00e-10*pow(tov300, -0.50);
      rate[Z243] = 3.10e-10*pow(tov300, -0.50);
      rate[Z244] = 6.90e-10*pow(tov300, -0.50);
      rate[Z245] = 6.20e-10*pow(tov300, -0.50);
      rate[Z246] = 6.10e-10*pow(tov300, -0.50);
      rate[Z247] = 7.51e-08*pow(tov300, -0.50);
      rate[Z248] = 7.51e-08*pow(tov300, -0.50);
      rate[Z249] = 3.76e-08*pow(tov300, -0.50);
      rate[Z250] = 7.51e-08*pow(tov300, -0.50);
      rate[Z251] = 2.69e-12* exp(-2.36e+04 * tinv);
      rate[Z252] = 2.25e-11*pow(tov300, 0.50)* exp(-1.48e+04 * tinv);
      rate[Z253] = 4.00e-10* exp(-5.00e+03 * tinv);
      rate[Z254] = 7.13e-12* exp(-5.05e+03 * tinv);
      rate[Z255] = 2.48e-10*pow(tov300, -3.30)* exp(-1.44e+03 * tinv);
      rate[Z256] = 4.98e-10* exp(-6.00e+03 * tinv);
      rate[Z257] = 1.44e-11*pow(tov300, 0.50)* exp(-3.00e+03 * tinv);
      rate[Z258] = 1.44e-11*pow(tov300, 0.50)* exp(-3.00e+03 * tinv);
      rate[Z259] = 7.13e-12* exp(-5.05e+03 * tinv);
      rate[Z260] = 2.30e-15*pow(tov300, 3.47)* exp(-6.68e+03 * tinv);
      rate[Z261] = 3.60e-11* exp(-2.02e+02 * tinv);
      rate[Z262] = 3.27e-14*pow(tov300, 2.20)* exp(-2.24e+03 * tinv);
      rate[Z263] = 1.20e-10* exp(-1.40e+03 * tinv);
      rate[Z264] = 3.77e-13*pow(tov300, 2.42)* exp(-1.16e+03 * tinv);
      rate[Z265] = 1.14e-11*pow(tov300, -0.48);
      rate[Z266] = 7.60e-12*pow(tov300, -0.48);
      rate[Z267] = 2.52e-11* exp(-2.38e+03 * tinv);
      rate[Z268] = 3.16e-10* exp(-2.19e+04 * tinv);
      rate[Z269] = 2.29e-12*pow(tov300, 2.20)* exp(-3.82e+03 * tinv);
      rate[Z270] = 1.85e-11*pow(tov300, 0.95)* exp(-8.57e+03 * tinv);
      rate[Z271] = 3.14e-18*pow(tov300, -0.15)* exp(-6.80e+01 * tinv);
      rate[Z272] = 5.00e-10*pow(tov300, -3.70)* exp(-8.00e+02 * tinv);
      rate[Z273] = 4.69e-19*pow(tov300, 1.52)* exp(5.05e+01 * tinv);
      rate[Z274] = 3.92e-16*pow(tov300, -2.29)* exp(-2.13e+01 * tinv);
      rate[Z275] = 5.09e-18*pow(tov300, -0.71)* exp(-1.16e+01 * tinv);
      rate[Z276] = 1.70e-17;
      rate[Z277] = 5.26e-18*pow(tov300, -5.22)* exp(-9.00e+01 * tinv);
      rate[Z278] = 1.10e-10*pow(tov300, -0.50);
      rate[Z279] = 5.36e-12*pow(tov300, -0.50);
      rate[Z280] = 3.24e-12*pow(tov300, -0.66);
      rate[Z281] = 2.50e-10* exp(-2.12e+04 * tinv);

      if (CRX > 0){
      rate[CRX1] = 1.02e+03 * zeta * (2*xH2) ;
      rate[CRX2] = 2.00e+01 * zeta * (2*xH2) ;
      rate[CRX3] = 3.52e+02 * zeta * (2*xH2) ;
      rate[CRX4] = 1.46e+03 * zeta * (2*xH2) ;
      rate[CRX5] = 1.94e+03 * zeta * (2*xH2) ;
      rate[CRX6] = 1.02e+03 * zeta * (2*xH2) ;
      rate[CRX7] = 4.68e+03 * zeta * (2*xH2) ;
      rate[CRX8] = 2.34e+02 * zeta * (2*xH2) ;
      rate[CRX9] = 1.50e+03 * zeta * (2*xH2) ;
      rate[CRX13] = 1.00e+03  * zeta * (2*xH2) ;
      rate[CRX14] = 1.00e+03 * zeta * (2*xH2) ;
      rate[CRX15] = 1.00e+03 * zeta * (2*xH2) ;
      rate[CRX16] = 1.00e+03 * zeta * (2*xH2) ;

      rate[CRX17] = 1.09e-16 * zeta_16;
      rate[CRX18] = 1.00e-16 * zeta_16;
      rate[CRX19] = 2.02e-16 * zeta_16;
      rate[CRX20] = 4.35e-18 * zeta_16;
      
      rate[CRX21] = 1.00e+03 * zeta * (2*xH2) ;
      }

    } else{

      double klow = 1.18e-10 * exp(-6.95e4 * tinv);
      double khigh = 8.125e-8 * pow( temp, -0.5) * exp(-5.2e4*tinv) * ( (double) 1.0 - exp(-6.0e3 *tinv));
      double ep2 = 4.845 - 1.3* log10( temp * 1.0e-4) + 1.62* pow( log10( temp * 1.0e-4) , 2.0);
      double ncrit = pow(10, ep2);
      double a = (double) 1.0 / ( (double) 1.0 + n / ncrit);
      rate[H7] = pow( khigh, 1.0 - a) * pow( klow, a);

      rate[D1] = 3.7e-10*pow(temp,0.25)* exp(-43.0*tinv);
      rate[D2] = 3.7e-10*pow(temp,0.28);
      rate[D3] = 9.00e-11*exp(-3876.0*tinv);
      rate[D4] = 2.1e-9;
      rate[D5] = 3.2e-11 * exp(-3624.0*tinv);
      rate[D6] = 1.0e-9*exp(-464.0*tinv);

      //Metallicity Reactions
      rate[Z1]  = 6.80e-10;
      rate[Z2]  = 7.00e-10 * exp( -232.0 * tinv);
      rate[Z3]  = 1.70e-9;
      rate[Z4]  = 1.01e-9;
      rate[Z5]  = 8.30e-10;
      rate[Z6]  = 2.00e-7   * pow ( tov300, -0.50 );
      rate[Z7]  = 1.60e-7   * pow ( tov300, -0.50 );
      rate[Z8]  = 6.50e-7   * pow ( tov300, -0.50 );
      rate[Z9]  = 3.50e-7   * pow ( tov300, -0.50 );
      rate[Z10] = 9.90e-19 * pow ( tov300, -0.38 );
      rate[Z11] = 3.43e-13 * pow ( tov300, 2.67  )   * exp ( -3160.0 * tinv );
      rate[Z12] = 1.55e-12 * pow ( tov300, 1.60  )   * exp ( -1660.0 * tinv );
      rate[Z13] = 1.65e-12 * pow ( tov300, 1.14  )   * exp ( -50.0   * tinv );
      rate[Z14] = 2.10e-9;
      rate[Z15] = 6.90e-9;
      rate[Z16] = 7.00e-14 * pow ( tov300, 2.80  )   * exp ( -1950.0  * tinv );
      rate[Z17] = 6.83e-12 * pow ( tov300, 1.60  )   * exp ( -9720.0 * tinv );
      rate[Z18] = 4.90e-20 * pow ( tov300, 1.58  );
      rate[Z19] = 4.34e-11 * pow ( tov300, -0.50 )   * exp ( -30.0    * tinv );
      rate[Z20] = 3.30e-10 * exp ( -8460.0 * tinv );
      rate[Z21] = 2.00e-9;
      rate[Z22] = 1.95e-7  * pow(  tov300, -0.70);
      rate[Z23] = 6.60e-11;
      rate[Z24] = 1.10e-10 * pow(  tov300, 0.50);
      rate[Z25] = 3.30e-11;
      rate[Z26] = 6.20e-10;
      rate[Z27] = 1.00e-13;
      rate[Z28] = 4.40e-12 * pow(  tov300, -0.61);
      rate[Z29] = 7.70e-10;
      rate[Z30] = 7.50e-10;
      rate[Z31] = 1.00e-17;
      rate[Z32] = 6.64e-10 * exp( -11700.0 * tinv);
      rate[Z33] = 4.98e-11;
      rate[Z34] = 2.38e-10 * exp( -1760.0 * tinv);
      rate[Z35] = 2.70e-10;
      rate[Z36] = 5.18e-11 * pow(  tov300, 0.17) * exp( -6400.0 * tinv);
      rate[Z37] = 1.00e-10* exp( -7600.0 * tinv);
      rate[Z38] = 6.86e-14 * pow(  tov300, 2.74) * exp( -4740.0 * tinv);
      rate[Z39] = 5.82e-13 * pow(  tov300, 3.00) * exp( -4045.0 * tinv);
      rate[Z40] = 1.00e-17;
    }
   return 0;
  }
