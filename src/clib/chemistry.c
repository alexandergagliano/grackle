#include <stdlib.h>
#include <stdio.h>
#include "grackle_chemistry_data.h"
#include "jac_rhs_builder.h"

/****************** Declaration of all chemical species ****************************/
// Multispecies = 1, withWater = 1
int H, Hplus, el, O, OH, H2O, O2, Oplus, OHplus, H2Oplus, H3Oplus, O2plus, Cplus, C, CH, CH2, CH3, CH4, CO, COplus, CO2; 

// Multispecies = 2, withWater = 1
int H2m, Hmin; 

// Multispecies = 3, withWater = 1
int D, Dplus, HD;

// withWater = 1, water_rates = 3 (Using Bialy's expanded network)
int CHplus, CH2plus, He, Heplus, H3plus, HCOplus, CH3plus, H2plus, HeHplus,
    CH4plus, CH5plus, O2Hplus;

// Needed for building reactions
int NRXN = -1;

// Global variables for nSpecies and nReactions
int nSpecies, nReactions;

/***************** Declaration of all chemical reactions ***************************/  
int H1, H2, H3, H4, H5, H6, H7, H8, H9, H10, H11, H12, H13, H14, H15, H16, H17, H18, H19, H20, 
    H21, H22, H23, H24, H25;

int D1, D2, D3, D4, D5, D6;

int Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8, Z9, Z10, Z11, Z12, Z13, Z14, Z15, Z16, 
    Z17, Z18, Z19, Z20, Z21, Z22, Z23, Z24, Z25, Z26, Z27, Z28, Z29, Z30, 
    Z31, Z32, Z33, Z34, Z35, Z36, Z37, Z38, Z39, Z40, Z43, Z44, 
    Z45, Z46, Z48, Z49, Z50, Z51, Z52, Z53, Z54, Z55, Z56, Z57, 
    Z58, Z59, Z60, Z61, Z62, Z63, Z64, Z65, Z66, Z67, Z68, Z69, 
    Z70, Z71, Z72, Z73, 
    Z74, Z75, Z76, Z77, Z78, Z79, Z80, Z81,
    Z82, Z83, Z84, Z85, Z85, Z86, Z87, Z88, Z89,
    Z90, Z91, Z92, Z93, Z94, Z95, Z96, Z97,
    Z98, Z99, Z100, Z101, Z102, Z103, Z104, Z105,
    Z106, Z107, Z108, Z109, Z110, Z111, Z112, Z113,
    Z114, Z115, Z116, Z117, Z118, Z119, Z120, Z121,
    Z122, Z123, Z124, Z125, Z126, Z127, Z128, Z129,
    Z130, Z131, Z132, Z133, Z134, Z135, Z136, Z137,
    Z138, Z139, Z140, Z141, Z142, Z143, Z144, Z145,
    Z146, Z147, Z148, Z149, Z150, Z151, Z152, Z153,
    Z154, Z155, Z156, Z157, Z158, Z159, Z160, Z161,
    Z162, Z163, Z164, Z165, Z166, Z167, Z168, Z169,
    Z170, Z171, Z172, Z173, Z174, Z175, Z176, Z177,
    Z178, Z179, Z180, Z181, Z182, Z183, Z184, Z185,
    Z186, Z187, Z188, Z189, Z190, Z191, Z192, Z193,
    Z194, Z195, Z196, Z197, Z198, Z199, Z200, Z201,
    Z202, Z203, Z204, Z205, Z206, Z207, Z208, Z209,
    Z210, Z211, Z212, Z213, Z214, Z215, Z216, Z217,
    Z218, Z219, Z220, Z221, Z222, Z223, Z224, Z225,
    Z226, Z227, Z228, Z229, Z230, Z231, Z232, Z233,
    Z234, Z235, Z236, Z237, Z238, Z239, Z240, Z241,
    Z242, Z243, Z244, Z245, Z246, Z247, Z248, Z249,
    Z250, Z251, Z252, Z253, Z254, Z255, Z256, Z257,
    Z258, Z259, Z260, Z261, Z262, Z263, Z264, Z265,
    Z266, Z267, Z268, Z269, Z270, Z271, Z272, Z273,
    Z274, Z275, Z276, Z277, Z278, Z279, Z280, Z281;

int UV1, UV2, UV3, UV4, UV5, UV6, UV7, UV8, UV9, UV10, UV11, UV12, UV13, UV14, UV15, UV16,
    UV17, UV18, UV19, UV20, UV21, UV22, UV23, UV24, UV25, UV26, UV27, UV28, UV29, UV30, UV31,
    UV32, UV33, UV34, UV35, UV36, UV37, UV38, UV39, UV40;

int CRX1, CRX2, CRX3, CRX4, CRX5, CRX6, CRX7, CRX8, CRX9, CRX13, CRX14, CRX15, CRX16, CRX17,
    CRX18, CRX19, CRX20, CRX21;

/*************** Initializing variables for reactions and species *****************/ 
void setup_rxns(int ispecies, int UV, int CRX, int water_rates)
{
   int rxn_counter = 0;

   H1 = rxn_counter;
   Z1 = ++rxn_counter;
   Z2 = ++rxn_counter;
   Z7 = ++rxn_counter;
   Z8 = ++rxn_counter; 
   Z9 = ++rxn_counter;
   Z10 = ++rxn_counter;
   Z13 = ++rxn_counter;
   Z14 = ++rxn_counter;
   Z15 = ++rxn_counter;
   Z18 = ++rxn_counter;
   Z19 = ++rxn_counter;
   Z20 = ++rxn_counter;
   Z21 = ++rxn_counter;
   Z22 = ++rxn_counter;
   Z23 = ++rxn_counter;
   Z24 = ++rxn_counter;
   Z25 = ++rxn_counter;
   Z26 = ++rxn_counter;

   if (water_rates != 3){
   	Z27 = ++rxn_counter;
   }
   Z28 = ++rxn_counter;
   Z29 = ++rxn_counter;
   Z30 = ++rxn_counter;
   Z31 = ++rxn_counter;
   
   if (ispecies > 1){
      H2 = ++rxn_counter;
      H3 = ++rxn_counter;
      H4 = ++rxn_counter;
      H5 = ++rxn_counter;
      H6 = ++rxn_counter;
      H7 = ++rxn_counter;
      H8 = ++rxn_counter;
      Z3 = ++rxn_counter;
      Z4 = ++rxn_counter;
      Z5 = ++rxn_counter;
      Z6 = ++rxn_counter;
      Z11 = ++rxn_counter;
      Z12 = ++rxn_counter;
      Z16 = ++rxn_counter;
      Z17 = ++rxn_counter;
      Z32 = ++rxn_counter;
      Z33 = ++rxn_counter;
      Z34 = ++rxn_counter;
      Z35 = ++rxn_counter;
      Z36 = ++rxn_counter;
      Z37 = ++rxn_counter;
      Z38 = ++rxn_counter;
      Z39 = ++rxn_counter;
      Z40 = ++rxn_counter;

      if (ispecies > 2){
         D1  = ++rxn_counter;
         D2  = ++rxn_counter;
         D3  = ++rxn_counter;
         D4  = ++rxn_counter;
         D5  = ++rxn_counter;
         D6  = ++rxn_counter;
      }
   }
     
   if (UV) {
      UV1 = ++rxn_counter;
      UV2 = ++rxn_counter;
      UV3 = ++rxn_counter;
      UV4 = ++rxn_counter;
      UV5 = ++rxn_counter;
      UV6 = ++rxn_counter;
      UV7 = ++rxn_counter;
      UV8 = ++rxn_counter;
      UV9 = ++rxn_counter;
      UV10 = ++rxn_counter;
      UV11 = ++rxn_counter;
      UV12 = ++rxn_counter;
      UV13 = ++rxn_counter;
      UV14 = ++rxn_counter;
      UV15 = ++rxn_counter;
      UV16 = ++rxn_counter;
      UV17 = ++rxn_counter;
      UV18 = ++rxn_counter;
      UV19 = ++rxn_counter;
      UV20 = ++rxn_counter;
      UV21 = ++rxn_counter;
      UV22 = ++rxn_counter;
      UV23 = ++rxn_counter;
      UV24 = ++rxn_counter;
      UV25 = ++rxn_counter;
      UV26 = ++rxn_counter;
      UV27 = ++rxn_counter;
      UV28 = ++rxn_counter;
      UV29 = ++rxn_counter;
      UV30 = ++rxn_counter;
      UV31 = ++rxn_counter;
      UV32 = ++rxn_counter;
      UV33 = ++rxn_counter;
      UV34 = ++rxn_counter;
      UV35 = ++rxn_counter;
      UV36 = ++rxn_counter;
      UV37 = ++rxn_counter;
      UV38 = ++rxn_counter;
      UV39 = ++rxn_counter;
      UV40 = ++rxn_counter;
   }
    
   if (water_rates == 3){
      H9  = ++rxn_counter;
      H10 = ++rxn_counter;
      H11 = ++rxn_counter;
      H12 = ++rxn_counter;
      H13 = ++rxn_counter;
      H14 = ++rxn_counter;
      H15 = ++rxn_counter;
      H16 = ++rxn_counter;
      H17 = ++rxn_counter;
      H18 = ++rxn_counter;
      H19 = ++rxn_counter;
      H20 = ++rxn_counter;
      H21 = ++rxn_counter;
      H22 = ++rxn_counter;
      H23 = ++rxn_counter;
      H24 = ++rxn_counter;
      H25 = ++rxn_counter;

      Z43 = ++rxn_counter;
      Z44 = ++rxn_counter;
      Z45 = ++rxn_counter;
      Z46 = ++rxn_counter;
      Z48 = ++rxn_counter;
      Z49 = ++rxn_counter;
      Z50 = ++rxn_counter;
      Z51 = ++rxn_counter;
      Z52 = ++rxn_counter;
      Z53 = ++rxn_counter;

      // DOUBLE-COUNTED RATES! 
      Z54 = ++rxn_counter;
      Z55 = ++rxn_counter;
      Z56 = ++rxn_counter;
      Z57 = ++rxn_counter;
      //
      Z58 = ++rxn_counter;
      Z59 = ++rxn_counter;
      Z60 = ++rxn_counter;
      Z61 = ++rxn_counter;
      Z62 = ++rxn_counter;
      Z63 = ++rxn_counter;
      Z64 = ++rxn_counter;
      Z65 = ++rxn_counter;
      Z66 = ++rxn_counter;
      Z67 = ++rxn_counter;
      Z68 = ++rxn_counter;
      Z69 = ++rxn_counter;
      Z70 = ++rxn_counter;
      Z71 = ++rxn_counter;
      Z72 = ++rxn_counter;
      // Double-counted
      Z73 = ++rxn_counter;
      //
      Z74 = ++rxn_counter;
      Z75 = ++rxn_counter;
      Z76 = ++rxn_counter;
      Z77 = ++rxn_counter;
      Z78 = ++rxn_counter;
      Z79 = ++rxn_counter;
      Z80 = ++rxn_counter;
      Z81 = ++rxn_counter;
      Z82 = ++rxn_counter;
      Z83 = ++rxn_counter;
      Z84 = ++rxn_counter;
      Z85 = ++rxn_counter;
      Z86 = ++rxn_counter;
      Z87 = ++rxn_counter;
      Z88 = ++rxn_counter;
      Z89 = ++rxn_counter;
      Z90 = ++rxn_counter;
      Z91 = ++rxn_counter;
      Z92 = ++rxn_counter;
      Z93 = ++rxn_counter;
      Z94 = ++rxn_counter;
      Z95 = ++rxn_counter;
      Z96 = ++rxn_counter;
      Z97 = ++rxn_counter;
      Z98 = ++rxn_counter;
      Z99 = ++rxn_counter;
      Z100 = ++rxn_counter;
      Z101 = ++rxn_counter;
      Z102 = ++rxn_counter;
      Z103 = ++rxn_counter;
      Z104 = ++rxn_counter;
      Z105 = ++rxn_counter;
      Z106 = ++rxn_counter;
      Z107 = ++rxn_counter;
      Z108 = ++rxn_counter;
      Z109 = ++rxn_counter;
      Z110 = ++rxn_counter;
      Z111 = ++rxn_counter;
      Z112 = ++rxn_counter;
      Z113 = ++rxn_counter;
      Z114 = ++rxn_counter;
      Z115 = ++rxn_counter;
      Z116 = ++rxn_counter;
      Z117 = ++rxn_counter;
      Z118 = ++rxn_counter;
      Z119 = ++rxn_counter;
      Z120 = ++rxn_counter;
      Z121 = ++rxn_counter;
      Z122 = ++rxn_counter;
      Z123 = ++rxn_counter;
      Z124 = ++rxn_counter;
      Z125 = ++rxn_counter;
      Z126 = ++rxn_counter;
      Z127 = ++rxn_counter;
      Z128 = ++rxn_counter;
      Z129 = ++rxn_counter;
      Z130 = ++rxn_counter;
      Z131 = ++rxn_counter;
      Z132 = ++rxn_counter;
      Z133 = ++rxn_counter;
      Z134 = ++rxn_counter;
      Z135 = ++rxn_counter;
      Z136 = ++rxn_counter;
      Z137 = ++rxn_counter;
      Z138 = ++rxn_counter;
      Z139 = ++rxn_counter;
      Z140 = ++rxn_counter;
      Z141 = ++rxn_counter;
      Z142 = ++rxn_counter;
      Z143 = ++rxn_counter;
      Z144 = ++rxn_counter;
      Z145 = ++rxn_counter;
      Z146 = ++rxn_counter;
      Z147 = ++rxn_counter;
      Z148 = ++rxn_counter;
      Z149 = ++rxn_counter;
      Z150 = ++rxn_counter;
      Z151 = ++rxn_counter;
      Z152 = ++rxn_counter;
      Z153 = ++rxn_counter;
      Z154 = ++rxn_counter;
      Z155 = ++rxn_counter;
      Z156 = ++rxn_counter;
      Z157 = ++rxn_counter;
      Z158 = ++rxn_counter;
      Z159 = ++rxn_counter;
      Z160 = ++rxn_counter;
      Z161 = ++rxn_counter;
      Z162 = ++rxn_counter;
      Z163 = ++rxn_counter;
      Z164 = ++rxn_counter;
      Z165 = ++rxn_counter;
      Z166 = ++rxn_counter;
      Z167 = ++rxn_counter;
      Z168 = ++rxn_counter;
      Z169 = ++rxn_counter;
      Z170 = ++rxn_counter;
      Z171 = ++rxn_counter;
      Z172 = ++rxn_counter;
      Z173 = ++rxn_counter;
      Z174 = ++rxn_counter;
      Z175 = ++rxn_counter;
      Z176 = ++rxn_counter;
      Z177 = ++rxn_counter;
      Z178 = ++rxn_counter;
      Z179 = ++rxn_counter;
      Z180 = ++rxn_counter;
      Z181 = ++rxn_counter;
      Z182 = ++rxn_counter;
      Z183 = ++rxn_counter;
      Z184 = ++rxn_counter;
      Z185 = ++rxn_counter;
      Z186 = ++rxn_counter;
      Z187 = ++rxn_counter;
      Z188 = ++rxn_counter;
      Z189 = ++rxn_counter;
      Z190 = ++rxn_counter;
      Z191 = ++rxn_counter;
      Z192 = ++rxn_counter;
      Z193 = ++rxn_counter;
      Z194 = ++rxn_counter;
      Z195 = ++rxn_counter;
      Z196 = ++rxn_counter;
      Z197 = ++rxn_counter;
      Z198 = ++rxn_counter;
      Z199 = ++rxn_counter;
      Z200 = ++rxn_counter;
      Z201 = ++rxn_counter;
      Z202 = ++rxn_counter;
      Z203 = ++rxn_counter;
      Z204 = ++rxn_counter;
      Z205 = ++rxn_counter;
      Z206 = ++rxn_counter;
      Z207 = ++rxn_counter;
      Z208 = ++rxn_counter;
      Z209 = ++rxn_counter;
      Z210 = ++rxn_counter;
      Z211 = ++rxn_counter;
      Z212 = ++rxn_counter;
      Z213 = ++rxn_counter;
      Z214 = ++rxn_counter;
      Z215 = ++rxn_counter;
      Z216 = ++rxn_counter;
      Z217 = ++rxn_counter;
      Z218 = ++rxn_counter;
      Z219 = ++rxn_counter;
      Z220 = ++rxn_counter;
      Z221 = ++rxn_counter;
      Z222 = ++rxn_counter;
      Z223 = ++rxn_counter;
      Z224 = ++rxn_counter;
      Z225 = ++rxn_counter;
      Z226 = ++rxn_counter;
      Z227 = ++rxn_counter;
      Z228 = ++rxn_counter;
      Z229 = ++rxn_counter;
      Z230 = ++rxn_counter;
      Z231 = ++rxn_counter;
      Z232 = ++rxn_counter;
      Z233 = ++rxn_counter;
      Z234 = ++rxn_counter;
      Z235 = ++rxn_counter;
      Z236 = ++rxn_counter;
      Z237 = ++rxn_counter;
      Z238 = ++rxn_counter;
      Z239 = ++rxn_counter;
      Z240 = ++rxn_counter;
      Z241 = ++rxn_counter;
      Z242 = ++rxn_counter;
      Z243 = ++rxn_counter;
      Z244 = ++rxn_counter;
      Z245 = ++rxn_counter;
      Z246 = ++rxn_counter;
      Z247 = ++rxn_counter;
      Z248 = ++rxn_counter;
      Z249 = ++rxn_counter;
      Z250 = ++rxn_counter;
      Z251 = ++rxn_counter;
      Z252 = ++rxn_counter;
      Z253 = ++rxn_counter;
      Z254 = ++rxn_counter;
      Z255 = ++rxn_counter;
      Z256 = ++rxn_counter;
      Z257 = ++rxn_counter;
      Z258 = ++rxn_counter;
      Z259 = ++rxn_counter;
      Z260 = ++rxn_counter;
      Z261 = ++rxn_counter;
      Z262 = ++rxn_counter;
      Z263 = ++rxn_counter;
      Z264 = ++rxn_counter;
      Z265 = ++rxn_counter;
      Z266 = ++rxn_counter;
      Z267 = ++rxn_counter;
      Z268 = ++rxn_counter;
      Z269 = ++rxn_counter;
      Z270 = ++rxn_counter;
      Z271 = ++rxn_counter;
      Z272 = ++rxn_counter;
      Z273 = ++rxn_counter;
      Z274 = ++rxn_counter;
      Z275 = ++rxn_counter;
      Z276 = ++rxn_counter;
      Z277 = ++rxn_counter;
      Z278 = ++rxn_counter;
      Z279 = ++rxn_counter;
      Z280 = ++rxn_counter;
      Z281 = ++rxn_counter;

      // cosmic-ray reactions for metals!
      if (CRX > 0){
      CRX1 = ++rxn_counter;
      CRX2 = ++rxn_counter;
      CRX3 = ++rxn_counter;
      CRX4 = ++rxn_counter;
      CRX5 = ++rxn_counter;
      CRX6 = ++rxn_counter;
      CRX7 = ++rxn_counter;
      CRX8 = ++rxn_counter;
      CRX9 = ++rxn_counter;
      CRX13 = ++rxn_counter;
      CRX14 = ++rxn_counter;
      CRX15 = ++rxn_counter;
      CRX16 = ++rxn_counter;
      CRX17 = ++rxn_counter;
      CRX18 = ++rxn_counter;
      CRX19 = ++rxn_counter;
      CRX20 = ++rxn_counter;
      CRX21 = ++rxn_counter;
      }
   }
   
   nReactions = ++rxn_counter;
}

void setup_species(int ispecies, int UV, int CRX, int water_rates)
{
  int sp_counter = 0;

  H       = sp_counter; 
  Hplus   = ++sp_counter; 
  el      = ++sp_counter; 
  O       = ++sp_counter; 
  OH      = ++sp_counter; 
  H2O     = ++sp_counter; 
  O2      = ++sp_counter; 
  Oplus   = ++sp_counter;
  OHplus  = ++sp_counter; 
  H2Oplus = ++sp_counter; 
  H3Oplus = ++sp_counter; 
  O2plus  = ++sp_counter; 
  Cplus   = ++sp_counter; 
  C       = ++sp_counter; 
  CH      = ++sp_counter;  
  CH2     = ++sp_counter;
  CH3     = ++sp_counter; 
  CH4     = ++sp_counter; 
  CO      = ++sp_counter; 
  COplus  = ++sp_counter; 
  CO2     = ++sp_counter;

  if (ispecies > 1){
     H2m  = ++sp_counter;  
     Hmin = ++sp_counter;

     if (ispecies > 2){ 
        D     = ++sp_counter; 
        Dplus = ++sp_counter; 
        HD    = ++sp_counter;
     }
  }

  if (water_rates == 3){
     CHplus  = ++sp_counter;
     CH2plus = ++sp_counter;
     He      = ++sp_counter;
     Heplus  = ++sp_counter;
     H3plus  = ++sp_counter;
     HCOplus = ++sp_counter;
     CH3plus = ++sp_counter;
     H2plus  = ++sp_counter; 
     HeHplus = ++sp_counter;
     CH4plus = ++sp_counter;
     CH5plus = ++sp_counter;
     O2Hplus = ++sp_counter;
  } 
  
  nSpecies = ++sp_counter;
}

void build_reactions(reaction_t *reactions,int ispecies, int UV, int CRX, int water_rates){
  // water_rates = 1, multispecies = 1
  build_reaction(&reactions[H1],H1,
                 1,1,Hplus,el,
                 1,0,0,0,H,NRXN,NRXN,NRXN);

  build_reaction(&reactions[Z1],Z1,
                 1,1,Oplus,H,
                 1,1,0,0,Hplus,O,NRXN,NRXN);

  build_reaction(&reactions[Z2],Z2,
                 1,1,Hplus,O,
                 1,1,0,0,Oplus,H,NRXN,NRXN);

  build_reaction(&reactions[Z7],Z7,
                 1,1,H2Oplus,el,
                 1,1,0,0,OH,H,NRXN,NRXN);

  build_reaction(&reactions[Z8],Z8,
                 1,1,H3Oplus,el,
                 1,2,0,0,OH,H,NRXN,NRXN);

  build_reaction(&reactions[Z9],Z9,
                 1,1,H3Oplus,el,
                 1,1,0,0,H2O,H,NRXN,NRXN);

  build_reaction(&reactions[Z10],Z10,
                 1,1,O,H,
                 1,0,0,0,OH,NRXN,NRXN,NRXN);

  build_reaction(&reactions[Z13],Z13,
                 2,0,OH,NRXN,
                 1,1,0,0,H2O,O,NRXN,NRXN);

  build_reaction(&reactions[Z14],Z14,
                 1,1,Hplus,OH,
                 1,1,0,0,OHplus,H,NRXN,NRXN);

  build_reaction(&reactions[Z15],Z15,
                 1,1,Hplus,H2O,
                 1,1,0,0,H2Oplus,H,NRXN,NRXN);

  build_reaction(&reactions[Z18],Z18,
                 2,0,O,NRXN,
                 1,0,0,0,O2,NRXN,NRXN,NRXN);

  build_reaction(&reactions[Z19],Z19,
                 1,1,O,OH,
                 1,1,0,0,O2,H,NRXN,NRXN);

  build_reaction(&reactions[Z20],Z20,
                 1,1,H,O2,
                 1,1,0,0,OH,O,NRXN,NRXN);

  build_reaction(&reactions[Z21],Z21,
                 1,1,Hplus,O2,
                 1,1,0,0,O2plus,H,NRXN,NRXN);

  build_reaction(&reactions[Z22],Z22,
                 1,1,O2plus,el,
                 2,0,0,0,O,NRXN,NRXN,NRXN);

  build_reaction(&reactions[Z23],Z23,
                 1,1,O,CH,
                 1,1,0,0,CO,H,NRXN,NRXN);

  build_reaction(&reactions[Z24],Z24,
                 1,1,C,OH,
                 1,1,0,0,CO,H,NRXN,NRXN);

  build_reaction(&reactions[Z25],Z25,
                 1,1,C,O2,
                 1,1,0,0,CO,O,NRXN,NRXN);
  
  build_reaction(&reactions[Z26],Z26,
                 1,1,Cplus,O2,
                 1,1,0,0,Oplus,CO,NRXN,NRXN);

  if (water_rates != 3){
     build_reaction(&reactions[Z27],Z27,
                    1,1,OH,CO,
                    1,1,0,0,CO2,H,NRXN,NRXN);
  }

  build_reaction(&reactions[Z28],Z28,
                 1,1,Cplus,el,
                 1,0,0,0,C,NRXN,NRXN,NRXN);

  build_reaction(&reactions[Z29],Z29,
                 1,1,Cplus,OH,
                 1,1,0,0,COplus,H,NRXN,NRXN);

  build_reaction(&reactions[Z30],Z30,
                 1,1,COplus,H,
                 1,1,0,0,Hplus,CO,NRXN,NRXN);

  build_reaction(&reactions[Z31],Z31,
                 1,1,C,H,
                 1,0,0,0,CH,NRXN,NRXN,NRXN);

  if (ispecies > 1){

    build_reaction(&reactions[H2],H2,
                   1,1,H,el,
                   1,0,0,0,Hmin,NRXN,NRXN,NRXN);

    build_reaction(&reactions[H3],H3,
                   1,1,Hmin,H,
                   1,1,0,0,H2m,el,NRXN,NRXN);
    
    build_reaction(&reactions[H4],H4,
                   1,1,H2m,H,
                   3,0,0,0,H,NRXN,NRXN,NRXN);

    build_reaction(&reactions[H5],H5,
                   3,0,H,NRXN,
                   1,1,0,0,H2m,H,NRXN,NRXN);

    build_reaction(&reactions[H6],H6,
                   2,1,H,H2m,
                   2,0,0,0,H2m,NRXN,NRXN,NRXN);

    build_reaction(&reactions[H7],H7,
                   2,0,H2m,NRXN,
                   2,1,0,0,H,H2m,NRXN,NRXN);

    // HACKY FIX FOR NOW
    build_reaction(&reactions[H8],H8,
                   1,1,H,H,
                   1,0,0,0,H2m,NRXN,NRXN,NRXN);
		   
    /*
    build_reaction(&reactions[H8],H8,
                   1,0,H,NRXN,
                   1,0,0,0,H2m,NRXN,NRXN,NRXN);
		   */
		   

    build_reaction(&reactions[Z3],Z3,
                   1,1,Oplus,H2m,
                   1,1,0,0,OHplus,H,NRXN,NRXN);

    build_reaction(&reactions[Z4],Z4,
                   1,1,OHplus,H2m,
                   1,1,0,0,H2Oplus,H,NRXN,NRXN);

    build_reaction(&reactions[Z5],Z5,
                   1,1,H2Oplus,H2m,
                   1,1,0,0,H3Oplus,H,NRXN,NRXN);

    build_reaction(&reactions[Z6],Z6,
                   1,1,H2Oplus,el,
                   1,1,0,0,O,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z11],Z11,
                   1,1,O,H2m,
                   1,1,0,0,OH,H,NRXN,NRXN);
    
    build_reaction(&reactions[Z12],Z12,
                   1,1,H2m,OH,
                   1,1,0,0,H2O,H,NRXN,NRXN);

    build_reaction(&reactions[Z16],Z16,
                   1,1,H,OH,
                   1,1,0,0,H2m,O,NRXN,NRXN);

    build_reaction(&reactions[Z17],Z17,
                   1,1,H,H2O,
                   1,1,0,0,OH,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z32],Z32,
                   1,1,C,H2m,
                   1,1,0,0,CH,H,NRXN,NRXN);

    build_reaction(&reactions[Z33],Z33,
                   1,1,H,CH,
                   1,1,0,0,C,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z34],Z34,
                   1,1,H2m,CH,
                   1,1,0,0,CH2,H,NRXN,NRXN);

    build_reaction(&reactions[Z35],Z35,
                   1,1,H,CH2,
                   1,1,0,0,CH,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z36],Z36,
                   1,1,H2m,CH2,
                   1,1,0,0,CH3,H,NRXN,NRXN);
    
    build_reaction(&reactions[Z37],Z37,
                   1,1,H,CH3,
                   1,1,0,0,CH2,H2m,NRXN,NRXN);
    
    build_reaction(&reactions[Z38],Z38,
                   1,1,H2m,CH3,
                   1,1,0,0,CH4,H,NRXN,NRXN);

    build_reaction(&reactions[Z39],Z39,
                   1,1,H,CH4,
                   1,1,0,0,H2m,CH3,NRXN,NRXN);

    build_reaction(&reactions[Z40],Z40,
                   1,1,H2m,C,
                   1,1,0,0,CH2,NRXN,NRXN,NRXN);

    if (ispecies > 2){
      build_reaction(&reactions[D1],D1,
                     1,1,D,Hplus,
                     1,1,0,0,Dplus,H,NRXN,NRXN);
      
      build_reaction(&reactions[D2],D2,
                     1,1,Dplus,H,
                     1,1,0,0,D,Hplus,NRXN,NRXN);

      build_reaction(&reactions[D3],D3,
                     1,1,D,H2m,
                     1,1,0,0,H,HD,NRXN,NRXN);

      build_reaction(&reactions[D4],D4,
                     1,1,Dplus,H2m,
                     1,1,0,0,Hplus,HD,NRXN,NRXN);

      build_reaction(&reactions[D5],D5,
                     1,1,HD,H,
                     1,1,0,0,H2m,D,NRXN,NRXN);
      
      build_reaction(&reactions[D6],D6,
                     1,1,HD,Hplus,
                     1,1,0,0,H2m,Dplus,NRXN,NRXN);
    }
  }

  // UV Reactions
  if (UV){
  build_reaction(&reactions[UV1],UV1,
                 1,0,OH,NRXN,
                 1,1,0,0,O,H,NRXN,NRXN);

  build_reaction(&reactions[UV2],UV2,
                 1,0,H2O,NRXN,
                 1,1,0,0,O,H2m,NRXN,NRXN);

  build_reaction(&reactions[UV3],UV3,
                 1,0,H2O,NRXN,
                 1,1,0,0,OH,H,NRXN,NRXN);

  build_reaction(&reactions[UV4],UV4,
                 1,0,O2,NRXN,
                 2,0,0,0,O,NRXN,NRXN,NRXN);

  build_reaction(&reactions[UV5],UV5,
                 1,0,CH,NRXN,
                 1,1,0,0,C,H,NRXN,NRXN);

  build_reaction(&reactions[UV6],UV6,
                 1,0,CO,NRXN,
                 1,1,0,0,C,O,NRXN,NRXN);

  build_reaction(&reactions[UV7],UV7,
                 1,0,C,NRXN,
                 1,1,0,0,Cplus,el,NRXN,NRXN);

  build_reaction(&reactions[UV8],UV8,
                 1,0,O,NRXN,
                 1,1,0,0,Oplus,el,NRXN,NRXN);

  build_reaction(&reactions[UV9],UV9,
                 1,0,Hmin,NRXN,
                 1,1,0,0,H,el,NRXN,NRXN);
  
  build_reaction(&reactions[UV10],UV10,
                 1,0,CH,NRXN,
                 1,1,0,0,CHplus,el,NRXN,NRXN);

  build_reaction(&reactions[UV11],UV11,
		  1,0,CH2,NRXN,
		  1,1,0,0,CH,H,NRXN,NRXN);

  build_reaction(&reactions[UV12],UV12,
		  1,0,CH2plus,NRXN,
		  1,1,0,0,CH,Hplus,NRXN,NRXN);

  build_reaction(&reactions[UV13],UV13,
		  1,0,CH3,NRXN,
		  1,1,0,0,CH2,H,NRXN,NRXN);

  build_reaction(&reactions[UV14],UV14,
		  1,0,CH4,NRXN,
		  1,1,0,0,CH2,H2m,NRXN,NRXN);
  
  build_reaction(&reactions[UV15],UV15,
		  1,0,CH4plus,NRXN,
		  1,1,0,0,CH2plus,H2m,NRXN,NRXN);
  
  build_reaction(&reactions[UV16],UV16,
		  1,0,CHplus,NRXN,
		  1,1,0,0,C,Hplus,NRXN,NRXN);
  
  build_reaction(&reactions[UV17],UV17,
		  1,0,COplus,NRXN,
		  1,1,0,0,Cplus,O,NRXN,NRXN);
  
  build_reaction(&reactions[UV18],UV18,
		  1,0,H2plus,NRXN,
		  1,1,0,0,Hplus,H,NRXN,NRXN);
  
  build_reaction(&reactions[UV19],UV19,
		  1,0,HCOplus,NRXN,
		  1,1,0,0,COplus,H,NRXN,NRXN);
  
  build_reaction(&reactions[UV20],UV20,
		  1,0,O2plus,NRXN,
		  1,1,0,0,Oplus,O,NRXN,NRXN);
  
  build_reaction(&reactions[UV21],UV21,
		  1,0,OHplus,NRXN,
		  1,1,0,0,Oplus,H,NRXN,NRXN);
  
  build_reaction(&reactions[UV22],UV22,
		  1,0,CH4,NRXN,
		  1,1,0,0,CH4plus,el,NRXN,NRXN);
  
  build_reaction(&reactions[UV23],UV23,
		  1,0,H2O,NRXN,
		  1,1,0,0,H2Oplus,el,NRXN,NRXN);
  
  build_reaction(&reactions[UV24],UV24,
		  1,0,O2,NRXN,
		  1,1,0,0,O2plus,el,NRXN,NRXN);
  
  build_reaction(&reactions[UV25],UV25,
		  1,0,CH2,NRXN,
		  1,1,0,0,CH2plus,el,NRXN,NRXN);
  
  build_reaction(&reactions[UV26],UV26,
		  1,0,CH3,NRXN,
		  1,1,0,0,CH3plus,el,NRXN,NRXN);
  
  build_reaction(&reactions[UV27],UV27,
		  1,0,H2Oplus,NRXN,
		  1,1,0,0,H2plus,O,NRXN,NRXN);
  
  build_reaction(&reactions[UV28],UV28,
		  1,0,H2Oplus,NRXN,
		  1,1,0,0,Oplus,H2m,NRXN,NRXN);
  
  build_reaction(&reactions[UV29],UV29,
		  1,0,H2Oplus,NRXN,
		  1,1,0,0,OHplus,H,NRXN,NRXN);
  
  build_reaction(&reactions[UV30],UV30,
		  1,0,H3Oplus,NRXN,
		  1,1,0,0,Hplus,H2O,NRXN,NRXN);
  
  build_reaction(&reactions[UV31],UV31,
		  1,0,H3Oplus,NRXN,
		  1,1,0,0,H2plus,OH,NRXN,NRXN);
  
  build_reaction(&reactions[UV32],UV32,
		  1,0,H3Oplus,NRXN,
		  1,1,0,0,H2Oplus,H,NRXN,NRXN);
  
  build_reaction(&reactions[UV33],UV33,
		  1,0,H3Oplus,NRXN,
		  1,1,0,0,OHplus,H2m,NRXN,NRXN);
  
  build_reaction(&reactions[UV34],UV34,
		  1,0,CH2plus,NRXN,
		  1,1,0,0,CHplus,H,NRXN,NRXN);
  
  build_reaction(&reactions[UV35],UV35,
		  1,0,CH3plus,NRXN,
		  1,1,0,0,CH2plus,H,NRXN,NRXN);
  
  build_reaction(&reactions[UV36],UV36,
		  1,0,CH3plus,NRXN,
		  1,1,0,0,CHplus,H2m,NRXN,NRXN);
  
  build_reaction(&reactions[UV37],UV37,
		  1,0,CH4plus,NRXN,
		  1,1,0,0,CH2plus,H2m,NRXN,NRXN);
  
  build_reaction(&reactions[UV38],UV38,
		  1,0,CH4plus,NRXN,
		  1,1,0,0,CH3plus,H,NRXN,NRXN);
  
  build_reaction(&reactions[UV39],UV39,
		  1,0,CH5plus,NRXN,
		  1,1,0,0,CH4plus,H,NRXN,NRXN);
  
  build_reaction(&reactions[UV40],UV40,
		  1,0,CH5plus,NRXN,
		  1,1,0,0,CH3plus,H2m,NRXN,NRXN);
  }

  // Bialy reactions
  if(water_rates == 3){
    //Hydrogen reactions
    build_reaction(&reactions[H9],H9,
                   1,1,H2plus,H2m,
                   1,1,0,0,H3plus,H,NRXN,NRXN);
    
    build_reaction(&reactions[H10],H10,
                   1,1,H3plus,el,
                   3,0,0,0,H,NRXN,NRXN,NRXN);

    build_reaction(&reactions[H11],H11,
                   1,1,H3plus,el,
                   1,1,0,0,H2m,H,NRXN,NRXN);

    build_reaction(&reactions[H12],H12,
                   1,1,Heplus,H,
                   1,1,0,0,He,Hplus,NRXN,NRXN);

    build_reaction(&reactions[H13],H13,
                   1,1,Heplus,H2m,
                   1,1,0,0,He,H2plus,NRXN,NRXN);
   
    build_reaction(&reactions[H14],H14,
                   1,1,Heplus,H2m,
                   1,1,1,0,He,Hplus,H,NRXN);

    build_reaction(&reactions[H15],H15,
                   1,1,H2plus,He,
                   1,1,0,0,HeHplus,H,NRXN,NRXN);

    build_reaction(&reactions[H16],H16,
                   1,1,HeHplus,H2m,
                   1,1,0,0,H3plus,He,NRXN,NRXN);

    build_reaction(&reactions[H17],H17,
                   1,1,HeHplus,H,
                   1,1,0,0,H2plus,He,NRXN,NRXN);

    build_reaction(&reactions[H18],H18,
                   1,1,He,Hplus,
                   1,0,0,0,HeHplus,NRXN,NRXN,NRXN);

    //Added H reactions (that are in grackle already but we 
    //want to test just the water network)
    build_reaction(&reactions[H19],H19,
		   1,1,Hmin,Hplus,
		   2,0,0,0,H,NRXN,NRXN,NRXN);

    build_reaction(&reactions[H20],H20,
		    1,1,H,Hplus,
		    1,0,0,0,H2plus,NRXN,NRXN,NRXN);

     build_reaction(&reactions[H21],H21,
                 1,1,H2plus,H,
                 1,1,0,0,H2m,Hplus,NRXN,NRXN);

    build_reaction(&reactions[H22],H22,
		    1,1,H2m,el,
		    2,1,0,0,H,el,NRXN,NRXN);

    build_reaction(&reactions[H23],H23,
		    1,1,Hmin,H2plus,
		    1,1,0,0,H2m,H,NRXN,NRXN);

    build_reaction(&reactions[H24],H24,
		    1,1,Hmin,H3plus,
		    2,1,0,0,H,H2m,NRXN,NRXN);

    build_reaction(&reactions[H25],H25,
		     1,1,Hmin,Heplus,
		     1,1,0,0,H,He,NRXN,NRXN);

    // Metal reactions
    build_reaction(&reactions[Z43],Z43,
                   1,1,CH2,O,
                   1,1,0,0,CO,H2m,NRXN,NRXN);
    
    build_reaction(&reactions[Z44],Z44,
                   1,1,CH2,O,
                   1,2,0,0,CO,H,NRXN,NRXN);

    build_reaction(&reactions[Z45],Z45,
                   1,1,Heplus,CO,
                   1,1,1,0,He,Cplus,O,NRXN);
		   
    build_reaction(&reactions[Z46],Z46,
                   1,1,O2,Cplus,
                   1,1,0,0,COplus,O,NRXN,NRXN);

    build_reaction(&reactions[Z48],Z48,
                   1,1,H3plus,O,
                   1,1,0,0,OHplus,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z49],Z49,
                   1,1,H3Oplus,el,
                   1,1,0,0,OH,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z50],Z50,
                   1,1,H3Oplus,el,
                   1,1,1,0,O,H2m,H,NRXN);

    build_reaction(&reactions[Z51],Z51,
                   1,1,H2Oplus,el,
                   1,2,0,0,O,H,NRXN,NRXN);

    build_reaction(&reactions[Z52],Z52,
                   1,1,OHplus,el,
                   1,1,0,0,O,H,NRXN,NRXN);

    build_reaction(&reactions[Z53],Z53,
                   1,1,H3plus,CO,
                   1,1,0,0,HCOplus,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z54],Z54,
                   1,1,H3Oplus,el,
                   1,2,0,0,OH,H,NRXN,NRXN);

    build_reaction(&reactions[Z55],Z55,
                   1,1,OH,C,
                   1,1,0,0,CO,H,NRXN,NRXN);
                   
    build_reaction(&reactions[Z56],Z56,
                   1,1,OH,Cplus,
                   1,1,0,0,COplus,H,NRXN,NRXN);

    build_reaction(&reactions[Z57],Z57,
                   1,1,OH,Hplus,
                   1,1,0,0,OHplus,H,NRXN,NRXN);
    // END OF DOUBLE COUNTED REACTIONS


    build_reaction(&reactions[Z58],Z58,
                   1,1,H2O,Cplus,
                   1,1,0,0,HCOplus,H,NRXN,NRXN);

    build_reaction(&reactions[Z59],Z59,
                   1,1,HCOplus,H2O,
                   1,1,0,0,H3Oplus,CO,NRXN,NRXN);

    build_reaction(&reactions[Z60],Z60,
                   1,1,H2O,H3plus,
                   1,1,0,0,H3Oplus,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z61],Z61,
                   1,1,H3plus,C,
                   1,1,0,0,CHplus,H2m,NRXN,NRXN);
    
    build_reaction(&reactions[Z62],Z62,
                   1,1,CHplus,H2m,
                   1,1,0,0,CH2plus,H,NRXN,NRXN);

    build_reaction(&reactions[Z63],Z63,
                   1,1,Cplus,H2m,
                   1,0,0,0,CH2plus,NRXN,NRXN,NRXN);

    build_reaction(&reactions[Z64],Z64,
                   1,1,CH2plus,H2m,
                   1,1,0,0,CH3plus,H,NRXN,NRXN);

    build_reaction(&reactions[Z65],Z65,
                   1,1,CH3plus,el,
                   1,1,0,0,CH2,H,NRXN,NRXN);

    build_reaction(&reactions[Z66],Z66,
                   1,1,CH3plus,el,
                   1,1,0,0,CH,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z67],Z67,
                   1,1,CH3plus,el,
                   1,2,0,0,CH,H,NRXN,NRXN);
    
    build_reaction(&reactions[Z68],Z68,
                   1,1,CHplus,el,
                   1,1,0,0,C,H,NRXN,NRXN);

    build_reaction(&reactions[Z69],Z69,
                   1,1,CH2plus,el,
                   1,1,0,0,CH,H,NRXN,NRXN);

    build_reaction(&reactions[Z70],Z70,
                   1,1,CH2plus,el,
                   1,2,0,0,C,H,NRXN,NRXN);

    build_reaction(&reactions[Z71],Z71,
                   1,1,CH2plus,el,
                   1,1,0,0,C,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z72],Z72,
                   1,1,CHplus,H,
                   1,1,0,0,Cplus,H2m,NRXN,NRXN);

    // DOUBLE COUNTED
    build_reaction(&reactions[Z73],Z73,
                   1,1,CH,H,
                   1,1,0,0,C,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z74],Z74,
                   1,1,COplus,H2m,
                   1,1,0,0,HCOplus,H,NRXN,NRXN);

    build_reaction(&reactions[Z75],Z75,
                   1,1,HCOplus,el,
                   1,1,0,0,CO,H,NRXN,NRXN);

    build_reaction(&reactions[Z76],Z76,
                   1,1,O,CH,
                   1,1,0,0,HCOplus,el,NRXN,NRXN);

    build_reaction(&reactions[Z77],Z77,
		   1,1,CO,H2plus,
		   1,1,0,0,COplus,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z78],Z78,
		    1,1,Oplus,CO,
		    1,1,0,0,COplus,O,NRXN,NRXN);

    build_reaction(&reactions[Z79],Z79,
		    1,1,H2plus,CO,
		    1,1,0,0,HCOplus,H,NRXN,NRXN);

    build_reaction(&reactions[Z80],Z80,
		    1,1,H2Oplus,CO,
		    1,1,0,0,HCOplus,OH,NRXN,NRXN);

    build_reaction(&reactions[Z81],Z81,
		    1,1,OHplus,CO,
		    1,1,0,0,HCOplus,O,NRXN,NRXN);

    build_reaction(&reactions[Z82],Z82,
		    1,1,H,CO,
		    1,1,0,0,OH,C,NRXN,NRXN);

    build_reaction(&reactions[Z83],Z83,
		    1,1,Hmin,Cplus,
		    1,1,0,0,C,H,NRXN,NRXN);

    build_reaction(&reactions[Z84],Z84,
		    1,1,Hmin,C,
		    1,1,0,0,CH,el,NRXN,NRXN);

    build_reaction(&reactions[Z85],Z85,
		    1,1,Hmin,CH2,
		    1,1,0,0,CH3,el,NRXN,NRXN);
 
    build_reaction(&reactions[Z86],Z86,
		    1,1,Hmin,CH3,
		    1,1,0,0,CH4,el,NRXN,NRXN);
 
    build_reaction(&reactions[Z87],Z87,
		    1,1,Hmin,CH,
		    1,1,0,0,CH2,el,NRXN,NRXN);
 
    build_reaction(&reactions[Z88],Z88,
		    1,1,Hmin,O,
		    1,1,0,0,OH,el,NRXN,NRXN);
 
    build_reaction(&reactions[Z89],Z89,
		    1,1,Hmin,OH,
		    1,1,0,0,H2O,el,NRXN,NRXN);
  
    build_reaction(&reactions[Z90],Z90,
		    1,1,H2m,CH,
		    1,1,1,0,C,H2m,H,NRXN);
 
    build_reaction(&reactions[Z91],Z91,
		    1,1,H2m,H2O,
		    1,1,1,0,OH,H2m,H,NRXN);
 
    build_reaction(&reactions[Z92],Z92,
		    1,1,H2m,O2,
		    2,1,0,0,O,H2m,NRXN,NRXN);
 
    build_reaction(&reactions[Z93],Z93,
		    1,1,H2m,OH,
		    1,1,1,0,O,H2m,H,NRXN);
 
    build_reaction(&reactions[Z94],Z94,
		    1,1,H,CH,
		    2,1,0,0,H,C,NRXN,NRXN);

    build_reaction(&reactions[Z95],Z95,
		    1,1,H,H2O,
		    2,1,0,0,H,OH,NRXN,NRXN);

    build_reaction(&reactions[Z96],Z96,
		    1,1,H,O2,
		    2,1,0,0,O,H,NRXN,NRXN);

    build_reaction(&reactions[Z97],Z97,
		    1,1,H,OH,
		    2,1,0,0,H,O,NRXN,NRXN);

    build_reaction(&reactions[Z98],Z98,
		    1,1,Cplus,CH2,
		    1,1,0,0,CH2plus,C,NRXN,NRXN);

    build_reaction(&reactions[Z99],Z99,
		    1,1,Cplus,CH,
		    1,1,0,0,CHplus,C,NRXN,NRXN);

    build_reaction(&reactions[Z100],Z100,
		    1,1,C,COplus,
		    1,1,0,0,CO,Cplus,NRXN,NRXN);

    build_reaction(&reactions[Z101],Z101,
		    1,1,C,O2plus,
		    1,1,0,0,O2,Cplus,NRXN,NRXN);

    build_reaction(&reactions[Z102],Z102,
		    1,1,CH2,COplus,
		    1,1,0,0,CO,CH2plus,NRXN,NRXN);

    build_reaction(&reactions[Z103],Z103,
		    1,1,CH2,H2Oplus,
		    1,1,0,0,H2O,CH2plus,NRXN,NRXN);

    build_reaction(&reactions[Z104],Z104,
		    1,1,CH2,Oplus,
		    1,1,0,0,O,CH2plus,NRXN,NRXN);

    build_reaction(&reactions[Z105],Z105,
		    1,1,CH2,O2plus,
		    1,1,0,0,O2,CH2plus,NRXN,NRXN);

    build_reaction(&reactions[Z106],Z106,
		    1,1,CH2,OHplus,
		    1,1,0,0,OH,CH2plus,NRXN,NRXN);

    build_reaction(&reactions[Z107],Z107,
		    1,1,CH4plus,O2,
		    1,1,0,0,O2plus,CH4,NRXN,NRXN);

    build_reaction(&reactions[Z108],Z108,
		    1,1,CH4,COplus,
		    1,1,0,0,CO,CH4plus,NRXN,NRXN);

    build_reaction(&reactions[Z109],Z109,
		    1,1,CH,COplus,
		    1,1,0,0,CO,CHplus,NRXN,NRXN);

    build_reaction(&reactions[Z110],Z110,
		    1,1,CH,H2Oplus,
		    1,1,0,0,H2O,CHplus,NRXN,NRXN);

    build_reaction(&reactions[Z111],Z111,
		    1,1,CH,Oplus,
		    1,1,0,0,O,CHplus,NRXN,NRXN);

    build_reaction(&reactions[Z112],Z112,
		    1,1,CH,O2plus,
		    1,1,0,0,O2,CHplus,NRXN,NRXN);

    build_reaction(&reactions[Z113],Z113,
		    1,1,CH,OHplus,
		    1,1,0,0,OH,CHplus,NRXN,NRXN);

    build_reaction(&reactions[Z114],Z114,
		    1,1,COplus,O2,
		    1,1,0,0,O2plus,CO,NRXN,NRXN);

    build_reaction(&reactions[Z115],Z115,
		    1,1,Hplus,CH2,
		    1,1,0,0,CH2plus,H,NRXN,NRXN);

    build_reaction(&reactions[Z116],Z116,
		    1,1,Hplus,CH3,
		    1,1,0,0,CH3plus,H,NRXN,NRXN);

    build_reaction(&reactions[Z117],Z117,
		    1,1,Hplus,CH4,
		    1,1,0,0,CH4plus,H,NRXN,NRXN);

    build_reaction(&reactions[Z118],Z118,
		    1,1,Hplus,CH,
		    1,1,0,0,CHplus,H,NRXN,NRXN);

    build_reaction(&reactions[Z119],Z119,
		    1,1,H2plus,CH2,
		    1,1,0,0,CH2plus,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z120],Z120,
		    1,1,H2plus,CH4,
		    1,1,0,0,CH4plus,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z121],Z121,
		    1,1,H2plus,CH,
		    1,1,0,0,CHplus,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z122],Z122,
		    1,1,H2plus,H2O,
		    1,1,0,0,H2Oplus,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z123],Z123,
		    1,1,H2plus,O2,
		    1,1,0,0,O2plus,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z124],Z124,
		    1,1,H2plus,OH,
		    1,1,0,0,OHplus,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z125],Z125,
		    1,1,H2Oplus,O2,
		    1,1,0,0,O2plus,H2O,NRXN,NRXN);

    build_reaction(&reactions[Z126],Z126,
		    1,1,H2O,COplus,
		    1,1,0,0,CO,H2Oplus,NRXN,NRXN);

    build_reaction(&reactions[Z127],Z127,
		    1,1,Heplus,C,
		    1,1,0,0,Cplus,He,NRXN,NRXN);

    build_reaction(&reactions[Z128],Z128,
		    1,1,Heplus,CH4,
		    1,1,0,0,CH4plus,He,NRXN,NRXN);

    build_reaction(&reactions[Z129],Z129,
		    1,1,Heplus,CH,
		    1,1,0,0,CHplus,He,NRXN,NRXN);

    build_reaction(&reactions[Z130],Z130,
		    1,1,Heplus,H2O,
		    1,1,0,0,H2Oplus,He,NRXN,NRXN);

    build_reaction(&reactions[Z131],Z131,
		    1,1,Heplus,O2,
		    1,1,0,0,O2plus,He,NRXN,NRXN);

    build_reaction(&reactions[Z132],Z132,
		    1,1,Oplus,CH4,
		    1,1,0,0,CH4plus,O,NRXN,NRXN);

    build_reaction(&reactions[Z133],Z133,
		    1,1,Oplus,H2O,
		    1,1,0,0,H2Oplus,O,NRXN,NRXN);

    build_reaction(&reactions[Z134],Z134,
		    1,1,Oplus,O2,
		    1,1,0,0,O2plus,O,NRXN,NRXN);

    build_reaction(&reactions[Z135],Z135,
		    1,1,Oplus,OH,
		    1,1,0,0,OHplus,O,NRXN,NRXN);

    build_reaction(&reactions[Z136],Z136,
		    1,1,O,COplus,
		    1,1,0,0,CO,Oplus,NRXN,NRXN);

    build_reaction(&reactions[Z137],Z137,
		    1,1,OHplus,H2O,
		    1,1,0,0,H2Oplus,OH,NRXN,NRXN);

    build_reaction(&reactions[Z138],Z138,
		    1,1,OHplus,O2,
		    1,1,0,0,O2plus,OH,NRXN,NRXN);

    build_reaction(&reactions[Z139],Z139,
		    1,1,OH,COplus,
		    1,1,0,0,CO,OHplus,NRXN,NRXN);

    build_reaction(&reactions[Z140],Z140,
		    1,1,CH4plus,el,
		    2,1,0,0,H,CH2,NRXN,NRXN);

    build_reaction(&reactions[Z141],Z141,
		    1,1,CH4plus,el,
		    1,1,0,0,CH3,H,NRXN,NRXN);

    build_reaction(&reactions[Z142],Z142,
		    1,1,CH5plus,el,
		    1,1,1,0,CH2,H2m,H,NRXN);

    build_reaction(&reactions[Z143],Z143,
		    1,1,CH5plus,el,
		    1,1,0,0,CH3,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z144],Z144,
		    1,1,CH5plus,el,
		    2,1,0,0,H,CH3,NRXN,NRXN);

    build_reaction(&reactions[Z145],Z145,
		    1,1,CH5plus,el,
		    1,1,0,0,CH4,H,NRXN,NRXN);

    build_reaction(&reactions[Z146],Z146,
		    1,1,CH5plus,el,
		    2,1,0,0,H2m,CH,NRXN,NRXN);

    build_reaction(&reactions[Z147],Z147,
		    1,1,COplus,el,
		    1,1,0,0,O,C,NRXN,NRXN);

    build_reaction(&reactions[Z148],Z148,
		    1,1,H2plus,el,
		    2,0,0,0,H,NRXN,NRXN,NRXN);

    build_reaction(&reactions[Z149],Z149,
		    1,1,HeHplus,el,
		    1,1,0,0,He,H,NRXN,NRXN);

    build_reaction(&reactions[Z150],Z150,
		    1,1,O2Hplus,el,
		    1,1,0,0,O2,H,NRXN,NRXN);

    build_reaction(&reactions[Z151],Z151,
		    1,1,C,CH5plus,
		    1,1,0,0,CH4,CHplus,NRXN,NRXN);

    build_reaction(&reactions[Z152],Z152,
		    1,1,C,H2Oplus,
		    1,1,0,0,OH,CHplus,NRXN,NRXN);

    build_reaction(&reactions[Z153],Z153,
		    1,1,C,H3Oplus,
		    1,1,0,0,HCOplus,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z154],Z154,
		    1,1,C,HCOplus,
		    1,1,0,0,CO,CHplus,NRXN,NRXN);

    build_reaction(&reactions[Z155],Z155,
		    1,1,C,O2plus,
		    1,1,0,0,COplus,O,NRXN,NRXN);

    build_reaction(&reactions[Z156],Z156,
		    1,1,C,O2Hplus,
		    1,1,0,0,O2,CHplus,NRXN,NRXN);

    build_reaction(&reactions[Z157],Z157,
		    1,1,C,OHplus,
		    1,1,0,0,O,CHplus,NRXN,NRXN);

    build_reaction(&reactions[Z158],Z158,
		    1,1,CHplus,H2O,
		    1,1,0,0,H3Oplus,C,NRXN,NRXN);

    build_reaction(&reactions[Z159],Z159,
		    1,1,CHplus,H2O,
		    1,1,0,0,HCOplus,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z160],Z160,
		    1,1,CHplus,O2,
		    1,1,0,0,COplus,OH,NRXN,NRXN);

    build_reaction(&reactions[Z161],Z161,
		    1,1,CHplus,O2,
		    1,1,0,0,HCOplus,O,NRXN,NRXN);

    build_reaction(&reactions[Z162],Z162,
		    1,1,CHplus,O,
		    1,1,0,0,COplus,H,NRXN,NRXN);

    build_reaction(&reactions[Z163],Z163,
		    1,1,CHplus,OH,
		    1,1,0,0,COplus,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z164],Z164,
		    1,1,CH2plus,O2,
		    1,1,0,0,HCOplus,OH,NRXN,NRXN);

    build_reaction(&reactions[Z165],Z165,
		    1,1,CH2plus,O,
		    1,1,0,0,HCOplus,H,NRXN,NRXN);

    build_reaction(&reactions[Z166],Z166,
		    1,1,CH2,CH5plus,
		    1,1,0,0,CH4,CH3plus,NRXN,NRXN);

    build_reaction(&reactions[Z167],Z167,
		    1,1,CH2,COplus,
		    1,1,0,0,HCOplus,CH,NRXN,NRXN);

    build_reaction(&reactions[Z168],Z168,
		    1,1,CH2,H2Oplus,
		    1,1,0,0,OH,CH3plus,NRXN,NRXN);

    build_reaction(&reactions[Z169],Z169,
		    1,1,CH2,H3Oplus,
		    1,1,0,0,H2O,CH3plus,NRXN,NRXN);

    build_reaction(&reactions[Z170],Z170,
		    1,1,CH2,HCOplus,
		    1,1,0,0,CO,CH3plus,NRXN,NRXN);

    build_reaction(&reactions[Z171],Z171,
		    1,1,CH2,O2Hplus,
		    1,1,0,0,O2,CH3plus,NRXN,NRXN);

    build_reaction(&reactions[Z172],Z172,
		    1,1,CH2,OHplus,
		    1,1,0,0,O,CH3plus,NRXN,NRXN);

    build_reaction(&reactions[Z173],Z173,
		    1,1,CH3plus,O,
		    1,1,0,0,HCOplus,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z174],Z174,
		    1,1,CH4plus,CH4,
		    1,1,0,0,CH5plus,CH3,NRXN,NRXN);

    build_reaction(&reactions[Z175],Z175,
		    1,1,CH4plus,CO,
		    1,1,0,0,HCOplus,CH3,NRXN,NRXN);

    build_reaction(&reactions[Z176],Z176,
		    1,1,CH4plus,H2O,
		    1,1,0,0,H3Oplus,CH3,NRXN,NRXN);

    build_reaction(&reactions[Z177],Z177,
		    1,1,CH4,COplus,
		    1,1,0,0,HCOplus,CH3,NRXN,NRXN);

    build_reaction(&reactions[Z178],Z178,
		    1,1,CH4,H2Oplus,
		    1,1,0,0,H3Oplus,CH3,NRXN,NRXN);

    build_reaction(&reactions[Z179],Z179,
		    1,1,CH4,OHplus,
		    1,1,0,0,CH5plus,O,NRXN,NRXN);

    build_reaction(&reactions[Z180],Z180,
		    1,1,CH4,OHplus,
		    1,1,0,0,H3Oplus,CH2,NRXN,NRXN);

    build_reaction(&reactions[Z181],Z181,
		    1,1,CH5plus,CO,
		    1,1,0,0,HCOplus,CH4,NRXN,NRXN);

    build_reaction(&reactions[Z182],Z182,
		    1,1,CH5plus,H2O,
		    1,1,0,0,H3Oplus,CH4,NRXN,NRXN);

    build_reaction(&reactions[Z183],Z183,
		    1,1,CH,CH5plus,
		    1,1,0,0,CH4,CH2plus,NRXN,NRXN);

    build_reaction(&reactions[Z184],Z184,
		    1,1,CH,COplus,
		    1,1,0,0,HCOplus,C,NRXN,NRXN);

    build_reaction(&reactions[Z185],Z185,
		    1,1,CH,H2Oplus,
		    1,1,0,0,OH,CH2plus,NRXN,NRXN);

    build_reaction(&reactions[Z186],Z186,
		    1,1,CH,H3Oplus,
		    1,1,0,0,H2O,CH2plus,NRXN,NRXN);

    build_reaction(&reactions[Z187],Z187,
		    1,1,CH,HCOplus,
		    1,1,0,0,CO,CH2plus,NRXN,NRXN);

    build_reaction(&reactions[Z188],Z188,
		    1,1,CH,Oplus,
		    1,1,0,0,COplus,H,NRXN,NRXN);

    build_reaction(&reactions[Z189],Z189,
		    1,1,CH,O2plus,
		    1,1,0,0,HCOplus,O,NRXN,NRXN);

    build_reaction(&reactions[Z190],Z190,
		    1,1,CH,O2Hplus,
		    1,1,0,0,O2,CH2plus,NRXN,NRXN);

    build_reaction(&reactions[Z191],Z191,
		    1,1,CH,OHplus,
		    1,1,0,0,O,CH2plus,NRXN,NRXN);

    build_reaction(&reactions[Z192],Z192,
		    1,1,CO,O2Hplus,
		    1,1,0,0,O2,HCOplus,NRXN,NRXN);

    build_reaction(&reactions[Z193],Z193,
		    1,1,Hplus,CH2,
		    1,1,0,0,CHplus,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z194],Z194,
		    1,1,Hplus,CH4,
		    1,1,0,0,CH3plus,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z195],Z195,
		    1,1,H2plus,C,
		    1,1,0,0,CHplus,H,NRXN,NRXN);

    build_reaction(&reactions[Z196],Z196,
		    1,1,H2plus,CH2,
		    1,1,0,0,CH3plus,H,NRXN,NRXN);

    build_reaction(&reactions[Z197],Z197,
		    1,1,H2plus,CH4,
		    1,1,1,0,CH3plus,H2m,H,NRXN);

    build_reaction(&reactions[Z198],Z198,
		    1,1,H2plus,CH4,
		    1,1,0,0,CH5plus,H,NRXN,NRXN);

    build_reaction(&reactions[Z199],Z199,
		    1,1,H2plus,CH,
		    1,1,0,0,CH2plus,H,NRXN,NRXN);

    build_reaction(&reactions[Z200],Z200,
		    1,1,H2plus,H2O,
		    1,1,0,0,H3Oplus,H,NRXN,NRXN);

    build_reaction(&reactions[Z201],Z201,
		    1,1,H2plus,O2,
		    1,1,0,0,O2Hplus,H,NRXN,NRXN);

    build_reaction(&reactions[Z202],Z202,
		    1,1,H2plus,O,
		    1,1,0,0,OHplus,H,NRXN,NRXN);

    build_reaction(&reactions[Z203],Z203,
		    1,1,H2plus,OH,
		    1,1,0,0,H2Oplus,H,NRXN,NRXN);

    build_reaction(&reactions[Z204],Z204,
		    1,1,H2m,Cplus,
		    1,1,0,0,CHplus,H,NRXN,NRXN);

    build_reaction(&reactions[Z205],Z205,
		    1,1,H2m,CH4plus,
		    1,1,0,0,CH5plus,H,NRXN,NRXN);

    build_reaction(&reactions[Z206],Z206,
		    1,1,H2m,O2Hplus,
		    1,1,0,0,O2,H3plus,NRXN,NRXN);

    build_reaction(&reactions[Z207],Z207,
		    1,1,H2Oplus,H2O,
		    1,1,0,0,H3Oplus,OH,NRXN,NRXN);

    build_reaction(&reactions[Z208],Z208,
		    1,1,H2O,COplus,
		    1,1,0,0,HCOplus,OH,NRXN,NRXN);

    build_reaction(&reactions[Z209],Z209,
		    1,1,H2O,O2Hplus,
		    1,1,0,0,O2,H3Oplus,NRXN,NRXN);

    build_reaction(&reactions[Z210],Z210,
		    1,1,H3plus,CH2,
		    1,1,0,0,CH3plus,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z211],Z211,
		    1,1,H3plus,CH3,
		    1,1,0,0,CH4plus,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z212],Z212,
		    1,1,H3plus,CH4,
		    1,1,0,0,CH5plus,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z213],Z213,
		   1,1,H3plus,CH,
		   1,1,0,0,CH2plus,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z214],Z214,
		    1,1,H3plus,O2,
		    1,1,0,0,O2Hplus,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z215],Z215,
		    1,1,H3plus,O,
		    1,1,0,0,H2Oplus,H,NRXN,NRXN);

    build_reaction(&reactions[Z216],Z216,
		    1,1,H3plus,OH,
		    1,1,0,0,H2Oplus,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z217],Z217,
		    1,1,H,CH2plus,
		    1,1,0,0,CHplus,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z218],Z218,
		    1,1,H,CH3plus,
		    1,1,0,0,CH2plus,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z219],Z219,
		    1,1,H,CH4plus,
		    1,1,0,0,CH3plus,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z220],Z220,
		    1,1,H,CH5plus,
		    1,1,0,0,CH4plus,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z221],Z221,
		    1,1,Heplus,CH2,
		    1,1,1,0,Cplus,He,H2m,NRXN);

    build_reaction(&reactions[Z222],Z222,
		    1,1,Heplus,CH2,
		    1,1,1,0,CHplus,He,H,NRXN);

    build_reaction(&reactions[Z223],Z223,
		    1,1,Heplus,CH3,
		    1,1,1,0,CHplus,He,H2m,NRXN);

    build_reaction(&reactions[Z224],Z224,
                    1,1,Heplus,CH4,
                    1,1,1,1,CHplus,He,H2m,H);

    build_reaction(&reactions[Z225],Z225,
		    1,1,Heplus,CH4,
		    1,1,1,0,CH2plus,He,H2m,NRXN);

    build_reaction(&reactions[Z226],Z226,
		    1,1,Heplus,CH4,
		    1,1,1,0,CH3plus,He,H,NRXN);

    build_reaction(&reactions[Z227],Z227,
		    1,1,Heplus,CH4,
		    1,1,1,0,CH3,He,H,NRXN);

    build_reaction(&reactions[Z228],Z228,
		    1,1,Heplus,CH,
		    1,1,1,0,Cplus,He,H,NRXN);

    build_reaction(&reactions[Z229],Z229,
		    1,1,Heplus,H2O,
		    1,1,1,0,OHplus,He,H,NRXN);

    build_reaction(&reactions[Z230],Z230,
		    1,1,Heplus,H2O,
		    1,1,1,0,OH,He,Hplus,NRXN);

    build_reaction(&reactions[Z231],Z231,
		    1,1,Heplus,O2,
		    1,1,1,0,Oplus,O,He,NRXN);

    build_reaction(&reactions[Z232],Z232,
		    1,1,Heplus,OH,
		    1,1,1,0,Oplus,He,H,NRXN);

    build_reaction(&reactions[Z233],Z233,
		    1,1,Oplus,CH4,
		    1,1,0,0,OH,CH3plus,NRXN,NRXN);

    build_reaction(&reactions[Z234],Z234,
		    1,1,Oplus,OH,
		    1,1,0,0,O2plus,H,NRXN,NRXN);

    build_reaction(&reactions[Z235],Z235,
		    1,1,O,CH4plus,
		    1,1,0,0,OH,CH3plus,NRXN,NRXN);

    build_reaction(&reactions[Z236],Z236,
		    1,1,O,CH5plus,
		    1,1,0,0,H3Oplus,CH2,NRXN,NRXN);

    build_reaction(&reactions[Z237],Z237,
		    1,1,O,H2Oplus,
		    1,1,0,0,O2plus,H2m,NRXN,NRXN);

    build_reaction(&reactions[Z238],Z238,
		    1,1,O,O2Hplus,
		    1,1,0,0,O2,OHplus,NRXN,NRXN);

    build_reaction(&reactions[Z239],Z239,
		    1,1,O,OHplus,
		    1,1,0,0,O2plus,H,NRXN,NRXN);

    build_reaction(&reactions[Z240],Z240,
		    1,1,OHplus,H2O,
		    1,1,0,0,H3Oplus,O,NRXN,NRXN);

    build_reaction(&reactions[Z241],Z241,
		    1,1,OHplus,OH,
		    1,1,0,0,H2Oplus,O,NRXN,NRXN);

    build_reaction(&reactions[Z242],Z242,
		    1,1,OH,CH5plus,
		    1,1,0,0,H2Oplus,CH4,NRXN,NRXN);

    build_reaction(&reactions[Z243],Z243,
		    1,1,OH,COplus,
		    1,1,0,0,HCOplus,O,NRXN,NRXN);

    build_reaction(&reactions[Z244],Z244,
		    1,1,OH,H2Oplus,
		    1,1,0,0,H3Oplus,O,NRXN,NRXN);

    build_reaction(&reactions[Z245],Z245,
		    1,1,OH,HCOplus,
		    1,1,0,0,CO,H2Oplus,NRXN,NRXN);

    build_reaction(&reactions[Z246],Z246,
		    1,1,OH,O2Hplus,
		    1,1,0,0,O2,H2Oplus,NRXN,NRXN);

    build_reaction(&reactions[Z247],Z247,
		    1,1,Hmin,CH3plus,
		    1,1,0,0,H,CH3,NRXN,NRXN);

    build_reaction(&reactions[Z248],Z248,
		    1,1,Hmin,H3Oplus,
		    2,1,0,0,H,H2O,NRXN,NRXN);

    build_reaction(&reactions[Z249],Z249,
		    1,1,Hmin,HCOplus,
		    2,1,0,0,H,CO,NRXN,NRXN);

    build_reaction(&reactions[Z250],Z250,
		    1,1,Hmin,Oplus,
		    1,1,0,0,H,O,NRXN,NRXN);

    build_reaction(&reactions[Z251],Z251,
		    1,1,C,CH2,
		    1,1,0,0,CH,CH,NRXN,NRXN);

    build_reaction(&reactions[Z252],Z252,
		    1,1,C,OH,
		    1,1,0,0,O,CH,NRXN,NRXN);

    build_reaction(&reactions[Z253],Z253,
		    2,0,CH2,NRXN,
		    1,1,0,0,CH3,CH,NRXN,NRXN);

    build_reaction(&reactions[Z254],Z254,
		    1,1,CH2,CH4,
		    2,0,0,0,CH3,NRXN,NRXN,NRXN);

    build_reaction(&reactions[Z255],Z255,
		    1,1,CH2,O2,
		    1,1,0,0,CO,H2O,NRXN,NRXN);

    build_reaction(&reactions[Z256],Z256,
		    1,1,CH2,O,
		    1,1,0,0,OH,CH,NRXN,NRXN);

    build_reaction(&reactions[Z257],Z257,
		    1,1,CH2,OH,
		    1,1,0,0,H2O,CH,NRXN,NRXN);

    build_reaction(&reactions[Z258],Z258,
		    1,1,CH2,OH,
		    1,1,0,0,O,CH3,NRXN,NRXN);

    build_reaction(&reactions[Z259],Z259,
		    2,0,CH3,NRXN,
		    1,1,0,0,CH4,CH2,NRXN,NRXN);

    build_reaction(&reactions[Z260],Z260,
		    1,1,CH3,H2O,
		    1,1,0,0,OH,CH4,NRXN,NRXN);

    build_reaction(&reactions[Z261],Z261,
		    1,1,CH3,O,
		    1,1,1,0,CO,H2m,H,NRXN);

    build_reaction(&reactions[Z262],Z262,
		    1,1,CH3,OH,
		    1,1,0,0,CH4,O,NRXN,NRXN);

    build_reaction(&reactions[Z263],Z263,
		    1,1,CH3,OH,
		    1,1,0,0,H2O,CH2,NRXN,NRXN);

    build_reaction(&reactions[Z264],Z264,
		    1,1,CH4,OH,
		    1,1,0,0,H2O,CH3,NRXN,NRXN);

    build_reaction(&reactions[Z265],Z265,
		    1,1,CH,O2,
		    1,1,1,0,CO,O,H,NRXN);

    build_reaction(&reactions[Z266],Z266,
		    1,1,CH,O2,
		    1,1,0,0,CO,OH,NRXN,NRXN);

    build_reaction(&reactions[Z267],Z267,
		    1,1,CH,O,
		    1,1,0,0,OH,C,NRXN,NRXN);

    build_reaction(&reactions[Z268],Z268,
		    1,1,H2m,O2,
		    2,0,0,0,OH,NRXN,NRXN,NRXN);

    build_reaction(&reactions[Z269],Z269,
		    1,1,O,CH4,
		    1,1,0,0,OH,CH3,NRXN,NRXN);

    build_reaction(&reactions[Z270],Z270,
		    1,1,O,H2O,
		    2,0,0,0,OH,NRXN,NRXN,NRXN);

    build_reaction(&reactions[Z271],Z271,
		    1,1,Cplus,O,
		    1,0,0,0,COplus,NRXN,NRXN,NRXN);

    build_reaction(&reactions[Z272],Z272,
		    1,1,C,Oplus,
		    1,0,0,0,COplus,NRXN,NRXN,NRXN);

    build_reaction(&reactions[Z273],Z273,
		    1,1,C,O,
		    1,0,0,0,CO,NRXN,NRXN,NRXN);

    build_reaction(&reactions[Z274],Z274,
		    1,1,H2m,CH3plus,
		    1,0,0,0,CH5plus,NRXN,NRXN,NRXN);

    build_reaction(&reactions[Z275],Z275,
		    1,1,H2m,CH,
		    1,0,0,0,CH3,NRXN,NRXN,NRXN);

    build_reaction(&reactions[Z276],Z276,
		    1,1,H,Cplus,
		    1,0,0,0,CHplus,NRXN,NRXN,NRXN);

    build_reaction(&reactions[Z277],Z277,
		    1,1,H,OH,
		    1,0,0,0,H2O,NRXN,NRXN,NRXN);

    build_reaction(&reactions[Z278],Z278,
		    1,1,CH3plus,el,
		    1,0,0,0,CH3,NRXN,NRXN,NRXN);

    build_reaction(&reactions[Z279],Z279,
		    1,1,Heplus,el,
		    1,0,0,0,He,NRXN,NRXN,NRXN);

    build_reaction(&reactions[Z280],Z280,
		    1,1,Oplus,el,
		    1,0,0,0,O,NRXN,NRXN,NRXN);

    build_reaction(&reactions[Z281],Z281,
		    1,1,H2m,Hplus,
		    1,1,0,0,H2plus,H,NRXN,NRXN);

    // CRX Reactions
    if (CRX > 0){
    build_reaction(&reactions[CRX1],CRX1,
                   1,0,C,NRXN,
                   1,1,0,0,Cplus,el,NRXN,NRXN);

    build_reaction(&reactions[CRX2],CRX2,
                   1,0,CO,NRXN,
                   1,1,0,0,C,O,NRXN,NRXN);

    build_reaction(&reactions[CRX3],CRX3,
                   1,0,CHplus,NRXN,
                   1,1,0,0,C,Hplus,NRXN,NRXN);

    build_reaction(&reactions[CRX4],CRX4,
                   1,0,CH,NRXN,
                   1,1,0,0,C,H,NRXN,NRXN);

    build_reaction(&reactions[CRX5],CRX5,
                   1,0,H2O,NRXN,
                   1,1,0,0,OH,H,NRXN,NRXN);

    build_reaction(&reactions[CRX6],CRX6,
                   1,0,OH,NRXN,
                   1,1,0,0,O,H,NRXN,NRXN);

    build_reaction(&reactions[CRX7],CRX7,
                   1,0,CH4,NRXN,
                   1,1,0,0,CH2,H2m,NRXN,NRXN);

    build_reaction(&reactions[CRX8],CRX8,
                   1,0,O2,NRXN,
                   1,1,0,0,O2plus,el,NRXN,NRXN);

    build_reaction(&reactions[CRX9],CRX9,
                   1,0,O2,NRXN,
                   1,1,0,0,O,O,NRXN,NRXN);

    build_reaction(&reactions[CRX13],CRX13,
                   1,0,CH2,NRXN,
                   1,1,0,0,CH2plus,el,NRXN,NRXN);

    build_reaction(&reactions[CRX14],CRX14,
                   1,0,CH2,NRXN,
                   1,1,0,0,CH,H,NRXN,NRXN);

    build_reaction(&reactions[CRX15],CRX15,
                   1,0,CH3,NRXN,
                   1,1,0,0,CH2,H,NRXN,NRXN);

    build_reaction(&reactions[CRX16],CRX16,
                   1,0,CH3,NRXN,
                   1,1,0,0,CH,H2m,NRXN,NRXN);

    build_reaction(&reactions[CRX17],CRX17,
                   1,0,He,NRXN,
                   1,1,0,0,Heplus,el,NRXN,NRXN);

    build_reaction(&reactions[CRX18],CRX18,
		  1,0,H,NRXN, 
		  1,1,0,0,Hplus,el,NRXN,NRXN);

    build_reaction(&reactions[CRX19],CRX19,
		   1,0,H2m,NRXN,
		   1,1,0,0,H2plus,el,NRXN,NRXN);

    build_reaction(&reactions[CRX20],CRX20,
		   1,0,H2m,NRXN,
		   1,1,1,0,Hplus,H,el,NRXN);

    build_reaction(&reactions[CRX21],CRX21,
		    1,0,CH3,NRXN,
		    1,1,0,0,CH3plus,el,NRXN,NRXN);
    }
  }
}
