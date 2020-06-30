#include <stdlib.h>
#include <stdio.h>
#include "jac_rhs_builder.h"

void setup_rxns(int ispecies, int UV, int CRX, int water_rates);
void setup_species(int ispecies, int UV, int CRX, int water_rates);
void build_reactions(reaction_t *reactions, int ispecies, int UV, int CRX, int water_rates);

extern reaction_t *my_reactions;

extern int nSpecies   ;
extern int nReactions ;

extern int H       ;
extern int Hplus   ; 
extern int Hmin    ; 
extern int H2m     ; 
extern int el      ; 
extern int D       ; 
extern int Dplus   ; 
extern int HD      ; 
extern int O       ;
extern int OH      ; 
extern int H2O     ; 
extern int O2      ;
extern int Oplus   ; 
extern int OHplus  ; 
extern int H2Oplus ; 
extern int H3Oplus ; 
extern int O2plus  ; 
extern int Cplus   ; 
extern int C       ; 
extern int CH      ; 
extern int CH2     ; 
extern int CH3     ; 
extern int CH4     ; 
extern int CH4plus ;
extern int CH5plus ;
extern int CO      ; 
extern int COplus  ; 
extern int CO2     ; 
extern int CHplus  ;
extern int CH2plus ;
extern int He      ;
extern int Heplus  ;
extern int H3plus  ;
extern int HCOplus ;
extern int CH3plus ;
extern int H2plus  ;
extern int HeHplus ;
extern int O2Hplus ;

extern int H1  ;
extern int H2  ;
extern int H3  ;
extern int H4  ;
extern int H5  ;
extern int H6  ;
extern int H7  ;
extern int H8  ;
extern int H9  ;
extern int H10 ; 
extern int H11 ;
extern int H12 ; 
extern int H13 ;
extern int H14 ;
extern int H15 ;
extern int H16 ;
extern int H17 ;
extern int H18 ;
extern int H19 ;
extern int H20 ;
extern int H21 ;
extern int H22 ;
extern int H23 ;
extern int H24 ;
extern int H25 ;

extern int D1  ;
extern int D2  ;
extern int D3  ;
extern int D4  ;
extern int D5  ;
extern int D6  ;

extern int Z1  ;
extern int Z2  ;
extern int Z3  ;
extern int Z4  ;
extern int Z5  ;
extern int Z6  ;
extern int Z7  ;
extern int Z8  ;
extern int Z9  ;
extern int Z10 ;
extern int Z11 ;
extern int Z12 ;
extern int Z13 ;
extern int Z14 ;
extern int Z15 ;
extern int Z16 ;
extern int Z17 ;
extern int Z18 ;
extern int Z19 ;
extern int Z20 ;
extern int Z21 ;
extern int Z22 ;
extern int Z23 ;
extern int Z24 ;
extern int Z25 ;
extern int Z26 ;
extern int Z27 ;
extern int Z28 ;
extern int Z29 ;
extern int Z30 ;
extern int Z31 ;
extern int Z32 ;
extern int Z33 ;
extern int Z34 ;
extern int Z35 ;
extern int Z36 ;
extern int Z37 ;
extern int Z38 ;
extern int Z39 ;
extern int Z40 ;
extern int Z43 ;
extern int Z44 ;
extern int Z45 ;
extern int Z46 ;
extern int Z48 ;
extern int Z49 ;
extern int Z50 ;
extern int Z51 ;
extern int Z52 ;
extern int Z53 ;
extern int Z54 ;
extern int Z55 ;
extern int Z56 ;
extern int Z57 ;
extern int Z58 ;
extern int Z59 ;
extern int Z60 ;
extern int Z61 ;
extern int Z62 ;
extern int Z63 ;
extern int Z64 ;
extern int Z65 ;
extern int Z66 ;
extern int Z67 ;
extern int Z68 ;
extern int Z69 ;
extern int Z70 ;
extern int Z71 ;
extern int Z72 ;
extern int Z73 ;
extern int Z74 ;
extern int Z75 ;
extern int Z76 ;
extern int Z77 ;
extern int Z78 ;
extern int Z79 ;
extern int Z80 ;
extern int Z81 ;
extern int Z82 ;
extern int Z83 ;
extern int Z84 ;
extern int Z85 ;
extern int Z86 ;
extern int Z87 ;
extern int Z88 ;
extern int Z89 ;
extern int Z90 ;
extern int Z91 ;
extern int Z92 ;
extern int Z93 ;
extern int Z94 ;
extern int Z95 ;
extern int Z96 ;
extern int Z97 ;
extern int Z98 ;
extern int Z99 ;
extern int Z100 ;
extern int Z101 ;
extern int Z102 ;
extern int Z103 ;
extern int Z104 ;
extern int Z105 ;
extern int Z106 ;
extern int Z107 ;
extern int Z108 ;
extern int Z109 ;
extern int Z110 ;
extern int Z111 ;
extern int Z112 ;
extern int Z113 ;
extern int Z114 ;
extern int Z115 ;
extern int Z116 ;
extern int Z117 ;
extern int Z118 ;
extern int Z119 ;
extern int Z120 ;
extern int Z121 ;
extern int Z122 ;
extern int Z123 ;
extern int Z124 ;
extern int Z125 ;
extern int Z126 ;
extern int Z127 ;
extern int Z128 ;
extern int Z129 ;
extern int Z130 ;
extern int Z131 ;
extern int Z132 ;
extern int Z133 ;
extern int Z134 ;
extern int Z135 ;
extern int Z136 ;
extern int Z137 ;
extern int Z138 ;
extern int Z139 ;
extern int Z140 ;
extern int Z141 ;
extern int Z142 ;
extern int Z143 ;
extern int Z144 ;
extern int Z145 ;
extern int Z146 ;
extern int Z147 ;
extern int Z148 ;
extern int Z149 ;
extern int Z150 ;
extern int Z151 ;
extern int Z152 ;
extern int Z153 ;
extern int Z154 ;
extern int Z155 ;
extern int Z156 ;
extern int Z157 ;
extern int Z158 ;
extern int Z159 ;
extern int Z160 ;
extern int Z161 ;
extern int Z162 ;
extern int Z163 ;
extern int Z164 ;
extern int Z165 ;
extern int Z166 ;
extern int Z167 ;
extern int Z168 ;
extern int Z169 ;
extern int Z170 ;
extern int Z171 ;
extern int Z172 ;
extern int Z173 ;
extern int Z174 ;
extern int Z175 ;
extern int Z176 ;
extern int Z177 ;
extern int Z178 ;
extern int Z179 ;
extern int Z180 ;
extern int Z181 ;
extern int Z182 ;
extern int Z183 ;
extern int Z184 ;
extern int Z185 ;
extern int Z186 ;
extern int Z187 ;
extern int Z188 ;
extern int Z189 ;
extern int Z190 ;
extern int Z191 ;
extern int Z192 ;
extern int Z193 ;
extern int Z194 ;
extern int Z195 ;
extern int Z196 ;
extern int Z197 ;
extern int Z198 ;
extern int Z199 ;
extern int Z200 ;
extern int Z201 ;
extern int Z202 ;
extern int Z203 ;
extern int Z204 ;
extern int Z205 ;
extern int Z206 ;
extern int Z207 ;
extern int Z208 ;
extern int Z209 ;
extern int Z210 ;
extern int Z211 ;
extern int Z212 ;
extern int Z213 ;
extern int Z214 ;
extern int Z215 ;
extern int Z216 ;
extern int Z217 ;
extern int Z218 ;
extern int Z219 ;
extern int Z220 ;
extern int Z221 ;
extern int Z222 ;
extern int Z223 ;
extern int Z224 ;
extern int Z225 ;
extern int Z226 ;
extern int Z227 ;
extern int Z228 ;
extern int Z229 ;
extern int Z230 ;
extern int Z231 ;
extern int Z232 ;
extern int Z233 ;
extern int Z234 ;
extern int Z235 ;
extern int Z236 ;
extern int Z237 ;
extern int Z238 ;
extern int Z239 ;
extern int Z240 ;
extern int Z241 ;
extern int Z242 ;
extern int Z243 ;
extern int Z244 ;
extern int Z245 ;
extern int Z246 ;
extern int Z247 ;
extern int Z248 ;
extern int Z249 ;
extern int Z250 ;
extern int Z251 ;
extern int Z252 ;
extern int Z253 ;
extern int Z254 ;
extern int Z255 ;
extern int Z256 ;
extern int Z257 ;
extern int Z258 ;
extern int Z259 ;
extern int Z260 ;
extern int Z261 ;
extern int Z262 ;
extern int Z263 ;
extern int Z264 ;
extern int Z265 ;
extern int Z266 ;
extern int Z267 ;
extern int Z268 ;
extern int Z269 ;
extern int Z270 ;
extern int Z271 ;
extern int Z272 ;
extern int Z273 ;
extern int Z274 ;
extern int Z275 ;
extern int Z276 ;
extern int Z277 ;
extern int Z278 ;
extern int Z279 ;
extern int Z280 ;
extern int Z281 ;

extern int UV1 ;
extern int UV2 ;
extern int UV3 ;
extern int UV4 ;
extern int UV5 ;
extern int UV6 ;
extern int UV7 ; 
extern int UV8 ;
extern int UV9 ;
extern int UV10;
extern int UV11;
extern int UV12;
extern int UV13;
extern int UV14;
extern int UV15;
extern int UV16;
extern int UV17;
extern int UV18;
extern int UV19;
extern int UV20;
extern int UV21;
extern int UV22;
extern int UV23;
extern int UV24;
extern int UV25;
extern int UV26;
extern int UV27;
extern int UV28;
extern int UV29;
extern int UV30;
extern int UV31;
extern int UV32;
extern int UV33;
extern int UV34;
extern int UV35;
extern int UV36;
extern int UV37;
extern int UV38;
extern int UV39;
extern int UV40;

extern int CRX1;
extern int CRX2;
extern int CRX3;
extern int CRX4;
extern int CRX5;
extern int CRX6;
extern int CRX7;
extern int CRX8;
extern int CRX9;
extern int CRX13;
extern int CRX14;
extern int CRX15;
extern int CRX16;
extern int CRX17;
extern int CRX18;
extern int CRX19;
extern int CRX20;
extern int CRX21;
