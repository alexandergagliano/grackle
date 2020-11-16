/***********************************************************************
/
/ Solve the chemistry and cooling
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"
#include "phys_constants.h"
#include "chemistry.h"
#include <solver.h>
#include <calc_mass.h>
#include "err.h"
#include "jac_rhs_builder.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <const.h>
#include <stdbool.h>
#ifdef _OPENMP
#include <omp.h>
#endif

extern chemistry_data *grackle_data;
extern chemistry_data_storage grackle_rates;
reaction_t *my_reactions;

#define tiny 1.e-50

/* function prototypes */

int update_UVbackground_rates(chemistry_data *my_chemistry,
                              chemistry_data_storage *my_rates,
                              photo_rate_storage *my_uvb_rates,
                              code_units *my_units);

extern void FORTRAN_NAME(solve_rate_cool_g)(
        int *icool,
	gr_float *d, gr_float *e, gr_float *u, gr_float *v, gr_float *w, gr_float *de,
	gr_float *HI, gr_float *HII, gr_float *HeI, gr_float *HeII, gr_float *HeIII,
	int *in, int *jn, int *kn, int *nratec, int *iexpand,
        int *ispecies, int *imetal, int *imcool, int *idust, int *idustall,
        int *idustfield, int *idim, int *iwater, int *water_rates,
	int *is, int *js, int *ks, int *ie, int *je, int *ke,
        int *ih2co, int *ipiht, int *igammah,
	double *dx, double *dt, double *aye, double *temstart, double *temend,
	double *utem, double *uxyz, double *uaye, double *urho, double *utim,
	double *gamma, double *fh, double *dtoh, double *z_solar, double *fgr,
	double *k1a, double *k2a, double *k3a, double *k4a, double *k5a,
	double *k6a, double *k7a, double *k8a, double *k9a, double *k10a,
	double *k11a, double *k12a, double *k13a, double *k13dda, double *k14a,
	double *k15a, double *k16a, double *k17a, double *k18a, double *k19a,
        double *k22a,	double *k24, double *k25, double *k26, double *k27,
        double *k28, double *k29, double *k30, double *k31,
	double *k50a, double *k51a, double *k52a, double *k53a, double *k54a,
	double *k55a, double *k56a, double *k57a, double *k58a,
	int *ndratec, double *dtemstart, double *dtemend, double *h2dusta,
	double *ncrna, double *ncrd1a, double *ncrd2a,
	double *ceHIa, double *ceHeIa, double *ceHeIIa, double *ciHIa,
	double *ciHeIa, double *ciHeISa, double *ciHeIIa,
        double *reHIIa, double *reHeII1a, double *reHeII2a, double *reHeIIIa,
        double *brema, double *compa, double *gammaha, double *isrf,
        double *regra, double *gamma_isrfa, double *comp_xraya, double *comp_temp,
	double *piHI, double *piHeI, double *piHeII,
	gr_float *HM, gr_float *H2I, gr_float *H2II,
        gr_float *DI, gr_float *DII, gr_float *HDI,
        gr_float *Water_density, gr_float *O_density, gr_float *OH_density,
        gr_float *O2_density, gr_float *Oplus_density, gr_float *OHplus_density,
        gr_float *H2Oplus_density, gr_float *H3Oplus_density, gr_float *O2plus_density,
        gr_float *Cplus_density, gr_float *C_density, gr_float *CH_density,
        gr_float *CH2_density, gr_float *CH3_density, gr_float *CH4_density,
        gr_float *CO_density, gr_float *COplus_density, gr_float *CO2_density,
        gr_float *CHplus_density, gr_float *CH2plus_density, gr_float *H3plus_density,
        gr_float *HCOplus_density, gr_float *HeHplus_density, gr_float *CH3plus_density,
        gr_float *CH4plus_density, gr_float *CH5plus_density, gr_float *O2Hplus_density,
        gr_float *metal, gr_float *dust,
	double *hyd01ka, double *h2k01a, double *vibha,
        double *rotha, double *rotla,
	double *gpldl, double *gphdl, double *HDltea, double *HDlowa,
	double *gaHIa, double *gaH2a, double *gaHea, double *gaHpa, double *gaela,
	double *h2ltea, double *gasgra, int *iH2shield,
        int *iradshield, double *avgsighi, double *avgsighei, double *avgsigheii,
        int *iradtrans, int *iradcoupled, int *iradstep, int *irt_honly,
        gr_float *kphHI, gr_float *kphHeI, gr_float *kphHeII, gr_float *kdissH2I,
        gr_float *photogamma, gr_float *xH2shield,
	int *ierr,
	int *ih2optical, int *iciecool, int *ithreebody, double *ciecoa,
 	int *icmbTfloor, int *iClHeat, double *clEleFra,
        long long *priGridRank, long long *priGridDim,
        double *priPar1, double *priPar2, double *priPar3,
        double *priPar4, double *priPar5,
 	long long *priDataSize, double *priCooling,
        double *priHeating, double *priMMW,
        long long *metGridRank, long long *metGridDim,
 	double *metPar1, double *metPar2, double *metPar3,
        double *metPar4, double *metPar5,
 	long long *metDataSize, double *metCooling,
        double *metHeating, int *clnew,
        int *iVheat, int *iMheat, gr_float *Vheat, gr_float *Mheat,
        int *iisrffield, gr_float* isrf_habing);

int local_solve_chemistry(chemistry_data *my_chemistry,
                          chemistry_data_storage *my_rates,
                          code_units *my_units,
                          grackle_field_data *my_fields,
                          double dt_value)
{

  /* Return if this doesn't concern us. */

  if (!my_chemistry->use_grackle)
    return SUCCESS;

  /* Update UV background rates. */
  photo_rate_storage my_uvb_rates;

  my_uvb_rates.k24 = my_uvb_rates.k25 = my_uvb_rates.k26 =
    my_uvb_rates.k27 = my_uvb_rates.k28 = my_uvb_rates.k29 =
    my_uvb_rates.k30 = my_uvb_rates.k31 = my_uvb_rates.piHI =
    my_uvb_rates.piHeI = my_uvb_rates.piHeII = my_uvb_rates.crsHI =
    my_uvb_rates.crsHeI = my_uvb_rates.crsHeII =
    my_uvb_rates.comp_xray = my_uvb_rates.temp_xray = 0.;

  if (my_chemistry->UVbackground == 1) {
    if (update_UVbackground_rates(my_chemistry, my_rates,
                                  &my_uvb_rates, my_units) == FAIL) {
      fprintf(stderr, "Error in update_UVbackground_rates.\n");
      return FAIL;
    }
  }
  else {
    my_uvb_rates.k24       = my_rates->k24;
    my_uvb_rates.k25       = my_rates->k25;
    my_uvb_rates.k26       = my_rates->k26;
    my_uvb_rates.k27       = my_rates->k27;
    my_uvb_rates.k28       = my_rates->k28;
    my_uvb_rates.k29       = my_rates->k29;
    my_uvb_rates.k30       = my_rates->k30;
    my_uvb_rates.k31       = my_rates->k31;
    my_uvb_rates.piHI      = my_rates->piHI;
    my_uvb_rates.piHeI     = my_rates->piHeI;
    my_uvb_rates.piHeII    = my_rates->piHeII;
    my_uvb_rates.crsHI     = my_rates->crsHI;
    my_uvb_rates.crsHeI    = my_rates->crsHeI;
    my_uvb_rates.crsHeII   = my_rates->crsHeII;
    my_uvb_rates.comp_xray = my_rates->comp_xray;
    my_uvb_rates.temp_xray = my_rates->temp_xray;
  }

  /* Check for a metal field. */

  int metal_field_present = TRUE;
  if (my_fields->metal_density == NULL)
    metal_field_present = FALSE;

  double co_length_units, co_density_units;
  if (my_units->comoving_coordinates == TRUE) {
    co_length_units = my_units->length_units;
    co_density_units = my_units->density_units;
  }
  else {
    co_length_units = my_units->length_units *
      my_units->a_value * my_units->a_units;
    co_density_units = my_units->density_units /
      POW(my_units->a_value * my_units->a_units, 3);
  }

  /* Error checking for H2 shielding approximation */
  if (my_chemistry->H2_self_shielding == 1 && my_fields->grid_rank != 3){
    fprintf(stderr, "Error in solve_chemistry: H2 self-shielding option 1 "
                    "will only work for 3D Cartesian grids. Use option 2 "
                    "to provide an array of shielding lengths with "
                    "H2_self_shielding_length or option 3 to use the "
                    "local Jeans length.");
    return FAIL;
  }


  if (my_chemistry->water_rates ==3 && my_chemistry->primordial_chemistry == 1){

    fprintf(stderr, "Error in solve_chemistry: Bialy (2019) water network "
                    "will only work with H2 species included. Set"
                    "multispecies to 2 or 3 or change water chemistry "
                    "option.");
    return FAIL;
   }

  /* Calculate temperature units. */

  double temperature_units =  mh * POW(my_units->velocity_units, 2) / kboltz;

  /* Call the fortran routine to solve cooling equations. */

  int ierr = 0;

  //if (my_chemistry->dust_chemistry){
  //  printf("gas_grain = %.2e\n", my_rates->gas_grain);
  //  printf("regr = %.2e\n", my_rates->regr);
  //}

  if(!my_chemistry->water_only){
     FORTRAN_NAME(solve_rate_cool_g)(
       &my_chemistry->with_radiative_cooling,
       my_fields->density,
       my_fields->internal_energy,
       my_fields->x_velocity,
       my_fields->y_velocity,
       my_fields->z_velocity,
       my_fields->e_density,
       my_fields->HI_density,
       my_fields->HII_density,
       my_fields->HeI_density,
       my_fields->HeII_density,
       my_fields->HeIII_density,
       my_fields->grid_dimension,
       my_fields->grid_dimension+1,
       my_fields->grid_dimension+2,
       &my_chemistry->NumberOfTemperatureBins,
       &my_units->comoving_coordinates,
       &my_chemistry->primordial_chemistry,
       &metal_field_present,
       &my_chemistry->metal_cooling,
       &my_chemistry->h2_on_dust,
       &my_chemistry->dust_chemistry,
       &my_chemistry->use_dust_density_field,
       &(my_fields->grid_rank),
       &my_chemistry->withWater, 
       &my_chemistry->water_rates,
       my_fields->grid_start,
       my_fields->grid_start+1,
       my_fields->grid_start+2,
       my_fields->grid_end,
       my_fields->grid_end+1,
       my_fields->grid_end+2,
       &my_chemistry->ih2co,
       &my_chemistry->ipiht,
       &my_chemistry->photoelectric_heating,
       &(my_fields->grid_dx),
       &dt_value,
       &my_units->a_value,
       &my_chemistry->TemperatureStart,
       &my_chemistry->TemperatureEnd,
       &temperature_units,
       &co_length_units,
       &my_units->a_units,
       &co_density_units,
       &my_units->time_units,
       &my_chemistry->Gamma,
       &my_chemistry->HydrogenFractionByMass,
       &my_chemistry->DeuteriumToHydrogenRatio,
       &my_chemistry->SolarMetalFractionByMass,
       &my_chemistry->local_dust_to_gas_ratio,
       my_rates->k1,
       my_rates->k2,
       my_rates->k3,
       my_rates->k4,
       my_rates->k5,
       my_rates->k6,
       my_rates->k7,
       my_rates->k8,
       my_rates->k9,
       my_rates->k10,
       my_rates->k11,
       my_rates->k12,
       my_rates->k13,
       my_rates->k13dd,
       my_rates->k14,
       my_rates->k15,
       my_rates->k16,
       my_rates->k17,
       my_rates->k18,
       my_rates->k19,
       my_rates->k22,
       &my_uvb_rates.k24,
       &my_uvb_rates.k25,
       &my_uvb_rates.k26,
       &my_uvb_rates.k27,
       &my_uvb_rates.k28,
       &my_uvb_rates.k29,
       &my_uvb_rates.k30,
       &my_uvb_rates.k31,
       my_rates->k50,
       my_rates->k51,
       my_rates->k52,
       my_rates->k53,
       my_rates->k54,
       my_rates->k55,
       my_rates->k56,
       my_rates->k57,
       my_rates->k58,
       &my_chemistry->NumberOfDustTemperatureBins,
       &my_chemistry->DustTemperatureStart,
       &my_chemistry->DustTemperatureEnd,
       my_rates->h2dust,
       my_rates->n_cr_n,
       my_rates->n_cr_d1,
       my_rates->n_cr_d2,
       my_rates->ceHI,
       my_rates->ceHeI,
       my_rates->ceHeII,
       my_rates->ciHI,
       my_rates->ciHeI,
       my_rates->ciHeIS,
       my_rates->ciHeII,
       my_rates->reHII,
       my_rates->reHeII1,
       my_rates->reHeII2,
       my_rates->reHeIII,
       my_rates->brem,
       &my_rates->comp,
       &my_rates->gammah,
       &my_chemistry->interstellar_radiation_field,
       my_rates->regr,
       &my_rates->gamma_isrf,
       &my_uvb_rates.comp_xray,
       &my_uvb_rates.temp_xray,
       &my_uvb_rates.piHI,
       &my_uvb_rates.piHeI,
       &my_uvb_rates.piHeII,
       my_fields->HM_density,
       my_fields->H2I_density,
       my_fields->H2II_density,
       my_fields->DI_density,
       my_fields->DII_density,
       my_fields->HDI_density,
       my_fields->Water_density,
       my_fields->O_density,
       my_fields->OH_density,
       my_fields->O2_density,
       my_fields->Oplus_density,
       my_fields->OHplus_density,
       my_fields->H2Oplus_density,
       my_fields->H3Oplus_density,
       my_fields->O2plus_density,
       my_fields->Cplus_density,
       my_fields->C_density,
       my_fields->CH_density,
       my_fields->CH2_density,
       my_fields->CH3_density,
       my_fields->CH4_density,
       my_fields->CO_density,
       my_fields->COplus_density,
       my_fields->CO2_density,
       my_fields->CHplus_density,
       my_fields->CH2plus_density,
       my_fields->H3plus_density,
       my_fields->HCOplus_density,
       my_fields->HeHplus_density,
       my_fields->CH3plus_density,
       my_fields->CH4plus_density,
       my_fields->CH5plus_density,
       my_fields->O2Hplus_density,
       my_fields->metal_density,
       my_fields->dust_density,
       my_rates->hyd01k,
       my_rates->h2k01,
       my_rates->vibh,
       my_rates->roth,
       my_rates->rotl,
       my_rates->GP99LowDensityLimit,
       my_rates->GP99HighDensityLimit,
       my_rates->HDlte,
       my_rates->HDlow,
       my_rates->GAHI,
       my_rates->GAH2,
       my_rates->GAHe,
       my_rates->GAHp,
       my_rates->GAel,
       my_rates->H2LTE,
       my_rates->gas_grain,
       &my_chemistry->H2_self_shielding,
       &my_chemistry->self_shielding_method,
       &my_uvb_rates.crsHI,
       &my_uvb_rates.crsHeI,
       &my_uvb_rates.crsHeII,
       &my_chemistry->use_radiative_transfer,
       &my_chemistry->radiative_transfer_coupled_rate_solver,
       &my_chemistry->radiative_transfer_intermediate_step,
       &my_chemistry->radiative_transfer_hydrogen_only,
       my_fields->RT_HI_ionization_rate,
       my_fields->RT_HeI_ionization_rate,
       my_fields->RT_HeII_ionization_rate,
       my_fields->RT_H2_dissociation_rate,
       my_fields->RT_heating_rate,
       my_fields-> H2_self_shielding_length,
       &ierr,
       &my_chemistry->h2_optical_depth_approximation,
       &my_chemistry->cie_cooling,
       &my_chemistry->three_body_rate,
       my_rates->cieco,
       &my_chemistry->cmb_temperature_floor,
       &my_chemistry->UVbackground,
       &my_chemistry->cloudy_electron_fraction_factor,
       &my_rates->cloudy_primordial.grid_rank,
       my_rates->cloudy_primordial.grid_dimension,
       my_rates->cloudy_primordial.grid_parameters[0],
       my_rates->cloudy_primordial.grid_parameters[1],
       my_rates->cloudy_primordial.grid_parameters[2],
       my_rates->cloudy_primordial.grid_parameters[3],
       my_rates->cloudy_primordial.grid_parameters[4],
       &my_rates->cloudy_primordial.data_size,
       my_rates->cloudy_primordial.cooling_data,
       my_rates->cloudy_primordial.heating_data,
       my_rates->cloudy_primordial.mmw_data,
       &my_rates->cloudy_metal.grid_rank,
       my_rates->cloudy_metal.grid_dimension,
       my_rates->cloudy_metal.grid_parameters[0],
       my_rates->cloudy_metal.grid_parameters[1],
       my_rates->cloudy_metal.grid_parameters[2],
       my_rates->cloudy_metal.grid_parameters[3],
       my_rates->cloudy_metal.grid_parameters[4],
       &my_rates->cloudy_metal.data_size,
       my_rates->cloudy_metal.cooling_data,
       my_rates->cloudy_metal.heating_data,
       &my_rates->cloudy_data_new,
       &my_chemistry->use_volumetric_heating_rate,
       &my_chemistry->use_specific_heating_rate,
       my_fields->volumetric_heating_rate,
       my_fields->specific_heating_rate,
       &my_chemistry->use_isrf_field,
       my_fields->isrf_habing);
}

//Decide whether to turn to the water network
  if (!my_chemistry->withWater || my_fields->metal_density == NULL)
     return SUCCESS;

// Iteration over all cells of grid
  int i, j, k, nstp;

  int i_start = *my_fields->grid_start;
  int i_end = *my_fields->grid_end;
  int j_start = *(my_fields->grid_start+1);
  int j_end = *(my_fields->grid_end+1);
  int k_start = *(my_fields->grid_start+2);
  int k_end   = *(my_fields->grid_end+2);

  ierr = 0;

  int di = i_end - i_start + 1;
  int dj = j_end - j_start + 1;
  int dk = k_end - k_start + 1;

  int nghost = i_start;

  int index_start = nghost * ((2*nghost + di) * (2*nghost + dj)) + nghost * (2*nghost + di) + nghost;

  /* flag to turn on UV rates in water network! */
  double UV_water = (double) my_chemistry->UVbackground;
  double d_to_n = co_density_units/mh;

  setup_rxns(my_chemistry->primordial_chemistry, UV_water, my_chemistry->crx_ionization, my_chemistry->water_rates);
  setup_species(my_chemistry->primordial_chemistry, UV_water, my_chemistry->crx_ionization, my_chemistry->water_rates);

  // Building the reactions is time consuming. We should only do it once.
  //static int first = 1;
   //static double *Y;
  //if(first){
    my_reactions = (reaction_t*) malloc(nReactions*sizeof(reaction_t));
    build_reactions(my_reactions,my_chemistry->primordial_chemistry,UV_water,my_chemistry->crx_ionization, my_chemistry->water_rates);
    double *Y = (double *) malloc(nSpecies * sizeof(double));
  //  first = 0;
  //}

  //# ifdef _OPENMP
  //# pragma omp parallel for schedule( runtime ) collapse(3) private(Y)
  //# endif
  for (k = 0; k < dk; k++){
    for (j = 0; j < dj; j++){
        for (i = 0; i < di; i++){

          //flatten 3D cube of space into 1d array
           int index = index_start + i + j*(di + 2*nghost) + k*(di*dj + di*2*nghost + dj*2*nghost + 4*nghost*nghost);

           double Z_solar = my_chemistry->SolarMetalFractionByMass;
           double metallicity = my_fields->metal_density[index] / my_fields->density[index] / Z_solar;

               //double temperature = (double) my_fields->internal_energy[index]* (double) temperature_units;
               double temperature;
               if(my_chemistry->water_only){
                temperature = 100.0;  //setting high for now
               }
               else{
                 temperature = (double) my_fields->internal_energy[index]* (double) temperature_units;
                 if( metallicity < 1.e-8){
                   continue;
                 }
               }

               Y[H] = my_fields->HI_density[index];
               Y[Hplus] = my_fields->HII_density[index];

               // convert to number density
              // (grackle solves its species
              // in code density)
               Y[H]     *= d_to_n;
               Y[Hplus] *= d_to_n;

               Y[el] = my_fields->e_density[index];
               Y[O] = my_fields->O_density[index];
               Y[OH] = my_fields->OH_density[index];
               Y[H2O] = my_fields->Water_density[index];
               Y[O2] = my_fields->O2_density[index];
               Y[Oplus] = my_fields->Oplus_density[index];
               Y[OHplus] = my_fields->OHplus_density[index];
               Y[H2Oplus] = my_fields->H2Oplus_density[index];
               Y[H3Oplus] = my_fields->H3Oplus_density[index];
               Y[O2plus] = my_fields->O2plus_density[index];
               Y[Cplus] = my_fields->Cplus_density[index];
               Y[C] = my_fields->C_density[index];
               Y[CH] = my_fields->CH_density[index];
               Y[CH2] = my_fields->CH2_density[index];
               Y[CH3] = my_fields->CH3_density[index];
               Y[CH4] = my_fields->CH4_density[index];
               Y[CO] = my_fields->CO_density[index];
               Y[COplus] = my_fields->COplus_density[index];
               Y[CO2] = my_fields->CO2_density[index];

               if (my_chemistry->primordial_chemistry > 1)
               {
                   Y[Hmin] = my_fields->HM_density[index];
                   Y[H2m] = my_fields->H2I_density[index];

                   Y[Hmin] *= d_to_n;
                   Y[H2m]  *= d_to_n;

                   if (my_chemistry->primordial_chemistry > 2)
                   {
                       Y[D] = my_fields->DI_density[index];
                       Y[Dplus] = my_fields->DII_density[index];
                       Y[HD] = my_fields->HDI_density[index];

                       // convert to number density!
                       Y[D]     *= d_to_n;
                       Y[Dplus] *= d_to_n;
                       Y[HD]    *= d_to_n;
                   }
               }

              //Bialy (2019) network
               if (my_chemistry->water_rates == 3){
                  Y[CHplus]  = my_fields->CHplus_density[index];
                  Y[CH2plus] = my_fields->CH2plus_density[index];
                  Y[He]      = my_fields->HeI_density[index];
                  Y[Heplus]  = my_fields->HeII_density[index];
                  Y[H3plus]  = my_fields->H3plus_density[index];
                  Y[HCOplus] = my_fields->HCOplus_density[index];
                  Y[H2plus]  = my_fields->H2II_density[index];
                  Y[HeHplus] = my_fields->HeHplus_density[index];
                  Y[CH3plus] = my_fields->CH3plus_density[index];
                  Y[CH4plus] = my_fields->CH4plus_density[index];
                  Y[CH5plus] = my_fields->CH5plus_density[index];
                  Y[O2Hplus] = my_fields->O2Hplus_density[index];

                  Y[H2plus] *= d_to_n;
                  Y[He]     *= d_to_n;
                  Y[Heplus] *= d_to_n;
               }


               //debugging printout 
              //for (int j = 0; j < nSpecies; j++){
              //    printf("%i, %.2e\n", j, Y[j]);
              //}

              double C_num_pre, O_num_pre, C_num_post, O_num_post, scale, sum_metl, metal_cgs, delta_mass, metl_frac, metal_frac;

              //if (my_chemistry->water_rates == 1){
                  // keep the number density of C species the same
                  C_num_pre = Y[C] + Y[Cplus] + Y[CH] + Y[CH2] + Y[CH3] + Y[CH4] + Y[CO] + Y[COplus] + Y[CO2];

                 if (my_chemistry->water_rates ==3){
                         C_num_pre += Y[CHplus] + Y[CH2plus] + Y[HCOplus] + Y[CH3plus] + Y[CH4plus] + Y[CH5plus];
                 }

                  // keep the number density of O species the same
                  O_num_pre = Y[O] + Y[OH] + Y[H2O] + 2.0*Y[O2] + Y[Oplus] + Y[OHplus] + Y[H2Oplus] + Y[H3Oplus] + 2.0*Y[O2plus] + Y[CO] + Y[COplus] + 2.0*Y[CO2];

                 if (my_chemistry->water_rates ==3){
                         O_num_pre += Y[HCOplus] + 2.0*Y[O2Hplus];
                 }
              //} //how much is it differing before or after?? Err on the side of not correcting?
              // Maybe indicate that you need to subcycle more?

               // complete one iteration of the water network
               ierr = integrate_network(my_chemistry->water_rates, Y, temperature, temperature, my_fields->density[index], metallicity*Z_solar, UV_water, my_chemistry->UVbackground_molec_redshift_on, my_chemistry->crx_ionization, dt_value * my_units->time_units, &nstp, my_units, my_chemistry->primordial_chemistry, my_chemistry->H2_self_shielding, my_uvb_rates.crsHI, my_uvb_rates.k24,my_chemistry->water_only, my_rates);

               if (ierr != 0 && ierr != MXSTP)
               {
                   printf("Error in network integration: %s\n", errmsg[ierr]);
                   exit(99);
               }
              //if (my_chemistry->water_rates == 1){
                  // keep the number density of C species the same
                  C_num_post = Y[C] + Y[Cplus] + Y[CH] + Y[CH2] + Y[CH3] + Y[CH4] + Y[CO] + Y[COplus] + Y[CO2];
                   if (my_chemistry->water_rates ==3){
                          C_num_post += Y[CHplus] + Y[CH2plus] + Y[HCOplus] + Y[CH3plus] + Y[CH4plus] + Y[CH5plus];
                  }

                  scale = C_num_pre/C_num_post;
                 if (abs(scale - 1.0) > 0.1){

                  Y[C]      *= scale;
                  Y[Cplus]  *= scale;
                  Y[CH]     *= scale;
                  Y[CH2]    *= scale;
                  Y[CH3]    *= scale;
                  Y[CH4]    *= scale;
                  Y[CO]     *= scale;
                  Y[COplus] *= scale;
                  Y[CO2]    *= scale;

                  if (my_chemistry->water_rates ==3){
                        Y[CHplus]  *= scale;
                        Y[CH2plus] *= scale;
                        Y[CH3plus] *= scale;
                        Y[CH4plus] *= scale;
                        Y[CH5plus] *= scale;
                  }
                 }
                 // keep the number density of O species the same

                  O_num_post = Y[O] + Y[OH] + Y[H2O] + 2.0*Y[O2] + Y[Oplus] + Y[OHplus] + Y[H2Oplus] + Y[H3Oplus] + 2.0*Y[O2plus] + Y[CO] + Y[COplus] + 2.0*Y[CO2];

                if (my_chemistry->water_rates ==3){
                         O_num_post += Y[HCOplus] + 2.0*Y[O2Hplus];
                  }

                  scale = O_num_pre/O_num_post;
                 if (abs(scale - 1.0) > 0.1){

                  Y[O]       *= scale;
                  Y[OH]      *= scale;
                  Y[H2O]     *= scale;
                  Y[O2]      *= scale;
                  Y[Oplus]   *= scale;
                  Y[OHplus]  *= scale;
                  Y[H2Oplus] *= scale;
                  Y[H3Oplus] *= scale;
                  Y[O2plus]  *= scale;
                  Y[CO]      *= scale;
                  Y[COplus]  *= scale;

                 if (my_chemistry->water_rates ==3){
                        Y[HCOplus] *= scale;
                        Y[O2Hplus] *= scale;
                     }
                 }

                  // Calculate mass densities of metal field and summed water metals
                  metal_cgs = my_fields->metal_density[index] * co_density_units;
                  sum_metl = calculate_metl_mass(Y, my_chemistry->primordial_chemistry, my_chemistry->water_rates);

                  // calculate deviation of network mass from all metals that should be in the network
                  //  NOTE: The fractions of metals enriched from supernovae that are included in the
                  //  water network are 0.553 for Oxygen and 0.228 for Carbon, summing to 0.781
                 if (my_units->comoving_coordinates == TRUE)
                 {
                    metal_frac = 0.781;
                 }
                 else{
                     metal_frac = 1.00;
                 }

                  double metal_exp = metal_frac*metal_cgs;
                  delta_mass = fabs((sum_metl - metal_exp)/(metal_exp));

                  // Ensure water metals don't sum to greater than metal field //
                  if ((sum_metl > metal_cgs) || (delta_mass > tiny)) {
                     // scale the water species back to appropriate values //
                     metl_frac = min(metal_cgs/sum_metl * metal_frac, 1.e50); //don't let us scale up by more than a certain amount in one iteration

                     Y[O]       *= metl_frac;
                     Y[OH]      *= metl_frac;
                     Y[H2O]     *= metl_frac;
                     Y[O2]      *= metl_frac;
                     Y[Oplus]   *= metl_frac;
                     Y[OHplus]  *= metl_frac;
                     Y[H2Oplus] *= metl_frac;
                     Y[H3Oplus] *= metl_frac;
                     Y[O2plus]  *= metl_frac;
                     Y[Cplus]   *= metl_frac;
                     Y[C]       *= metl_frac;
                     Y[CH]      *= metl_frac;
                     Y[CH2]     *= metl_frac;
                     Y[CH3]     *= metl_frac;
                     Y[CH4]     *= metl_frac;
                     Y[CO]      *= metl_frac;
                     Y[COplus]  *= metl_frac;
                     Y[CO2]     *= metl_frac;

                    if (my_chemistry->water_rates ==3){
                       Y[CHplus]  *= metl_frac;
                        Y[CH2plus] *= metl_frac;
                        Y[HCOplus] *= metl_frac;
                        Y[CH3plus] *= metl_frac;
                        Y[CH4plus] *= metl_frac;
                        Y[CH5plus] *= metl_frac;
                        Y[O2Hplus] *= metl_frac;
                    }
                   }

              // Set tiny floor for metal species - everything past
              for (int j = 0; j < nSpecies; j++){
                  Y[j] = max(Y[j], tiny);
              }

               //account for ionization at high temperature
               if (temperature >= 1.e5){
                  Y[Cplus] += Y[C];
                  Y[Oplus] += Y[O];
                  Y[C] = tiny;
                  Y[O] = tiny;
               }


               /* Write updated number densities back to metal fields */
               Y[H]   /= d_to_n;
               Y[Hplus] /= d_to_n;

               my_fields->HI_density[index] = Y[H];
               my_fields->HII_density[index] = Y[Hplus];
               my_fields->e_density[index] = Y[el];
               my_fields->O_density[index] = Y[O];
               my_fields->OH_density[index] = Y[OH];
               my_fields->Water_density[index] = Y[H2O];
               my_fields->O2_density[index] = Y[O2];
               my_fields->Oplus_density[index] = Y[Oplus];
               my_fields->OHplus_density[index] = Y[OHplus];
               my_fields->H2Oplus_density[index] = Y[H2Oplus];
               my_fields->H3Oplus_density[index] = Y[H3Oplus];
               my_fields->O2plus_density[index] = Y[O2plus];
               my_fields->Cplus_density[index] = Y[Cplus];
               my_fields->C_density[index] = Y[C];
               my_fields->CH_density[index] = Y[CH];
               my_fields->CH2_density[index] = Y[CH2];
               my_fields->CH3_density[index] = Y[CH3];
               my_fields->CH4_density[index] = Y[CH4];
               my_fields->CO_density[index] = Y[CO];
               my_fields->COplus_density[index] = Y[COplus];
               my_fields->CO2_density[index] = Y[CO2];

               if (my_chemistry->primordial_chemistry > 1)
               {
                   Y[Hmin] /= d_to_n;
                   Y[H2m] /= d_to_n;

                   my_fields->HM_density[index] = Y[Hmin];
                   my_fields->H2I_density[index] = Y[H2m];
                   if (my_chemistry->primordial_chemistry > 2)
                   {
                       Y[D]     /= d_to_n;
                       Y[Dplus] /= d_to_n;
                       Y[HD]    /= d_to_n;
                       my_fields->DI_density[index] = Y[D];
                       my_fields->DII_density[index] = Y[Dplus];
                       my_fields->HDI_density[index] = Y[HD];
                   }
               }
               if (my_chemistry->water_rates == 3){
                  Y[H2plus] /= d_to_n;
                  Y[He] /= d_to_n;
                  Y[Heplus] /= d_to_n;
                  my_fields->CHplus_density[index] = Y[CHplus];
                  my_fields->CH2plus_density[index] = Y[CH2plus];
                  my_fields->HeI_density[index] = Y[He];
                  my_fields->HeII_density[index] = Y[Heplus];
                  my_fields->H3plus_density[index] = Y[H3plus];
                  my_fields->HCOplus_density[index] = Y[HCOplus];
                  my_fields->H2II_density[index] = Y[H2plus];
                  my_fields->HeHplus_density[index] = Y[HeHplus];
                  my_fields->CH3plus_density[index] = Y[CH3plus];
                  my_fields->CH4plus_density[index] = Y[CH4plus];
                  my_fields->CH5plus_density[index] = Y[CH5plus];
                  my_fields->O2Hplus_density[index] = Y[O2Hplus];
               }
        }
     }
  }
        free(Y);
        free(my_reactions);
	return SUCCESS;
    }

int _solve_chemistry(chemistry_data *my_chemistry,
                     chemistry_data_storage *my_rates,
                     code_units *my_units, double dt_value, double dx_value,
                     int grid_rank, int *grid_dimension,
                     int *grid_start, int *grid_end,
                     gr_float *density, gr_float *internal_energy,
                     gr_float *x_velocity, gr_float *y_velocity, gr_float *z_velocity,
                     gr_float *HI_density, gr_float *HII_density, gr_float *HM_density,
                     gr_float *HeI_density, gr_float *HeII_density, gr_float *HeIII_density,
                     gr_float *H2I_density, gr_float *H2II_density,
                     gr_float *DI_density, gr_float *DII_density, gr_float *HDI_density,
                     gr_float *Water_density,
                     gr_float *O_density, gr_float *OH_density,
                     gr_float *O2_density, gr_float *Oplus_density, gr_float *OHplus_density,
                     gr_float *H2Oplus_density, gr_float *H3Oplus_density, gr_float *O2plus_density,
                     gr_float *Cplus_density, gr_float *C_density, gr_float *CH_density,
                     gr_float *CH2_density, gr_float *CH3_density, gr_float *CH4_density,
                     gr_float *CO_density, gr_float *COplus_density,
                     gr_float *CO2_density, gr_float *CHplus_density, gr_float *CH2plus_density,
                     gr_float *H3plus_density, gr_float *HCOplus_density, gr_float *HeHplus_density,
                     gr_float *CH3plus_density, gr_float *CH4plus_density,
                     gr_float *CH5plus_density, gr_float *O2Hplus_density,
                     gr_float *e_density, gr_float *metal_density, gr_float *dust_density,
                     gr_float *volumetric_heating_rate, gr_float *specific_heating_rate,
                     gr_float *RT_heating_rate, gr_float *RT_HI_ionization_rate, gr_float *RT_HeI_ionization_rate,
                     gr_float *RT_HeII_ionization_rate, gr_float *RT_H2_dissociation_rate,
                     gr_float *H2_self_shielding_length)
{

  grackle_field_data my_fields;
  my_fields.grid_dx                  = dx_value;
  my_fields.grid_rank                = grid_rank;
  my_fields.grid_dimension           = grid_dimension;
  my_fields.grid_start               = grid_start;
  my_fields.grid_end                 = grid_end;
  my_fields.density                  = density;
  my_fields.internal_energy          = internal_energy;
  my_fields.x_velocity               = x_velocity;
  my_fields.y_velocity               = y_velocity;
  my_fields.z_velocity               = z_velocity;
  my_fields.HI_density               = HI_density;
  my_fields.HII_density              = HII_density;
  my_fields.HM_density               = HM_density;
  my_fields.HeI_density              = HeI_density;
  my_fields.HeII_density             = HeII_density;
  my_fields.HeIII_density            = HeIII_density;
  my_fields.H2I_density              = H2I_density;
  my_fields.H2II_density             = H2II_density;
  my_fields.DI_density               = DI_density;
  my_fields.DII_density              = DII_density;
  my_fields.HDI_density              = HDI_density;
  my_fields.Water_density            = Water_density;
  my_fields.O_density                = O_density; 
  my_fields.OH_density               = OH_density;
  my_fields.O2_density               = O2_density; 
  my_fields.Oplus_density            = Oplus_density;
  my_fields.OHplus_density           = OHplus_density; 
  my_fields.H2Oplus_density          = H2Oplus_density;
  my_fields.H3Oplus_density          = H3Oplus_density; 
  my_fields.O2plus_density           = O2plus_density;
  my_fields.Cplus_density            = Cplus_density;
  my_fields.C_density                = C_density;
  my_fields.CH_density               = CH_density; 
  my_fields.CH2_density              = CH2_density;
  my_fields.CH3_density              = CH3_density;  
  my_fields.CH4_density              = CH4_density;
  my_fields.CO_density               = CO_density;
  my_fields.COplus_density           = COplus_density;
  my_fields.CO2_density              = CO2_density;
  my_fields.CHplus_density           = CHplus_density;
  my_fields.CH2plus_density          = CH2plus_density; 
  my_fields.H3plus_density           = H3plus_density;
  my_fields.COplus_density           = COplus_density;
  my_fields.HeHplus_density          = HeHplus_density;
  my_fields.CH3plus_density          = CH3plus_density;
  my_fields.CH4plus_density          = CH4plus_density;
  my_fields.CH5plus_density          = CH5plus_density;  
  my_fields.O2Hplus_density          = O2Hplus_density;
  my_fields.e_density                = e_density;
  my_fields.metal_density            = metal_density;
  my_fields.dust_density             = dust_density;
  my_fields.volumetric_heating_rate  = volumetric_heating_rate;
  my_fields.specific_heating_rate    = specific_heating_rate;
  my_fields.RT_heating_rate          = RT_heating_rate;
  my_fields.RT_HI_ionization_rate    = RT_HI_ionization_rate;
  my_fields.RT_HeI_ionization_rate   = RT_HeI_ionization_rate;
  my_fields.RT_HeII_ionization_rate  = RT_HeII_ionization_rate;
  my_fields.RT_H2_dissociation_rate  = RT_H2_dissociation_rate;
  my_fields.H2_self_shielding_length = H2_self_shielding_length;

  if (local_solve_chemistry(my_chemistry, my_rates,
                            my_units, &my_fields, dt_value) == FAIL) {
    fprintf(stderr, "Error in local_solve_chemistry.\n");
    return FAIL;
  }
  return SUCCESS;
}

int solve_chemistry(code_units *my_units,
                    grackle_field_data *my_fields,
                    double dt_value)
{
  if (local_solve_chemistry(grackle_data, &grackle_rates,
                            my_units, my_fields, dt_value) == FAIL) {
    fprintf(stderr, "Error in local_solve_chemistry.\n");
    return FAIL;
  }
  return SUCCESS;
}
