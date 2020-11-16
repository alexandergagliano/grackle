/***********************************************************************
/
/ Example executable using libgrackle
/
/
/ Copyright (c) 2013, Enzo/Grackle Development Team.
/
/ Distributed under the terms of the Enzo Public Licence.
/
/ The full license is in the file LICENSE, distributed with this 
/ software.
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>

extern "C" {
#include <grackle.h>
}

#define mh     1.67262171e-24   
#define kboltz 1.3806504e-16

int main(int argc, char *argv[])
{

  /*********************************************************************
  / Initial setup of units and chemistry objects.
  / This should be done at simulation start.
  *********************************************************************/

  // Set initial redshift (for internal units).
  double initial_redshift = 3.;

  // Enable output
  grackle_verbose = 1;

  // First, set up the units system.
  // These are conversions from code units to cgs.
  code_units my_units;
  my_units.comoving_coordinates = 0; // 1 if cosmological sim, 0 if not
  my_units.density_units = 1.67e-24;
  my_units.length_units = 1.0;
  my_units.time_units = 1.0e12;
  my_units.velocity_units = my_units.length_units / my_units.time_units;
  my_units.a_units = 1.0; // units for the expansion factor
  // Set expansion factor to 1 for non-cosmological simulation.
  my_units.a_value = 1. / (1. + initial_redshift) / my_units.a_units;

  // Second, create a chemistry object for parameters.  This needs to be a pointer.
  chemistry_data *my_grackle_data;
  my_grackle_data = new chemistry_data;
  if (set_default_chemistry_parameters(my_grackle_data) == 0) {
    fprintf(stderr, "Error in set_default_chemistry_parameters.\n");
    return EXIT_FAILURE;
  }
  // Set parameter values for chemistry.
  // Access the parameter storage with the struct you've created
  // or with the grackle_data pointer declared in grackle.h (see further below).
  grackle_data->use_grackle = 1;            // chemistry on
  grackle_data->use_isrf_field = 1;
  grackle_data->with_radiative_cooling = 1; // cooling on
  grackle_data->primordial_chemistry = 2;   // molecular network with H, He, D
  grackle_data->dust_chemistry = 1;
  grackle_data->metal_cooling = 1;          // metal cooling on
  grackle_data->UVbackground = 1;           // UV background on
  grackle_data->h2_on_dust = 1;
  grackle_data->grackle_data_file = "../../input/CloudyData_UVB=HM2012_shielded.h5"; // data file
  grackle_data->cmb_temperature_floor = 1;
  grackle_data->Gamma = 1.6667;
  grackle_data->use_dust_density_field            = 0;

  grackle_data->use_isrf_field                    = 0;
  grackle_data->interstellar_radiation_field      = 1.7;
  grackle_data->use_volumetric_heating_rate       = 0;
  grackle_data->use_specific_heating_rate         = 0;
  grackle_data->three_body_rate                   = 0;
  grackle_data->cie_cooling                       = 0;
  grackle_data->h2_optical_depth_approximation    = 0;
  grackle_data->ih2co                             = 1;
  grackle_data->ipiht                             = 1;
  grackle_data->HydrogenFractionByMass            = 0.76;
  grackle_data->DeuteriumToHydrogenRatio          = 6.8e-05;
  grackle_data->SolarMetalFractionByMass          = 0.01295;
  grackle_data->local_dust_to_gas_ratio           = 0.009387;
  grackle_data->LWbackground_sawtooth_suppression = 0;
  grackle_data->LWbackground_intensity            = 0;
  grackle_data->cloudy_electron_fraction_factor   = 0.00915396;
  grackle_data->use_radiative_transfer            = 1;
  grackle_data->radiative_transfer_coupled_rate_solver = 0;
  grackle_data->radiative_transfer_intermediate_step = 0;
  grackle_data->radiative_transfer_hydrogen_only  = 0;

  //testing waternet
  grackle_data->withWater = 1;
  grackle_data->water_only = 0;
  grackle_data->water_rates = 3;
  grackle_data->crx_ionization = 1;
  grackle_data->grackle_molecular_data = "/lustre/scratch3/turquoise/agagliano/WATER/summer2020/fromSky/grackle/input/UVB_HM2012_waterNetwork.h5";

  grackle_data->photoelectric_heating             = 0;
  grackle_data->photoelectric_heating_rate        = 8.5e-26;
  grackle_data->use_volumetric_heating_rate       = 0;
  grackle_data->use_specific_heating_rate         = 0;
  grackle_data->three_body_rate                   = 0;
  grackle_data->cie_cooling                       = 0;
  grackle_data->h2_optical_depth_approximation    = 0;
  grackle_data->ih2co                             = 1;
  grackle_data->ipiht                             = 1;
  grackle_data->HydrogenFractionByMass            = 0.76;
  grackle_data->DeuteriumToHydrogenRatio          = 6.8e-05;
  grackle_data->SolarMetalFractionByMass          = 0.01295;
  grackle_data->local_dust_to_gas_ratio           = 0.009387;
  grackle_data->NumberOfTemperatureBins           = 600;
  grackle_data->CaseBRecombination                = 0;
  grackle_data->TemperatureStart                  = 1;
  grackle_data->TemperatureEnd                    = 1e+09;
  grackle_data->NumberOfDustTemperatureBins       = 250;
  grackle_data->DustTemperatureStart              = 1;
  grackle_data->DustTemperatureEnd                = 1500;
  grackle_data->Compton_xray_heating              = 0;
  grackle_data->UVbackground_redshift_on          = 7;
  grackle_data->UVbackground_redshift_off         = 0;
  grackle_data->UVbackground_redshift_fullon      = 6;
  grackle_data->UVbackground_redshift_drop        = 0;
  grackle_data->cloudy_electron_fraction_factor   = 0.00915396;
  grackle_data->use_radiative_transfer            = 1;

  grackle_data->H2_self_shielding = 1;

  // Finally, initialize the chemistry object.
  if (initialize_chemistry_data(&my_units) == 0) {
    fprintf(stderr, "Error in initialize_chemistry_data.\n");
    return EXIT_FAILURE;
  }

  gr_float tiny_number = 1.e-20;

  // Create struct for storing grackle field data
  grackle_field_data my_fields;

  // Set grid dimension and size.
  // grid_start and grid_end are used to ignore ghost zones.
  int field_size = 1;
  my_fields.grid_rank = 3;
  my_fields.grid_dimension = new int[3];
  my_fields.grid_start = new int[3];
  my_fields.grid_end = new int[3];
  for (int i = 0;i < 3;i++) {
    my_fields.grid_dimension[i] = 1; // the active dimension not including ghost zones.
    my_fields.grid_start[i] = 0;
    my_fields.grid_end[i] = 0;
  }
  my_fields.grid_dimension[0] = field_size;
  my_fields.grid_end[0] = field_size - 1;
  my_fields.grid_dx = 0.0; // used only for H2 self-shielding approximation

  my_fields.density         = new gr_float[field_size];
  my_fields.internal_energy = new gr_float[field_size];
  my_fields.x_velocity      = new gr_float[field_size];
  my_fields.y_velocity      = new gr_float[field_size];
  my_fields.z_velocity      = new gr_float[field_size];
  // for primordial_chemistry >= 1
  my_fields.HI_density      = new gr_float[field_size];
  my_fields.HII_density     = new gr_float[field_size];
  my_fields.HeI_density     = new gr_float[field_size];
  my_fields.HeII_density    = new gr_float[field_size];
  my_fields.HeIII_density   = new gr_float[field_size];
  my_fields.e_density       = new gr_float[field_size];
  // for primordial_chemistry >= 2
  my_fields.HM_density      = new gr_float[field_size];
  my_fields.H2I_density     = new gr_float[field_size];
  my_fields.H2II_density    = new gr_float[field_size];
  // for primordial_chemistry >= 3
  my_fields.DI_density      = new gr_float[field_size];
  my_fields.DII_density     = new gr_float[field_size];
  my_fields.HDI_density     = new gr_float[field_size];
  // for metal_cooling = 1
  my_fields.metal_density   = new gr_float[field_size];

  my_fields.Water_density     = new gr_float[field_size];
  my_fields.O_density         = new gr_float[field_size];
  my_fields.OH_density        = new gr_float[field_size];
  my_fields.O2_density        = new gr_float[field_size];
  my_fields.Oplus_density     = new gr_float[field_size];
  my_fields.OHplus_density    = new gr_float[field_size];
  my_fields.H2Oplus_density   = new gr_float[field_size];
  my_fields.H3Oplus_density   = new gr_float[field_size];
  my_fields.O2plus_density    = new gr_float[field_size];
  my_fields.Cplus_density     = new gr_float[field_size];
  my_fields.C_density         = new gr_float[field_size];
  my_fields.CH_density        = new gr_float[field_size];
  my_fields.CH2_density       = new gr_float[field_size];
  my_fields.CH3_density       = new gr_float[field_size];
  my_fields.CH4_density       = new gr_float[field_size];
  my_fields.CO_density        = new gr_float[field_size];
  my_fields.COplus_density    = new gr_float[field_size];
  my_fields.CO2_density       = new gr_float[field_size];
  my_fields.CHplus_density    = new gr_float[field_size];
  my_fields.CH2plus_density   = new gr_float[field_size];
  my_fields.H3plus_density    = new gr_float[field_size];
  my_fields.HCOplus_density   = new gr_float[field_size];
  my_fields.HeHplus_density   = new gr_float[field_size];
  my_fields.CH3plus_density   = new gr_float[field_size];
  my_fields.CH4plus_density   = new gr_float[field_size];
  my_fields.CH5plus_density   = new gr_float[field_size];
  my_fields.O2Hplus_density   = new gr_float[field_size];

  // volumetric heating rate (provide in units [erg s^-1 cm^-3])
  my_fields.volumetric_heating_rate = new gr_float[field_size];
  // specific heating rate (provide in units [egs s^-1 g^-1]
  my_fields.specific_heating_rate = new gr_float[field_size];

  // radiative transfer ionization / dissociation rate fields (provided in units of [1/s])
  my_fields.RT_HI_ionization_rate = new gr_float[field_size];
  my_fields.RT_HeI_ionization_rate = new gr_float[field_size];
  my_fields.RT_HeII_ionization_rate = new gr_float[field_size];
  my_fields.RT_H2_dissociation_rate = new gr_float[field_size];
  // radiative transfer heating rate (provide in units [erg s^-3 cm^-3])
  my_fields.RT_heating_rate = new gr_float[field_size];

  // interstellar radiation field strength
  my_fields.isrf_habing = new gr_float[field_size];

  // set temperature units
  double temperature_units = mh * pow(my_units.a_units * 
                                      my_units.length_units /
                                      my_units.time_units, 2) / kboltz;

  int i;
  for (i = 0;i < field_size;i++) {
    my_fields.density[i] = 1.0;
    my_fields.HI_density[i] = grackle_data->HydrogenFractionByMass *
      my_fields.density[i];
    my_fields.HII_density[i] = tiny_number * my_fields.density[i];
    my_fields.HM_density[i] = tiny_number * my_fields.density[i];
    my_fields.HeI_density[i] = (1.0 - grackle_data->HydrogenFractionByMass) *
      my_fields.density[i];
    my_fields.HeII_density[i] = tiny_number * my_fields.density[i];
    my_fields.HeIII_density[i] = tiny_number * my_fields.density[i];
    my_fields.H2I_density[i] = tiny_number * my_fields.density[i];
    my_fields.H2II_density[i] = tiny_number * my_fields.density[i];
    my_fields.DI_density[i] = 2.0 * 3.4e-5 * my_fields.density[i];
    my_fields.DII_density[i] = tiny_number * my_fields.density[i];
    my_fields.HDI_density[i] = tiny_number * my_fields.density[i];
    my_fields.e_density[i] = tiny_number * my_fields.density[i];
    // solar metallicity
    my_fields.metal_density[i] = grackle_data->SolarMetalFractionByMass *
      my_fields.density[i];

    my_fields.x_velocity[i] = 0.0;
    my_fields.y_velocity[i] = 0.0;
    my_fields.z_velocity[i] = 0.0;

    // initilize internal energy (here 1000 K for no reason)
    my_fields.internal_energy[i] = 1000. / temperature_units;

    my_fields.volumetric_heating_rate[i] = 0.0;
    my_fields.specific_heating_rate[i] = 0.0;

    my_fields.RT_HI_ionization_rate[i] = 0.0;
    my_fields.RT_HeI_ionization_rate[i] = 0.0;
    my_fields.RT_HeII_ionization_rate[i] = 0.0;
    my_fields.RT_H2_dissociation_rate[i] = 0.0;
    my_fields.RT_heating_rate[i] = 0.0;

    my_fields.isrf_habing[i] = grackle_data->interstellar_radiation_field;

    my_fields.Water_density[i] = tiny_number * my_fields.density[i];
    my_fields.O_density[i] = tiny_number * my_fields.density[i];
    my_fields.OH_density[i] = tiny_number * my_fields.density[i];
    my_fields.O2_density[i] = tiny_number * my_fields.density[i];
    my_fields.Oplus_density[i] = tiny_number * my_fields.density[i];
    my_fields.OHplus_density[i] = tiny_number * my_fields.density[i];
    my_fields.H2Oplus_density[i] = tiny_number * my_fields.density[i];
    my_fields.H3Oplus_density[i] = tiny_number * my_fields.density[i];
    my_fields.O2plus_density[i] = tiny_number * my_fields.density[i];
    my_fields.Cplus_density[i] = tiny_number * my_fields.density[i];
    my_fields.C_density[i] = tiny_number * my_fields.density[i];
    my_fields.CH_density[i] = tiny_number * my_fields.density[i];
    my_fields.CH2_density[i] = tiny_number * my_fields.density[i];
    my_fields.CH3_density[i] = tiny_number * my_fields.density[i];
    my_fields.CH4_density[i] = tiny_number * my_fields.density[i];
    my_fields.CO_density[i] = tiny_number * my_fields.density[i];
    my_fields.COplus_density[i] = tiny_number * my_fields.density[i];
    my_fields.CO2_density[i] = tiny_number * my_fields.density[i];
    my_fields.CHplus_density[i] = tiny_number * my_fields.density[i];
    my_fields.CH2plus_density[i] = tiny_number * my_fields.density[i];
    my_fields.H3plus_density[i] = tiny_number * my_fields.density[i];
    my_fields.HCOplus_density[i] = tiny_number * my_fields.density[i];
    my_fields.HeHplus_density[i] = tiny_number * my_fields.density[i];
    my_fields.CH3plus_density[i] = tiny_number * my_fields.density[i];
    my_fields.CH4plus_density[i] = tiny_number * my_fields.density[i];
    my_fields.CH5plus_density[i] = tiny_number * my_fields.density[i];
    my_fields.O2Hplus_density[i] = tiny_number * my_fields.density[i];
  }

  /*********************************************************************
  / Calling the chemistry solver
  / These routines can now be called during the simulation.
  *********************************************************************/

  // Evolving the chemistry.
  // some timestep
  double dt = 3.15e7 * 1e6 / my_units.time_units;

  if (solve_chemistry(&my_units, &my_fields, dt) == 0) {
    fprintf(stderr, "Error in solve_chemistry.\n");
    return EXIT_FAILURE;
  }

  // Calculate cooling time.
  gr_float *cooling_time;
  cooling_time = new gr_float[field_size];
  if (calculate_cooling_time(&my_units, &my_fields,
                             cooling_time) == 0) {
    fprintf(stderr, "Error in calculate_cooling_time.\n");
    return EXIT_FAILURE;
  }

  fprintf(stderr, "Cooling time = %le s.\n", cooling_time[0] *
          my_units.time_units);

  // Calculate temperature.
  gr_float *temperature;
  temperature = new gr_float[field_size];
  if (calculate_temperature(&my_units, &my_fields,
                            temperature) == 0) {
    fprintf(stderr, "Error in calculate_temperature.\n");
    return EXIT_FAILURE;
  }

  fprintf(stderr, "Temperature = %le K.\n", temperature[0]);

  // Calculate pressure.
  gr_float *pressure;
  double pressure_units = my_units.density_units *
    pow(my_units.velocity_units, 2);
  pressure = new gr_float[field_size];
  if (calculate_pressure(&my_units, &my_fields,
                         pressure) == 0) {
    fprintf(stderr, "Error in calculate_pressure.\n");
    return EXIT_FAILURE;
  }

  fprintf(stderr, "Pressure = %le dyne/cm^2.\n", pressure[0]*pressure_units);

  // Calculate gamma.
  gr_float *gamma;
  gamma = new gr_float[field_size];
  if (calculate_gamma(&my_units, &my_fields,
                      gamma) == 0) {
    fprintf(stderr, "Error in calculate_gamma.\n");
    return EXIT_FAILURE;
  }

  fprintf(stderr, "gamma = %le.\n", gamma[0]);

    // Calculate dust temperature.
  gr_float *dust_temperature;
  dust_temperature = new gr_float[field_size];
  if (calculate_dust_temperature(&my_units, &my_fields,
                                 dust_temperature) == 0) {
    fprintf(stderr, "Error in calculate_dust_temperature.\n");
    return EXIT_FAILURE;
  }

  fprintf(stderr, "dust_temperature = %g K.\n", dust_temperature[0]);

  _free_chemistry_data(my_grackle_data, &grackle_rates);

  return EXIT_SUCCESS;
}
