/***********************************************************************
/
/ Initialize UV background data
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
#include <math.h>
#include "hdf5.h"
#include "grackle_macros.h"
#include "grackle_types.h"
#include "grackle_chemistry_data.h"

extern int grackle_verbose;

// function prototypes
int read_dataset(hid_t file_id, char *dset_name, double *buffer);

// Initialize UV Background data
int initialize_UVbackground_data(chemistry_data *my_chemistry,
                                 chemistry_data_storage *my_rates)
{
  long long Nz;

  // Return if no UV background selected or using fully tabulated cooling.
  if (my_chemistry->UVbackground == 0 ||
      my_chemistry->primordial_chemistry == 0)
    return SUCCESS;


  if (grackle_verbose)
    fprintf(stdout, "Initializing UV background.\n");


  // Read in UV background data from hdf5 file.

  hid_t       file_id, dset_id, dspace_id;
  herr_t      status;
  herr_t      h5_error = -1;

  if (grackle_verbose)
    fprintf(stdout, "Reading UV background data from %s.\n",
            my_chemistry->grackle_data_file);
  file_id = H5Fopen(my_chemistry->grackle_data_file, 
                    H5F_ACC_RDONLY, H5P_DEFAULT);


  // Read Info dataset

  dset_id =  H5Dopen(file_id, "/UVBRates/Info");
  if (dset_id == h5_error) {
    fprintf(stderr, "Can't open 'Info' dataset in %s.\n",
            my_chemistry->grackle_data_file);
    return FAIL;
  }

  int strlen = (int)(H5Dget_storage_size(dset_id));
  char info_string[strlen+1];

  hid_t memtype = H5Tcopy(H5T_C_S1);
  H5Tset_size(memtype, strlen+1);

  status = H5Dread(dset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, info_string);
  if (status == h5_error) {
    fprintf(stderr, "Failed to read dataset 'Info'.\n");
    return FAIL;
  }

  H5Tclose(memtype); 
  H5Dclose(dset_id);



  // Open redshift dataset and get number of elements

  dset_id =  H5Dopen(file_id, "/UVBRates/z");
  if (dset_id == h5_error) {
    fprintf(stderr, "Can't open redshift dataset ('z') in %s.\n",
            my_chemistry->grackle_data_file);
    return FAIL;
  }

  dspace_id = H5Dget_space(dset_id);
  if (dspace_id == h5_error) {
    fprintf(stderr, "Error opening dataspace for dataset 'z' in %s.\n",
            my_chemistry->grackle_data_file);
    return FAIL;
  }

  Nz = H5Sget_simple_extent_npoints(dspace_id);
  if(Nz <= 0) {
    fprintf(stderr, "Redshift dataset ('z') has inappropriate size = %lld in %s.\n",
            Nz, my_chemistry->grackle_data_file);
    return FAIL;
  }

  H5Sclose(dspace_id);
  H5Dclose(dset_id);



  // Now allocate memory for UV background table.
  my_rates->UVbackground_table.Nz = Nz;

  my_rates->UVbackground_table.z = malloc(Nz * sizeof(double));
  my_rates->UVbackground_table.k24 = malloc(Nz * sizeof(double));
  my_rates->UVbackground_table.k25 = malloc(Nz * sizeof(double));
  my_rates->UVbackground_table.k26 = malloc(Nz * sizeof(double));

  if (my_chemistry->primordial_chemistry > 1) {
    my_rates->UVbackground_table.k27 = malloc(Nz * sizeof(double));
    my_rates->UVbackground_table.k28 = malloc(Nz * sizeof(double));
    my_rates->UVbackground_table.k29 = malloc(Nz * sizeof(double));
    my_rates->UVbackground_table.k30 = malloc(Nz * sizeof(double));
    my_rates->UVbackground_table.k31 = malloc(Nz * sizeof(double));
  }    

  my_rates->UVbackground_table.piHI = malloc(Nz * sizeof(double));
  my_rates->UVbackground_table.piHeII = malloc(Nz * sizeof(double));
  my_rates->UVbackground_table.piHeI = malloc(Nz * sizeof(double));

  if (my_chemistry->self_shielding_method > 0){
    my_rates->UVbackground_table.crsHI   = malloc(Nz * sizeof(double));
    my_rates->UVbackground_table.crsHeII = malloc(Nz * sizeof(double));
    my_rates->UVbackground_table.crsHeI  = malloc(Nz * sizeof(double));
  }


  // Now read everything.


  // *** Redshift ***
  if(! read_dataset(file_id, "/UVBRates/z",
                    my_rates->UVbackground_table.z) ) {
    fprintf(stderr, "Error reading dataset 'z' in %s.\n",
            my_chemistry->grackle_data_file);
    return FAIL;
  }

  // *** k24 ***
  if(! read_dataset(file_id, "/UVBRates/Chemistry/k24",
                    my_rates->UVbackground_table.k24) ) {
    fprintf(stderr, "Error reading dataset '/UVBRates/Chemistry/k24' in %s.\n",
            my_chemistry->grackle_data_file);
    return FAIL;
  }

  // *** k25 ***
  if(! read_dataset(file_id, "/UVBRates/Chemistry/k25",
                    my_rates->UVbackground_table.k25) ) {
    fprintf(stderr, "Error reading dataset '/UVBRates/Chemistry/k25' in %s.\n",
            my_chemistry->grackle_data_file);
    return FAIL;
  }

  // *** k26 ***
  if(! read_dataset(file_id, "/UVBRates/Chemistry/k26",
                    my_rates->UVbackground_table.k26) ) {
    fprintf(stderr, "Error reading dataset '/UVBRates/Chemistry/k26' in %s.\n",
            my_chemistry->grackle_data_file);
    return FAIL;
  }

  if (my_chemistry->primordial_chemistry > 1) {

    // *** k27 ***
    if(! read_dataset(file_id, "/UVBRates/Chemistry/k27",
                      my_rates->UVbackground_table.k27) ) {
      fprintf(stderr, "Error reading dataset '/UVBRates/Chemistry/k27' in %s.\n",
              my_chemistry->grackle_data_file);
      return FAIL;      
    }

    // *** k28 ***
    if(! read_dataset(file_id, "/UVBRates/Chemistry/k28",
                      my_rates->UVbackground_table.k28) ) {
      fprintf(stderr, "Error reading dataset '/UVBRates/Chemistry/k28' in %s.\n",
              my_chemistry->grackle_data_file);
      return FAIL;      
    }

    // *** k29 ***
    if(! read_dataset(file_id, "/UVBRates/Chemistry/k29",
                      my_rates->UVbackground_table.k29) ) {
      fprintf(stderr, "Error reading dataset '/UVBRates/Chemistry/k29' in %s.\n",
              my_chemistry->grackle_data_file);
      return FAIL;      
    }

    // *** k30 ***
    if(! read_dataset(file_id, "/UVBRates/Chemistry/k30",
                      my_rates->UVbackground_table.k30) ) {
      fprintf(stderr, "Error reading dataset '/UVBRates/Chemistry/k30' in %s.\n",
              my_chemistry->grackle_data_file);
      return FAIL;      
    }

    // *** k31 ***
    if(! read_dataset(file_id, "/UVBRates/Chemistry/k31",
                      my_rates->UVbackground_table.k31) ) {
      fprintf(stderr, "Error reading dataset '/UVBRates/Chemistry/k31' in %s.\n",
              my_chemistry->grackle_data_file);
      return FAIL;      
    }
    
  }

  // *** piHI ***
  if(! read_dataset(file_id, "/UVBRates/Photoheating/piHI",
                    my_rates->UVbackground_table.piHI) ) {
    fprintf(stderr, "Error reading dataset '/UVBRates/Photoheating/piHI' in %s.\n",
            my_chemistry->grackle_data_file);
    return FAIL;
  }

  // *** piHeII ***
  if(! read_dataset(file_id, "/UVBRates/Photoheating/piHeII",
                    my_rates->UVbackground_table.piHeII) ) {
    fprintf(stderr, "Error reading dataset '/UVBRates/Photoheating/piHeII' in %s.\n",
            my_chemistry->grackle_data_file);
    return FAIL;
  }

  // *** piHeI ***
  if(! read_dataset(file_id, "/UVBRates/Photoheating/piHeI",
                    my_rates->UVbackground_table.piHeI) ) {
    fprintf(stderr, "Error reading dataset '/UVBRates/Photoheating/piHeI' in %s.\n",
            my_chemistry->grackle_data_file);
    return FAIL;
  }


  if (my_chemistry->self_shielding_method > 0) {
    // *** crsHI ***
    if(! read_dataset(file_id, "/UVBRates/CrossSections/hi_avg_crs",
                      my_rates->UVbackground_table.crsHI) ) {
      fprintf(stderr, "Error reading dataset '/UVBRates/CrossSections/hi_avg_crs' in %s.\n",
              my_chemistry->grackle_data_file);
      fprintf(stderr, "In order to use self-shielding, you must use the shielding datasets\n");
      return FAIL;
    }

    // *** crsHeII ***
    if(! read_dataset(file_id, "/UVBRates/CrossSections/heii_avg_crs",
                    my_rates->UVbackground_table.crsHeII) ) {
      fprintf(stderr, "Error reading dataset '/UVBRates/CrossSections/heii_avg_crs' in %s.\n",
              my_chemistry->grackle_data_file);
      fprintf(stderr, "In order to use self-shielding, you must use the shielding datasets\n");
      return FAIL;
    }

    // *** crsHeI ***
    if(! read_dataset(file_id, "/UVBRates/CrossSections/hei_avg_crs",
                      my_rates->UVbackground_table.crsHeI) ) {
      fprintf(stderr, "Error reading dataset '/UVBRates/CrossSections/hei_avg_crs' in %s.\n",
              my_chemistry->grackle_data_file);
      fprintf(stderr, "In order to use self-shielding, you must use the shielding datasets\n");
      return FAIL;
    }
  }


  H5Fclose(file_id);


  // Get min/max of redshift vector
  my_rates->UVbackground_table.zmin = my_rates->UVbackground_table.z[0];
  my_rates->UVbackground_table.zmax = my_rates->UVbackground_table.z[Nz-1];

  // Print out some information about the dataset just read in.
  if (grackle_verbose) {
    fprintf(stdout, "UV background information:\n");
    fprintf(stdout, "  %s\n",info_string);
    fprintf(stdout, "  z_min = %6.3f\n  z_max = %6.3f\n",
            my_rates->UVbackground_table.zmin,
            my_rates->UVbackground_table.zmax);
  }

  // Set redshift on/off flags from data if not set.
  if (my_chemistry->UVbackground_redshift_on <= FLOAT_UNDEFINED) {
    my_chemistry->UVbackground_redshift_on =
      my_rates->UVbackground_table.zmax;
    if (grackle_verbose)
      fprintf(stdout, "Setting UVbackground_redshift_on to %f.\n",
              my_chemistry->UVbackground_redshift_on);
  }
  if (my_chemistry->UVbackground_redshift_fullon <= FLOAT_UNDEFINED) {
    my_chemistry->UVbackground_redshift_fullon =
      my_rates->UVbackground_table.zmax;
    if (grackle_verbose)
      fprintf(stdout, "Setting UVbackground_redshift_fullon to %f.\n",
              my_chemistry->UVbackground_redshift_fullon);
  }
  if (my_chemistry->UVbackground_redshift_drop <= FLOAT_UNDEFINED) {
    my_chemistry->UVbackground_redshift_drop =
      my_rates->UVbackground_table.zmin;
    if (grackle_verbose)
      fprintf(stdout, "Setting UVbackground_redshift_drop to %f.\n",
              my_chemistry->UVbackground_redshift_drop);
  }
  if (my_chemistry->UVbackground_redshift_off <= FLOAT_UNDEFINED) {
    my_chemistry->UVbackground_redshift_off =
      my_rates->UVbackground_table.zmin;
    if (grackle_verbose)
      fprintf(stdout, "Setting UVbackground_redshift_off to %f.\n",
              my_chemistry->UVbackground_redshift_off);
  }

  /**********************************************************************/
  /* Adding skeleton for molecular chemistry! We'll use this in rates.c */
  /**********************************************************************/
    long long Nz_molec;

  // Return if no UV background selected or if water network is off
  if (my_chemistry->UVbackground == 0 ||
      my_chemistry->withWater == 0 || my_chemistry->water_only == 1)
    return SUCCESS;


  if (grackle_verbose)
    fprintf(stdout, "Initializing UV molecular rates.\n");


  // Read in UV molecular rates from hdf5 file.

  if (grackle_verbose)
    fprintf(stdout, "Reading UV molecular rates from %s.\n",
            my_chemistry->grackle_molecular_data);
  file_id = H5Fopen(my_chemistry->grackle_molecular_data,
                    H5F_ACC_RDONLY, H5P_DEFAULT);


  // Open redshift dataset and get number of elements

  dset_id =  H5Dopen(file_id, "z");
  if (dset_id == h5_error) {
    fprintf(stderr, "Can't open redshift dataset ('z') in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  dspace_id = H5Dget_space(dset_id);
  if (dspace_id == h5_error) {
    fprintf(stderr, "Error opening dataspace for dataset 'z' in %s.\n",
            my_chemistry->grackle_data_file);
    return FAIL;
  }

  Nz_molec = 60;

  H5Sclose(dspace_id);
  H5Dclose(dset_id);

  // Now allocate memory for UV background table.
  my_rates->UVbackground_table.Nz_molec = Nz_molec;

  my_rates->UVbackground_table.z_molec = malloc(Nz_molec * sizeof(double));
  my_rates->UVbackground_table.UV1 = malloc(Nz_molec * sizeof(double));
  my_rates->UVbackground_table.UV2 = malloc(Nz_molec * sizeof(double));
  my_rates->UVbackground_table.UV3 = malloc(Nz_molec * sizeof(double));
  my_rates->UVbackground_table.UV4 = malloc(Nz_molec * sizeof(double));
  my_rates->UVbackground_table.UV5 = malloc(Nz_molec * sizeof(double));
  my_rates->UVbackground_table.UV6 = malloc(Nz_molec * sizeof(double));
  my_rates->UVbackground_table.UV7 = malloc(Nz_molec * sizeof(double));
  my_rates->UVbackground_table.UV8 = malloc(Nz_molec * sizeof(double));
  my_rates->UVbackground_table.UV9 = malloc(Nz_molec * sizeof(double));
  my_rates->UVbackground_table.UV10 = malloc(Nz_molec * sizeof(double));
  my_rates->UVbackground_table.UV11 = malloc(Nz_molec * sizeof(double));
  my_rates->UVbackground_table.UV12 = malloc(Nz_molec * sizeof(double));
  my_rates->UVbackground_table.UV13 = malloc(Nz_molec * sizeof(double));
  my_rates->UVbackground_table.UV14 = malloc(Nz_molec * sizeof(double));
  my_rates->UVbackground_table.UV15 = malloc(Nz_molec * sizeof(double));
  my_rates->UVbackground_table.UV16 = malloc(Nz_molec * sizeof(double));
  my_rates->UVbackground_table.UV17 = malloc(Nz_molec * sizeof(double));
  my_rates->UVbackground_table.UV18 = malloc(Nz_molec * sizeof(double));
  my_rates->UVbackground_table.UV19 = malloc(Nz_molec * sizeof(double));
  my_rates->UVbackground_table.UV20 = malloc(Nz_molec * sizeof(double));
  my_rates->UVbackground_table.UV21 = malloc(Nz_molec * sizeof(double));
  my_rates->UVbackground_table.UV22 = malloc(Nz_molec * sizeof(double));
  my_rates->UVbackground_table.UV23 = malloc(Nz_molec * sizeof(double));
  my_rates->UVbackground_table.UV24 = malloc(Nz_molec * sizeof(double));
  my_rates->UVbackground_table.UV25 = malloc(Nz_molec * sizeof(double));
  my_rates->UVbackground_table.UV26 = malloc(Nz_molec * sizeof(double));
  my_rates->UVbackground_table.UV34 = malloc(Nz_molec * sizeof(double));
  my_rates->UVbackground_table.UV37 = malloc(Nz_molec * sizeof(double));
  my_rates->UVbackground_table.UV38 = malloc(Nz_molec * sizeof(double));

  // Now read everything.

  // *** Redshift ***
  if(! read_dataset(file_id, "z",
                    my_rates->UVbackground_table.z_molec) ) {
    fprintf(stderr, "Error reading dataset 'z' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  // *** UV1 ***
  if(! read_dataset(file_id, "/UV1",
                    my_rates->UVbackground_table.UV1) ) {
    fprintf(stderr, "Error reading dataset '/UV1' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

    // *** UV2 ***
  if(! read_dataset(file_id, "/UV2",
                    my_rates->UVbackground_table.UV2) ) {
    fprintf(stderr, "Error reading dataset '/UV2' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  // *** UV3 ***
  if(! read_dataset(file_id, "/UV3",
                    my_rates->UVbackground_table.UV3) ) {
    fprintf(stderr, "Error reading dataset '/UV3' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  // *** UV4 ***
  if(! read_dataset(file_id, "/UV4",
                    my_rates->UVbackground_table.UV4) ) {
    fprintf(stderr, "Error reading dataset '/UV4' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  // *** UV5 ***
  if(! read_dataset(file_id, "/UV5",
                    my_rates->UVbackground_table.UV5) ) {
    fprintf(stderr, "Error reading dataset '/UV5' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  // *** UV6 ***
  if(! read_dataset(file_id, "/UV6",
                    my_rates->UVbackground_table.UV6) ) {
    fprintf(stderr, "Error reading dataset '/UV6' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  // *** UV7 ***
  if(! read_dataset(file_id, "/UV7",
                    my_rates->UVbackground_table.UV7) ) {
    fprintf(stderr, "Error reading dataset '/UV7' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  // *** UV8 ***
  if(! read_dataset(file_id, "/UV8",
                    my_rates->UVbackground_table.UV8) ) {
    fprintf(stderr, "Error reading dataset '/UV8' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  // *** UV9 ***
  if(! read_dataset(file_id, "/UV9",
                    my_rates->UVbackground_table.UV9) ) {
    fprintf(stderr, "Error reading dataset '/UV9' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  // *** UV10 ***
  if(! read_dataset(file_id, "/UV10",
                    my_rates->UVbackground_table.UV10) ) {
    fprintf(stderr, "Error reading dataset '/UV10' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  // *** UV11 ***
  if(! read_dataset(file_id, "/UV11",
                    my_rates->UVbackground_table.UV11) ) {
    fprintf(stderr, "Error reading dataset '/UV11' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  // *** UV12 ***
  if(! read_dataset(file_id, "/UV12",
                    my_rates->UVbackground_table.UV12) ) {
    fprintf(stderr, "Error reading dataset '/UV12' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  // *** UV13 ***
  if(! read_dataset(file_id, "/UV13",
                    my_rates->UVbackground_table.UV13) ) {
    fprintf(stderr, "Error reading dataset '/UV13' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  // *** UV14 ***
  if(! read_dataset(file_id, "/UV14",
                    my_rates->UVbackground_table.UV14) ) {
    fprintf(stderr, "Error reading dataset '/UV14' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  // *** UV15 ***
  if(! read_dataset(file_id, "/UV15",
                    my_rates->UVbackground_table.UV15) ) {
    fprintf(stderr, "Error reading dataset '/UV15' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  // *** UV16 ***
  if(! read_dataset(file_id, "/UV16",
                    my_rates->UVbackground_table.UV16) ) {
    fprintf(stderr, "Error reading dataset '/UV16' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  // *** UV17 ***
  if(! read_dataset(file_id, "/UV17",
                    my_rates->UVbackground_table.UV17) ) {
    fprintf(stderr, "Error reading dataset '/UV17' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  // *** UV18 ***
  if(! read_dataset(file_id, "/UV18",
                    my_rates->UVbackground_table.UV18) ) {
    fprintf(stderr, "Error reading dataset '/UV18' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  // *** UV19 ***
  if(! read_dataset(file_id, "/UV19",
                    my_rates->UVbackground_table.UV19) ) {
    fprintf(stderr, "Error reading dataset '/UV19' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  // *** UV20 ***
  if(! read_dataset(file_id, "/UV20",
                    my_rates->UVbackground_table.UV20) ) {
    fprintf(stderr, "Error reading dataset '/UV20' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  // *** UV21 ***
  if(! read_dataset(file_id, "/UV21",
                    my_rates->UVbackground_table.UV21) ) {
    fprintf(stderr, "Error reading dataset '/UV21' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  // *** UV22 ***
  if(! read_dataset(file_id, "/UV22",
                    my_rates->UVbackground_table.UV22) ) {
    fprintf(stderr, "Error reading dataset '/UV22' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  // *** UV23 ***
  if(! read_dataset(file_id, "/UV23",
                    my_rates->UVbackground_table.UV23) ) {
    fprintf(stderr, "Error reading dataset '/UV23' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  // *** UV24 ***
  if(! read_dataset(file_id, "/UV24",
                    my_rates->UVbackground_table.UV24) ) {
    fprintf(stderr, "Error reading dataset '/UV24' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  // *** UV25 ***
  if(! read_dataset(file_id, "/UV25",
                    my_rates->UVbackground_table.UV25) ) {
    fprintf(stderr, "Error reading dataset '/UV25' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  // *** UV26 ***
  if(! read_dataset(file_id, "/UV26",
                    my_rates->UVbackground_table.UV26) ) {
    fprintf(stderr, "Error reading dataset '/UV26' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  // *** UV34 ***
  if(! read_dataset(file_id, "/UV34",
                    my_rates->UVbackground_table.UV34) ) {
    fprintf(stderr, "Error reading dataset '/UV34' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  // *** UV37 ***
  if(! read_dataset(file_id, "/UV37",
                    my_rates->UVbackground_table.UV37) ) {
    fprintf(stderr, "Error reading dataset '/UV37' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  // *** UV38 ***
  if(! read_dataset(file_id, "/UV38",
                    my_rates->UVbackground_table.UV38) ) {
    fprintf(stderr, "Error reading dataset '/UV38' in %s.\n",
            my_chemistry->grackle_molecular_data);
    return FAIL;
  }

  H5Fclose(file_id);

  // Get min/max of redshift vector
  my_rates->UVbackground_table.zmin_molec = my_rates->UVbackground_table.z_molec[0];
  my_rates->UVbackground_table.zmax_molec = my_rates->UVbackground_table.z_molec[Nz_molec-1];

  // Print out some information about the dataset just read in.
  if (grackle_verbose) {
    fprintf(stdout, "UV molecular rates information:\n");
    fprintf(stdout, "  z_min = %6.3f\n  z_max = %6.3f\n",
            my_rates->UVbackground_table.zmin_molec,
            my_rates->UVbackground_table.zmax_molec);
  }

  return SUCCESS;
}



int read_dataset(hid_t file_id, char *dset_name, double *buffer) {
  hid_t dset_id;
  herr_t status;
  herr_t h5_error = -1;

  dset_id =  H5Dopen(file_id, dset_name);
  if (dset_id == h5_error) {
    fprintf(stderr, "Failed to open dataset 'z'.\n");
    return FAIL;
  }

  status = H5Dread(dset_id, HDF5_R8, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
  if (status == h5_error) {
    fprintf(stderr, "Failed to read dataset 'z'.\n");
    return FAIL;
  }
 
  H5Dclose(dset_id);

  return SUCCESS;
}
