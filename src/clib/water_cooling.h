// Authors: Shmuel Bialy,
// converted to C++ by Alex Gagliano.
// Date: 03/17/2020
// descript:

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

int find_nearest(double tempx, double tempy, int N, double x[], double y[], double z[]);

double* interp1(int N1, double logT[], double L0[], int N2, double newlogT[]);

double *interp2(int Nx, int Ny, double x[Nx][Ny], double y[Nx][Ny], double z[Nx][Ny], int N2, double tx_arr[], double ty_arr[]);

double *OI_cooling_coef(int N, double n[],double T[],double x_H2,double x_e);

double * LyA_cooling_coef(int N, double T[], double xe);

double *L_Ortho_H2O_rot(int N, double n[], double T[], double Nt[]);

double *L_COrot_NK93(int N, double n[], double T[], double Nt[]);

double *HD_cooling_coef(int N, double n[], double T[]);

double *H2_cooling_coef_Hp_GA08(int N, double T[],int O2P);

double *H2_cooling_coef_Hp_G15(int N, double T[]);

double *H2_cooling_coef_He_GA08(int N, double T[],int O2P);

double *H2_cooling_coef_H2_GA08(int N, double T[],int O2P);

double *H2_cooling_coef_e_GA08(int N, double T[],double O2P);

double *H2_cooling_coef_GA08(int N, double n[],double T[],double xH,double xHe,double xH2,double xHp, double xe, double O2P);

double *CII_cooling_coef(int N, double n[],double T[],double x_H2, double x_e);

double *dust_cooling_coef(int N, double T[], float IUV);


double *CI_cooling_coef(int N, double n[],double T[],double x_H2);
