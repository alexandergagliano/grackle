// Authors: Shmuel Bialy,
// Impemented in C++ by Alex Gagliano.
// Date: 04/27/2020
// descript: A set of routines to calculate the cooling rate
//           contributed by individual metal species at
//           temperature T, density n, and column Nt.
//           For use in the water network, with withWater = 1.
//

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

int find_nearest(double tempx, double tempy, int N, double x[], double y[], double z[]){
    double *idxs= (double *) malloc(2 * sizeof(double));
    double dist[N];
    double minDist = 1.e99;
    int idx = 0;
    for (int i = 0; i < N; i++){
      dist[i] = sqrt(pow(tempx - x[i], 2) + pow(tempy  - y[i], 2));
      if (dist[i] < minDist){
        minDist = dist[i];
        idx = i;
      }
    }
    return idx;
}

// routine to do 1d linear interpolation on an array, based on the values from a previous array
// input:
// N1 is the length of the original arrays (logT and L0)
// N2 is the length of the new array (newlogT)
// outputs:
double* interp1(int N1, double logT[], double L0[], int N2, double newlogT[]){
 		double *newL0 = (double *) malloc(N2 * sizeof(double)); //dynamically allocate memory so we can return it
		double slope, tempL0;
		int j;
		for (int i = 0; i < N2; i++){
			if (newlogT[i] < logT[0]){
				//do linear extrapolation with the first and second point
				slope = (L0[1] - L0[0])/(logT[1] - logT[0]);
				tempL0 = slope*(newlogT[i] - logT[0]) + L0[0];
			} else if (newlogT[i] > logT[N1-1]){
				//do linear extrapolation with the last two points
				slope = (L0[N1-2] - L0[N1-1])/(logT[N1-2] - logT[N1-1]);
				tempL0 = slope*(newlogT[i] - logT[N1-1]) + L0[N1];
			} else {
				//not pretending this is the most efficient way to do things
				j = 0;
				while (newlogT[i] > logT[j]){
					//do linear interpolation with the two nearest points
					//assume the array is sorted
					j++;
				}
				slope = (L0[j-1] - L0[j])/(logT[j-1] - logT[j]);
				tempL0 = slope*(newlogT[i] - logT[j-1]) + L0[j-1];
			}
			newL0[i] = tempL0;
		}
    return newL0;
}

// routine to do 2d linear interpolation on an array, based on the values from a previous array
// input:
// N1 is the length of the original arrays (x,y,z)
// N2 is the length of the new arrays (x,y)
// outputs:
double *interp2(int Nx, int Ny, double x[Nx][Ny], double y[Nx][Ny], double z[Nx][Ny], int N2, double tx_arr[], double ty_arr[]){
    //computes a bilinear interpolation of the data
    //inputs - Nx and Ny, the dimensions of the 2d arrays
    //         x,y,z, the 2d arrays of grid points
    //         N2, the dimension of the 1d array to interpolate
    //         tx_arr, the x vals of the 1d array to interpolate
    //         ty_arr, the y vals of the 1d array
    // returns - z_arr, the bilinear interpolated z vals of the 1d array
    double *z_arr = (double *) malloc(N2 * sizeof(double)); //dynamically allocate memory for 4 integers
    double *tz1;
    double *tz2;
    double *tz;
    int N = Nx * Ny;
    double x_flat[N];
    double y_flat[N];
    double z_flat[N];
    for (int j = 0;j < N; j++){
      int a = floor(j/Ny);
      int b = j%Ny;
      x_flat[j] = x[a][b];
      y_flat[j] = y[a][b];
      z_flat[j] = z[a][b];
    }

    for (int i = 0; i < N2; i++){
        double tx = tx_arr[i];
        double ty = ty_arr[i];
        if ((tx < x[0][0]) || (tx > x[Nx-1][Ny-1]) || (ty < y[0][0]) || (ty > y[Nx-1][Ny-1])){
            z_arr[i] = NAN;
        } else{
            int idx;
            idx = find_nearest(tx, ty, N, x_flat,y_flat,z_flat);
            int idx_x = floor(idx/Nx);
            int idx_y = idx%Ny;
            double p1[3] = {x[idx_x][idx_y], y[idx_x][idx_y], z[idx_x][idx_y]};
            int changex = pow(-1, 1+(p1[0] < tx));
            int changey = pow(-1, 1+(p1[1] < ty));
            int p2x = idx_x;
            int p2y = (idx_y + changex+ Ny)%Ny;
            int p3x = (idx_x + changey + Nx)%Nx;
            int p3y = idx_y;
            int p4x = p3x;
            int p4y = p2y;

            double p2[3] = {x[p2x][p2y], y[p2x][p2y], z[p2x][p2y]};
            double p3[3] = {x[p3x][p3y], y[p3x][p3y], z[p3x][p3y]};
            double p4[3] = {x[p4x][p4y], y[p4x][p4y], z[p4x][p4y]};
            double interpPoint1a[2] = {p1[0], p1[0]};
            double interpPoint1b[2] = {p1[2], p1[2]};
            //hacky fix to ensure our independent variables are strictly increasing
            //otherwise our 1d interpolation won't work!
            //we set our full array as one val and then fill the other one
            //such that they're ordered
            interpPoint1a[(p1[0] < p2[0])] = p2[0];
            interpPoint1b[(p1[0] < p2[0])] = p2[2];
            double interpPoint2a[2] = {p3[0], p3[0]};
            double interpPoint2b[2] = {p3[2], p3[2]};
            interpPoint2a[(p3[0] < p4[0])] = p4[0];
            interpPoint2b[(p3[0] < p4[0])] = p4[2];
            double tx_1val[1] = {tx};
            double ty_1val[1] = {ty};
            tz1 = interp1(2, interpPoint1a, interpPoint1b, 1, tx_1val);
            tz2 = interp1(2, interpPoint2a, interpPoint2b, 1, tx_1val);
            double interpPoint3a[2] = {p1[1], p1[1]};
            double interpPoint3b[2] = {tz1[0], tz1[0]};
            interpPoint3a[(p1[1] < p3[1])] = p3[1];
            interpPoint3b[(p1[1] < p3[1])] = tz2[0];
            tz = interp1(2, interpPoint3a, interpPoint3b, 1, ty_1val);
            z_arr[i] = tz[0];
        }
    }
    return z_arr;
}

//calculates the cooling coefficient of [OI]
//input:
//T is the temperature, valid in the range of 20-10e10. acurrate only for the range 20-1000K
//ni is a row array with the densities of the collisional ellements
//[n-pH2,n-oH2,n-H,n-e]
//nOI is the OI density
//N is the number of density and temperature bins
//output: the cooling coeff in [erg*sec^-1*cm^+3]
double *OI_cooling_coef(int N, double n[],double T[],double x_H2,double x_e){

  double *L = (double *) malloc(N * sizeof(double)); //dynamically allocate memory for the rate

	double OtoP=3;
	double x_para=(1/(OtoP+1))*x_H2;
	double x_ortho=(OtoP/(OtoP+1))*x_H2;
	double x_H=1-2*x_H2;

	//===========we neglect here inelastic collisions with H+==========
	//data:
	double h=6.62606957e-27;
	double k=1.3806488e-16;
	double v10=4744.77749e9;
	double v20=6804.84658e9;
	double v21=2060.06909e9;
	double E10=v10*h;
	double E20=v20*h;
	double E21=v21*h;

	double A10=8.910e-05;
	double A20=1.340e-10;
	double A21=1.750e-05;

	double g0=5.;
	double g1=3.;
	double g2=1.;

	// Using Draine's fitts (p.501)
	double T2[N],q10_H[N],q20_H[N],q21_H[N],
	q10_para[N],q20_para[N],q21_para[N],
	q10_ortho[N],q20_ortho[N],q21_ortho[N],
	T4[N], q10_e[N],q20_e[N],q21_e[N],
	R10[N],R20[N],R21[N],R01[N],R02[N],R12[N];
	double x1_to_x0[N],x2_to_x0[N];
	double x0[N];
	double x1[N];
	double x2[N];

	for (int i = 0; i < N; i++){
		T2[i]=T[i]/1.e2;

		q10_H[i]=3.57e-10*pow(T2[i],(0.419-0.003*log(T2[i])));
		q20_H[i]=3.19e-10*pow(T2[i],(0.369-0.006*log(T2[i])));
		q21_H[i]=4.34e-10*pow(T2[i],(0.755-0.160*log(T2[i])));

		q10_para[i]=1.49e-10*pow(T2[i],(0.264+0.025*log(T2[i])));
		q20_para[i]=1.90e-10*pow(T2[i],(0.203+0.041*log(T2[i])));
		q21_para[i]=2.10e-12*pow(T2[i],(0.889+0.043*log(T2[i])));

		q10_ortho[i]=1.37e-10*pow(T2[i],(0.296+0.043*log(T2[i])));
		q20_ortho[i]=2.23e-10*pow(T2[i],(0.237+0.058*log(T2[i])));
		q21_ortho[i]=3.00e-12*pow(T2[i],(1.198+0.525*log(T2[i])));

		T4[i]=T[i]/1e4;
		q10_e[i]=8.629e-8*pow(T4[i],-0.5)*(0.041*pow(T4[i],0.69)+0.064*pow(T4[i],1.72))/g1;
		q20_e[i]=8.629e-8*pow(T4[i],-0.5)*(0.0136*pow(T4[i],0.61)+0.0186*pow(T4[i],1.49))/g2;
		q21_e[i]=8.629e-8*pow(T4[i],-0.5)*(0.00166*pow(T4[i],0.71)+0.0288*pow(T4[i],1.97))/g2;

		R10[i]=n[i]*(q10_H[i]*x_H+q10_para[i]*x_para+x_ortho*q10_ortho[i]+x_e*q10_e[i])+A10;
		R20[i]=n[i]*(q20_H[i]*x_H+q20_para[i]*x_para+x_ortho*q20_ortho[i]+x_e*q20_e[i])+A20;
		R21[i]=n[i]*(q21_H[i]*x_H+q21_para[i]*x_para+x_ortho*q21_ortho[i]+x_e*q21_e[i])+A21;

	  R01[i]=n[i]*(q10_H[i]*x_H+q10_para[i]*x_para+x_ortho*q10_ortho[i]+x_e*q10_e[i])*(exp(-E10/(k*T[i]))*(g1/g0));
		R02[i]=n[i]*(q20_H[i]*x_H+q20_para[i]*x_para+x_ortho*q20_ortho[i]+x_e*q20_e[i])*(exp(-E20/(k*T[i]))*(g2/g0));
		R12[i]=n[i]*(q21_H[i]*x_H+q21_para[i]*x_para+x_ortho*q21_ortho[i]+x_e*q21_e[i])*(exp(-E21/(k*T[i]))*(g2/g1));

		x1_to_x0[i]=(R01[i]*R20[i]+R01[i]*R21[i]+R21[i]*R02[i])/(R10[i]*R20[i]+R10[i]*R21[i]+R12[i]*R20[i]);
		x2_to_x0[i]=(R02[i]*R10[i]+R02[i]*R12[i]+R12[i]*R01[i])/(R10[i]*R20[i]+R10[i]*R21[i]+R12[i]*R20[i]);
		x0[i]=1./(1+x1_to_x0[i]+x2_to_x0[i]);

		x1[i]=x1_to_x0[i]*x0[i];
		x2[i]=x2_to_x0[i]*x0[i];

		L[i]=(x1[i]*A10*E10+x2[i]*(A21*E21+A20*E20))/n[i];
	}
  return L;
}

//calculates the cooling coefficient of LymanAlfa line assuming n<n_crit (~1e17 cm^-3)
//input:
// T is the tempreture - in the range of 10 to 1e10.
// n - the total H density
//output: LL = the cooling function coefficient [erg*sec^-1*cm^-3]

// The cooling rate per unit volume is L = n^2 * xH * LL
// assuming n << n_crit
double * LyA_cooling_coef(int N, double T[], double xe){

  double *LL = (double *) malloc(N * sizeof(double)); //dynamically allocate memory for the rate

	// constants
	double hv = 1.6366e-11;  //in erg
	double k = 1.3806488e-16; // in erg K^-1

	// 2-1 rate coefficient
	//  rate coefficient  is a fit to both 2s-1, and 2p-1 excitations (Spitzer
	//  book), see the fitting in "Ly_alpha_col_rate_coef.m"

	double q21[N];
	for (int i = 0; i < N; i++){
		q21[i] = 1.70e-6/sqrt(T[i]);
	}

	// 1-2 rate coefficient
	double g2_to_g1=4;
	double q12[N];
	for (int i = 0; i < N; i++){
		q12[i] = g2_to_g1 * q21[i] * exp(-hv/(k*T[i]));
		LL[i] = xe * q12[i] * hv;
	}
  return LL;
}



double *L_Ortho_H2O_rot(int N, double n[], double T[], double Nt[]){
  double *L = (double *) malloc(N * sizeof(double)); //dynamically allocate memory for 4 integers
  int Nx = 10;
  int Ny = 11;
  double logNt[N];
  double logn[N];
  double logT[N];
  double *logL0;
  double L0[N];

  double logT_[] = {1.0000, 1.3010, 1.4771, 1.6990, 1.9031, 2.0000,
      2.3010, 2.6021, 3.0000, 3.3010, 3.6021};
  double logNt_[] = {10.0000, 11.0000, 12.0000, 13.0000, 14.0000, 15.0000,
      16.0000, 17.0000, 18.0000, 19.0000};
  double logL0_[] =
      {-26.8133,-25.8779,-25.4274,-24.9556,-24.5820,-24.4095,
      -23.9075,-23.4457,-22.8830,-22.4955,-22.1327};

  for (int i = 0; i < N; i++){
    logn[i] = log10(n[i]);
    logT[i] = log10(T[i]);
    logNt[i] = log10(Nt[i]);
  }

  logL0 = interp1(Ny, logT_, logL0_, N, logT);

  for (int i = 0; i < N; i++){
    L0[i]  = pow(10, logL0[i]);
  }


  double logT_mesh[10][11];
  double logNt_mesh[10][11];

  for (int i = 0; i < 11; i++){
    for (int j = 0; j < 10; j++){
      logT_mesh[j][i] = logT_[i];
      logNt_mesh[j][i] = logNt_[j];
    }
  }

    double logLTE_[10][11] = {{-17.943 , -16.7079, -16.0839, -15.413 , -14.8546, -14.601 ,
        -13.855 , -13.1585, -12.3099, -11.8499, -11.6265},
       {-17.9647, -16.724 , -16.0947, -15.4183, -14.8569, -14.6026,
        -13.8555, -13.1587, -12.3099, -11.8499, -11.6265},
       {-18.1372, -16.8552, -16.1854, -15.466 , -14.8784, -14.6174,
        -13.8604, -13.1604, -12.3103, -11.8501, -11.6267},
       {-18.7671, -17.3607, -16.5781, -15.7181, -15.025 , -14.7287,
        -13.9041, -13.1763, -12.3145, -11.8522, -11.6281},
       {-19.6995, -18.1087, -17.2472, -16.2714, -15.468 , -15.1142,
        -14.1295, -13.2889, -12.3515, -11.8723, -11.6415},
       {-20.6659, -18.9338, -18.0505, -17.0116, -16.1183, -15.7238,
        -14.6107, -13.6399, -12.5498, -12.0183, -11.7515},
       {-21.5792, -19.8109, -18.9124, -17.8186, -16.877 , -16.4531,
        -15.2414, -14.1789, -13.0024, -12.486 , -12.1796},
       {-22.5254, -20.7085, -19.8006, -18.6861, -17.6987, -17.2497,
        -15.957 , -14.8195, -13.6346, -13.1447, -12.8265},
       {-23.5005, -21.641 , -20.7233, -19.5788, -18.5638, -18.0994,
        -16.7349, -15.5296, -14.3564, -13.8718, -13.5436},
       {-24.409 , -22.5793, -21.6545, -20.4941, -19.4566, -18.9882,
        -17.5593, -16.2934, -15.1454, -14.6656, -14.3246}};

    double lognc_[10][11] = {{8.8146, 8.9046, 9.0329, 9.1395, 9.1941, 9.2045, 9.2142, 9.2875,
        9.5279, 9.641 , 9.5598},
       {8.786 , 8.8834, 9.0142, 9.1261, 9.2018, 9.1987, 9.2115, 9.2861,
        9.5275, 9.6407, 9.5597},
       {8.5971, 8.728 , 8.879 , 9.0306, 9.1277, 9.1546, 9.194 , 9.2728,
        9.5233, 9.6388, 9.5583},
       {7.9562, 8.1364, 8.313 , 8.5537, 8.7818, 8.8484, 9.0048, 9.1722,
        9.4851, 9.6198, 9.545 },
       {7.0116, 7.2052, 7.3961, 7.6862, 7.9962, 8.1062, 8.3807, 8.6956,
        9.2521, 9.4831, 9.4483},
       {6.0247, 6.2048, 6.4102, 6.7077, 7.0394, 7.1694, 7.4782, 7.8482,
        8.5825, 8.9688, 9.0183},
       {5.034 , 5.2088, 5.4174, 5.7028, 6.0567, 6.19  , 6.4966, 6.884 ,
        7.6846, 8.1364, 8.2468},
       {4.0185, 4.2137, 4.4058, 4.7068, 5.0519, 5.1894, 5.5405, 5.931 ,
        6.7523, 7.234 , 7.3716},
       {3.0239, 3.2201, 3.41  , 3.7112, 4.0584, 4.1951, 4.5393, 4.9793,
        5.8272, 6.343 , 6.4944},
       {2.0317, 2.208 , 2.4158, 2.704 , 3.0633, 3.2001, 3.5415, 3.9836,
        4.906 , 5.4459, 5.6036}};


    double a_[10][11] = {{0.7050, 0.4860, 0.4760, 0.4550, 0.4520, 0.4540, 0.4410, 0.4050,
          0.3780, 0.3500, 0.3400},
        {0.6430, 0.4870, 0.4750, 0.4540, 0.4500, 0.4530, 0.4400, 0.4050,
          0.3780, 0.3500, 0.3400},
        {0.5660, 0.4970, 0.4740, 0.4460, 0.4390, 0.4410, 0.4330, 0.4000,
          0.3750, 0.3480, 0.3390},
        {0.5590, 0.5330, 0.4780, 0.4410, 0.4180, 0.4140, 0.4040, 0.3770,
          0.3590, 0.3380, 0.3330},
        {0.5890, 0.5970, 0.5270, 0.4670,0.4470, 0.4290, 0.3870, 0.3480,
          0.3340, 0.3240, 0.3110},
        {0.7160, 0.6430, 0.5750, 0.5170, 0.4880,0.4700, 0.4010, 0.3530,
          0.3240, 0.3030, 0.2950},
        {0.8450, 0.6790, 0.6110, 0.5580, 0.5230, 0.5000, 0.4220, 0.3570,
          0.3250, 0.2890, 0.2910},
        {0.8570,0.7000, 0.6320, 0.5800, 0.5500, 0.5280, 0.4430, 0.3710,
          0.3170, 0.2760, 0.2770},
        {0.9270, 0.7170,0.6550, 0.6050, 0.5720, 0.5480, 0.4580, 0.3860,
          0.3100, 0.2640, 0.2670},
        {0.8700, 0.7270, 0.6720, 0.6220, 0.5920, 0.5640, 0.4750, 0.3980,
          0.3090, 0.2570, 0.2580}};

  double *lognc;
  double *a;
  double nc[N];
  double L_LTE[N];
  lognc = interp2(Nx, Ny, logT_mesh, logNt_mesh, lognc_, N, logT, logNt);
  for (int i = 0; i < N; i++){
    nc[i] = pow(10, lognc[i]);
  }
  a = interp2(Nx, Ny, logT_mesh, logNt_mesh, a_, N, logT, logNt);
  double *logL_LTE;
  logL_LTE = interp2(Nx, Ny, logT_mesh, logNt_mesh, logLTE_, N, logT, logNt);

  for (int i = 0; i < N; i++){
    L_LTE[i]   = pow(10., logL_LTE[i]);
    L[i] = pow(pow(L0[i], -1) + n[i]/L_LTE[i] + pow(L0[i], -1)*pow(n[i]/nc[i], a[i])*(1 - nc[i]*L0[i]/L_LTE[i]), -1);
  }

  free(lognc);
  free(a);
  free(logL_LTE);

	return L;
}

double *L_COrot_NK93(int N, double n[], double T[], double Nt[]){
  double *L = (double *) malloc(N * sizeof(double)); //dynamically allocate memory for 4 integers
  int Nx = 10;
  int Ny = 11;

  // based on David Neufeld and Michael Kaufmann (1993)
  // INPUT
  // n is the collider volume density (cm^-3)
  // T is the gas temperature (K)
  // N_tilde is the column density (cm^-2) per km/s of Hydrogen nuclei
  // T should be either a vector or a scalar.
  // n and N_tilde should either be scalars or vectors the size of T

  //if length(n) == 1
  //    n = n * ones(size(T));
  //end
  //if length(Nt) == 1
  //    Nt = Nt * ones(size(T));
  //end

  // OUTPUT
  // the cooling coefficient (erg cm^3 s^-1)

  // data
  double logNt[N];
  double logn[N];
  double logT[N];
  double *logL0;
  double L0[N];

  double logT_[] = {1.0000, 1.3010, 1.4771, 1.6990, 1.9031, 2.0000,
          2.4771, 2.7782, 3.0000, 3.1761, 3.3010};
  double logNt_[] =
         {14.5000,15.0000,15.5000,16.0000,16.5000,17.0000,
         17.5000,18.0000,18.5000,19.0000};
  double logL0_[] = {-24.9553, -24.5569, -24.356 , -24.1165, -23.8846, -23.8013,
         -23.4032, -23.0729, -22.8126, -22.6073, -22.4738};

  for (int i = 0; i < N; i++){
     logn[i] = log10(n[i]);
     logT[i] = log10(T[i]);
     logNt[i] = log10(Nt[i]);
  }

  logL0 = interp1(Ny, logT_, logL0_, N, logT);

  for (int i = 0; i < N; i++){
     L0[i]  = pow(10, logL0[i]);
  }


  double logT_mesh[10][11];
  double logNt_mesh[10][11];

   for (int i = 0; i < 11; i++){
       for (int j = 0; j < 10; j++){
           logT_mesh[j][i] = logT_[i];
           logNt_mesh[j][i] = logNt_[j];
         }
     }

  double logLTE_[10][11] = {{-21.1462, -20.3891, -19.9734, -19.4701, -19.0225, -18.8142, -17.8243, -17.2328,
           -16.822 , -16.5203, -16.3264},
         {-21.25, -20.4489, -20.0154, -19.4963, -19.0392, -18.8277, -17.8287, -17.235 ,
           -16.8233, -16.5211, -16.327},
         {-21.4471, -20.5821, -20.117 , -19.5655, -19.0861, -18.8663, -17.8425, -17.2419,
           -16.8273, -16.5237, -16.3289},
         {-21.7345, -20.8077, -20.3065, -19.7114, -19.1956, -18.9604, -17.8815, -17.2624,
           -16.8397, -16.5318, -16.3349},
         {-22.0842, -21.1097, -20.5793, -19.9452, -19.3916, -19.1382, -17.9762, -17.318 ,
          -16.8751, -16.5558, -16.3529},
         {-22.473 , -21.4622, -20.9103, -20.2476, -19.6653, -19.3972, -18.1535, -17.4413,
           -16.9626, -16.6196, -16.4028},
         {-22.8873, -21.847 , -21.2788, -20.5953, -19.9926, -19.7142, -18.41, -17.6494,
           -17.131 , -16.7558, -16.5173},
         {-23.3163, -22.2534, -21.6722, -20.9727, -20.355 , -20.0691, -18.7232, -17.9287,
           -17.3799, -16.9773, -16.7188},
         {-23.7607, -22.6752, -22.083 , -21.3704, -20.741 , -20.4496, -19.0737, -18.2564,
           -17.6877, -17.2677, -16.9979},
         {-24.2121, -23.1083, -22.5066, -21.783 , -21.1441, -20.8482, -19.4497, -18.6163,
           -18.0346, -17.6052, -17.3338}};

  double lognc_[10][11] = {{3.3116, 3.7168, 3.9577, 4.2410, 4.3926, 4.5149, 4.9462, 5.2445,
          5.4266, 5.5730, 5.6778},
          {3.1147, 3.5864, 3.8652, 4.1653, 4.3542, 4.4671, 4.9263, 5.2340,
          5.4206, 5.5692, 5.6750},
          {2.7964, 3.3351, 3.6439, 4.0065, 4.2329, 4.3634, 4.8831, 5.2030,
          5.4168, 5.5572, 5.6815},
          {2.3659, 2.9467, 3.3111, 3.7333, 3.9902, 4.1308, 4.7536, 5.1330,
          5.3607, 5.5354, 5.6526},
          {1.8982, 2.5016, 2.8742, 3.3324, 3.6089, 3.7829, 4.5023, 4.9539,
          5.2388, 5.4361, 5.5755},
          {1.3965, 2.0233, 2.4039, 2.8620, 3.1680, 3.3420, 4.1179, 4.6470,
          4.9971, 5.2489, 5.4130},
          {0.9008, 1.5183, 1.9023, 2.3783, 2.6913, 2.8668, 3.6742, 4.2478,
          4.6324, 4.9184, 5.1323},
          {0.4036, 1.0215, 1.4063, 1.8856, 2.1867, 2.3765, 3.1962, 3.7738,
          4.1802, 4.5055, 4.7296},
          {-0.0941, 0.5236, 0.9087, 1.3900, 1.6902, 1.8811, 2.6930, 3.2881,
          3.7032, 4.0255, 4.2613},
          {-0.5915, 0.0258, 0.4109, 0.8786, 1.1926, 1.3843, 2.1962, 2.7944,
          3.2014, 3.5375, 3.7765}};

  double a_[10][11] =
         {{0.4220, 0.3880, 0.3780, 0.3890, 0.4180, 0.4310, 0.3780, 0.3750,
           0.3690, 0.3590, 0.3460},
          {0.4130, 0.3770, 0.3660, 0.3790, 0.4080, 0.4220, 0.3720, 0.3710,
            0.3670, 0.3570, 0.3450},
          {0.4180, 0.3750, 0.3620, 0.3710, 0.3970, 0.4100, 0.3610, 0.3630,
            0.3600, 0.3520, 0.3410},
          {0.4360, 0.3860, 0.3720, 0.3750, 0.3940, 0.4050, 0.3480, 0.3500,
            0.3490, 0.3420, 0.3320},
          {0.4600, 0.4260, 0.4090, 0.4050, 0.4060, 0.4120, 0.3510, 0.3390,
            0.3350, 0.3290, 0.3190},
          {0.4900, 0.4650, 0.4500, 0.4440, 0.4480, 0.4520, 0.3680, 0.3450,
            0.3290, 0.3190, 0.3080},
          {0.5210, 0.4950, 0.4830, 0.4840, 0.4930, 0.4900, 0.3980, 0.3730,
            0.3500, 0.3220, 0.3070},
          {0.5480, 0.5250, 0.5180, 0.5190, 0.5250, 0.5210, 0.4300, 0.4020,
            0.3780, 0.3540, 0.3260},
          {0.5660, 0.5480, 0.5450, 0.5540, 0.5480, 0.5440, 0.4580, 0.4330,
            0.4100, 0.3830, 0.3570},
          {0.5830, 0.5680, 0.5650, 0.5750, 0.5660, 0.5610, 0.4810, 0.4620,
            0.4380, 0.4140, 0.3870}};

  double *lognc;
  double *a;
  double nc[N];
  double L_LTE[N];
  lognc = interp2(Nx, Ny, logT_mesh, logNt_mesh, lognc_, N, logT, logNt);
  for (int i = 0; i < N; i++){
    nc[i] = pow(10, lognc[i]);
  }
  a = interp2(Nx, Ny, logT_mesh, logNt_mesh, a_, N, logT, logNt);
  double *logL_LTE;
  logL_LTE = interp2(Nx, Ny, logT_mesh, logNt_mesh, logLTE_, N, logT, logNt);

  for (int i = 0; i < N; i++){
    L_LTE[i]   = pow(10., logL_LTE[i]);
    L[i] = pow(pow(L0[i], -1) + n[i]/L_LTE[i] + pow(L0[i], -1)*pow(n[i]/nc[i], a[i])*(1 - nc[i]*L0[i]/L_LTE[i]), -1);
  }

  free(lognc);
  free(a);
  free(logL_LTE);

  return L;
}


//Following Anton Lipovka, Ramona Nunez-Lopez, and Vladimir Avila-Reese 2005
//http://arxiv.org/pdf/astro-ph/0503682.pdf
//the cooling function is writen there but I had to add a correction to the
//normalization: log(W)=log10(1/2.5)+...

//n_ and T_ must be vectors of the same lengths.
double *HD_cooling_coef(int N, double n[], double T[]) {
  // N - number of temperature and density bins
  // M - number of cooling lines
  int M = 5;
  double *LL = (double *) malloc(N * sizeof(double)); //dynamically allocate memory for the rate

  double a;
  double log_n, log_T;

  double D[5][5] = {{-42.57688, 0.92433, 0.54962, -0.07676, 0.00275},
    {21.93385, 0.77952, -1.06447, 0.11864, -0.00366},
    {-10.19097, -0.54263, 0.62343, -0.07366, 0.002514},
    {2.19906, 0.11711, -0.13768, 0.01759, -0.000666317},
    {-0.17334, -0.00835, 0.0106, -0.001482, 0.000061926}};

  for (int i=0; i<N; i++){
      a=1.;
      if (n[i]<1.) {// %under n=1 the cooling coefficient L(erg s^-1 cm^3) is constant with n
          a=1.;
          n[i]=1.;
      } else if (n[i]>1.e8) { // %above n=1e8 the cooling coefficient L(erg s^-1 cm^3) is proportional to 1/n
          a=1.e8/n[i];
          n[i]=1.e8;
      }
      log_n =log10(n[i]);
      log_T =log10(T[i]);

      double logT_arr[]={1.,log_T,pow(log_T,2),pow(log_T,3),pow(log_T,4)};
      double logn_arr[]={1.,log_n,pow(log_n,2),pow(log_n,3),pow(log_n,4)};

      double val = 1.;
      double row[M];
      double sum = 0;
      for (int k=0; k<M;k++){
        row[k] = 0;
        for (int j=0; j<M;j++){
          row[k] += D[k][j]*logn_arr[j];
        }
        sum+= logT_arr[k]*row[k];
      }
      LL[i]=a*pow(10.,sum)/n[i];
  }
  return LL;
}


double *H2_cooling_coef_Hp_GA08(int N, double T[],int O2P){
  double *Lamb = (double *) malloc(N * sizeof(double)); //dynamically allocate memory for the rate
	//calculate the ortho-H2 cooling coefficient (erg cm^3/s) from collisions by H only.
  int p[6];
  for (int i=0; i<6;i++){
    p[i] = i;
  }

  double logT3[N];
  for (int i = 0; i < N; i++){
      logT3[i]=log10(T[i]/1.e3);
  }

  double logT3_mesh[6][N];
  double p_mesh[6][N];
  for (int i = 0; i < 6; i++){
      for (int j = 0; j < N; j++){
        logT3_mesh[i][j] = logT3[j];
        p_mesh[i][j] = p[i];
    }
  }

  // para
  double a_p[]={-21.757160,	1.3998367,
  -0.37209530,	0.061554519,
  -0.37238286,	0.23314157};

  // ortho
  double a_o[]={-21.706641,	1.3901283,
  -0.34993699,	0.075402398,
  -0.23170723,	0.068938876};

 double Wp[N];
 double Wo[N];
 // total
 double xo=(double) O2P/(1+ (double) O2P);
 double xp=1./(1+ (double) O2P);

 //double loop maybe?
 for (int j = 0; j < N; j++){
   Wp[j] = 0;
   Wo[j] = 0;
   for (int i = 0; i < 6; i++){
     //funky indexing to deal with the fact that we're doing 1d vectors only
     Wp[j] +=  a_p[i]*pow(logT3_mesh[i][j],p_mesh[i][j]);
     Wo[j] +=  a_o[i]*pow(logT3_mesh[i][j],p_mesh[i][j]);
    }
    Wp[j] = pow(10.,Wp[j]);
    Wo[j] = pow(10.,Wo[j]);
    Lamb[j] = xo*Wo[j]+xp*Wp[j];
    if (Lamb[j] < 0){
      Lamb[j] = 0;
    }
  }

  return Lamb;
}


double *H2_cooling_coef_Hp_G15(int N, double T[]){
//updated as for Glover 2015
//%calculate the ortho-H2 cooling coefficient (erg cm^3/s) from collisions by H only.
// assuming O2P=3

	double T3[N];
  double logT3[N];
  double *Lamb = (double *) malloc(N * sizeof(double)); //dynamically allocate memory so we can return it

  for (int i = 0; i < N; i++){
    T3[i] = T[i]/1.e3;
    logT3[i] = log10(T3[i]);
    Lamb[i] =pow(10., -22.089523
	    + 1.571471100 * logT3[i]
	    + 0.015391166 * pow(logT3[i], 2)
	    -0.2361998500 * pow(logT3[i], 3)
	    -0.5100222100 * pow(logT3[i], 4)
	    +0.3216873000 * pow(logT3[i], 5));
    }
    return Lamb;
}

double *H2_cooling_coef_H_GA08(int N, double T[],int O2P){
	//calculate the ortho-H2 cooling coefficient (erg cm^3/s) from collisions by H only.
	double xp=1./((double) O2P+1);
	double xo= (double) O2P/( (double) O2P+1);

  double *Lamb = (double *) malloc(N * sizeof(double)); //dynamically allocate memory for the rate


	double W_oH2[N];
  double T3[N];
  for (int i = 0; i < N; i++){
	   T3[i]=T[i]/1.e3;
  }


  //double W_oH2[N];
  for (int i = 0; i < N; i++){
    W_oH2[i] = 0.0;
  }

  for (int i = 0; i < N; i++){
    //I1=T<100;
    if (T[i] < 100){
	     W_oH2[i]=5.09e-27*pow(T3[i],0.5)*exp(-825.5/T[i]);
   }
  }

	double a[8][2]={{-24.330855,  -24.329086},
	    {4.4404496,   4.6105087},
	   {-4.0460989,  -3.9505350},
	   {-1.1390725,   12.363818},
	    {9.8094223,  -32.403165},
	    {8.6273872,   48.853562},
	    {0.0,        -38.542008},
	    {0.0,         12.066770}};


  //consolidate into far fewer for loops later!
  //explore a potential indexing problem here - why do we only go to eight?
  //logW is a vector here - does it need to be?
  //redo this tomorros
	double logW[N];
  for (int j =0; j<N; j++){
    logW[j] = 0;
  	for (int i=0; i<8; i++){
      //	I2=T>=100 & T<1000;
      if ((T[i] >= 1.e2) && (T[i] < 1.e3)){
  	     logW[j]+= a[i][0]*pow(log10(T3[j]),i);
        }
  	}
}

  for (int i=0; i<N; i++){
    //	I2=T>=100 & T<1000;
    if ((T[i] >= 1.e2) && (T[i] < 1.e3)){
	     W_oH2[i]=pow(10.,logW[i]);
     }
  }

  for (int j =0; j<N; j++){
  	logW[j] = 0;
  	for (int i = 0; i <8; i++){
      //	I3=T>=1000 & T<=1e4;
        if ((T[i]>=1.e3) && (T[i] <= 1.e4)){
  	       logW[j] += a[i][1]*pow(log10(T3[j]),i);
        }
  	}
  }

  for (int i = 0; i <N; i++){
    //	I3=T>=1000 & T<=1e4;
      if ((T[i]>=1.e3) && (T[i] <= 1.e4)){
         W_oH2[i]=pow(10.,logW[i]);
      }
  }


//	I4=T>1e4;
  for (int i = 0; i <N; i++){
    if (T[i] > 1.e4){
	     W_oH2[i]=4.6759e-22;
     }
   }
	//calculate the para-H2 cooling function (erg/s) from collisions by H only.
	double W_pH2[N];
  for (int i = 0; i < N; i++){
	   T3[i]=T[i]/1.e3;
   }

  for(int i = 0; i < N; i++){
    if (T[i] < 1.e2){
	     W_pH2[i]=8.16e-26*pow(T3[i],0.5)*exp(-509.85/T[i]);
     }
   }

	double a2[8][2]={{-24.216387, -24.216387},
	    {3.3237480,  4.2046488},
	   {-11.642384, -1.3155285},
	   {-35.553366, -1.6552763},
	   {-35.105689,  4.1780102},
	   {-10.922078, -0.56949697},
	   {0.0,       -3.3824407},
	   {0.0,        1.0904027}};

  for (int j =0; j<N; j++){
    logW[j] = 0;
  	for (int i=0; i<8; i++){
        if ((T[i] >= 1.e2) && (T[i] < 1.e3)){
  	       logW[j] += a2[i][0]*pow(log10(T3[j]),i);
         }
  	}
  }

  for (int i=0; i<N; i++){
      if ((T[i] >= 1.e2) && (T[i] < 1.e3)){
        W_pH2[i]=pow(10.,logW[i]);
      }
  }

  for (int j =0; j<N; j++){
    logW[j] = 0;
  	for (int i=0; i<8; i++){
      	//I3=T>=1000 & T<=1e4;
        if ((T[i] >= 1.e3) && (T[i] <= 1.e4)){
  	       logW[j] += a2[i][1]*pow(log10(T3[j]),i);
        }
      }
  }

  for (int i=0; i<N; i++){
    	//I3=T>=1000 & T<=1e4;
      if ((T[i] >= 1.e3) && (T[i] <= 1.e4)){
	       W_pH2[i]=pow(10.,logW[i]);
       }
  }


  //	I4=T>1e4;
  for (int i=0; i<N; i++){
    if (T[i]> 1.e4){
	     W_pH2[i]=2.1574e-22;
    }
  }

  for (int i=0; i<N; i++){
	   Lamb[i]=xo*W_oH2[i]+xp*W_pH2[i];
   }

   return Lamb;
}

double *H2_cooling_coef_He_GA08(int N, double T[],int O2P){
  double *Lamb = (double *) malloc(N * sizeof(double)); //dynamically allocate memory for the rate


	//calculate the ortho-H2 cooling coefficient (erg cm^3/s) from collisions by H only.
  int p[6];
	for (int i=0; i<6;i++){
		p[i] = i;
	}

	double logT3[N];

  for (int i = 0; i < N; i++){
    logT3[i] = log10(T[i]/1.e3);
  }

	//[logT3,p]=meshgrid(logT3,p);
  double logT3_mesh[6][N];
  double p_mesh[6][N];
  for (int i = 0; i < 6; i++){
      for (int j = 0; j < N; j++){
        logT3_mesh[i][j] = logT3[j];
        p_mesh[i][j] = p[i];
    }
  }

	// para
	double a_p[6]={-23.489029,
	1.8210825,
	-0.59110559,
	0.42280623,
	-0.30171138,
	0.12872839};

	// ortho
	double a_o[6]={-23.7749,
	2.40654,
	-1.23449,
	0.739874,
	-0.258940,
	0.120573};

  double Wp[N];
  double Wo[N];
  //double loop maybe?
  for (int j = 0; j < N; j++){
    Wp[j] = 0;
    Wo[j] = 0;
    for (int i = 0; i < 6; i++){
    //funky indexing to deal with the fact that we're doing 1d vectors only
    Wp[j] +=  a_p[i]*pow(logT3_mesh[i][j],p_mesh[i][j]);
    Wo[j] +=  a_o[i]*pow(logT3_mesh[i][j],p_mesh[i][j]);
   }
   Wp[j] = pow(10,Wp[j]);
   Wo[j] = pow(10,Wo[j]);
 }
  // total
  double xo=(double) O2P/(1+(double) O2P);
  double xp=1./(1+ (double) O2P);

  for (int i = 0; i < N; i++){
    Lamb[i]=xo*Wo[i]+xp*Wp[i];
  }

  return Lamb;

}

double *H2_cooling_coef_H2_GA08(int N, double T[],int O2P){
	//calculate the H2 cooling coefficient (erg cm^3/s) from collisions with H2 (para and ortho, ground state) only.
  double *Lamb = (double *) malloc(N * sizeof(double)); //dynamically allocate memory for the rate

  double T3[N];
  double logT3[N];

  for (int i = 0; i < N; i++){
    T3[i] = T[i]/1.e3;
    logT3[i] = log10(T3[i]);
  }

	//if nargin<2
	//    O2P=3;
	//end
	// para-H2 excited by para-H2
	double a_pp[6]={-23.889798,
	1.8550774,
	-0.55593388,
	0.28429361,
	-0.20581113,
	0.13112378};

	double logW[N];
  for (int j = 0; j < N; j++){
    logW[j] = 0;
    for (int i=0; i<6; i++){
  	    logW[j] += a_pp[i]*pow(logT3[j],i);
  	}
}

	double W_pH2_pH2[N];
  for (int i = 0; i < N; i++){
    W_pH2_pH2[i]=pow(10.,logW[i]);
  }

	// para-H2 excited by ortho-H2
	double a_po[6]={-23.748534,
	1.76676480,
	-0.58634325,
	0.31074159,
	-0.17455629,
	0.18530758};

  for (int j = 0; j < N; j++){
    logW[j] = 0;
   for (int i=0; i<6; i++){
     //double check that this is getting indexed right
  	    logW[j] += a_po[i]*pow(logT3[j], i);
  	}
  }

  double W_pH2_oH2[N];
  for (int i = 0; i < N; i++){
	   W_pH2_oH2[i]=pow(10.,logW[i]);
  }

	// ortho-H2 excited by para-H2
	double a_op[6]={-24.126177,
	2.3258217,
	-1.0082491,
	0.54823768,
	-0.33679759,
	0.20771406};

  for (int j = 0; j < N; j++){
    logW[j] = 0;
  	for (int i=0; i < 6; i++){
  	    logW[j] += a_op[i] * pow(logT3[j], (i));
  	}
  }

  double W_oH2_pH2[N];
  for (int i = 0; i < N; i++){
	   W_oH2_pH2[i]=pow(10.,logW[i]);
  }

	// ortho-H2 excited by ortho-H2
	double a_oo[6]={-24.020047,
	2.2687566,
	-1.0200304,
	0.83561432,
	-0.40772247,
	0.096025713};

  for (int j = 0; j < N; j++){
    logW[j] = 0;
    for (int i = 0; i < 6; i++){
  	    logW[j] += a_oo[i]*pow(logT3[j], (i));
  	}
  }
  double W_oH2_oH2[N];
  for (int i = 0; i < N; i++){
	   W_oH2_oH2[i]=pow(10., logW[i]);
  }

	// total H2 cooling (p and o) from collisions with p/o-H2
	double xp=1./((double)O2P+1);
	double xo=(double) O2P/( (double) O2P+1);


  for (int i = 0; i < N; i++){
    	Lamb[i]= pow(xp,2.) *W_pH2_pH2[i]  +xo*xp *W_oH2_pH2[i]+
	     xo*xp*W_pH2_oH2[i] +pow(xo,2.)*W_oH2_oH2[i];
     }

  return Lamb;
}

double *H2_cooling_coef_e_GA08(int N, double T[],double O2P){
	//calculate the ortho-H2 cooling coefficient (erg cm^3/s) from collisions by H only.
  double *Lamb = (double *) malloc(N * sizeof(double)); //dynamically allocate memory for the rate

  double Wp[N];
 //Wp=zeros(size(T));

	double p[6];
	for (int i = 0; i<6; i++){
			p[i]=i;
	}

  double logT3[N];
  for (int i = 0; i < N; i++){
	   logT3[i]=log10(T[i]/1e3);
  }

  double logT3_mesh[6][N];
  double p_mesh[6][N];
  for (int i = 0; i < 6; i++){
      for (int j = 0; j < N; j++){
        logT3_mesh[i][j] = logT3[j];
        p_mesh[i][j] = p[i];
    }
  }

	//[logT3,p]=meshgrid(logT3,p);
	//double a[];

	// para below 1e3K
	//I=(T<=1e3);
	double a_pL[6]={-22.817869,
	    0.95653474,
	    0.79283462,
	    0.56811779,
	    0.27895033,
	    0.056049813};

  double Ep=509.85;  // in K

  // para above 1e3K
  double a_pH[6]={-22.817869,
          0.66916141,
          7.1191428,
          -11.176835,
          7.0467275,
          -1.6471816};

      for (int i = 0; i < N; i++){
        Wp[i] = 0;
        if (T[i] <= 1.e3){
          for (int j = 0; j < 6; j++){
            Wp[i]+=a_pL[j]*pow(logT3_mesh[j][i],p_mesh[j][i]);
        }
      } else{
        for (int j = 0; j < 6; j++){
          Wp[i]+=a_pH[j]*pow(logT3_mesh[j][i],p_mesh[j][i]);
        }
      }
      Wp[i] = exp(-Ep/T[i])*pow(10.,Wp[i]);
    }

	// ortho
	double a_o[6]={-21.703215,
	    0.76059565,
	    0.50644890,
	    0.050371349,
	    -0.10372467,
	    -0.035709409};

  double Eo=845; // % in K
  double Wo[N];
  for (int i = 0; i < N; i++){
    Wo[i] = 0;
    for (int j = 0; j < 6; j++){
      Wo[i]+=a_o[j]*pow(logT3_mesh[j][i],p_mesh[j][i]);
    }
    Wo[i] = exp(-Eo/T[i])*pow(10.,Wo[i]);
  }

	// total
	double xo=(double) O2P/(1+ (double) O2P);
	double xp=1./(1+ (double) O2P);
  for (int i = 0; i < N; i++){
	   Lamb[i]=xo*Wo[i]+xp*Wp[i];
  }

  return Lamb;
}


//Convert to C!
double *H2_LTE(int N, double T[]){
  // returns the H2 cooling function per particle (erg s^-1)
  // and the level populations at LTE
    double *L = (double *) malloc(N * sizeof(double)); //dynamically allocate memory for the rate

    //first one is 24, second one is 27
    double O2P=3;
    double xO=O2P/(1+O2P);
    double xP=1./(1+O2P);
    FILE *E_fn;
    FILE *k_fn;
    FILE *A_fn;
    int twoCases[2] = {24, 27};
    double X;
    for (int i=0; i < 2; i++){
        int itot = twoCases[i];
        int jtot = twoCases[i];
        float fValues[36450];
        float E[itot][jtot];
        float k[itot][jtot*50];
        float A[itot][jtot];
        double g[twoCases[i]];
        double Erow[twoCases[i]];
        int n = -1;
        if (i==0){
            E_fn = fopen("deltaE_24.txt", "r");
            k_fn = fopen("k_24.txt", "r");
            A_fn = fopen("A_24.txt", "r");
            double v[24] = {0,0,0,0,1,1,0,1,0,1,2,2,1,0,2,1,2,3,0,3,2,1,3,3};
            double J[24] = {1,3,5,7,1,3,9,5,11,7,1,3,9,13,5,11,7,1,15,3,9,13,5,7};
            X=xO;
            for (int j = 0; j < 24; j++){
              g[j] = 2*J[j] + 1;
            }
        } else {
            E_fn = fopen("deltaE_27.txt", "r");
            k_fn = fopen("k_27.txt", "r");
            A_fn = fopen("A_27.txt", "r");
            double v[27] = {0,0,0,0,0,1,1,1,0,1,1,2,0,2,2,1,2,0,2,3,1,3,3,0,2,3,3};
            double J[27] = {0,2,4,6,8,0,2,4,10,6,8,0,12,2,4,10,6,14,8,0,12,2,4,16,10,6,8};
            X=xP;
            for (int j = 0; j < 27; j++){
              g[j] = 2*J[j] + 1;
            }
        }

        if (E_fn == NULL) {
          printf("failed to open file\n");
        }
        while (fscanf(E_fn, "%f", &fValues[++n]) == 1) {
          fscanf(E_fn, ",");
          E[n/itot][n%jtot] = fValues[n];
        }
        fclose(E_fn);
        n = -1;
        if (A_fn == NULL) {
          printf("failed to open file\n");
        }
        while (fscanf(A_fn, "%e", &fValues[++n]) == 1) {
          fscanf(A_fn, ",");
          A[n/itot][n%jtot] = fValues[n];
        }
        fclose(A_fn);
        n = -1;
        if (k_fn == NULL) {
          printf("failed to open file\n");
        }
        while (fscanf(k_fn, "%e", &fValues[++n]) == 1) {
          fscanf(k_fn, ",");
          int i = n/(itot*50);
          int j = n%(jtot*50);
          k[i][j] = fValues[n];
        }
        fclose(k_fn);

          //loop through 24 and fill g = 2*J+1
          for (int j = 0; j < twoCases[i]; j++){
            Erow[j] = E[0][j];
          }

          //[E_,T_]=meshgrid(E,T);
          double E_mesh[N][twoCases[i]];
          double T_mesh[N][twoCases[i]];
          for (int j = 0; j < N; j++){
            for (int k = 0; k < twoCases[i]; k++){
                T_mesh[j][k] = T[j];
                E_mesh[j][k] = Erow[k];
            }
          }

          //gg=meshgrid(g,T);
          double gg_mesh[N][twoCases[i]];
          for (int j = 0; j < N; j++){
              for (int k = 0; k < twoCases[i]; k++){
                gg_mesh[j][k] = g[k];
            }
          }
          double S[N];
          double x_LTE[N][twoCases[i]];
          for (int j = 0; j < N; j++){
          S[j] = 0;
            for (int k = 0; k < twoCases[i]; k++){
              S[j]+=gg_mesh[j][k]*exp(-E_mesh[j][k]/T_mesh[j][k]);
            }
          }

          for (int j = 0; j < N; j++){
            for (int k = 0; k < twoCases[i]; k++){
              x_LTE[j][k] = X*gg_mesh[j][k]*exp(-E_mesh[j][k]/T_mesh[j][k])/S[j];
            }
          }
          double kB=1.381e-16; // ergs per kelvin
        for (int r = 0; r < N; r++){
          if (i == 0){L[r] = 0;}
        }

        double temp2[twoCases[i]][N];

        for (int j = 0; j < twoCases[i]; j++){
          for (int r = 0; r < N; r++){
            temp2[j][r] = 0.;
          }
        }


        for (int j = 0; j < twoCases[i]; j++){
          for (int k = 0; k < twoCases[i]; k++){
            for (int r = 0; r < N; r++){
              temp2[j][r] += A[j][k]*E[j][k]*kB*x_LTE[r][k];
            }
          }
        }

        for (int r = 0; r < N; r++){
          for (int j = 0; j < twoCases[i]; j++){
            L[r] += temp2[j][r];
          }
        }

      } //twoCases end loop
      return L;
}

double *H2_cooling_coef_GA08(int N, double n[],double T[],double xH,double xHe,double xH2,double xHp, double xe, double O2P){
  double *L = (double *) malloc(N * sizeof(double)); //dynamically allocate memory for the rate

  double *H;
  double *He;
  double *Hp;
  double *H2;
  double *e;
  double *H2_LTE_arr;

  H = H2_cooling_coef_H_GA08(N, T,O2P);
  He = H2_cooling_coef_He_GA08(N, T,O2P);
  H2 =  H2_cooling_coef_H2_GA08(N, T,O2P);
  Hp =  H2_cooling_coef_Hp_GA08(N, T,O2P);
  e = H2_cooling_coef_e_GA08(N, T,O2P);

	double L_low_n[N];
  double L_LTE[N];
  H2_LTE_arr = H2_LTE(N, T);

  for (int i = 0; i < N; i++){
    L_low_n[i] = xH  * H[i] + xHe * He[i] +
    xH2 * H2[i] + xHp * Hp[i] + xe  * e[i];
    L_LTE[i]=H2_LTE_arr[i]/n[i];
    L[i]=L_LTE[i]/(1+L_LTE[i]/L_low_n[i]);
  }
  return L;
}


double *CII_cooling_coef(int N, double n[],double T[],double x_H2, double x_e){
  double *L = (double *) malloc(N * sizeof(double)); //dynamically allocate memory for the rate
	double OtoP=3;
	double x_para=(1/(OtoP+1))*x_H2;
	double x_ortho=(OtoP/(OtoP+1))*x_H2;
	double x_H=1-2*x_H2;

  //these are going to be integers, but we multiply them by floats later
  //so why risk it?
	double g2=4;
	double g1=2;

	double hv=1.2593e-14; //in erg
	double A=2.300e-6; //in s^-1
	double k=1.3806488e-16; //in K^-1 * erg

	// Using Draine's fitts (p.501)
  double T2[N], T4[N], logT2[N];
  for (int i = 0; i < N; i++){
	   T2[i] = T[i]/1.e2;
	   T4[i] = T[i]/1.e4;
     logT2[i] = log(T2[i]);
  }

	double Q_H[N], Q_pH2[N], Q_oH2[N], Q_e[N], Qni_21[N], Qni_12[N], n1_to_n2[N], x2[N];

  for (int i = 0; i < N; i++){
  	Q_H[i] = 7.58e-10*pow(T2[i], (0.128+0.009*logT2[i]));
  	Q_pH2[i] = 4.25e-10*pow(T2[i], (0.124-0.018*logT2[i]));
  	Q_oH2[i] = 5.14e-10*pow(T2[i], (0.095+0.023*logT2[i]));
  	Q_e[i] = 8.629e-8*pow(T4[i], -0.5*(1.55+1.25*T4[i])/(1+0.35*pow(T4[i],1.25))/g2);
  	Qni_21[i] = n[i]*(Q_pH2[i]*x_para+Q_oH2[i]*x_ortho+Q_H[i]*x_H+Q_e[i]*x_e);
  	Qni_12[i] = Qni_21[i]*g2/g1*exp(-hv/(k*T[i]));
  	n1_to_n2[i] = (Qni_21[i]+A)/Qni_12[i];
  	x2[i] = 1./(1+n1_to_n2[i]);
  	L[i] = x2[i]*hv*A/n[i];
  }
  return L;
}

double *dust_cooling_coef(int N, double T[], float IUV)
{
//  if (nargin<2) {
//      double Td = 0;
//  } else {
      double T0 = 16.4;
      double a_0_1_mu_m = 1;
      double Td = T0 * (a_0_1_mu_m) * pow(IUV, 1./6.); // Draine  pp.290
//  }
// let's assume that we always feed it an IUV rate
   // glover & Jappsen (2007, table 7)
  for (int i = 0; i < N; i++){
    L[i] = 3.8e-33 * pow(T[i], 0.5) * (T[i] - Td) * (1 - 0.8 * exp(-75./T[i]));
  }
}

//calculates the cooling coefficient of [OI]
//input:
//T is the temperature. valid in the range of 20-10e10. acurate only for the range 20-1000K
//ni is a row array with the densities of the collisional ellements
//[n-pH2,n-oH2,n-H,n-e]
//nOI is the OI density
//output: the cooling coeff in [erg*sec^-1*cm^+3]
double *CI_cooling_coef(int N, double n[],double T[],double x_H2){
	double OtoP=3;
  double *L = (double *) malloc(N * sizeof(double)); //dynamically allocate memory for the rate
	double  x_para=(1/(OtoP+1))*x_H2;
	double x_ortho=(OtoP/(OtoP+1))*x_H2;
	double x_H=1-2*x_H2;

	//===========we neglect here unelastic collisions with H+==========
	//data:
	//from Lamda database
	double h=6.62606957e-27;
	double k=1.3806488e-16;
	double v10=492.160651e9;
	double v20=1301.50262e9;
	double v21=809.34197e9;

	double E10=v10*h;
	double E20=v20*h;
	double E21=v21*h;

	double A10=7.88e-8;
	double A20=1.810e-14;
	double A21=2.65e-7;

	double g0=1;
	double g1=3;
	double g2=5;

	// Using Draine's fitts (p.501)
  double T2[N];
  double q10_H[N], q20_H[N], q21_H[N], q10_para[N], q20_para[N], q21_para[N],
      q10_ortho[N], q20_ortho[N], q21_ortho[N];
  for (int i = 0; i < N; i++){
	  T2[i]=T[i]/100.;
  	q10_H[i] = 1.26e-11*pow(T2[i], 0.115-0.057*log(T2[i]));
  	q20_H[i] = 2.64e-11*pow(T2[i], 0.231-0.046*log(T2[i]));
  	q21_H[i] = 8.90e-11*pow(T2[i], 0.228-0.046*log(T2[i]));

  	q10_para[i] = 0.67e-11*pow(T2[i], -0.085+0.102*log(T2[i]));
  	q20_para[i] = 0.86e-10*pow(T2[i], -0.010+0.048*log(T2[i]));
  	q21_para[i] = 1.75e-10*pow(T2[i], 0.072+0.064*log(T2[i]));

  	q10_ortho[i] = 0.71e-10*pow(T2[i], -0.004+0.049*log(T2[i]));
  	q20_ortho[i] = 0.69e-10*pow(T2[i], 0.169+0.038*log(T2[i]));
  	q21_ortho[i] = 1.48e-10*pow(T2[i], 0.263+0.031*log(T2[i]));
  }

	// rates for H+ and e- are of the same order of mag
	// or at most factor of 10 higher, and thus not important
	// in gas with xe, xHp < 0.1
	// rates for helium are similar as for H, and AHe=0.1, thus not important

	double R10[N], R20[N], R21[N], R01[N], R02[N], R12[N], x1_to_x0[N], x2_to_x0[N], x0[N], x1[N], x2[N];
	for (int i = 0; i < N; i++){
		R10[i]=n[i]*(q10_H[i]*x_H+q10_para[i]*x_para+x_ortho*q10_ortho[i])+A10;
		R20[i]=n[i]*(q20_H[i]*x_H+q20_para[i]*x_para+x_ortho*q20_ortho[i])+A20;
		R21[i]=n[i]*(q21_H[i]*x_H+q21_para[i]*x_para+x_ortho*q21_ortho[i])+A21;
		R01[i]=n[i]*(q10_H[i]*x_H+q10_para[i]*x_para+x_ortho*q10_ortho[i])*(exp(-E10/(k*T[i]))*(g1/g0));
		R02[i]=n[i]*(q20_H[i]*x_H+q20_para[i]*x_para+x_ortho*q20_ortho[i])*(exp(-E20/(k*T[i]))*(g2/g0));
		R12[i]=n[i]*(q21_H[i]*x_H+q21_para[i]*x_para+x_ortho*q21_ortho[i])*(exp(-E21/(k*T[i]))*(g2/g1));

		x1_to_x0[i]=(R01[i]*R20[i]+R01[i]*R21[i]+R21[i]*R02[i])/(R10[i]*R20[i]+R10[i]*R21[i]+R12[i]*R20[i]);
		x2_to_x0[i]=(R02[i]*R10[i]+R02[i]*R12[i]+R12[i]*R01[i])/(R10[i]*R20[i]+R10[i]*R21[i]+R12[i]*R20[i]);
		x0[i]=1./(1+x1_to_x0[i]+x2_to_x0[i]);

		x1[i]=x1_to_x0[i]*x0[i];
		x2[i]=x2_to_x0[i]*x0[i];

		L[i]=(x1[i]*A10*E10+x2[i]*(A21*E21+A20*E20))/n[i];
	}
  return L;
}
