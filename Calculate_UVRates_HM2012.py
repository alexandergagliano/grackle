import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy import interpolate
os.chdir('/home/alexgagliano/Documents/Research/water_network/Fresh_1126_crxFlag/input')

all_species = ['C', 'CH', 'CH2', 'CH3', 'CH4', 'CO', 'H2O', 'H-', 'O2', 'O', 'OH']
reaction_types = ['photoabsorption', 'photodissociation', 'photoionization']
h = 6.626196e-27
df_UVB = pd.read_csv("UVB_HM12.out", delim_whitespace=True) # currently in units of ergs/s/cm^2/Hz/sr
# get the redshifts for which we have a UV Background
z = np.array([float(x) for x in df_UVB.columns.values[1:]])
z_cols = df_UVB.columns.values[1:]
N_z = len(z)


for i in np.arange(N_z):
    df_UVB[z_cols[i]] *= 4*np.pi # assuming isotropic; now in units of ergs/s/cm^2/Hz
    df_UVB[z_cols[i]] *= 1/(h*df_UVB['Lambda(Angstroms)'])  # now in units of phot/s/cm^2/A

os.chdir('/home/alexgagliano/Documents/Research/water_network/Fresh_1126_crxFlag/input/metal_UVB_HM12_crs')
for species in all_species:
    for rxn in reaction_types:
        #try:
        df_crs = pd.read_csv("%s_0.1nm.txt" % species, delim_whitespace=True)
            # convert from nm to Angstroms
        df_crs['wavelength'] *= 10

        f = interpolate.interp1d(df_UVB['Lambda(Angstroms)'], df_UVB['0.0000E+000'])
        lambda_UVB_new = f(df_crs['wavelength'])
            #plt.plot(df_crs['wavelength'], lambda_UVB_new, 'o')
            #plt.plot(df_crs['wavelength'], df_crs['photoionisation'],'o')

        dlambda = 1.0 # 1 angstrom separation between each data point
        k_z = []
        for i in np.arange(N_z):
            #trapezoidal estimate to the integration
            tempK = np.sum(0.5*dlambda*(df_crs['%s'%rxn]*df_UVB[z_cols[i]]))
            k_z.append(tempK)
        df_calcRates = pd.DataFrame({'z': z_cols, 'k(s^-1)' : k_z})
        df_calcRates.to_csv("%s_%s.txt"%(species, rxn), index=False)
        #except:
        #    continue
