import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
from scipy import interpolate
os.chdir('/home/alexgagliano/Documents/Research/water_network/Fresh_1126_crxFlag/input')

species_set = ['C', 'CH', 'CH+', 'CH2', 'CH2+', 'CH3', 'CH3+', 'CH4', 'CH4+', 'CH5+', 'CO', 'CO+', 'H-', 'H2+', 'H2O', 'H2O+', 'H3O+', 'HCO+', 'O', 'O2', 'O2+', 'OH', 'OH+']
reaction_set = ['photoabsorption', 'photodissociation', 'photoionisation']

h = 6.626196e-27
df_UVB = pd.read_csv("UVB_HM12.out", delim_whitespace=True) # currently in units of ergs/s/cm^2/Hz/sr
# get the redshifts for which we have a UV Background
z = np.array([float(x) for x in df_UVB.columns.values[1:]])
z_cols = df_UVB.columns.values[1:]
N_z = len(z)
c = 3.e10 #cm/s
for i in np.arange(N_z):
    df_UVB[z_cols[i]] *= 4*np.pi # assuming isotropic; now in units of ergs/s/cm^2/Hz
    df_UVB[z_cols[i]] *= 1/(h*df_UVB['Lambda(Angstroms)'])  # now in units of phot/s/cm^2/A

df_UVB_Bialy = df_UVB
for i in np.arange(N_z):
    df_UVB_Bialy[z_cols[i]] /= (4*np.pi) # assuming isotropic; now in units of phot cm^-2 s^-1 str^-1 A^-1
    df_UVB_Bialy[z_cols[i]] *= df_UVB['Lambda(Angstroms)'] # now in units of phot cm^-2 s^-1 str^-1

    df_UVB_Bialy[z_cols[i]] *= (df_UVB['Lambda(Angstroms)']/(h*c)) # now in units of photons cm^-2 s^-1 erg^-1  str^-1
    
plt.semilogy((h*c*1.e8)/df_UVB['Lambda(Angstroms)']*6.2415e11, df_UVB_Bialy[z_cols[0]]*(1.602e-12)/1.e6) #to convert to eV
plt.xlim(xmin=6, xmax=14)

os.chdir('/home/alexgagliano/Documents/Research/water_network/Fresh_1126_crxFlag/input/metal_UVB_HM12_crs')
for species in species_set:
    for rxn in reaction_set:
        try:
            df_crs = pd.read_csv("%s_0.1nm.txt" % species, delim_whitespace=True)
            # convert from nm to Angstroms
            df_crs['wavelength'] *= 10

            f = interpolate.interp1d(df_UVB['Lambda(Angstroms)'], df_UVB['0.0000E+000'])
            lambda_UVB_new = f(df_crs['wavelength'])
            #plt.plot(df_crs['wavelength'], lambda_UVB_new, 'o')
            #plt.plot(df_crs['wavelength'], df_crs['%s'%rxn],'o')
            #plt.show()
            dlambda = 1.0 # 1 angstrom separation between each data point
            k_z = []
            for i in np.arange(N_z):
                #trapezoidal estimate to the integration
                tempK = np.sum(0.5*dlambda*(df_crs['%s' % rxn]*df_UVB[z_cols[i]]))
                k_z.append(tempK)
            df_rxns = pd.DataFrame({'z' : z_cols, 'k(s^-1)' : k_z})
            df_rxns.to_csv("%s_%s.txt" % (species, rxn), index=False)
        except:
            continue

rxn_dict = {}
rxn_dict[1] = ['OH', 'photodissociation']
rxn_dict[2] = ['H2O', 'photodissociation']
rxn_dict[3] = ['H2O', 'photodissociation']
rxn_dict[4] = ['O2', 'photodissociation']
rxn_dict[5] = ['CH', 'photodissociation']
rxn_dict[6] = ['CO', 'photodissociation']
rxn_dict[7] = ['C', 'photoionisation']
rxn_dict[8] = ['O', 'photoionisation']
rxn_dict[9] = ['H-', 'photoionisation']
rxn_dict[10] = ['CH', 'photoionisation']
rxn_dict[11] = ['CH2', 'photodissociation']
rxn_dict[12] = ['CH2+', 'photodissociation']
rxn_dict[13] = ['CH3', 'photodissociation']
rxn_dict[14] = ['CH4', 'photodissociation']
rxn_dict[15] = ['CH4+', 'photodissociation']
rxn_dict[16] = ['CH+', 'photodissociation']
rxn_dict[17] = ['CO+', 'photodissociation']
rxn_dict[18] = ['H2+', 'photodissociation']
rxn_dict[19] = ['HCO+', 'photodissociation']
rxn_dict[20] = ['O2+', 'photodissociation']
rxn_dict[21] = ['OH+', 'photodissociation']
rxn_dict[22] = ['CH4', 'photoionisation']
rxn_dict[23] = ['H2O', 'photoionisation']
rxn_dict[24] = ['O2', 'photoionisation']
rxn_dict[25] = ['CH2', 'photoionisation']
rxn_dict[26] = ['CH3', 'photoionisation']
rxn_dict[27] = ['H2O+', 'photodissociation'] #none below 911.6 angstroms!
rxn_dict[28] = ['H2O+', 'photodissociation']
rxn_dict[29] = ['H2O+', 'photodissociation']
rxn_dict[30] = ['H3O+', 'photodissociation']
rxn_dict[31] = ['H3O+', 'photodissociation']
rxn_dict[32] = ['H3O+', 'photodissociation']
rxn_dict[33] = ['H3O+', 'photodissociation']
rxn_dict[34] = ['CH2+', 'photodissociation']
rxn_dict[35] = ['CH3+', 'photodissociation']
rxn_dict[36] = ['CH3+', 'photodissociation']
rxn_dict[37] = ['CH4+', 'photodissociation']
rxn_dict[38] = ['CH4+', 'photodissociation']
rxn_dict[39] = ['CH5+', 'photodissociation']
rxn_dict[40] = ['CH5+', 'photodissociation']

pd_final_df = pd.DataFrame({'z':z_cols})

for i in np.arange(1, 41):
    try:
        [species, rxn] = rxn_dict[i]
        tempName = "UV" + str(i)
        fn = "%s_%s.txt" % (species, rxn)
        tempDF = pd.read_csv(fn)
        pd_final_df[tempName] = np.array(tempDF['k(s^-1)'])
    except:
        continue
os.chdir('/home/alexgagliano/Documents/Research/water_network/Fresh_1126_crxFlag/input/')
pd_final_df

pd_final_df.columns.values
pd_final_df.to_csv("UVB_HM2012_waterNetwork.txt", index=False)

# UV27-UV33, UV35, UV36, UV39, UV40
# CH3+ CH5+ H2O+ H3O+
