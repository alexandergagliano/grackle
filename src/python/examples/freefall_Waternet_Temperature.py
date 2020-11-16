########################################################################
#
# Free-fall example script
#
#
# Copyright (c) 2013-2016, Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

from matplotlib import pyplot
import os
import yt
import numpy as np
import pandas as pd

from pygrackle import \
    chemistry_data, \
    FluidContainer, \
    evolve_constant_density, \
    evolve_freefall

from pygrackle.utilities.physical_constants import \
    mass_hydrogen_cgs, \
    mass_electron_cgs, \
    sec_per_Myr, \
    cm_per_mpc

from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib import cm

cm = cm.get_cmap('coolwarm', 4)
cols = [cm(0.2), cm(0.4), cm(0.6), cm(0.8)]

tiny_number = 1e-60

if __name__=="__main__":
    current_redshift = 0.

    # Set solver parameters
    metallicity = [1.e-6, 1.e-4, 1.e-2, 1.e0]
    logZ = [-6., -4., -2., 0.]
    
    #pyplot.figure(dpi=1000)
    for i in range(4):
        my_chemistry = chemistry_data()
        my_chemistry.use_grackle = 1
        my_chemistry.with_radiative_cooling = 1
        my_chemistry.primordial_chemistry = 3
        my_chemistry.UVbackground = 0
        my_chemistry.self_shielding_method = 0
        my_chemistry.H2_self_shielding = 0
        my_chemistry.Gamma = 5. / 3.
        my_chemistry.three_body_rate = 1
        my_chemistry.water_rates = 1
        my_chemistry.CaseBRecombination = 0
        my_chemistry.cie_cooling = 1
        my_chemistry.h2_optical_depth_approximation = 1
        my_chemistry.cmb_temperature_floor = 1
        my_chemistry.withWater = 1 
        my_chemistry.metal_cooling = 1
        grackle_dir = os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.dirname(os.path.abspath(__file__)))))
        my_chemistry.grackle_data_file = os.sep.join(
            [grackle_dir, "input", "cloudy_metals_2008_3D.h5"])
        my_chemistry.photoelectric_heating = 0
        my_chemistry.h2_on_dust = 1
        my_chemistry.interstellar_radiation_field = 0
        my_chemistry.dust_chemistry = 1 
 
        # Set units
        my_chemistry.comoving_coordinates = 0 # proper units
        my_chemistry.a_units = 1.0
        my_chemistry.a_value = 1. / (1. + current_redshift) / \
            my_chemistry.a_units
        my_chemistry.density_units  = mass_hydrogen_cgs # rho = 1.0 is 1.67e-24 g
        my_chemistry.length_units   = cm_per_mpc        # 1 Mpc in cm
        my_chemistry.time_units     = sec_per_Myr       # 1 Myr in s
        my_chemistry.velocity_units = my_chemistry.a_units * \
            (my_chemistry.length_units / my_chemistry.a_value) / \
            my_chemistry.time_units

        # set initial density and temperature
        initial_temperature = 300. # start the gas at this temperature
        # then begin collapse
        initial_density     = 1.e-1 * mass_hydrogen_cgs # g / cm^3
        # stopping condition
        final_density       = 1.e12 * mass_hydrogen_cgs

        rval = my_chemistry.initialize()


        fc = FluidContainer(my_chemistry, 1)
        fc["density"][:] = initial_density / my_chemistry.density_units
        fc["HI"][:] =  0.758354919 * fc["density"]
        fc["HII"][:] = 0.0002785698815577569 * 0.76 * fc["density"]
        fc["HeI"][:] = 0.0972 * fc["density"]
        fc["HeII"][:] = 4.098176546009e-7 * fc["density"]
        fc["HeIII"][:] = 7.56797580911136e-22 * fc["density"]
        fc["de"][:] = 2.786722968998593e-4 * fc["density"]
        if my_chemistry.primordial_chemistry > 1:
            fc["H2I"][:] = 1.e-6 * fc["density"]
            #fc["H2I"][:] = 1.e-20 * fc["density"]
            fc["H2II"][:] = 1.195304067126978e-12 * fc["density"]
            fc["HM"][:] = 3.966920888635451e-11 * fc["density"]
        if my_chemistry.primordial_chemistry > 2:
            fc["DI"][:] = 3.e-5 * metallicity[i] * fc["density"]
            fc["DII"][:] = 1.676644194346e-8 * fc["density"]
            fc["HDI"][:] = 2.36163087768956e-7 * fc["density"]
        if my_chemistry.metal_cooling == 1:
            fc["metal"][:] = metallicity[i] * my_chemistry.SolarMetalFractionByMass * \
                fc["density"]
        if my_chemistry.withWater == 1:
            fc["Water_density"][:] = tiny_number * fc["density"]
            fc["O_density"][:] = 3.568e-4 * metallicity[i] * fc["density"]
            fc["OH_density"][:] = tiny_number * fc["density"]
            fc["O2_density"][:] = tiny_number * fc["density"]
            fc["Oplus_density"][:] = tiny_number * fc["density"]
            fc["OHplus_density"][:] = tiny_number * fc["density"]
            fc["H2Oplus_density"][:] = tiny_number * fc["density"]
            fc["H3Oplus_density"][:] = tiny_number * fc["density"]
            fc["O2plus_density"][:] = tiny_number * fc["density"]
            fc["Cplus_density"][:] = 0.927e-4 * metallicity[i] * fc["density"]
            fc["C_density"][:] = tiny_number * fc["density"]
            fc["CH_density"][:] = tiny_number * fc["density"]
            fc["CH2_density"][:] = tiny_number * fc["density"]
            fc["CH3_density"][:] = tiny_number * fc["density"]
            fc["CH4_density"][:] = tiny_number * fc["density"]
            fc["CO_density"][:] = tiny_number * fc["density"]
            fc["COplus_density"][:] = tiny_number * fc["density"]
            fc["CO2_density"][:] = tiny_number * fc["density"]

        fc["energy"][:] = initial_temperature / \
            fc.chemistry_data.temperature_units
        fc["x-velocity"][:] = 0.0
        fc["y-velocity"][:] = 0.0
        fc["z-velocity"][:] = 0.0

        # timestepping safety factor
        safety_factor = 0.01

        # then begin collapse 
        # evolve density and temperature according to free-fall collapse
        data = evolve_freefall(fc, final_density,
                               safety_factor=safety_factor)
      

        # make a plot of nH2 vs. T
        pyplot.plot(np.log10( (data["HI"]+data["HII"]+data["H2I"]+data["H2II"]) / mass_hydrogen_cgs), np.log10(data["temperature"]), label ='log(Z/Z$_{\odot}$) = %d, Sim' % logZ[i], color=cols[i])
        if i == 0:
            pyplot.xlabel("log(n$_H$ [cm$^{-3}$])")
            pyplot.ylabel("log(T [K])")

        ## OPEN OMUKAI 2005 FOR COMPARISON
        with open('waterData/T_vs_logn__Z=%d_WaterOmukai2005.csv' % logZ[i], 'r') as omukaiData:
                lines = omukaiData.readlines()
                Omukai_n = []
                Omukai_T = []
                for line in lines[6:]:
                        p = line.split(',')
                        Omukai_n.append(float(p[0]))
                        Omukai_T.append(float(p[1]))
                pyplot.plot(Omukai_n, np.log10(Omukai_T), label='log (Z/Z$_{\odot}$) = %d, Omukai' % logZ[i], color=cols[i], linestyle=':')

                Omukai_n = []
                Omukai_T = []

        #writing fH2 data!
        logfH2_df = pd.DataFrame({"log(n_H)" : np.log10( (data["HI"]+data["HII"]+data["H2I"]+data["H2II"]) / mass_hydrogen_cgs ), "log(fH2)" : np.log10( data["H2I"] / (data["HI"]+data["HII"]+data["H2I"]+data["H2II"]) )})
        logfH2_df.to_csv("logfH2_v_logN_[Z]=%d_2005.csv" % logZ[i], index=None)
 
    pyplot.xlim(xmin = -1.5, xmax=12)
    pyplot.ylim(ymax=4)

    pyplot.text(5, 2.78, '-6', fontsize=16)
    pyplot.text(7, 2.40, '-4', fontsize=16)
    pyplot.text(2.5, 1.7, '-2', fontsize=16)
    pyplot.text(6, 0.7, '0', fontsize=16)


    #leg = pyplot.legend(fancybox = True, labelspacing=0.0, loc='upper right')
    #leg.get_frame().set_alpha(0.5)
 
    if my_chemistry.metal_cooling:
        output = "freefall_metal"
    else:
        output = "freefall"

    pyplot.tight_layout()
    pyplot.savefig("%s_Water.pdf" % (output))

