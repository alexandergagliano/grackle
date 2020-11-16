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
import sys
import matplotlib as mpl
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib import cm

mpl.rcParams['agg.path.chunksize'] = 10000

from pygrackle import \
    chemistry_data, \
    FluidContainer, \
    evolve_constant_density, \
    evolve_freefall, \
    evolve_freefall_metal

from pygrackle.utilities.physical_constants import \
    mass_hydrogen_cgs, \
    mass_electron_cgs, \
    sec_per_Myr, \
    cm_per_mpc

tiny_number = 1.e-60

amu_cgs           = 1.660538921e-24  # g

dark2 = cm.get_cmap('Dark2', 7)
cols_metl = dark2(range(7))
brbg = cm.get_cmap('BrBG', 10)
cols_h2 = [brbg(0.25), brbg(0.75)]
set2 = cm.get_cmap('Set2', 3)
col_ion = set2(range(3))


if __name__=="__main__":
    current_redshift = 0.

    # Set solver parameters
    final_metallicity = 1.0
    N_pts = 1
    dtmax = 1.e12 # in seconds
    #metallicity = np.logspace(-2,0, N_pts)dd
    metallicity = np.array([1.0])
    n_tot = []
    XO = []
    XC = []
    XH2O = []
    XOH = []
    XO2 = []
    XCH = []
    XCO = []
    XO_tot = []
    XC_tot = []
    H = []
    H2 = []
    HM = []
    H2plus = []
    Hplus = []
    XH = []
    XH2 = []
    XHplus = []
    Xel = []
    Xions = []
    el = []
    ions = []

    for i in range(len(metallicity)):
        print("Running the point for metallicity %.2e\n" % metallicity[i])
        my_chemistry = chemistry_data()
        my_chemistry.use_grackle = 1
        my_chemistry.with_radiative_cooling = 1
        my_chemistry.primordial_chemistry = 2
        my_chemistry.UVbackground = 1
        my_chemistry.self_shielding_method = 0
        my_chemistry.H2_self_shielding = 0
        my_chemistry.Gamma = 5. / 3.
        my_chemistry.three_body_rate = 0
        my_chemistry.CaseBRecombination = 0 
        my_chemistry.cie_cooling = 0 
        my_chemistry.h2_optical_depth_approximation = 0
        my_chemistry.water_rates = 3
        my_chemistry.cmb_temperature_floor = 1
        my_chemistry.withWater = 1
        my_chemistry.water_only = 1
        my_chemistry.dust_chemistry = 1
        my_chemistry.crx_ionization = 4 #strongest possible ionization rate 
        my_chemistry.metal_cooling = 1
        grackle_dir = os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.dirname(os.path.abspath(__file__)))))
        my_chemistry.grackle_data_file = os.sep.join(
            [grackle_dir, "input", "CloudyData_UVB=HM2012_shielded.h5"])
        my_chemistry.h2_on_dust = 1

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
        initial_temperature = 1.e5
        #initial_density   = 1.0e3/0.9 * mass_hydrogen_cgs
        initial_density = 1.e3 * mass_hydrogen_cgs
        final_time = 3.2e3 # 0.05 #1.e3 #1. # 1.e3 # 1.e3 # in Myr

        rval = my_chemistry.initialize()

        fc = FluidContainer(my_chemistry, 1)
        fc["density"][:] = initial_density / my_chemistry.density_units
        fc["HII"][:] = tiny_number * fc["density"]
        fc["HI"][:] = (0.65 - 4.9e-4 * metallicity[i] - 2.9e-4 * metallicity[i])* fc["density"]
        fc["HeI"][:] = tiny_number * fc["density"]
        fc["HeII"][:] = tiny_number * fc["density"]
        fc["HeIII"][:] = tiny_number * fc["density"]
        if my_chemistry.primordial_chemistry > 1:
            fc["H2I"][:] = 0.30 * fc["density"]
            fc["H2II"][:] = tiny_number * fc["density"]
            fc["HM"][:] = tiny_number * fc["density"]
        if my_chemistry.primordial_chemistry > 2:
            fc["HDI"][:] = 1.e-9 * fc["density"]
            fc["DI"][:] = (3.e-5 - fc["HDI"][:]) * fc["density"]
            fc["DII"][:] = tiny_number * fc["density"]
        if my_chemistry.metal_cooling == 1:
            fc["metal"][:] = metallicity[i] * my_chemistry.SolarMetalFractionByMass * \
                fc["density"]
        if my_chemistry.withWater == 1:
            fc["de"][:] = tiny_number * fc["density"]
            fc["HeI"][:] = 0.1 * fc["density"]
            fc["Water_density"][:] = tiny_number * fc["density"]
            fc["O_density"][:] = 0.553 * metallicity[i] * fc["density"]
            fc["OH_density"][:] = tiny_number * fc["density"]
            fc["O2_density"][:] = tiny_number * fc["density"]
            fc["Oplus_density"][:] = tiny_number * fc["density"]
            fc["OHplus_density"][:] = tiny_number * fc["density"]
            fc["H2Oplus_density"][:] = tiny_number * fc["density"]
            fc["H3Oplus_density"][:] = tiny_number * fc["density"]
            fc["O2plus_density"][:] = tiny_number * fc["density"]
            fc["Cplus_density"][:] = tiny_number * fc["density"]
            fc["C_density"][:] = 0.228 * metallicity[i] * fc["density"]
            fc["CH_density"][:] = tiny_number * fc["density"]
            fc["CH2_density"][:] = tiny_number * fc["density"]
            fc["CH3_density"][:] = tiny_number * fc["density"]
            fc["CH4_density"][:] = tiny_number * fc["density"]
            fc["CO_density"][:] = tiny_number * fc["density"]
            fc["COplus_density"][:] = tiny_number  * fc["density"]
            fc["CO2_density"][:] = tiny_number * fc["density"]
            if (my_chemistry.water_rates == 3):
                fc["CHplus_density"][:] = tiny_number * fc["density"]
                fc["CH2plus_density"][:] = tiny_number * fc["density"]
                fc["H3plus_density"][:] = tiny_number * fc["density"]
                fc["HCOplus_density"][:] = tiny_number * fc["density"]
                fc["HeHplus_density"][:] = tiny_number * fc["density"]
                fc["CH3plus_density"][:] = tiny_number * fc["density"]
                fc["CH4plus_density"][:] = tiny_number * fc["density"]
                fc["CH5plus_density"][:] = tiny_number * fc["density"]
                fc["O2Hplus_density"][:] = tiny_number * fc["density"]

        fc["energy"][:] = initial_temperature / \
            fc.chemistry_data.temperature_units

        fc["x-velocity"][:] = 0.0
        fc["y-velocity"][:] = 0.0
        fc["z-velocity"][:] = 0.0

        fc.calculate_hydrogen_number_density()

        # then begin collapse
        # evolve density and temperature according to free-fall collapse
        data = evolve_freefall_metal(fc, final_metallicity, final_time, dtmax=dtmax)

######### METALLICITY CONVERGENCE PLOTS ###################
        pyplot.loglog(data["time"], data["O_density"]/(data['O_density'] + data['Oplus_density']), label ='O', color='maroon',markevery=.1)
        pyplot.loglog(data["time"], data["C_density"]/(data['C_density'] + data['Cplus_density']), label ='C', color='goldenrod',markevery=.1)
        pyplot.loglog(data["time"], data["Oplus_density"]/(data['O_density'] + data['Oplus_density']), label ='O+', color='navy',markevery=.1)
        pyplot.loglog(data["time"], data["Cplus_density"]/(data['C_density'] + data['Cplus_density']), label ='C+', color='magenta',markevery=.1)

        pyplot.xlabel("Time (s)",fontsize=22)
        pyplot.ylabel("%",fontsize=22)
        pyplot.title("Convergence plot for [Z/H] = %i" % np.log10(metallicity[i]));
        pyplot.tight_layout()
        #pyplot.xlim(left=1.e3)
        pyplot.ylim(ymin=1.e-17)
        pyplot.legend(loc='best',fontsize=20)
        pyplot.savefig("metal_convergence_t%iGyr_N%i_Z%i_UMIST.pdf" %(int(np.log10(final_time/1.e3)),N_pts, np.log10(metallicity[i])),dpi=300)
        pyplot.clf()
