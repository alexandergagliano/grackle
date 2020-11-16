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
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib import cm

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

dark2 = cm.get_cmap('Dark2', 4)
cols = dark2(range(4))

tiny_number = 1.e-60

amu_cgs           = 1.660538921e-24  # g
mass_oxygen_cgs = 15.9994*amu_cgs  # g
mass_carbon_cgs = 12.0107*amu_cgs  # g

if __name__=="__main__":
    current_redshift = 0.

    # Set solver parameters
    colors = ["royalblue", "peru", "mediumorchid", "mediumseagreen"]
    metallicity = [1.e-6, 1.e-4, 1.e-2]
    logZ = [-6, -4, -2]
    for i in range(len(metallicity)):
        my_chemistry = chemistry_data()
        my_chemistry.use_grackle = 1
        my_chemistry.with_radiative_cooling = 1
        my_chemistry.primordial_chemistry = 3
        my_chemistry.UVbackground = 0
        my_chemistry.self_shielding_method = 0
        my_chemistry.H2_self_shielding = 0 
        my_chemistry.water_rates = 1
        my_chemistry.water_only = 0
        my_chemistry.Gamma = 5. / 3.
        #my_chemistry.three_body_rate = 1
        my_chemistry.CaseBRecombination = 0
        my_chemistry.cie_cooling = 1
        my_chemistry.h2_optical_depth_approximation = 1
        my_chemistry.cmb_temperature_floor = 0
        my_chemistry.withWater = 1
        my_chemistry.metal_cooling = 1
        my_chemistry.h2_on_dust = 1
        my_chemistry.dust_chemistry = 1
        my_chemistry.photoelectric_heating = 0
        my_chemistry.interstellar_radiation_field = 0
        grackle_dir = os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.dirname(os.path.abspath(__file__)))))
        my_chemistry.grackle_data_file = os.sep.join(
            [grackle_dir, "input", "cloudy_metals_2008_3D.h5"])
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
        initial_temperature = 300.
        # then begin collapse
        initial_density     = 1.e-1 * mass_hydrogen_cgs # g / cm^3
        # stopping condition
        final_density       = 5.e10 * mass_hydrogen_cgs
         
        rval = my_chemistry.initialize()

        fc = FluidContainer(my_chemistry, 1)
        fc["density"][:] = initial_density / my_chemistry.density_units
        fc["HII"][:] = 1.e-4 * fc["density"]
        fc["HeI"][:] = 8.333e-2 * fc["density"]
        fc["HeII"][:] = tiny_number * fc["density"]
        fc["HeIII"][:] = tiny_number * fc["density"]
        if my_chemistry.primordial_chemistry > 1:
            fc["H2I"][:] = 1.e-6 * fc["density"]
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
            fc["HI"][:] = (1.e0 - fc["HII"][:] - 2.e0 * fc["H2I"][:]) * fc["density"]
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
            fc["COplus_density"][:] = tiny_number  * fc["density"]
            fc["CO2_density"][:] = tiny_number * fc["density"]
            fc["de"][:] = 1.e-4 * fc["density"]
            fc["HeI"][:] = 8.333e-2 * fc["density"]
            fc["HeII"][:] = tiny_number * fc["density"]
            fc["HeIII"][:] = tiny_number * fc["density"]

        fc["energy"][:] = initial_temperature / \
            fc.chemistry_data.temperature_units


        fc["x-velocity"][:] = 0.0
        fc["y-velocity"][:] = 0.0
        fc["z-velocity"][:] = 0.0
 
        fc.calculate_hydrogen_number_density()
        #fc["nH"][:] = fc.calculate_hydrogen_number_density()

        # timestepping safety factor
        #safety_factor = 1.0
        safety_factor = 0.05
        data = evolve_freefall(fc, final_density,
                               safety_factor=safety_factor)

#        # FIRST SET OF PLOTS: H2O, OH, O
#        pyplot.figure(dpi=1000)

        f, (ax1, ax2, ax3) = pyplot.subplots(3, sharex=True, sharey=False)
        f.set_figheight(12)
        f.set_figwidth(12)
        ax2.plot(np.log10( data["nH"]), np.log10((data["Water_density"]+data["H2Oplus_density"]) / (data["nH"])) , label ='H$_2$O', color=cols[0])
        ax2.plot(np.log10( data["nH"]), np.log10((data["OH_density"]+data["OHplus_density"]) / (data["nH"])) , label ='OH', color=cols[1])
        ax2.plot(np.log10( data["nH"]), np.log10((data["O_density"] + data["Oplus_density"]) / (data["nH"])) , label ='O', color=cols[2])


        # OPEN OMUKAI 2005 FOR COMPARISON
        with open('waterData/logy(H2O)_v_logn_Z=%d_Omukai2005.csv' % logZ[i], 'r') as omukaiData:
                lines = omukaiData.readlines()
                Omukai_logn = []
                Omukai_spec = []
                for line in lines[6:]:
                        p = line.split(',')
                        Omukai_logn.append(float(p[0]))
                        Omukai_spec.append(float(p[1]))
                ax2.plot(Omukai_logn, Omukai_spec, color=cols[0], linestyle=':')

                Omukai_logn = []
                Omukai_spec = []

        ## OPEN OMUKAI 2005 FOR COMPARISON
        with open('waterData/logy(OH)_v_logn_Z=%d_Omukai2005.csv' % logZ[i], 'r') as omukaiData:
                lines = omukaiData.readlines()
                Omukai_logn = []
                Omukai_spec = []
                for line in lines[6:]:
                        p = line.split(',')
                        Omukai_logn.append(float(p[0]))
                        Omukai_spec.append(float(p[1]))
                ax2.plot(Omukai_logn, Omukai_spec, color=cols[1], linestyle=':')

                Omukai_logn = []
                Omukai_spec = []
#        ## OPEN OMUKAI 2005 FOR COMPARISON
        with open('waterData/logy(O)_v_logn_Z=%d_Omukai2005.csv' % logZ[i], 'r') as omukaiData:
                lines = omukaiData.readlines()
                Omukai_logn = []
                Omukai_spec = []
                for line in lines[6:]:
                        p = line.split(',')
                        Omukai_logn.append(float(p[0]))
                        Omukai_spec.append(float(p[1]))
                ax2.plot(Omukai_logn, Omukai_spec, color=cols[2], linestyle=':')
                Omukai_logn = []
                Omukai_spec = []
 
        ax2.set_ylabel(r"log y(O, OH, H$_2$O)")

        #xlim = [-1, 10]
        if (logZ[i] == -6):
                ylim = [-15, -8]
        elif (logZ[i] == -4):
                ylim = [-13,-7]
        elif (logZ[i] == -2):
                ylim = [-11,-4]
        ax2.set_ylim(ylim)

        #text for O, OH, H2O for different metallicities
        xO = [2, 1, 4.75]
        xOH = [0, 1, 1]
        xH2O = [4, 5, 2.75]

        yO = [-10, -8, -6]
        yOH = [-11.75, -9.8, -8]
        yH2O = [-13.5, -10.5, -9.75]

        ax2.text(xO[i], yO[i], 'O', fontsize=16)
        ax2.text(xOH[i], yOH[i], 'OH', fontsize=16)
        ax2.text(xH2O[i], yH2O[i], 'H$_2$O', fontsize=16)

        #ax2.legend(fancybox = True, labelspacing=0.0, loc='lower right', framealpha=0.5)

        # SECOND SET OF PLOTS: H2, e, HD 
#        pyplot.figure(dpi=1000)

        ax1.plot(np.log10( data["nH"]), np.log10((data["H2I"] + data["H2II"]) / (data["nH"]*my_chemistry.density_units)) , label =r'H$_2$', color=cols[0])
        ax1.plot(np.log10( data["nH"]), np.log10(data["de"] / (data["nH"]*my_chemistry.density_units)), label ='e', color=cols[1])
        if (my_chemistry.primordial_chemistry == 3):
            ax1.plot(np.log10( data["nH"]), np.log10(data["HDI"] / (data["nH"]*my_chemistry.density_units)) , label ='HD', color=cols[2])

        
        ## OPEN OMUKAI 2005 FOR COMPARISON
        with open('waterData/logy(H2)_v_logn_Z=%d_Omukai2005.csv' % logZ[i], 'r') as omukaiData:
                lines = omukaiData.readlines()
                Omukai_logn = []
                Omukai_spec = []
                for line in lines[6:]:
                        p = line.split(',')
                        Omukai_logn.append(float(p[0]))
                        Omukai_spec.append(float(p[1]))
                ax1.plot(Omukai_logn, Omukai_spec, color=cols[0], linestyle=':')

                Omukai_logn = []
                Omukai_spec = []

        ## OPEN OMUKAI 2005 FOR COMPARISON
        with open('waterData/logy(e)_v_logn_Z=%d_Omukai2005.csv' % logZ[i], 'r') as omukaiData:
                lines = omukaiData.readlines()
                Omukai_logn = []
                Omukai_spec = []
                for line in lines[6:]:
                        p = line.split(',')
                        Omukai_logn.append(float(p[0]))
                        Omukai_spec.append(float(p[1]))
                ax1.plot(Omukai_logn, Omukai_spec, color=cols[1], linestyle=':')

                Omukai_logn = []
                Omukai_spec = []

        ## OPEN OMUKAI 2005 FOR COMPARISON
        with open('waterData/logy(HD)_v_logn_Z=%d_Omukai2005.csv' % logZ[i], 'r') as omukaiData:
                lines = omukaiData.readlines()
                Omukai_logn = []
                Omukai_spec = []
                for line in lines[6:]:
                        p = line.split(',')
                        Omukai_logn.append(float(p[0]))
                        Omukai_spec.append(float(p[1]))
                ax1.plot(Omukai_logn, Omukai_spec, color=cols[2], linestyle=':')

                Omukai_logn = []
                Omukai_spec = []


        ax1.set_ylabel(r"log y(H$_2$, e, HD)")
        ylim = [-13,0]
        ax1.set_ylim(ylim)
        ax1.text(-0.5, ylim[1] - 1.5, '[Z/H] = %i' % logZ[i], fontsize=18)

                #text for H2, e, HD for different metallicities
        xH2 = [2, 6, 5]
        xe = [4, 1, 3.5]
        xHD = [1, 2, 8]

        yH2 = [-2.5, -3.5, -2]
        ye = [-5, -5.2, -9]
        yHD = [-9, -8, -6]

        ax1.text(xH2[i], yH2[i], r'H$_2$', fontsize=16)
        ax1.text(xe[i], ye[i], 'e', fontsize=16)
        ax1.text(xHD[i], yHD[i], 'HD', fontsize=16)

       # pyplot.xlim(xmin = xlim[0], xmax = xlim[1])

        #ax1.legend(fancybox = True, labelspacing=0.0, loc='lower right', framealpha=0.5)
        #leg.get_frame().set_framealpha(0.5)

       # if my_chemistry.metal_cooling:
       #     output = "freefall_metal"
       # else:
       #     output = "freefall"

        #pyplot.savefig("%s_Water_abundances2_%d.png" % (output, logZ[i]))
        #pyplot.clf()

        # THIRD SET OF PLOTS: C+, C, CO
        #pyplot.figure(dpi=1000)

        ax3.plot(np.log10( data["nH"]), np.log10((data["Cplus_density"]) / (data["nH"])), label =r'C$^+$', color=cols[0])
        ax3.plot(np.log10( data["nH"]), np.log10(data["C_density"] / (data["nH"])) , label ='C', color=cols[1])
        ax3.plot(np.log10( data["nH"]), np.log10((data["CO_density"] + data["COplus_density"]) / (data["nH"])) , label ='CO', color=cols[2])

        
        ## OPEN OMUKAI 2005 FOR COMPARISON
        with open('waterData/logy(C+)_v_logn_Z=%d_Omukai2005.csv' % logZ[i], 'r') as omukaiData:
                lines = omukaiData.readlines()
                Omukai_logn = []
                Omukai_spec = []
                for line in lines[6:]:
                        p = line.split(',')
                        Omukai_logn.append(float(p[0]))
                        Omukai_spec.append(float(p[1]))
                ax3.plot(Omukai_logn, Omukai_spec, color=cols[0], linestyle=':')

                Omukai_logn = []
                Omukai_spec = []

        ## OPEN OMUKAI 2005 FOR COMPARISON
        with open('waterData/logy(C)_v_logn_Z=%d_Omukai2005.csv' % logZ[i], 'r') as omukaiData:
                lines = omukaiData.readlines()
                Omukai_logn = []
                Omukai_spec = []
                for line in lines[6:]:
                        p = line.split(',')
                        Omukai_logn.append(float(p[0]))
                        Omukai_spec.append(float(p[1]))
                ax3.plot(Omukai_logn, Omukai_spec, color=cols[1], linestyle=':')

                Omukai_logn = []
                Omukai_spec = []

        ## OPEN OMUKAI 2005 FOR COMPARISON
        with open('waterData/logy(CO)_v_logn_Z=%d_Omukai2005.csv' % logZ[i], 'r') as omukaiData:
                lines = omukaiData.readlines()
                Omukai_logn = []
                Omukai_spec = []
                for line in lines[6:]:
                        p = line.split(',')
                        Omukai_logn.append(float(p[0]))
                        Omukai_spec.append(float(p[1]))
                ax3.plot(Omukai_logn, Omukai_spec, color=cols[2], linestyle=':')

                Omukai_logn = []
                Omukai_spec = []


        ax3.set_xlabel("log(n$_H$ [cm$^{-3}$])")
        ax3.set_ylabel(r"log y(C$^+$, C, CO)")

        xlim = [-1, 10]
        if (logZ[i] == -6): 
                ylim = [-15, -8]
        elif (logZ[i] == -4):
                ylim = [-13,-7]
        elif (logZ[i] == -2):
                ylim = [-11,-5]
        ax3.set_ylim(ylim)
        ax3.set_xlim(xlim)

                #text for C+, C, CO for different metallicities
        xC = [5, 0, 4.1]
        xCplus = [4.5, 4, 2.70]
        xCO = [8, 2.5, 7]

        yC = [-9.38, -9, -8]
        yCplus = [-11.25, -9, -9.5]
        yCO = [-11.25, -10.8, -6.5]

        ax3.text(xCplus[i], yCplus[i], 'C+', fontsize=16)
        ax3.text(xC[i], yC[i], 'C', fontsize=16)
        ax3.text(xCO[i], yCO[i], 'CO', fontsize=16)

       # ax3.legend(fancybox = True, labelspacing=0.0, loc='lower right', framealpha=0.5)

        f.subplots_adjust(hspace=0)
        pyplot.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        if my_chemistry.metal_cooling:
            output = "freefall_metal"
        else:
            output = "freefall"

        pyplot.savefig("%s_Water_abundances_%d.pdf" % (output, logZ[i]), dpi=1000)

