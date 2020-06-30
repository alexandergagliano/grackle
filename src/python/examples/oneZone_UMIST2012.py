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
    N_pts = 10
    dtmax = 1.e12 # in seconds
    metallicity = np.logspace(-5,0, N_pts)
    #metallicity = [1.0]
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
        my_chemistry.crx_ionization = 1
        my_chemistry.metal_cooling = 1
        grackle_dir = os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.dirname(os.path.abspath(__file__)))))
        my_chemistry.grackle_data_file = os.sep.join(
            [grackle_dir, "input", "CloudyData_UVB=FG2011_shielded.h5"])
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
        initial_temperature = 100.
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
            fc["O_density"][:] = 4.9e-4 * metallicity[i] * fc["density"]
            fc["OH_density"][:] = tiny_number * fc["density"]
            fc["O2_density"][:] = tiny_number * fc["density"]
            fc["Oplus_density"][:] = tiny_number * fc["density"]
            fc["OHplus_density"][:] = tiny_number * fc["density"]
            fc["H2Oplus_density"][:] = tiny_number * fc["density"]
            fc["H3Oplus_density"][:] = tiny_number * fc["density"]
            fc["O2plus_density"][:] = tiny_number * fc["density"]
            fc["Cplus_density"][:] = tiny_number * fc["density"]
            fc["C_density"][:] = 2.9e-4 * metallicity[i] * fc["density"]
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
        pyplot.loglog(data["time"], (np.divide(data["O_density"], (data["O_density"] + data["Water_density"] + data["OH_density"] + data["O2_density"]))), label ='O', color='maroon',markevery=.1)
        pyplot.loglog(data["time"], (np.divide(data["C_density"],(data["Cplus_density"] + data["COplus_density"] + data["CHplus_density"] + data["C_density"] + data["CO_density"] + data["CH_density"]))), label ='C', color='black',markevery=.1)
        pyplot.loglog(data["time"], (np.divide(data["CO_density"],(data["Cplus_density"] + data["COplus_density"] + data["CHplus_density"] + data["C_density"] + data["CO_density"] + data["CH_density"]))), label ='CO', color='blue',markevery=.1)
        pyplot.loglog(data["time"], (np.divide(data["Water_density"],(data["O_density"] + data["Water_density"] + data["OH_density"] + data["O2_density"]))), label ='H2O', color='red',markevery=.1)
        pyplot.loglog(data["time"], (np.divide(data["OH_density"],(data["O_density"] + data["Water_density"] + data["OH_density"] + data["O2_density"]))), label ='OH', color='goldenrod',markevery=.1)
        pyplot.loglog(data["time"], (np.divide(data["O2_density"],(data["O_density"] + data["Water_density"] + data["OH_density"] + data["O2_density"]))), label ='O2', color='magenta',markevery=.1)
        pyplot.loglog(data["time"], (np.divide(data["CH_density"],(data["Cplus_density"] + data["COplus_density"] + data["CHplus_density"] + data["C_density"] + data["CO_density"] + data["CH_density"]))), label ='CH', color='teal',markevery=.1)
        pyplot.loglog(data["time"], (np.divide(data["Cplus_density"],(data["Cplus_density"] + data["COplus_density"] + data["CHplus_density"] + data["C_density"] + data["CO_density"] + data["CH_density"]))), label ='Cplus', ls="dashed",color='black',markevery=.1)
        pyplot.xlabel("Time (s)",fontsize=22)
        pyplot.ylabel("Xi",fontsize=22)
        pyplot.title("Convergence plot for [Z/H] = %i" % np.log10(metallicity[i]));
        pyplot.tight_layout()
        pyplot.xlim(left=1.e3)
        pyplot.ylim(ymin=1.e-17)
        pyplot.legend(loc='best',fontsize=4)
        pyplot.savefig("metal_convergence_t%iGyr_N%i_Z%i_UMIST.pdf" %(int(np.log10(final_time/1.e3)),N_pts, np.log10(metallicity[i])),dpi=300)
        pyplot.clf()

################ H2 CONVERGENCE PLOTS! ############################
        pyplot.loglog(data["time"], np.divide(data["H2I"]/my_chemistry.density_units, data["HI"]/my_chemistry.density_units + data["HM"]/my_chemistry.density_units + data["HII"]/my_chemistry.density_units + 2.0*(data["H2I"]/my_chemistry.density_units + data["H2II"]/my_chemistry.density_units)), label ='XH2', color='black')
        pyplot.loglog(data["time"], np.divide(data["HI"]/my_chemistry.density_units, data["HI"]/my_chemistry.density_units + data["HII"]/my_chemistry.density_units + data["HM"]/my_chemistry.density_units +  2.0*(data["H2I"]/my_chemistry.density_units + data["H2II"]/my_chemistry.density_units)), label ='XH', color='orange')
        pyplot.xlabel("Time (s)",fontsize=22)
        pyplot.ylabel("Xi",fontsize=22)
        pyplot.title("Convergence plot for [Z/H] = %i" % np.log10(metallicity[i]));
        pyplot.tight_layout()
        pyplot.xlim(left=1.e3)
        pyplot.ylim(ymin=1.e-2)
        pyplot.legend(loc='best',fontsize=4)
        pyplot.savefig("H2_convergence_t%iGyr_N%i_Z%i_UMIST.pdf" %(int(np.log10(final_time/1.e3)),N_pts, np.log10(metallicity[i])),dpi=300)
        pyplot.clf()

################ ION CONVERGENCE PLOTS! #########################
        pyplot.loglog(data["time"], np.divide(data["de"]/my_chemistry.density_units, data["HI"]/my_chemistry.density_units + data["HII"]/my_chemistry.density_units + 2.0*(data["H2I"]/my_chemistry.density_units + data["H2II"]/my_chemistry.density_units)), label ='e', color='red')

        pyplot.loglog(data["time"], np.divide(data["HII"]/my_chemistry.density_units, data["HI"]/my_chemistry.density_units + data["HII"]/my_chemistry.density_units + 2.0*(data["H2I"]/my_chemistry.density_units + data["H2II"]/my_chemistry.density_units)), label ='H+', color='blue')
        pyplot.xlabel("Time (s)",fontsize=22)
        pyplot.ylabel("Xi",fontsize=22)
        pyplot.title("Convergence plot for [Z/H] = %i" % np.log10(metallicity[i]));
        pyplot.tight_layout()
        pyplot.xlim(left=1.e3)
        pyplot.legend(loc='best',fontsize=4)
        pyplot.savefig("ion_convergence_t%iGyr_N%i_Z%i_UMIST.pdf" %(int(np.log10(final_time/1.e3)),N_pts, np.log10(metallicity[i])),dpi=300)
        pyplot.clf()

#####################################################################


        n = data["density"][-1]/my_chemistry.density_units
        n_tot.append(n)
        XO.append(data["O_density"][-1]/n)
        XC.append(data["C_density"][-1]/n)
        XCO.append(data["CO_density"][-1]/n)
        XH2O.append(data["Water_density"][-1]/n)
        XOH.append(data["OH_density"][-1]/n)
        XO2.append(data["O2_density"][-1]/n)
        XCH.append(data["CH_density"][-1]/n)

        XC_tot.append((data["C_density"][-1] + data["CO_density"][-1] + data["CH_density"][-1])/n)
        XO_tot.append((data["O_density"][-1] + data["Water_density"][-1] + data["OH_density"][-1] + data["O2_density"][-1])/n)

        H2.append(data["H2I"][-1]/my_chemistry.density_units)
        H2plus.append(data['H2II'][-1]/my_chemistry.density_units)
        H.append(data["HI"][-1]/my_chemistry.density_units)
        HM.append(data["HM"][-1]/my_chemistry.density_units)
        el.append(data["de"][-1]/my_chemistry.density_units)
        Hplus.append(data["HII"][-1]/my_chemistry.density_units)
        ions.append(data["Oplus_density"][-1] + data["Cplus_density"][-1] + data["H2Oplus_density"][-1] + data["H3Oplus_density"][-1] + data["O2plus_density"][-1] + data["COplus_density"][-1] + data["CHplus_density"][-1] + data["CH2plus_density"][-1] + data["H3plus_density"][-1] + data["HCOplus_density"][-1] + data["HeHplus_density"][-1] + data["CH3plus_density"][-1] + data["CH4plus_density"][-1] + data["CH5plus_density"][-1] + data["O2Hplus_density"][-1])

        XH2.append((data["H2I"][-1] + data["H2II"][-1])/n/my_chemistry.density_units)

        XH.append((data["HI"][-1] + data["HII"][-1])/n/my_chemistry.density_units)

        Xel.append(data["de"][-1]/n/my_chemistry.density_units)
        XHplus.append(data["HII"][-1]/n/my_chemistry.density_units)
        Xions.append((data["Oplus_density"][-1] + data["Cplus_density"][-1] + data["H2Oplus_density"][-1] + data["H3Oplus_density"][-1] + data["O2plus_density"][-1] + data["COplus_density"][-1] + data["CHplus_density"][-1] + data["CH2plus_density"][-1] + data["H3plus_density"][-1] + data["HCOplus_density"][-1] + data["HeHplus_density"][-1] + data["CH3plus_density"][-1] + data["CH4plus_density"][-1] + data["CH5plus_density"][-1] + data["O2Hplus_density"][-1])/n)

    XO = np.array(XO)
    XC = np.array(XC)
    XCO = np.array(XCO)
    XH2O = np.array(XH2O)
    XOH = np.array(XOH)
    XO2 = np.array(XO2)
    XCH = np.array(XCH)
    XO_tot = np.array(XO_tot)
    XC_tot = np.array(XC_tot)
    el = np.array(el)
    Hplus = np.array(Hplus)
    H = np.array(H)
    H2 = np.array(H2)
    H2plus = np.array(H2plus)
    HM = np.array(HM)
    ions = np.array(ions)
    Xel = np.array(Xel)
    XHplus = np.array(XHplus)
    XH = np.array(XH)
    XH2 = np.array(XH2)
    Xions = np.array(Xions)

        #pyplot.figure(dpi=500)

############################################################
#                                                          #
#   PLOT ONE: METAL MASS FRACTIONS VS METALLICITY          #
#                                                          #
############################################################
    if True:
        if my_chemistry.withWater == 1:
            ## OPEN FOR COMPARISON
            with open('fC_v_Z_UMIST2012.csv') as bialyData:
                lines = bialyData.readlines()
                bialy_Z = []
                bialy_C = []
                for line in lines[6:]:
                    p = line.split(',')
                    bialy_Z.append(float(p[0]))
                    bialy_C.append(float(p[1]))
                pyplot.plot(np.log10(bialy_Z), np.log10(bialy_C), color=cols_metl[1], linestyle=':')
                bialy_Z = []

            with open('fCH_v_Z_UMIST2012.csv') as bialyData:
                lines = bialyData.readlines()
                bialy_Z = []
                bialy_CH = []
                for line in lines[6:]:
                    p = line.split(',')
                    bialy_Z.append(float(p[0]))
                    bialy_CH.append(float(p[1]))
                pyplot.plot(np.log10(bialy_Z), np.log10(bialy_CH), color=cols_metl[6], linestyle=':')
                bialy_Z = []

            with open('fO_v_Z_UMIST2012.csv') as bialyData:
                lines = bialyData.readlines()
                bialy_Z = []
                bialy_O = []
                for line in lines[6:]:
                    p = line.split(',')
                    bialy_Z.append(float(p[0]))
                    bialy_O.append(float(p[1]))
                pyplot.plot(np.log10(bialy_Z), np.log10(bialy_O), color=cols_metl[0], linestyle=':')
                bialy_Z = []

                with open('fH2O_v_Z_UMIST2012.csv') as bialyData:
                    lines = bialyData.readlines()
                    bialy_Z = []
                    bialy_H2O = []
                    for line in lines[6:]:
                        p = line.split(',')
                        bialy_Z.append(float(p[0]))
                        bialy_H2O.append(float(p[1]))
                    pyplot.plot(np.log10(bialy_Z), np.log10(bialy_H2O), color=cols_metl[3], linestyle=':')
                    bialy_Z = []

                with open('fO2_v_Z_UMIST2012.csv') as bialyData:
                    lines = bialyData.readlines()
                    bialy_Z = []
                    bialy_O2 = []
                    for line in lines[6:]:
                        p = line.split(',')
                        bialy_Z.append(float(p[0]))
                        bialy_O2.append(float(p[1]))
                    pyplot.plot(np.log10(bialy_Z), np.log10(bialy_O2), color=cols_metl[5], linestyle=':')
                    bialy_Z = []

                with open('fCO_v_Z_UMIST2012.csv') as bialyData:
                    lines = bialyData.readlines()
                    bialy_Z = []
                    bialy_CO = []
                    for line in lines[6:]:
                        p = line.split(',')
                        bialy_Z.append(float(p[0]))
                        bialy_CO.append(float(p[1]))
                    pyplot.plot(np.log10(bialy_Z), np.log10(bialy_CO), color=cols_metl[2], linestyle=':')
                    bialy_Z = []

                with open('fOH_v_Z_UMIST2012.csv') as bialyData:
                    lines = bialyData.readlines()
                    bialy_Z = []
                    bialy_OH = []
                    for line in lines[6:]:
                        p = line.split(',')
                        bialy_Z.append(float(p[0]))
                        bialy_OH.append(float(p[1]))
                    pyplot.plot(np.log10(bialy_Z), np.log10(bialy_OH), color=cols_metl[4], linestyle=':')

                metals_final = pd.DataFrame({'Z':metallicity, 'O':np.divide(XO,XO_tot), 'C':np.divide(XC,XC_tot), 'CO':np.divide(XCO,XC_tot), 'H2O':np.divide(XH2O,XO_tot), 'OH':np.divide(XOH,XO_tot), 'O2':np.divide(XO2,XO_tot), 'CH':np.divide(XCH,XC_tot)})
                metals_final.to_csv("BialyTest_metals_%.2e.csv" % dtmax, index=False)

                pyplot.plot(np.log10(metallicity), np.log10(np.divide(XO,XO_tot)), label ='O', color=cols_metl[0])
                pyplot.plot(np.log10(metallicity), np.log10(np.divide(XC,XC_tot)), label ='C', color=cols_metl[1])
                pyplot.plot(np.log10(metallicity), np.log10(np.divide(XCO,XC_tot)), label ='CO', color=cols_metl[2])
                pyplot.plot(np.log10(metallicity), np.log10(np.divide(XH2O,XO_tot)),label ='H2O', color=cols_metl[3])
                pyplot.plot(np.log10(metallicity), np.log10(np.divide(XOH,XO_tot)), label ='OH', color=cols_metl[4])
                pyplot.plot(np.log10(metallicity), np.log10(np.divide(XO2,XO_tot)), label ='O2', color=cols_metl[5])
                pyplot.plot(np.log10(metallicity), np.log10(np.divide(XCH,XC_tot)), label ='CH', color=cols_metl[6])

                pyplot.xlabel("log Z'",fontsize=22)
                pyplot.ylabel(r"log f$_i$",fontsize=22)

                # Add pyplot.text for each of the species
                pyplot.text(-2, -0.7, 'C', fontsize=14)
                pyplot.text(-4.824, 0.176, 'O', fontsize=14)
                pyplot.text(-4, -8.7, r'O$_2$', fontsize=14)
                pyplot.text(-4.7, -4.52, r'H$_2$O', fontsize=14)
                pyplot.text(-4.82, -2.3, 'OH', fontsize=14)
                pyplot.text(-4.30, -7.3, 'CH', fontsize=14)
                pyplot.text(-0.82, -2.22, 'CO', fontsize=14)


                #leg = pyplot.legend(fancybox = True, labelspacing=0.0, loc='best')
                #leg.get_frame().set_alpha(0.5)
                pyplot.ylim(ymin=-10,top=1)
                pyplot.xlim(left=-5, right=0)
                pyplot.tight_layout()
                pyplot.savefig("p_oneZone_UMIST_metals_t%iGyr_N%i_UMIST.pdf" %(int(np.log10(final_time/1.e3)),N_pts),dpi=300)
                pyplot.clf()

############################################################
#                                                          #
#             PLOT TWO: H2, HI VS METALLICITY              #
#                                                          #
############################################################

            #pyplot.figure(dpi=1000)
            ## OPEN FOR COMPARISON
            with open('yH_v_Z_UMIST2012.csv') as bialyData:
                    lines = bialyData.readlines()
                    bialy_Z = []
                    bialy_H = []
                    for line in lines[6:]:
                            p = line.split(',')
                            bialy_Z.append(float(p[0]))
                            bialy_H.append(float(p[1]))
                    pyplot.plot(np.log10(bialy_Z), np.log10(bialy_H), label='H, Bialy', color=cols_h2[1], linestyle=':')
                    bialy_Z = []

            with open('yH2_v_Z_UMIST2012.csv') as bialyData:
                    lines = bialyData.readlines()
                    bialy_Z = []
                    bialy_H2 = []
                    for line in lines[6:]:
                            p = line.split(',')
                            bialy_Z.append(float(p[0]))
                            bialy_H2.append(float(p[1]))
                    pyplot.plot(np.log10(bialy_Z), np.log10(bialy_H2), label='H2, Bialy', color=cols_h2[0], linestyle=':')
                    bialy_Z = []

#            pyplot.loglog(metallicity, np.divide(H2 + H2plus, (2*(H2+H2plus) + H + HM + Hplus)), label ='H2', color='black')
#            pyplot.loglog(metallicity, np.divide(H, (2*(H2+H2plus) + H + HM + Hplus)), label ='H', color='coral')
            pyplot.plot(np.log10(metallicity), np.log10(np.divide(H2, 2*H2 + H)), label ='H2', color=cols_h2[0])
            pyplot.plot(np.log10(metallicity), np.log10(np.divide(H, 2*H2 + H)), label ='H',color=cols_h2[1])

            #create pandas dataframe
            H2_final = pd.DataFrame({'Z': metallicity, 'H2': np.divide(H2, 2*H2 + H), 'H': np.divide(H, 2*H2 + H)})
            H2_final.to_csv("BialyTest_H2_%.2e.csv"% dtmax, index=False)
            #pyplot.loglog(metallicity, XH, label='H', color='coral')
            #pyplot.loglog(metallicity, XH2, label='H2', color='black')

            pyplot.xlabel("log Z'",fontsize=22)
            pyplot.ylabel(r" log f$_i$",fontsize=22)

            # Add pyplot.text for each of the species 
            pyplot.text(-4.52, -1.15, r'H$_2$', fontsize=14)
            pyplot.text(-4.52, -0.22, 'H', fontsize=14)


           # leg = pyplot.legend(fancybox = True, labelspacing=0.0, loc='best')
           # leg.get_frame().set_alpha(0.5)
            pyplot.ylim(ymin=-2,top=0)
            pyplot.xlim(left=-5,right=0)
            pyplot.tight_layout()
            pyplot.savefig("p_oneZone_UMIST_H2_t%iGyr_N%i_UMIST.pdf" %(int(np.log10(final_time/1.e3)),N_pts),dpi=300)

            pyplot.clf()


############################################################
#                                                          #
#             PLOT THREE: H+, e VS METALLICITY             #
#                                                          #
############################################################

            #pyplot.figure(dpi=1000)
            pyplot.loglog(metallicity, np.divide(el, (2*(H2+H2plus) + H + Hplus + HM)), label ='e', color=col_ion[0])
            pyplot.loglog(metallicity, np.divide(Hplus, (2*(H2+H2plus) + H + Hplus + HM)), label ='H+', color=col_ion[1])
            pyplot.loglog(metallicity, np.divide(ions, (2*(H2+H2plus) + H + Hplus + HM)), label='ions', color=col_ion[2])
            pyplot.xlabel("Z'",fontsize=22)
            pyplot.ylabel(r"X$_i$",fontsize=22)
            pyplot.ylim(bottom=1.e-8)
            pyplot.xlim(left=1.e-3)

            #leg = pyplot.legend(fancybox = True, labelspacing=0.0, loc='best')
            #leg.get_frame().set_alpha(0.5)
            pyplot.tight_layout()
            pyplot.savefig("p_oneZone_UMIST_ions_t%iGyr_N%i_UMIST.pdf" %(int(np.log10(final_time/1.e3)),N_pts),dpi=300)
            ions_final = pd.DataFrame({'Z': metallicity, 'el': np.divide(el, (2*(H2+H2plus) + H + Hplus + HM)), 'Hplus': np.divide(Hplus, (2*(H2+H2plus) + H + Hplus + HM)), 'ions':np.divide(ions, (2*(H2+H2plus) + H + Hplus + HM))})
            ions_final.to_csv("BialyTest_ions_%.2e.csv"% dtmax, index=False)
