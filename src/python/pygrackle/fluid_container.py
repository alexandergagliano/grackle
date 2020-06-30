########################################################################
#
# Grackle python fluid container
#
#
# Copyright (c) Enzo/Grackle Development Team. All rights reserved.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this
# software.
########################################################################

import numpy as np

from pygrackle.grackle_wrapper import \
    calculate_cooling_time, \
    calculate_gamma, \
    calculate_pressure, \
    calculate_temperature, \
    calculate_dust_temperature, \
    solve_chemistry

from pygrackle.utilities.misc import \
    issue_deprecation_warning

from pygrackle.utilities.physical_constants import \
    mass_hydrogen_cgs

_base_fluids = ["density", "metal", "dust"]
_nd_fields   = ["energy",
                "x-velocity", "y-velocity", "z-velocity",
                "temperature", "dust_temperature", "pressure",
                "gamma", "cooling_time", "mu", "nH",
                "mean_molecular_weight"]

_fluid_names = {}
_fluid_names[0] = _base_fluids
_fluid_names[1] = _fluid_names[0] + \
  ["HI", "HII", "HeI", "HeII", "HeIII", "de"]
_fluid_names[2] = _fluid_names[1] + \
  ["H2I", "H2II", "HM"]
_fluid_names[3] = _fluid_names[2] + \
  ["DI", "DII", "HDI"]

_water_names = \
  ["Water_density", "O_density", "OH_density", "O2_density", "Oplus_density", "OHplus_density", "H2Oplus_density", "H3Oplus_density", "O2plus_density", "Cplus_density", "C_density", "CH_density", "CH2_density", "CH3_density", "CH4_density", "CO_density", "COplus_density", "CO2_density"]

_bialy_names = \
["CHplus_density", "CH2plus_density", "H3plus_density", "HCOplus_density", "HeHplus_density", "CH3plus_density", "CH4plus_density", "CH5plus_density", "O2Hplus_density"]


_rad_trans_names = ['RT_heating_rate', 'RT_HI_ionization_rate',
                    'RT_HeI_ionization_rate', 'RT_HeII_ionization_rate',
                    'RT_H2_dissociation_rate']

_extra_fields = {}
_extra_fields[2] = ["H2_self_shielding_length"]
_extra_fields[3] = _extra_fields[2] + []

class FluidContainer(dict):
    def __init__(self, chemistry_data, n_vals, dtype="float64",
                 itype="int64"):
        super(FluidContainer, self).__init__()
        self.dtype = dtype
        self.chemistry_data = chemistry_data
        self.n_vals = n_vals
        for fluid in _fluid_names[self.chemistry_data.primordial_chemistry] + \
        _extra_fields.get(self.chemistry_data.primordial_chemistry, []) + \
        _nd_fields:
            self._setup_fluid(fluid)
        if self.chemistry_data.use_radiative_transfer:
            for fluid in _rad_trans_names:
                self._setup_fluid(fluid)
#Set up number density of hydrogen as well!
        self._setup_fluid("n_H")
        if self.chemistry_data.withWater:
            if (self.chemistry_data.water_rates == 3):
                for fluid in (_water_names + _bialy_names):
                   self._setup_fluid(fluid)
            else:
                for fluid in _water_names:
                   self._setup_fluid(fluid)
        for htype in ["specific", "volumetric"]:
            if getattr(self.chemistry_data, "use_%s_heating_rate" % htype, 0):
                self._setup_fluid("%s_heating_rate" % htype)

    def _setup_fluid(self, fluid_name):
        self[fluid_name] = np.zeros(self.n_vals, self.dtype)

    @property
    def cooling_units(self):
        warn = 'The cooling_units attribute is deprecated.\n' + \
          'For example, instead of fc.cooling_units, ' + \
          'use fc.chemistry_data.cooling_units.'
        issue_deprecation_warning(warn)
        return self.chemistry_data.cooling_units

    @property
    def density_fields(self):
        return _fluid_names[self.chemistry_data.primordial_chemistry]

    @property
    def n_density_fields(self):
        if self.chemistry_data.withWater:
            if (self.chemistry_data.water_rates == 3):
               return _water_names + _bialy_names
            else:
               return _water_names
        else:
            return []

    def calculate_hydrogen_number_density(self):
        my_chemistry = self.chemistry_data
        if my_chemistry.primordial_chemistry == 0:
            self["nH"] = my_chemistry.HydrogenFractionByMass * \
              self["density"] * my_chemistry.density_units / mass_hydrogen_cgs
            return
        nH = self["HI"] + self["HII"]
        if my_chemistry.primordial_chemistry > 1:
            nH += self["HM"] + self["H2I"] + self["H2II"]
        if my_chemistry.primordial_chemistry > 2:
            nH += self["HDI"] / 2.
        self["nH"] = nH * my_chemistry.density_units / mass_hydrogen_cgs

    def calculate_metallicity(self):
        my_chemistry = self.chemistry_data
        if my_chemistry.withWater == 0:
            return 0
        metal_density = 16.0*self["O_density"] + 17.0*self["OH_density"] + 18.0*self["Water_density"] + \
                32.0*self["O2_density"] + 16.0*self["Oplus_density"] + 17.0*self["OHplus_density"] + \
                18.0*self["H2Oplus_density"] + 19.0*self["H3Oplus_density"] + 32.0*self["O2plus_density"] + \
                12.0*self["Cplus_density"] + 12.0*self["C_density"] + 13.0*self["CH_density"] + \
                14.0*self["CH2_density"] + 15.0*self["CH3_density"] + 16.0*self["CH4_density"] + \
                28.0*self["CO_density"] + 28.0*self["COplus_density"] + 44.0*self["CO2_density"]
        if my_chemistry.water_rates == 3:
                metal_density += 13.0*self["CHplus_density"] + 14.0*self["CH2plus_density"] + \
                29.0*self["HCOplus_density"] + 15.0*self["CH3plus_density"]

        return (metal_density / (self["density"] * my_chemistry.density_units))

    def calculate_mean_molecular_weight(self):
        # If energy has been set, calculate mu from the energy
        if not (self["energy"] == 0).all():
            self.calculate_temperature()
            self.calculate_gamma()
            self["mu"] = self["temperature"] / \
                (self["energy"] * (self["gamma"] - 1.) *
                 self.chemistry_data.temperature_units)
            self["mean_molecular_weight"] = self["mu"]
            return
            
        # Default to mu=1
        self["mu"] = np.ones(self["energy"].size)
        self["mean_molecular_weight"] = self["mu"]

        if self.chemistry_data.primordial_chemistry == 0:
            return # mu=1

        # Check that (chemistry) density fields have been set.
        # Allow metals to be 0
        for field in self.density_fields:
            if field == 'metal':
                continue
            if (self[field] == 0).all():
                return

        # Calculate mu from the species densities; ignore deuterium
        nden = self["metal"]/16.
        nden += self["HI"]+self["HII"]+self["de"] + \
            (self["HeI"]+self["HeII"]+self["HeIII"])/4.
            
        if self.chemistry_data.primordial_chemistry > 1:
            nden += self["HM"]+(self["H2I"]+self["H2II"])/2.
            
        self["mu"] = self["density"]/nden
        self["mean_molecular_weight"] = self["mu"]

    def calculate_cooling_time(self):
        calculate_cooling_time(self)

    def calculate_gamma(self):
        calculate_gamma(self)

    def calculate_pressure(self):
        calculate_pressure(self)

    def calculate_temperature(self):
        calculate_temperature(self)

    def calculate_dust_temperature(self):
        calculate_dust_temperature(self)

    def solve_chemistry(self, dt):
        solve_chemistry(self, dt)
