# Waternet in Grackle

*Waternet* is a molecular chemistry solver. It is written on top of Grackle, a chemistry and radiative cooling library, and can be run with cosmology simulations in Enzo (compatibility with Gizmo is pending). Waternet implements two different chemical networks:

   1. [Omukai et al. (2005)](https://iopscience.iop.org/article/10.1086/429955/meta) 

   2. [Bialy & Sternberg (2019)](https://ui.adsabs.harvard.edu/abs/2019arXiv190206764B/abstract)


Heating and cooling is done in base Grackle, which also allows for:

- photo-heating and photo-ionization from two UV backgrounds with optional
  self-shielding corrections:

   1. [Faucher-Giguere et al. (2009)](http://adsabs.harvard.edu/abs/2009ApJ...703.1416F).

   2. [Haardt & Madau (2012)](http://adsabs.harvard.edu/abs/2012ApJ...746..125H).

Installation of grackle's waternet is done in the same way as base Grackle, instructions for which can be found in the [online documentation](https://grackle.readthedocs.io/).

## Resources

- documentation: https://grackle.readthedocs.io/

- source code repository: https://github.com/grackle-project/grackle
