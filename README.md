### Model Overview
`vol_accretor` is a package for generating mass abundance predictions of major volatile species during the formmation and accretion of Earth-sized terrestrial planets. The model is resolved across three different profiles ("boxes") and in time. The goal of this tooklit is to provide an understanding of gas-phase species behavior for a given N-body dynamical simulation; namely, the first-order modeling of the volatile content evolution in the atmospheres, interiors, and cores of solar system planets and exoplanets. Version 1 of this package was first built and used in Chen & Jacboson Earth and Planetary Science Letters (2022).

To acheive the goal of balancing between simplicity and realism, the calculations for the volatile delivery, exchange, and loss are idealized in by assuming that:

1.
2.
3.

The main programs are named "run_.py". Differing solely in the characteristics of , these programs read in various N-body simulation output files including aorig.dat (), ABsizes.dat, MODEL_OUT_emb6.csv, MODEL_OUT_emb8.csv. It also calls other subprograms such as SF.py and henry.py that calculate the .  
The main- and sub-programs were written with flexbility in mind and it is relatively straightforward for users to introduce new mechanismms, formalisms, and equations into the main program. This can be done simply by.


Deciding which script to use will depend on the intent of the science goal. For instance, if the goal is to obtain a reasonable esimate of the an N2-CO2-H2O-rich atmosphere right after the accretion phase, then one could use "run_nbody.py" to generate a suite of different total atmospheric pressure and N2-CO2-H2O partial pressures. The only modification needed from the original code is providing a new planetgrowth.out data for the specific system of interest. 

Within the main script, the user can specify the desired output and format, as well as....


### Crediting This Work
Please reference this work by citing Chen & Jacobson () and stating the version used (`vol_accretor.version`).

### Getting Started



### More Examples
