### Model Overview
`vgs` is the official name for the original Volatile Growth Simulator branch and is package for generating mass abundance predictions of major volatile species during the formmation and accretion of Earth-sized terrestrial planets. The model is resolved across three different profiles ("boxes") and in time. The goal of this tooklit is to provide an understanding of gas-phase species behavior for a given N-body dynamical simulation; namely, the first-order modeling of the volatile content evolution in the atmospheres, interiors, and cores of solar system planets and exoplanets. Version 1 of this package was first built and used in Chen & Jacboson Earth and Planetary Science Letters (2022).

To acheive the goal of balancing between simplicity and realism, the calculations for the volatile delivery, exchange, and loss are idealized in by assuming that:

1. Only three major volatile species are included in the mantle and atmosphere reservoirs; only carbon is included in the core.

3. The magma ocean depth is constant and persists across the entire duration of the accretion phase. Changes in MO depth and solidification time are not included.

3. The parent bodies are rapidly oxidized and the redox conditions are fixed throughout time.

4. Steam condensation is not included and no oceans are assumed to be present during impact events.

The main programs are named "run_.py". Differing only in the nature of the input and output files, these programs contains the same physical proceses and mechanisms as described the main paper. They read in N-body integrator outputs from Jacoboson & Morbidelli (2014) including aorig.dat (initial semi-major axis of the planetesimals and embryos), ABsizes.dat (planetesimal size distribution), MODEL_OUT_emb6.csv, MODEL_OUT_emb8.csv. It also calls other subprograms such as SF.py and henry.py that prescribes the volatile fraction for different chondritic materials and Henry's solubility coefficients for different gases.
The main- and sub-programs were written with flexbility in mind and it is relatively straightforward for users to introduce new mechanismms, formalisms, and equations into the main program. This can be done simply by replacing the existng expressions with the appropriate sections of the main program.


Deciding which script to use will depend on the intent of the science goal. For instance, if the goal is to obtain a reasonable esimate of the an N2-CO2-H2O-rich atmosphere right after the accretion phase, then one could use "run_nbody.py" to generate a suite of different total atmospheric pressure and N2-CO2-H2O partial pressures. The only modification needed from the original code is providing a new planetgrowth.out data for the specific system of interest. 


### Crediting This Work
Please reference this work by citing Chen & Jacobson (2022) (https://doi.org/10.1016/j.epsl.2022.117741) and stating the version used (`vgs.version`).

### Getting Started



### More Examples




P.S. The documentation and files will be regularly updated to fix issues or improve upon the existing version. Please check for discrepancies or email the author(s) for previous model versions.
