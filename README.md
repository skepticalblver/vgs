`vol_accret` is a package for generating a reasonable mass abundance  oThe simulated radiation is resolved in both wavelength and time. The intent is to provide consistent input for applications requiring time-dependent stellar UV radiation fields that balances simplicity with realism, namely for simulations of exoplanet atmospheres. Version 1 of this package was first built and used in Chen & Jacboson Earth and Planetary Science Letters (2022).

For this balance of simplicity and realism, the flares generated are idealized in the spectral and temporal distribution of their energy through the following assumptions:

This model was written with flexbility in mind and users can easily introduce new mechanismms, formalisms, and equations into the main program.

The main programs are entitled with "run_". They read in various output data files including aorig.dat, ABsizes.dat, MODEL_OUT_emb6.csv, MODEL_OUT_emb8.csv. It also calls other subprograms such as SF.py and henry.py. 


Deciding which script to use will depend on the intent of the science goal. For instance, 


Within the main script, the user can specify the desired output and format, as well as....
