`vol_accret` is a package for generating a reasonable mass abundance of major volatile species during the formmation and accretion of Earth-sized terrestrial planets. The model is resolved across three different profiles ("boxes") and in time.   that balances simplicity with realism, namely for the first-order modeling of the volatile content evolution in the atmospheres, interiors, and cores of solar system planets and exoplanets. Version 1 of this package was first built and used in Chen & Jacboson Earth and Planetary Science Letters (2022).

To acheive the goal of balancing between simplicity and realism, the calculations for  are idealized in   by assuming that:

1.
2.
3.

The main programs are entitled with "run_". They read in various output data files including aorig.dat, ABsizes.dat, MODEL_OUT_emb6.csv, MODEL_OUT_emb8.csv. It also calls other subprograms such as SF.py and henry.py that calculate the .  
This model was written with flexbility in mind and users can easily introduce new mechanismms, formalisms, and equations into the main program. It is relatively straightforward to blah blah blah.


Deciding which script to use will depend on the intent of the science goal. For instance, 


Within the main script, the user can specify the desired output and format, as well as....
