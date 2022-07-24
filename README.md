`vol_accret` is a package for generating a reasonable approximation of the UV emission of M dwarf stars over a single flare or a series of them. The simulated radiation is resolved in both wavelength and time. The intent is to provide consistent input for applications requiring time-dependent stellar UV radiation fields that balances simplicity with realism, namely for simulations of exoplanet atmospheres. Version 1 of this package was first built and used in Chen & Jacboson Earth and Planetary Science Letters (2022).

For this balance of simplicity and realism, the flares generated are idealized in the spectral and temporal distribution of their energy through the following assumptions:

This model was written with flexbility in mind and users can easily introduce new mechanismms, formalisms, and equations into the main program.

The main programs are entitled with "run_" that executes the same calculation with different scenarios and assumptions. For instance.

Within the main script, the user can specify the desired output and format, as well as....
