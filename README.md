### Model Overview
`vgs` is the official repository for the original Volatile Growth Simulator, a modeling framework for simulating the delivery, partitioning, and loss of major volatile species (e.g., H₂O, CO₂, N₂) during the formation of terrestrial planets. This package is designed to generate mass abundance predictions for key volatiles across time, using N-body-derived accretion histories as input. VGS tracks these volatile species through three interconnected planetary reservoirs—atmosphere, mantle, and core—while accounting for both delivery (e.g., from chondritic impactors) and depletion mechanisms (e.g., impact erosion, ingassing/degassing, and hydrodynamic escape).

Version 1 of this toolkit was first implemented in Chen & Jacobson (2022, EPSL), and a modified version was recently employed in Chen et al. (2025, ApJL) to assess the volatile evolution of TRAPPIST-1 analog systems.


### Key Features

1. Time-resolved modeling of volatile transport throughout planetary accretion, including both stochastic impact events and secular processes.

2. Three-box architecture representing the atmosphere, mantle, and core, enabling explicit tracking of volatile exchange and retention across reservoirs.

3. Incorporates physics for:

      a. Giant impact–induced atmospheric loss (Schlichting et al. 2015 formalism)

      b. Energy-limited hydrodynamic escape

      c. Henry’s law partitioning between atmosphere and magma ocean

      d. Carbon sequestration during core formation (based on Deguen et al. 2011)

4. Compatible with Mercury6-based N-body simulations (Jacobson & Morbidelli 2014), including aorig.dat, planet_growth.out, and other collision tracking files.

5. Prescribes volatile input compositions based on cosmochemically informed chondritic gradients (e.g., E-type, S-type, C-type).

###Scientific Applications

VGS was designed to support a wide range of science goals related to early planetary evolution and habitability:

  1. Estimating volatile budgets of Earth-sized planets during and after the accretion phase.

  2. Simulating TRAPPIST-1 analog systems, including the influence of stellar luminosity evolution and disk composition on volatile outcomes.

  3. Exploring volatile formation pathways in compact M-dwarf systems.

  4. Producing initial conditions for atmosphere and climate models of rocky exoplanets.

As demonstrated in Chen et al. (2025), VGS has been used to reproduce a broad spectrum of volatile outcomes, from dry, airless worlds to ocean-rich planets, based on a suite of more than 600 N-body simulations.

### Assumptions & Limitations

To balance realism and tractability, the model adopts several idealizations:

  1. Only three volatiles (H₂O, CO₂, N₂) are tracked in the atmosphere and mantle; only carbon is tracked in the core.

  2. The magma ocean is assumed to be globally molten and constant in depth for the duration of accretion.

  3. Steam condensation and ocean formation are not explicitly modeled.

  4. Redox state and oxidation efficiency of impactors are held fixed throughout time.

  5. Surface temperature is approximated as isothermal for volatile solubility calculations.

### Crediting This Work
Please reference this work by citing Chen & Jacobson (2022) (https://doi.org/10.1016/j.epsl.2022.117741) and stating the version used (`vgs.version`).


the main program is run_model.py. It reads in various output data files including aorig.dat, ABsizes.dat, MODEL_OUT_emb6.csv, MODEL_OUT_emb8.csv. It also calls other subprograms such as SF.py and henry.py. To begin, simply change the desired parameters in the file and execute.
### Citing VGS

If you use this code in your work, please cite:

Chen & Jacobson (2022), Earth and Planetary Science Letters, 594, 117741 (https://doi.org/10.1016/j.epsl.2022.117741)

Chen et al. (2025), Astrophysical Journal Letters, 

Also include the version number (e.g., vgs.version) used in your analysis.

### Getting Started

Execute run_trappist.py to begin simulations tailored for initially TRAPPIST-1-like disks. This main script reads in N-body output files and calls supporting modules like SF.py and henry.py. It generates time-resolved predictions for surface and interior volatile content throughout accretion. To apply the model to a new planetary system (e.g., Kepler-186 or a custom compact analog), simply replace the input planet_growth.out file with one derived from your own N-body simulation, adjust parameters as needed, and run the script.

Documentation is under active development and may be updated periodically. For questions or prior versions of the model, please contact the repository maintainer.

P.S. The documentation and files will be regularly updated to fix issues or improve upon the existing version. Please check for discrepancies or email the author(s) for previous model versions. Contact hchen@fit.edu for any questions, concerns, or requests for additional data.
