# Molecular simulation scripts for bulk solutions

This dataset contains GROMACS files used in NMR Investigation of Water in Salt Crusts: Insights from Experiments and Molecular Simulations by Simon Gravelle, Sabina Haber-Pohlmeier, Carlos Mattea, Siegfried Stapf, Christian Holm, and Alexander Schlaich

The ff/ folder contains the force field data, the input/ folder the GROMACS mdp files, and the NaCl_1M/ and Na2SO4_1M folders are
two examples of configuration files for performing simulation at concentration 1M. 

Within both NaCl_1M/ and Na2SO4_1M folders, a Python file named ConfigurationGenerator allows to regenerate the initial conf.gro and topol.top files with a given salt concentration.

