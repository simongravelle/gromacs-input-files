# Diffusion coefficient of bulk water

<img align="right" width="30%" src="bulkwater.jpg">

### Description

The simulation consists in a cubic box filled with water molecule.
The force field for the water is [tip4p/epsilon](https://doi.org/10.1021/jp410865y). 

### How to

Execute generate_system.py using python to generate the initial configuration.

Then, run gromacs using
```
gmx grompp -f input/nvt.mdp -o nvt -pp nvt -po nvt
gmx mdrun -v -deffnm nvt
mv nvt.gro conf.gro

gmx grompp -f input/npt.mdp -o npt -pp npt -po npt
gmx mdrun -v -deffnm npt
mv npt.gro conf.gro

gmx grompp -f input/run.mdp -o run -pp run -po run
gmx mdrun -v -deffnm run
```

### Diffusion coefficient measurement

#### Method 1 - using GROMACS

Execute the following command in the terminal:
```
gmx msd -f run.xtc -s run.tpr
```

#### Method 2 - using MDAnalysis

Convert the trajectory to remove the periodic boundary conditions using:
```
gmx trjconv -f run.xtc -o trajectory.xtc -pbc nojump
```

Using python, run

```
python diff_coeff_mda.py
```

### Contact

Feel free to contact me by email if you have inquiries. You can find contact details on my [personal page](https://simongravelle.github.io/).



