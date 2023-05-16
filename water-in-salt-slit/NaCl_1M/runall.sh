#!/bin/bash

# path to gromacs
gmx=/work/sgravelle/Softwares/gromacs-install/bin/gmx

${gmx} grompp -f ../input/nvt.mdp -o nvt -pp nvt -po nvt -r conf.gro
${gmx} mdrun -v -deffnm nvt
mv nvt.gro conf.gro

${gmx} grompp -f ../input/npt.mdp -o npt -pp npt -po npt -r conf.gro
${gmx} mdrun -v -deffnm npt
mv npt.gro conf.gro

${gmx} grompp -f ../input/prodHR.mdp -o prodHR -pp prodHR -po prodHR -r conf.gro
${gmx} mdrun -v -deffnm prodHR

${gmx} grompp -f ../input/prodLR.mdp -o prodLR -pp prodLR -po prodLR -r conf.gro
${gmx} mdrun -v -deffnm prodLR

