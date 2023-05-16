#!/bin/bash

# path to gromacs
gmx=gmx

${gmx} grompp -f ../input/nvt.mdp -o nvt -pp nvt -po nvt -r conf.gro
${gmx} mdrun -v -deffnm nvt -nt 12
mv nvt.gro conf.gro

${gmx} grompp -f ../input/npt.mdp -o npt -pp npt -po npt -r conf.gro
${gmx} mdrun -v -deffnm npt -nt 12
mv npt.gro conf.gro

${gmx} grompp -f ../input/prodHR.mdp -o prodHR -pp prodHR -po prodHR -r conf.gro
${gmx} mdrun -v -deffnm prodHR -nt 12

${gmx} grompp -f ../input/prodLR.mdp -o prodLR -pp prodLR -po prodLR -r conf.gro
${gmx} mdrun -v -deffnm prodLR -nt 12

