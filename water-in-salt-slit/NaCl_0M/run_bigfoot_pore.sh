#!/bin/bash
#OAR -n mixture-nanopore
#OAR -l /nodes=1/gpu=1/cpu=1/core=8,walltime=48:00:00
#OAR -p gpumodel='A100'
#OAR --stdout log.out
#OAR --stderr log.err
#OAR --project tamtam

export GMX_MAXBACKUP=-1

gmx=/home/gravells/softwares/gromacs-2023/build-gpu/bin/gmx

#${gmx} grompp -f ../input/nvt.mdp -p topol.top -o nvt -pp nvt -po nvt -r conf.gro
#${gmx} mdrun -deffnm nvt -v -rdd 1 -nt 8 -pin on
#cp nvt.gro conf.gro

${gmx} grompp -f ../input/npt.mdp -p topol.top -o npt -pp npt -po npt -r conf.gro -maxwarn 2
${gmx} mdrun -deffnm npt -v -rdd 1 -nt 8 -pin on
cp npt.gro conf.gro

