import MDAnalysis as mda
import MDAnalysis.analysis.msd as msd
import numpy as np

u=mda.Universe('conf.gro', 'trajectory.xtc');

MSD=msd.EinsteinMSD(u, select='name OW1', msd_type='xyz', fft=True)
MSD.run();
ResMSDXY=MSD.timeseries;
nframes=MSD.n_frames;
timestep=np.round(u.trajectory.dt, 2); # actual time between frames, ps
lagtimes=np.arange(nframes)*timestep # make the lag-time axis
z = np.polyfit(lagtimes[int(nframes*0.1):int(nframes*0.9)], ResMSDXY[int(nframes*0.1):int(nframes*0.9)], 1)
d=3; # number of dimensions
D=z[0]/(2*d); # A^2/ps
print('The diffusion coefficient is '+str(D*1e-16/1e-12)+' cm^2/s')
