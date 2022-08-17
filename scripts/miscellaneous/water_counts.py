#!/usr/bin/env python

'''
This script will select & write solvent atoms (that are contained in
specific solvent shells) to a text file.
Likewise, the hydrogen bonds will be written to a seperate file.
Run the script in the terminal as an
 argument of the command 'chimera', for example:
chimera --nogui --script "script_1.py -args"
chimera --nogui --script "water_count.py"
** Note that runscript also works **
'''

# Loading Python Modules:{{{
print 'Importing Python Modules...'
import os,sys,commands,glob
import numpy as np
from chimera import runCommand as rc
import Midas as m
from chimera import NonChimeraError
from chimera.replyobj import reportException
import FindHBond
from Trajectory.formats.Gromacs import loadEnsemble
import Trajectory
import matplotlib ######################## Plotting
# Import Publication Style Fonts for Figures:
matplotlib.use('Agg')
fontfamily={'family':'sans-serif','sans-serif':['Arial']}
from matplotlib import pyplot as plt
plt.rc('font', **fontfamily)
print 'Done!'
# }}}

# Initial Conditions:{{{
# Directory for water selection the output files:
path_1 = '/volumes/WD_Passport_1TB/waters_per_shell'
# Directory for H-bond output files:
path_2 = '/volumes/WD_Passport_1TB/H_bonds'
# Setting the model number:
model = 0
# Setting the residue descriptor:
R = ['48.A','49.A','87.B','88.B','92.B','93.B']
#res = R[0]

# Shell Width & Shell location. Harrigan et al suggests 3 Angtrom width.
width = 3.0
inner = 0.001 # essentially the residue's location
shells = 4
loc = [inner,inner+width,inner+2.*width,inner+3.*width,inner+4.*width]
# List of trajectories:
trrs = ['RUN7_CLONE41','RUN9_CLONE70','RUN13_CLONE53','RUN15_CLONE32','RUN18_CLONE96','RUN24_CLONE68']
# Set the trajectory you are analyzing:
trr = trrs[1]

# }}}

# Directories{{{
# Root directory
root_dir = '/volumes/rmr_4tb/research/PROJ6381/chimera/'
# Directory for water selection the output files:
path_1 = str(root_dir)+'eigenvector_selection/waters_per_shell'

output_dir = '/Volumes/RMR_4TB/research/PROJ6381/chimera/Water_counts/%s'%trr
if os.path.exists(output_dir)==False:
    commands.getoutput('mkdir %s'%output_dir)

trr_dir='/Volumes/RMR_4TB/Research/PROJ6381/Movies/Traj_movie/trr/centered/TRRS/_Traj_%s.trr'%trr
# Setting directory for tpr
tpr_dir = '/Volumes/RMR_4TB/Research/PROJ6381/Movies/Traj_movie/tpr/frame0.tpr'
start = 'first'; end = 'last'

# Load in tpr and trr files:
gromacs_files = (tpr_dir,trr_dir)
def ensemble_callback(ensemble, keepLongBonds=True):
    	traj_ensemble = Trajectory.Ensemble(ensemble)
	traj_ensemble.CreateMolecule()
	for frame_num in range(1, len(ensemble)+1):
		if not traj_ensemble.findCoordSet(frame_num):
			traj_ensemble.LoadFrame(frame_num, makeCurrent=False)
	traj_ensemble.AddMolecule()
print 'Loading in tpr and trr files...'
loadEnsemble(gromacs_files, start, end, ensemble_callback)
print 'Done!'
# }}}

# Main:{{{

# Load trajectory with all frames:
print 'Taking frames from Models'
traj = m.chimera.openModels.list(id = model)[0]
traj.loadAllFrames()
frames = traj.coordSets # list of every frame

def counts():
    # append new shell to the list of shells
    for j in range(0,len(R)):
        res = R[j]
        print 'Computing Water Counts for %s...'%res
        SH,HB = [],[]
        for k in range(1, (shells+1)):
            sh,hb = [],[]
            for i in range(1, (len(frames)+1)):
                # Load trajectory frame
                traj.activeCoordSet = frames[i]
                # Select everything within the specified shell.
                rc('sel :%s za>%d,za<%d'%(res,loc[k-1],loc[k]))
                # Don't select any atoms other than water-oxygens
                rc('~sel :.A,.B')
                rc('~sel :.;~sel :.water@HW?')
                sel =  m.chimera.selection.currentAtoms()
                sh.append(len(sel))
                h_bonds = FindHBond.currentAtoms()
                hb.append(len(h_bonds))
                rc('sel:clear')
            SH.append(sh);HB.append(hb)
        Shells,H_bonds = np.array(SH),np.array(HB)
        np.save('%s/Shells_per_frame_%s.npy'%(output_dir,res),Shells)
        np.save('%s/H_Bonds_per_frame_%s.npy'%(output_dir,res),H_bonds)


# }}}

# Plot Water and Hydrogen Bond Counts:{{{
def plot():
    for j in range(0, len(R)):
        res = R[j]
        counts = np.load('%s/Shells_per_frame_%s.npy'%(output_dir,res))
        time = np.linspace(1,len(counts[0,:]),len(counts[0,:]))
        t = np.array(time)
        fig = plt.figure(figsize=(12,10))
        ax1 = fig.add_subplot(1,1,1)
        for i in range(0,len(counts[:,0])):
            ax1.plot(t, counts[i,:], color='k')
        #if res == "87.B":
        #    for k in range(0,len(counts[:,0])):
        #        if t[k] >= 125 and t[k] < 200:
        #            ax1.axvspan(t[k]-0.5, t[k]+0.5, fill=True, linewidth=0, color='y', alpha=0.3)
        #        if t[k] >= 200 and t[k] < 312:
        #            ax1.axvspan(t[k]-0.5, t[k]+0.5, fill=True, linewidth=0, color='b', alpha=0.3)
        ax1.set_xlabel('Time, t /(ns)', fontsize=22)
        ax1.set_ylabel('Water Counts', fontsize=22)
        #ax1.set_ylim(bottom=0, top=(np.max(counts[:,3])+25))
        ax1.set_xlim(left=0, right=len(t))
        ticks = [ax1.yaxis.get_minor_ticks(),
                 ax1.yaxis.get_major_ticks()]
        marks = [ax1.get_xticklabels(),
                 ax1.get_yticklabels()]
        for k in range(0,len(ticks)):
            for tick in ticks[k]:
                tick.label.set_fontsize(20)
        for k in range(0,len(marks)):
            for mark in marks[k]:
                mark.set_size(fontsize=20)
                mark.set_rotation(s=15)
        #fig.tight_layout()
        fig.savefig('%s/Water_counts_%s.pdf'%(output_dir,res))


# }}}

if __name__ == '__main__':
    counts()
    plot()

