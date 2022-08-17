#!/usr/bin/env python

"""
This script will use the previously imported volmaps (.dx files) created from
VMD and use them with tIC1 which was obtained by solvent features.
Overall, this script will switch volmaps according to the
position of the metastable state transitions in order to create a movie.

Unfortunately, you are unable to run the script in terminal as:
chimera --nogui --nostatus --script "scriptname.py"

You must run the script inside MD Movie's "perframe script"
"""

# Importing Modules:{{{
import os,sys,commands,glob
import numpy as np
from chimera import runCommand as rc
from chimera import replyobj
from chimera.tkgui import saveReplyLog
import Midas as m
#import FindHBond

# }}}

# Initial Conditions:{{{
# Setting the model number:
model = 0
# }}}

# Directories{{{
# Root directory
root_dir = '/volumes/rmr_4tb/research/PROJ6381/chimera/'
# Directory for water selection the output files:
path_1 = str(root_dir)+'eigenvector_selection/waters_per_shell'
# Grabbing all Trajectories
#sol_dir='/Volumes/RMR_4TB/Research/PROJ6381/Old_Data/owlsnest_solvent/'
sol_dir='/Volumes/RMR_4TB/Research/PROJ6381/scratch/'
tica_sol = np.load(sol_dir+'solv_tica_distances.npy')

# Make a directories for our outputs:{{{
#if os.path.exists(path_1)==False:
#    commands.getoutput('mkdir %s'%path_1)
#if os.path.exists(path_2)==False:
#    commands.getoutput('mkdir %s'%path_2)
# Pdb directory
#work_dir = str(root_dir)+'/600_states/'
#pdb_dir = str(work_dir)+'*.pdb'
#pdbs = sorted(glob.glob(pdb_dir))
# }}}

# }}}

# Commands:{{{
rc('windowsize 650 550') # windowsize 200 150
# As frames go by run these commands:
i = mdInfo['frame']

# Reset the volume, v for every frame before its selected:
rc('volume #all hide')

#rc('close #0.1')
rc('color dark green :.B')
rc('color blue :.A')
# selectiing waters around p53 within 3 angstroms
rc('disp :.water@OW')
rc('~disp :.B za>3')
# Display all p53 residues:
rc('disp :87.B,88.B,92.B,93.B')
# Display all MDM2 residues:
rc('disp :48.A,49.A')
# Color all p53 resdiues:
rc('color cyan :87.B; color cornflower blue :93.B; color green :88.B; color yellow :92.B')
# Color all MDM2 residues:
rc('color orange :49.A; color purple :48.A')
# Take away the hydrogens on these residues
rc('~disp :.B@H??,.B@H?,.B@H???')
rc('~disp :.A@H???,.A@H??,.A@H?')
rc('~show :.')
# Color the oxygens on the waters red:
rc('color red :.water@OW')
# Select the atoms that you would like to use for the molmap:
rc('sel :.water@OW')
# Set the molmap with specific settings:

# Conditional statements for surfaces:{{{
t = [3309,511,1118,710,3057,1868]
name = ['RUN9_CLONE70','RUN13_CLONE53','RUN18_CLONE96',
        'RUN15_CLONE32','RUN7_CLONE41','RUN24_CLONE68']
n_surf = [ '#1', '#2', '#3', '#4', '#5', '#0.1', '#all']

# level = 0.5 or 0.2

x = tica_sol[t[0]][i-1,0]
if x > -1.1 and x <= -0.6:
    print i,str(n_surf[0])
    v = n_surf[0] # surface p5
    rc('volume %s style surface level 0.070 color gray step 1'%v)
    rc('volume %s show'%v)
elif x > -0.6 and x <= 0.0:
    print i,str(n_surf[1])
    v = n_surf[1] # surface p4
    rc('volume %s style surface level 0.070 color gray step 1'%v)
    rc('volume %s show'%v)
elif x > 0.0 and x <= 0.45:
    print i,str(n_surf[2])
    v = n_surf[2] # surface p3
    rc('volume %s style surface level 0.070 color gray step 1'%v)
    rc('volume %s show'%v)
elif x > 0.45 and x <= 1.0:
    print i,str(n_surf[3])
    v = n_surf[3] # surface p2
    rc('volume %s style surface level 0.070 color gray step 1'%v)
    rc('volume %s show'%v)
elif x > 1.0:
    print i,str(n_surf[4])
    v = n_surf[4] # surface p0
    rc('volume %s style surface level 0.070 color gray step 1'%v)
    rc('volume %s show'%v)
else:
    print i,str(n_surf[6])
    v = n_surf[6]
    rc('volume %s hide'%v)
    #print i,str(n_surf[5])
    #v = n_surf[5]
    #rc('molmap sel 7 center 0,0,0 replace true gridSpacing 1.3 displayThreshold 0.98950 showDialog true')
    #rc('volume %s style surface level 0.0297 color gray step 1'%v)
    #rc('volume %s show'%v)

# }}}

# Surface cutoff
# There is a surface on the inside that we want to see:
rc('sel :.A')
rc('sop zone %s sel 8'%v)
rc('~sel')
# Make the surface transparent:
rc('transp 60,s')
rc('~sel')

# }}}


