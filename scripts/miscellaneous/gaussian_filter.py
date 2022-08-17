#!/usr/bin/env python

"""
This script will create a water "molmap" for every frame
and create a movie..
Run script in chimera GUI MD Movie " perframe script"
This script cannot be used with nogui because you need to make a selection and
hold it steady.
"""

# Importing Modules:{{{
import os,sys,commands,glob
import numpy as np
from chimera import runCommand as rc
from chimera import replyobj
from chimera.tkgui import saveReplyLog
import Midas as m
# }}}

# Initial Conditions:{{{

# What do you want performed?
job = 'pngs'
# Setting the model number:
model = 0
lv = 0.015
sample = 3
# }}}

# Directories{{{
# Root directory
root_dir = './'
# Directory for water selection the output files:
path_1 = str(root_dir)+'eigenvector_selection/waters_per_shell'
# List of trajectories:
trrs = ['RUN7_CLONE41','RUN9_CLONE70','RUN13_CLONE53','RUN15_CLONE32','RUN18_CLONE96','RUN24_CLONE68']
# Set the trajectory you are analyzing:
trr = trrs[2]
output_dir = 'per_frame_pngs/%s'%trr
if os.path.exists(output_dir)==False:
    commands.getoutput('mkdir %s'%output_dir)
# }}}

# Commands:{{{
# Set the GUI window size:
rc('windowsize 650 550') # windowsize 200 150
# As frames go by run these commands:
i = mdInfo['frame']
# Close the molmaps from the previous frame:
rc('close #0.1')
rc('close #1')
# Color your chains:
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
rc('molmap sel 7 center 0,0,0 replace true gridSpacing 1.3 displayThreshold 0.98950 showDialog true')
# Setting the level of the latter molmap:
rc('volume #0.1 style surface color gray level %s step 1'%lv)
# Perform Gaussian smoothing on the latter molmap
rc('vop gaussian #0.1')
# Hide the old molmap and only show the Gaussian:
rc('volume #0.1 hide')
# Set the Gaussian filtered surface to a step of 1 to get all the features:
rc('volume #1 style surface color gray step 1')
# Surface cutoff
# There is a surface on the inside that we want to see:
rc('sel :.A')
rc('sop zone #1 sel 8')
rc('~sel')
# Make the surface transparent:
rc('transp 60,s')
#rc('~sel')
# Set the hydrogen bond line width along with the color:
rc('hbonds lineWidth 1.0')
rc('hbonds color #0000ccccdfff')
# Hide the hbonds within MDM2 and p53
rc('hbonds intraMol false')
rc('~sel')
if job == 'pngs':
    rc('copy file %s/%s_frame_%d.png supersample %s'%(
        output_dir,trr,i,sample))




# }}}


