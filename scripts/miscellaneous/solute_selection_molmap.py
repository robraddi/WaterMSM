#!/usr/bin/env python

"""
This script will take the most significant normalized
eigenvector selected residues and created a water "molmap" around them
and create a movie..

Run script in terminal as:
chimera --nogui --nostatus --script "scriptname.py"
"""

# Importing Modules:{{{
import os,sys,commands,glob
import numpy as np
from chimera import runCommand as rc
from chimera import replyobj
from chimera.tkgui import saveReplyLog
import Midas as m
#from chimera import NonChimeraError
#from chimera.replyobj import reportException
#import FindHBond

# }}}

# Initial Conditions:{{{
# Setting the model number:
model = 0
# Shell Width
width = 3.0
inner = 0.01 # essentially the residue's location
shells = 1
loc = [inner,inner+width,inner+2.*width]
# Setting the model number:
model = 0

# Centered residues:
res=['33.A','37.A','48.A','49.A','51.A','67.A','86.B','87.B','89.B','92.B','93.B']
# real eigenvector selection:

#res=[
#1469
#1471
#1483
#1505
#1514
#1569
#1571
#1591

# }}}

# Directories{{{
# Root directory
root_dir = '/volumes/rmr_4tb/research/PROJ6381/chimera/'
# Directory for water selection the output files:
path_1 = str(root_dir)+'eigenvector_selection/waters_per_shell'
# Grabbing all Trajectories
trrs = sorted(glob.glob('/Volumes/RMR_4TB/Research/PROJ6381/Movies/Traj_movie/trr/'))
trr='/Volumes/RMR_4TB/Research/PROJ6381/Movies/Traj_movie/trr/Traj_RUN9_CLONE70.trr'
# Setting directory for tpr
tpr = '/Volumes/RMR_4TB/Research/PROJ6381/Movies/Traj_movie/tpr'

# Trajectory to be made into a movie:

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

# Load in tpr and trr files:
proteins=chimera.openModels.open('%s %s'%(tpr,trr))


# Commands:{{{
rc('windowsize 650 550') # windowsize 200 150


# Load trajectory with all frames:
traj = m.chimera.openModels.list(id = model)[0]
traj.loadAllFrames()
frames = traj.coordSets # list of every frame
# Load trajectory with all frames:
frame = mdInfo['frame']
rc('movie record supersample 3;')

for j in range(0,len(frames)):
    rc('close #0.1')
    rc('color dark green :.B')
    rc('color blue :.A')
    rc('sel :.B za<3;# selectiing waters around p53 within 3 angstroms')
    rc('disp sel')
    rc('~sel @H??;')
    rc('~sel :.;')
    rc('~sel :.A,.B;')
    #rc('# you can change the spacing between density levels')
    #rc('#molmap gridSpacing 0.25')
    rc('molmap sel 10 center 0,0,0 replace true gridSpacing 0.25 displayThreshold 0.10 showDialog true')
    rc('transp 60,s')
    rc('sel invert;')
    rc('~disp sel;')
    rc('sel :.B za<3;')
    rc('~sel @H??;')
    rc('~sel :.;')
    rc('~sel :.A,.B;')
    rc('~sel')
    rc('disp :92.B')
    rc('color gray :92.B')


#rc('movie encode coordset.mp4 coordset.webm coordset.ogv bitrate 200')
rc('movie encode output ~/Desktop/m2.mov bitrate 200')



# }}}

# Notes for movie recording in command script:{{{

#movie record ;
## movie record with higher quality will need:
#movie record supersample 3
## For extremely high quality images you will want to do the following:
#movie record raytrace true
#
#
#turn y 3 120 ;
#wait 120 ; movie stop
#movie encode output ~/Desktop/m2.mov
#
# }}}





