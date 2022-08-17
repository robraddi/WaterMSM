
#!/user/bin/env python
"""
For the movie making process requirements are the use of VMD or
Chimera, as well as matplotlib.
This is the script that is used only after you have acquired
png's for every frame of the trajectory.
"""

# Importing Modules:{{{
import os,sys,commands,glob
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
matplotlib.use('Agg')
fontfamily={'family':'sans-serif','sans-serif':['Arial']}
plt.rc('font', **fontfamily)
from moviepy.editor import VideoClip #########
from moviepy.video.io.bindings import mplfig_to_npimage #############
import matplotlib.image as mpimg #############
import moviepy as mp #########################
from pylab import *  #########################
from moviepy.editor import * #################
import matplotlib.gridspec as gridspec
import mdtraj as md
import msmbuilder.utils
from sklearn.externals import joblib

# }}}

# Initial Conditions:{{{

# Set the number of frames for your trajectory:
# t = index for trajectory in tica
t = [3309,511,1118,710,3057,1868]
# List of trajectories:
trrs = ['RUN7_CLONE41','RUN9_CLONE70','RUN13_CLONE53','RUN15_CLONE32','RUN18_CLONE96','RUN24_CLONE68']
trr = trrs[0]

# }}}

#Directories{{{
# Root directory
root_dir = '/volumes/rmr_4tb/research/PROJ6381/chimera/'
# Directory for water selection the output files:
path_1 = str(root_dir)+'eigenvector_selection/waters_per_shell'
# Set the trajectory you are analyzing:
output_dir = '/Volumes/RMR_4TB/research/PROJ6381/chimera/movies/per_frame_pngs/%s'%trr
tica_dir = '/Volumes/RMR_4TB/Research/PROJ6381/Movies/tICA_every_frame/%s'%trr
concat_dir = '/Volumes/RMR_4TB/research/PROJ6381/chimera/movies/concat_subplots/%s'%trr
movie_dir = '/Volumes/RMR_4TB/Research/PROJ6381/chimera/movies'
if os.path.exists(concat_dir)==False:
    commands.getoutput('mkdir %s'%concat_dir)
# }}}

# Create a Unified Image of 2x2 Subplots:{{{
def two_by_two():
    for i in range(1,2):#num+1):
        fig = plt.figure(figsize=(18,10))
        # Two plots on top and protein image on bottom.
        gs = gridspec.GridSpec(2, 2,
                               width_ratios=[1, 1.25],
                               height_ratios=[4, 2.5]
                               )
        ax1 = fig.add_subplot(2,2,1) # hbonds
        ax2 = fig.add_subplot(2,2,3) # water count
        ax3 = fig.add_subplot(gs[1]) # protein
        ax4 = fig.add_subplot(gs[3]) # tica
        # Path to images for every frame.
        hbonds=str('/volumes/WD_Passport_1TB/PROJ6381/TRAJ_%s/H_bonds/analysis/Traj_%s_HBonds_%s_%s_smoothed_frame_%d.png'%(traj,traj,label[j],res[j],i))
        waters=str('/volumes/WD_Passport_1TB/PROJ6381/TRAJ_%s/waters_per_shell/analysis/Traj_%s_Waters_%s_%s_smoothed_frame_%d.png'%(traj,traj,label[j],res[j],i))
        protein=str('/volumes/WD_Passport_1TB/PROJ6381/TRAJ_%s/movie_frames/%s/%s_frame_%d.png'%(traj,label[j],res[j],i))
        tica=str('/volumes/WD_Passport_1TB/PROJ6381/TRAJ_1071/tICA_frames/tICA_solvent_%s_frame_%d.png'%(traj_name[k],i))
        img1,img2=mpimg.imread(hbonds),mpimg.imread(waters)
        img3,img4=mpimg.imread(protein),mpimg.imread(tica)
        imgplot1,imgplot2 = ax1.imshow(img1),ax2.imshow(img2)
        imgplot3,imgplot4 = ax3.imshow(img3),ax4.imshow(img4)
        ax3.set_axis_bgcolor('none')
        ax4.set_axis_bgcolor('none')
        Subplots = [ax1,ax2,ax3,ax4]
        for sub in Subplots:
            sub.set_axis_off()
        fig.tight_layout()
        fig.savefig('/volumes/WD_Passport_1TB/PROJ6381/TRAJ_%d/unified_images/Traj_%d_%s_%s_frame_%d.png'%(traj,traj,label[j],res[j],i), edgecolor='black') #, transparent=True)


    ## 2x2 Subplot in Matplotlib:
    #for i in range(1,n_frames+1):
    ## You can probably try figsize=(12,10)...
    #    fig = plt.figure(figsize=(14,10))
    #    ax1 = fig.add_subplot(2,2,2)
    #    ax2 = fig.add_subplot(2,2,1)
    #    ax3 = fig.add_subplot(2,1,2)
    #    hbonds=str('/volumes/WD_Passport_1TB/PROJ6381/TRAJ_1071/H_bonds/analysis/Traj_1071_HBonds_THR18_87.B_raw.png')
    #    waters=str('/volumes/WD_Passport_1TB/PROJ6381/TRAJ_1071/waters_per_shell/analysis/Traj_1071_Waters_THR18_87.B_raw.png')
    #    protein=str('/volumes/WD_Passport_1TB/PROJ6381/TRAJ_1071/Traj_1072_red_87B/frame_%d.png'%i)
    #    img1=mpimg.imread(hbonds)
    #    img2=mpimg.imread(waters)
    #    img3=mpimg.imread(protein)
    #    imgplot1 = ax1.imshow(img1)
    #    imgplot2 = ax2.imshow(img2)
    #    imgplot3 = ax3.imshow(img3)
    #    ax1.set_alpha(alpha=0.0)
    #    ax2.set_alpha(alpha=0.0)
    #    ax3.set_alpha(alpha=0.0)
    #    ax1.patch.set_alpha(0.0)
    #    ax2.patch.set_alpha(0.0)
    #    ax3.patch.set_alpha(0.0)
    #    ax1.set_axis_off()
    #    ax2.set_axis_off()
    #    ax3.set_axis_off()
    #    fig.tight_layout()
    #    fig.savefig('/volumes/WD_Passport_1TB/PROJ6381/TRAJ_1071/unified_images/test__%d.png'%i, edgecolor='black')#, transparent=True)

# }}}

# Create a Unified Image of 2x1 Subplots:{{{
def two_by_one():
    n_frames = glob.glob(str(output_dir)+str('/*.png'))
    for i in range(1,len(n_frames)+1):#num+1):
        dpi = 80
        protein=str('%s/%s_frame_%d.png'%(output_dir,trr,i))
        tica=str('%s/tICA_%s_frame_%d.png'%(tica_dir,trr,i))
        img1,img2=mpimg.imread(protein),mpimg.imread(tica)
        height1, width1, nbands1 = img1.shape
        height2, width2, nbands2 = img1.shape
        figsize1 = width1 / float(dpi), height1 / float(dpi)
        figsize2 = width2 / float(dpi), height2 / float(dpi)
        FigSize = (figsize1[0] + figsize2[0],figsize1[1] + figsize2[1])
        # What size does the figure need to be in inches to fit the image?
        fig = plt.figure(figsize=(FigSize))
        ax1 = fig.add_subplot(211) # protein
        ax2 = fig.add_subplot(212) # tica
        ax1.text(-0.4, 4.5, fontsize=16,s='%s'%str(trr))
        # Path to images for every frame.
        imgplot1,imgplot2 = ax1.imshow(img1,interpolation='nearest'),ax2.imshow(img2,interpolation='nearest')
        ax1.set_axis_bgcolor('none')
        ax2.set_axis_bgcolor('none')
        ax1.set(xlim=[0, width1], ylim=[height1, 0], aspect=1)
        ax2.set(xlim=[0, width2*2], ylim=[height2*2, 0], aspect=1)
        Subplots = [ax1,ax2]
        for sub in Subplots:
            sub.set_axis_off()
        fig.tight_layout()
        fig.savefig('%s/tICA_%s_frame_%d.png'%(concat_dir,trr,i), edgecolor='black')

# }}}

# Concatenation of Unified Subplots --> Movie Creation:{{{
# Note that Hongbin suggests "matplotlib.animation.FFMpegWriter" for a nice module to use.
# Unfortunately, I couldn't get it to work.
def concat_frames():
    img=[]
    n_frames = glob.glob(str(concat_dir)+str('/*.png'))
    for i in range(0,len(n_frames)):
        l = '%s/tICA_%s_frame_%d.png'%(concat_dir,trr,i+1)
        #print l
        img.append(l)
    # Set duration and # of frames per second to match that
    # of most publication limits.
    # These setting are equivalent to a movie duration of 17-18 seconds.
    clips = [mp.editor.ImageClip(m).set_duration('00:00:0.0333') for m in img]
    concat_clip = mp.editor.concatenate_videoclips(clips, method="compose")
    # Save as mp4
    concat_clip.write_videofile('%s/tICA_%s.mp4'%(movie_dir,trr), fps=30)

# }}}

if __name__ == "__main__":
    #two_by_two()
    two_by_one()
    concat_frames()

