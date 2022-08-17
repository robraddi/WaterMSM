#!/usr/bin/env python


import os,sys,commands,glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
fontfamily={'family':'sans-serif','sans-serif':['Arial']}
from matplotlib import pyplot as plt
plt.rc('font', **fontfamily)
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy.interpolate import spline

# What is the selection?

selection='residue'
#selection='protein'


if selection=='protein':
    # Directory for water selection the output files:
    path_1 = './'

if selection=='residue':
    ## Directory for water selection the output files:
    path_1 = 'waters_per_shell'
    ## Directory for H-bond output files:
    path_2 = 'H_bonds'
    shell = np.loadtxt(str(path_1)+'/'+str(res[j])+str('/_Traj_1071_Waters_in_shell_%d_at_frame_%d.txt'%(i,k)))
    frames = glob.glob(str(path_1)+'/'+str(res[j])+str('/_Traj_1071_Waters_in_shell_1_at_frame_*.txt'))

    # Setting the residue descriptor:
res = ['1.A','38.A','47.A','67.A','69.A','70.A','72.A','75.A','93.B','94.B','95.B']
label = ['GLU25','MET62','GLN71','PHE91','VAL93','LYS94','HIS96','ILE99','LYS24','LEU25','LEU26']

def water_count():
    # For loop through each residue:
    for j in range(0,len(res)):
        # Counting the number of frames from the chimera outputs.
        frames = glob.glob(str(path_1)+'/'+str(res[j])+str('/_Traj_1071_Waters_in_shell_1_at_frame_*.txt'))
        num = len(frames)
        # Creating lists for each individual shell:
        traj_shell_1=[]
        traj_shell_2=[]
        traj_shell_3=[]
        traj_shell_4=[]
        for i in range(1,5):
            for k in range(1,num+1):
                shell = np.loadtxt(str(path_1)+'/'+str(res[j])+str('/_Traj_1071_Waters_in_shell_%d_at_frame_%d.txt'%(i,k)))

    frames = glob.glob(str(path_1)+str('_Traj_1071_Waters_in_shell_1_at_frame_*.txt'))
    num = len(frames)
    # Creating lists for each individual shell:
    traj_shell_1=[]
    traj_shell_2=[]
    traj_shell_3=[]
    traj_shell_4=[]
    for i in range(1,5):
        for k in range(1,num+1):
            shell = np.loadtxt(str(path_1)+str('_Traj_1071_Waters_in_shell_%d_at_frame_%d.txt'%(i,k)))
            print shell
            print (k,i)
            if i==1:
                if shell==[]:
                    traj_shell_1.append(len([]))
                else:
                    try:
                        traj_shell_1.append(len(shell))
                    except TypeError:
                        traj_shell_1.append(1)
            if i==2:
                if shell==[]:
                    traj_shell_2.append(len([]))
                else:
                    try:
                        traj_shell_2.append(len(shell))
                    except TypeError:
                        traj_shell_2.append(1)
            if i==3:
                if shell==[]:
                    traj_shell_3.append(len([]))
                else:
                    try:
                        traj_shell_3.append(len(shell))
                    except TypeError:
                        traj_shell_3.append(1)
            if i==4:
                if shell==[]:
                    traj_shell_4.append(len([]))
                else:
                    try:
                        traj_shell_4.append(len(shell))
                    except TypeError:
                        traj_shell_4.append(1)
    # Save the shell counts:
    shells_ = ['',traj_shell_1,traj_shell_2,traj_shell_3,traj_shell_4]
    for i in range(1,5):
        np.save('%s%sTraj_1071_shell_%d_%s.npy'%(path_1,'analysis/',i,'p53'), shells_[i])
    # Create a figure for each reside along the entire Trajectory:
    # Lists of the Max and Min of each of the shell counts:
    Max = []
    MaxPos = []
    Min = []
    MinPos = []
    for j in range(0,2):
        fig = plt.figure(figsize=(14,10))
        ax1 = fig.add_subplot(111)
        x  = [i for i in range(0,len(traj_shell_1))]
        y1 = [traj_shell_1[i] for i in range(0,len(traj_shell_1))]
        y2 = [traj_shell_2[i] for i in range(0,len(traj_shell_2))]
        y3 = [traj_shell_3[i] for i in range(0,len(traj_shell_3))]
        y4 = [traj_shell_4[i] for i in range(0,len(traj_shell_4))]
        a="$Shell_{1}=[0.01\AA-3.01\AA$]"
        b="$Shell_{2}=[3.01\AA-6.01\AA$]"
        c="$Shell_{3}=[6.01\AA-9.01\AA$]"
        d="$Shell_{4}=[9.01\AA-12.01\AA$]"
        if j==0:
            s = 'raw'
            ax1.plot(x, y1, 'k', label=a)#, marker='o')
            ax1.plot(x, y2, 'r', label=b)#, marker='o')
            ax1.plot(x, y3, 'b', label=c)#, marker='o')
            ax1.plot(x, y4, 'g', label=d)#, marker='o')
            y1_max = np.max(y1);Max.append(y1_max)
            y2_max = np.max(y2);Max.append(y2_max)
            y3_max = np.max(y3);Max.append(y3_max)
            y4_max = np.max(y4);Max.append(y4_max)
            y1_min = np.min(y1);Min.append(y1_min)
            y2_min = np.min(y2);Min.append(y2_min)
            y3_min = np.min(y3);Min.append(y3_min)
            y4_min = np.min(y4);Min.append(y4_min)
        if j==1:
            s = 'smoothed'
            x_smooth = np.linspace(np.min(x),np.max(x), 100)
            y1_smooth = spline(x, y1, x_smooth)
            y2_smooth = spline(x, y2, x_smooth)
            y3_smooth = spline(x, y3, x_smooth)
            y4_smooth = spline(x, y4, x_smooth)
            ax1.plot(x_smooth, y1_smooth, 'k', label=a)
            ax1.plot(x_smooth, y2_smooth, 'r', label=b)
            ax1.plot(x_smooth, y3_smooth, 'b', label=c)
            ax1.plot(x_smooth, y4_smooth, 'g', label=d)
        # Append the x positions of Minimum and Maximum to list
        y1_maxPos = y1.index(y1_max);MaxPos.append(y1_maxPos)
        y2_maxPos = y2.index(y2_max);MaxPos.append(y2_maxPos)
        y3_maxPos = y3.index(y3_max);MaxPos.append(y3_maxPos)
        y4_maxPos = y4.index(y4_max);MaxPos.append(y4_maxPos)
        y1_minPos = y1.index(y1_min);MinPos.append(y1_minPos)
        y2_minPos = y2.index(y2_min);MinPos.append(y2_minPos)
        y3_minPos = y3.index(y3_min);MinPos.append(y3_minPos)
        y4_minPos = y4.index(y4_min);MinPos.append(y4_minPos)
        plt.legend((a,b,c,d), loc=1, fontsize=18)#, facecolor="white",
        ax1.set_xlabel('Frames', fontsize=20)
        ax1.set_ylabel('Water Count', fontsize=20)
        ax1.set_xlim(left=-5, right=535)
        # Setting a y-limit
        y_lim = np.max(y4)
        ax1.set_ylim(bottom=0, top=(y_lim+75.))
        ax1.xaxis.get_minor_ticks()
        ax1.yaxis.get_minor_ticks()
        ax1.xaxis.get_major_ticks()
        ax1.yaxis.get_major_ticks()
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        minorLocator = MultipleLocator(100)
        # for the minor ticks, use no labels; default NullFormatter
        ax1.xaxis.set_minor_locator(minorLocator)
        # Annotations:
        if s=='raw':
            for k in range(0,4):
                # Maximum of each shell:
                ax1.annotate('(%d,%d)'%(MaxPos[k],Max[k]),
                        xy=(MaxPos[k],Max[k]),
                        fontsize=18,
                        textcoords='data',
                        horizontalalignment='left',
                        verticalalignment='center')
                # Minimum of each shell
                ax1.annotate('(%d,%d)'%(MinPos[k],Min[k]),
                        xy=(MinPos[k],Min[k]),
                        fontsize=18,
                        textcoords='data',
                        horizontalalignment='left',
                        verticalalignment='ceneter')
                print MaxPos[k],Max[k]
                print MinPos[k],Min[k]
        fig.tight_layout()
        fig.savefig('%s%sTraj_1071_%s_%s_%s.png'%(path_1,'analysis/','p53','shell_count',s))

def HBond():
    H_bond_shell_1=[]
    H_bond_shell_2=[]
    H_bond_shell_3=[]
    H_bond_shell_4=[]

#frames = glob.glob('/volumes/WD_Passport_1TB/H_bonds/Traj_1071_HBonds_in_shell_1_at_frame_*.txt')
    for i in range(1,shells+1):
        for k in range(1,len(frames)+1):
            shell = np.loadtxt('/Traj_1071_HBonds_in_shell_%d_at_frame_%d.txt'%(i,k))
            pass
        pass
    pass



if __name__ == '__main__':
    water_count()
    #HBond()




