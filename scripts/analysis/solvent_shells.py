#!/usr/bin/env python

"""
Uses Python 2

This python script contains code written by M.P.Harrigan.
https://github.com/mpharrigan/wetmsm/tree/master/wetmsm

I wanted custom plotting
Edits last made: May 2015
"""


# Importing Libraries:{{{
import warnings
warnings.filterwarnings("ignore")
import matplotlib
import mdtraj as md
import msmbuilder.utils
import numpy as np
from itertools import combinations
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from csv import reader
# }}}

# Translate component to two dimensions{{{
def translate_component(component, n_solute, n_shells, deleted=None):
    """Take 1d component from tICA/PCA and expand to (solute, shell) indexing.
    see: https://github.com/mpharrigan/wetmsm/blob/master/wetmsm/vmd_write.py

    Parameters
    ----------
    component : array_like, shape=(n_solute * n_shells - n_deleted)
        1-d loadings (from tICA/PCA) for transformation
    n_solute : int
        Number of solute atoms
    n_shells : int
        Number of shells
    deleted : array_like, shape=(n_deleted), default=None
        Indices (1-d) of features that were removed (likely due to
        low-variance) before performing tICA.
    Returns
    -------
    component2d : np.ndarray, shape=(n_solute, n_shells)
        2-d component.
    """
    component2d = np.zeros((n_solute, n_shells))
    if deleted is None:
        deleted = np.asarray([])

    assert n_solute * n_shells - len(deleted) == len(component)

    absi = 0
    pruni = 0
    for ute in range(n_solute):
        for sh in range(n_shells):
            if not np.in1d(absi, deleted):
                component2d[ute, sh] = component[pruni]
                pruni += 1
            else:
                component2d[ute, sh] = 0.0

            absi += 1

    return component2d
# }}}

# Unpack Component:{{{
# Taken from M.PHarrigans Github
def unpack_component(component):
    """Translate component to two dimensions
    see: https://github.com/mpharrigan/wetmsm/blob/master/wetmsm/vmd_write.py
    """
    if component.ndim == 1:
        n_solute = len(solute_indices)
        n_shells = len(component) // n_solute
        assert len(component) % n_solute == 0
        component2d = translate_component(component,n_solute,n_shells)
    else:
        component2d = component
    return component2d

# }}}

# Main:{{{
if __name__ == "__main__":

    # Set your working directory:
    dir='./'
    ############################################################
    solute_indices = np.loadtxt(str(dir)+'AtomIndices.dat').astype(int)
    solvent_indices = np.loadtxt(str(dir)+'solvent.txt').astype(int)
    atom_sel = str(dir)+'1ycr_indices.dat'

    component = np.load("tICA_Eigenvectors.npy")[:,0] # first Component
    n_solute = len(solute_indices)
    # python operator // is floor division
    n_shells = len(component) // n_solute  # 4 shells
    deleted=None
    cutoff=0.0
    xx = np.arange(len(component))
    xlim1 = (0, len(component))
    # NOTE: Same as MPH's `plot_component()` function. See:
    # https://github.com/mpharrigan/wetmsm/blob/master/wetmsm/vmd_write.py#L229-L230
    tic_sq = component ** 2
    tic_sq /= np.max(tic_sq)
    loading1d = np.copy(tic_sq)
    loading1d[tic_sq < cutoff] = 0.0
    loading2d = unpack_component(component=component) #loading2d
    xlim2 = (0, loading2d.shape[0])
    utebins, shbins = np.meshgrid(np.arange(loading2d.shape[0]),
                                  np.arange(loading2d.shape[1]))
    file = open(solute_dir, "r")
    lines = reader(file)
    dataset = list(lines)
    File = np.array(dataset)
    for i in range(0,len(File)):
        print File[i,0]

    # Determine if they are alphaCarbons or betaCabons:
    file = open(atom_sel, "r")
    lines = reader(file)
    dataset = list(lines)
    File1 = np.array(dataset)
    for i in range(0,len(File1)):
        print File1[i,0]
    number = [i.split(" ")[0] for i in File1[:,0]]
    chain = [i.split(" ")[3] for i in File1[:,0]]
    carbons = [i.split(" ")[1] for i in File1[:,0]]
    residue = [i.split("\t")[2] for i in File[:,0]]
    res_label = [i.split(" ")[3] for i in File[:,0]]
    protein = [i.split("\t")[0] for i in File[:,0]]

    fig = plt.figure(figsize=(12,10))
    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2)
    ax1.scatter(xx, loading1d, linewidth=0, s=50, c='k',edgecolor='none')
    x = 0
    for k in range(0,len(loading2d[0,:])):
        for i in range(0,len(loading2d[:,k])):
            if loading1d[x] >= 0.40:
                if utebins[k,i] % 2 == 0: # if even number
                    if carbons[i] == "CA":
                        C = r"$\alpha$"
                    else:
                        C = r"$\beta$"
                    ax1.annotate(
                            '%s \n %s'%(res_label[i],C),
                            xy=(xx[x],loading1d[x]),
                            xytext=(xx[x]-6,loading1d[x]-0.15),
                            textcoords='data',
                            arrowprops=dict(arrowstyle="->",
                                ),
                            fontsize=11)
                else:
                    if carbons[i] == "CA":
                        C = r"$\alpha$"
                    else:
                        C = r"$\beta$"
                    ax1.annotate(
                            '%s \n %s'%(res_label[i],C),
                            xy=(xx[x],loading1d[x]),
                            xytext=(xx[x]+4,loading1d[x]-0.15),
                            textcoords='data',
                            arrowprops=dict(arrowstyle="->",
                                ),
                            fontsize=11)
            x += 1

    # Partitioning the Normalized eigenvector squared plot:
    x = 0
    for k in range(0,len(loading2d[0,:])):
        for i in range(0,len(loading2d[:,k])):
            if chain[i] == "A":
                ax1.axvspan(xx[x]-0.5,xx[x]+0.5, fill=True, linewidth=0, color='b', alpha=0.3)
            if chain[i] == "B":
                ax1.axvspan(xx[x]-0.5,xx[x]+0.5, fill=True, linewidth=0, color='g', alpha=0.3)
            x += 1

    ax1.set_xlabel('Solvent Features', fontsize=18)
    ax1.set_ylabel('Normalized $1^{St}$ Eigenvector$^{2}$', fontsize=18)
    ax1.set_ylim(bottom=cutoff, top=1.0)
    ax1.set_xlim(left=0, right=len(component))
    ax1.set_xticks(np.linspace(0,len(loading1d),len(loading1d)/25))
    ax1.get_xminorticklabels()
    ax1.xaxis.get_minorticklines()


    # Partitioning solutes by p53 and MDM2:
    x = 0
    for k in range(0,len(loading2d[0,:])):
        for i in range(0,len(loading2d[:,k])):
            if chain[i] == "A":
                ax2.axvspan(utebins[k,i]-0.5,utebins[k,i]+0.5, fill=True, linewidth=0, color='b', alpha=0.07)
            if chain[i] == "B":
                ax2.axvspan(utebins[k,i]-0.5,utebins[k,i]+0.5, fill=True, linewidth=0, color='g', alpha=0.07)
            x += 1
    # added +1 to shbins to make shells > 0
    plt.scatter(utebins, shbins+1, edgecolor='none',
                s=250*loading1d, # scaling might need to be tweaked
                c=loading2d.T, # can be positive of negative, ...
                # depends on the sign of your tICA component/sign of eigenvector
                cmap=plt.get_cmap('spring'))#'RdBu'))
    ax2.set_ylabel('Shells', fontsize=18)
    ax2.set_xlabel('Solute Atoms', fontsize=18)
    ax2.set_xlim(left=0, right=len(utebins[0,:])-1)
    ax2.set_yticks([1,2,3,4])
    ax2.get_xminorticklabels()
    ax2.xaxis.get_minorticklines()

    # Annotate the important features from the Logistic Eigenvector
    x = 0
    for k in range(0,len(loading2d[0,:])):
        for i in range(0,len(loading2d[:,k])):
            if loading1d[x] >= 0.20:
                if utebins[k,i] % 2 == 0: # if even number
                    if carbons[i] == "CA":
                        C = r"$\alpha$"
                    else:
                        C = r"$\beta$"
                    ax2.annotate(
                            '%s \n %s'%(res_label[i],C),
                            xy=(utebins[k,i],shbins[k,i]+1),
                            xytext=(utebins[k,i]-2,shbins[k,i]+1.25),
                            textcoords='data',
                            arrowprops=dict(arrowstyle="->",
                                ),
                            fontsize=11)
                else:
                    if carbons[i] == "CA":
                        C = r"$\alpha$"
                    else:
                        C = r"$\beta$"
                    ax2.annotate(
                            '%s \n %s'%(res_label[i],C),
                            xy=(utebins[k,i],shbins[k,i]+1),
                            xytext=(utebins[k,i]-2,shbins[k,i]+0.55),
                            textcoords='data',
                            arrowprops=dict(arrowstyle="->",
                                ),
                            fontsize=11)
            x += 1

    x_labelsize = ax1.get_xticklabels()[0].get_size()
    y_labelsize = ax1.get_yticklabels()[0].get_size()
    ax1.xaxis.set_tick_params(labelsize=16)
    ax1.yaxis.set_tick_params(labelsize=16)

    ticks = [ax2.xaxis.get_minor_ticks(),ax1.yaxis.get_minor_ticks(),
             ax2.xaxis.get_major_ticks(),ax1.yaxis.get_major_ticks()]
    marks = [ax1.get_xticklabels(),ax2.get_xticklabels(),
            ax1.get_yticklabels(),ax2.get_yticklabels()]
    for k in range(0,len(ticks)):
        for tick in ticks[k]:
            tick.label.set_fontsize(16)
    for k in range(0,len(marks)):
        for mark in marks[k]:
            mark.set_size(fontsize=16)
            mark.set_rotation(s=15)
    cbar = plt.colorbar()
    cbar.set_label('$1^{st}$ Eigenvector', fontsize=18)
    cbar.ax.tick_params(labelsize=16)

    fig.tight_layout()
    fig.savefig('Solvent_Shell_Features.pdf')


# }}}






