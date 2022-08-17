#!/usr/bin/env python

# Libraries:{{{
import glob, os, sys, subprocess, commands
import matplotlib
import mdtraj as md
import msmbuilder.utils
import numpy as np
from msmbuilder import *
from itertools import combinations
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from msmbuilder.cluster import KCenters, KMedoids
from msmbuilder.decomposition import tICA
from msmbuilder.featurizer import AtomPairsFeaturizer
from msmbuilder.msm import ContinuousTimeMSM, implied_timescales, MarkovStateModel
from sklearn.externals import joblib
from sklearn.pipeline import Pipeline
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
from matplotlib.colors import LogNorm
matplotlib.use('Agg')
fontfamily={'family':'sans-serif','sans-serif':['Arial']}
from matplotlib import pyplot as plt
plt.rc('font', **fontfamily)
# }}}


def tICA_analysis(data, lagtime, nComponents):
    print('performing tICA analysis...')
    tica_model = tICA(
        lag_time=lagtime,
        n_components=nComponents)
    model = tica_model.fit_transform(data)
    save_object(model, 'tICA_model.pkl')
    eigenvec = tica_model.eigenvectors_
    np.save('eigenvectors.npy',eigenvec)
    np.save('tica_components_.npy',tica_model.components_)
    tica = np.concatenate(model)
    np.save('tICA.npy',tica)
    return tica


def plot_tICA(data, tICs=[0,1], output='tIC_1_and_2'):
    print('Plotting...')
    fig = plt.figure(figsize=(12,10))
    ax = fig.add_subplot(111)
    # Sometimes the sign of the tICA components need to be flipped
    ax.hist2d(data[:,tICs[0]], -data[:,tICs[1]], bins=300, norm=LogNorm())
    ax.set_xlabel(r'$tIC_{1}$', fontsize=20)
    ax.set_ylabel(r'$tIC_{2}$', fontsize=20)
    ticks = [ax.xaxis.get_minor_ticks(),ax.yaxis.get_minor_ticks(),
             ax.xaxis.get_major_ticks(),ax.yaxis.get_major_ticks()]
    for k in range(0,len(ticks)):
        if k == 0 or 1:
            for tick in ticks[k]:
                tick.label.set_fontsize(16)
        else:
            for tick in ticks[k]:
                tick.label.set_fontsize(18)
    plt.tight_layout()
    plt.savefig('%s'%output)



if __name__ == '__main__':


    frame0 ='template_files/system0.pdb'
    Atom = 'template_files/AtomIndices.dat'
    solvent = 'TRAJECTORIES/centered/solvent/solvent.txt'
    trrs = 'TRAJECTORIES/centered/TRRS/Traj_RUN*_CLONE*.trr'
    waterfeats = sorted(glob.glob("transformed/*.npy"))
    n_clusters=600
    cluster_method = 'kcenters'
    tica_lagtime = 50   # lag time used in Zhou et al. which equals 5ns
    tica_components=10
    indices = np.loadtxt(solvent).astype(int)
    AtomIndices = np.loadtxt(Atom).astype(int)
    pdb = md.load(frame0)   # load protien data bank template
    trrfiles = sorted(glob.glob(trrs))

    traj =[]
    for i in range(0,len(waterfeats)):
        features = np.load(waterfeats[i])
        print 'Reading trajectory', i, 'of', (len(waterfeats)-1) , '...'
        xtc = waterfeats[i]
        if not (os.path.exists(xtc)):
            print '    File not found!'
            continue
        traj.append(features)

    tica = tICA_analysis(data, lagtime, nComponents):
    plot_tICA(data, tICs=[0,1], output='tIC_1_and_2'):


    # Cluster and get ready to build MSM
    tICA_distances = np.load("/home/tuc41004/Projects/water-binding-rxn/PROJ6381/TRAJECTORIES/centered/solvent/solv_tICA_distances_sorted.npy")
    tica_distances = joblib.load("/home/tuc41004/Projects/water-binding-rxn/PROJ6381/TRAJECTORIES/centered/solvent/solv_tica_distances_sorted.pkl")
    Tica = joblib.load('/home/tuc41004/Projects/water-binding-rxn/PROJ6381/TRAJECTORIES/centered/solvent/solv_tica_model.pkl')

    if cluster_method == 'kcenters':
        print("Clustering via KCenters...")
        clusters = KCenters(n_clusters)
    elif cluster_method == 'kmeans':
        print("Clustering via KMeans...")
        clusters = KMeans(n_clusters)
    else:
        sys.exit("Invalid cluster_method. Use kmeans or kcenters.")
    sequences = clusters.fit_transform(tica_distances)
    np.save('%s_lag_%d_clusters_%d_sequences.npy'%('solvent',tica_lagtime, n_clusters), sequences)
    np.save('%s_lag_%d_clusters_%d_center.npy'%('solvent',tica_lagtime, n_clusters),clusters.cluster_centers_)
    joblib.dump(sequences,'%s_lag_%d_clusters_%d.pkl'%('solvent',tica_lagtime, n_clusters))

    # Determining the cluster populations
    print("Determining cluster populations...")
    counts = np.array([len(np.where(np.concatenate(sequences)==i)[0]) for i in range(n_clusters)]) # how many frames are in each cluster
    normalized_counts =  counts/float(counts.sum())
    print 'Normalized Counts',normalized_counts
#    np.savetxt('%spopulations.dat'%('fullMSM_output/'), normalized_counts)
    print 'saving populations.dat'
    percentages = [ i*100 for i in normalized_counts ]

    # Plotting the tICA components
    print("Plotting tICA components with cluster centers...")
    fig = plt.figure(figsize=(12,10))
    ax = fig.add_subplot(111)
    ax.hexbin(-tICA_distances[:,0], tICA_distances[:,1], bins='log',
            mincnt=1, cmap='viridis')
    ax.set_xlabel('$tIC_{1}$', fontsize=20)
    ax.set_ylabel('$tIC_{2}$', fontsize=20)
    ticks = [ax.xaxis.get_minor_ticks(),ax.yaxis.get_minor_ticks(),
             ax.xaxis.get_major_ticks(),ax.yaxis.get_major_ticks()]
    for k in range(0,len(ticks)):
        if k == 0 or 1:
            for tick in ticks[k]:
                tick.label.set_fontsize(16)
        else:
            for tick in ticks[k]:
                tick.label.set_fontsize(18)
    fig.tight_layout()
    fig.savefig('solvent_lag_%d_clusters_%d.png' %(tica_lagtime, n_clusters))

    plt.tight_layout()
    plt.savefig('solvent_tIC_1_tIC_2.png')

    # tic1 and tic2
    # Plot the hexbin heat map
    txx=tICA_distances
    clustered_trajs = sequences


    # Plot tICA with clustered trajectories to build MSM
    fig = plt.figure(figsize=(12,10))
    ax = fig.add_subplot(111)
    from matplotlib import pyplot as plt
    ax.hexbin(-txx[:,0], txx[:,1], bins='log', mincnt=1, cmap='viridis')
    ax.plot(-tica_distances[3309][:,0],tica_distances[3309][:,1],'k')
    from msmbuilder.msm import MarkovStateModel
    msm = MarkovStateModel(lag_time=tica_lagtime, n_timescales=10)
    msm.fit(clustered_trajs)
    assignments = clusters.partial_transform(txx)
    assignments = msm.partial_transform(assignments)
    x_centers = clusters.cluster_centers_[msm.state_labels_, 0]
    y_centers = clusters.cluster_centers_[msm.state_labels_, 1]
    data = ax.scatter(-x_centers,
                y_centers,
                c=msm.left_eigenvectors_[:, 0], # color by eigenvector
                cmap="coolwarm",
                zorder=3)
    fig.gca().set_xlabel('$tIC_{1}$', size=20)
    fig.gca().set_ylabel('$tIC_{2}$', size=20)
    for stateid in range(n_clusters):
    	#ax.plot(-clusters.cluster_centers_[stateid,0],clusters.cluster_centers_[stateid,1],'ro')
    	ax.text(-clusters.cluster_centers_[stateid,0],clusters.cluster_centers_[stateid,1],'%d'%stateid)
    ticks = [ax.xaxis.get_minor_ticks(),ax.yaxis.get_minor_ticks(),
             ax.xaxis.get_major_ticks(),ax.yaxis.get_major_ticks()]
    for k in range(0,len(ticks)):
        if k == 0 or 1:
            for tick in ticks[k]:
                tick.label.set_fontsize(16)
        else:
            for tick in ticks[k]:
                tick.label.set_fontsize(18)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    cbar = fig.colorbar(data, cax=cax, orientation='vertical' )
    cbar.ax.set_ylabel('$Eigenvector_{1}$',fontsize=20)
    fig.tight_layout()
    fig.savefig('solvent_overlay_RUN9CLONE70_cluster_centers_1_2.png')




    # 6. Calculating Implied Timescales and plotting with Lagtimes
    time_step = 1.0
    n_timescales = 10
    lagtimes = np.array([1,5,20,50,100,150])
    print("\nCalculating Implied Timescales...")
    timescales = implied_timescales(sequences, lagtimes, n_timescales=n_timescales,
                                    msm=MarkovStateModel(verbose=False))

    implied_timescale_data = 'ipt_lag_%d_clusters_%d.pkl' %(tica_lagtime, n_clusters)
    joblib.dump(timescales, implied_timescale_data)
    numpy_timescale_data = 'lag_%d_clusters_%d_timescales.npy' %(tica_lagtime, n_clusters)
    np.savetxt('lagtimes.txt', lagtimes)
    np.save(numpy_timescale_data, timescales)
    print("Plotting Implied Timescales...")
    fig = plt.figure(figsize=(12,10))
    ax = fig.add_subplot(111)
    for i in range(n_timescales):
        ax.plot(lagtimes * time_step, timescales[:, i] * time_step, 'o-')
    ax.set_title('MSM Relaxation Timescales')
    ax.set_yscale('log')
    ax.set_xlabel('Lagtime (ns)', fontsize=20)
    ax.set_ylabel('Implied Timescales (ns)', fontsize=20)
    ticks = [ax.xaxis.get_minor_ticks(),ax.yaxis.get_minor_ticks(),
             ax.xaxis.get_major_ticks(),ax.yaxis.get_major_ticks()]
    for k in range(0,len(ticks)):
        if k == 0 or 1:
            for tick in ticks[k]:
                tick.label.set_fontsize(16)
        else:
            for tick in ticks[k]:
                tick.label.set_fontsize(18)
    fig.tight_layout()
    fig.savefig('solvent_lag_%d_clusters_%d.png' %(tica_lagtime, n_clusters))

