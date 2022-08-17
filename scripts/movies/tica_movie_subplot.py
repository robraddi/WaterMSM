#!/usr/bin/env python

# Load Python Modules:{{{
import matplotlib
import numpy as np
from matplotlib import pyplot as plt
from sklearn.externals import joblib
import matplotlib
matplotlib.use('Agg')
fontfamily={'family':'sans-serif','sans-serif':['Arial']}
from matplotlib import pyplot as plt
plt.rc('font', **fontfamily)
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
# }}}

# Directories:{{{
sol_dir='solvent_data/'
tica_sol = joblib.load(sol_dir+'solv_tica_distances_sorted.pkl')
tICA_sol = np.concatenate(tica_sol)
# set your output directory for pngs
output = 'Movies/tICA_every_frame'
# }}}

# Code:{{{

# tICA sequence index for each trajectory labeled below:
t = [3309,511,1118,710,3057,1868]
# trajectory names:
name = ['RUN9_CLONE70','RUN13_CLONE53','RUN18_CLONE96',
        'RUN15_CLONE32','RUN7_CLONE41','RUN24_CLONE68']

for i in range(0,len(t)):
    print name[i],':'
    start_x,start_y = tica_sol[t[i]][0,0],tica_sol[t[i]][0,1]
    end_x,end_y = tica_sol[t[i]][len(tica_sol[t[i]][:,0])-1,0],tica_sol[t[i]][len(tica_sol[t[i]][:,0])-1,1]
    print '(%s,%s)'%(start_x,start_y),' --> ','(%s,%s)'%(end_x,end_y)
    x,y = [],[]
    for j in range(1,len(tica_sol[t[i]][:,0])+1):
        x.append(tica_sol[t[i]][j-1,0]),y.append(tica_sol[t[i]][j-1,1])
        fig = plt.figure(figsize=(12,10))
        ax = fig.add_subplot(111)
        ax.hexbin(tICA_sol[:,0], tICA_sol[:,1], bins='log',
                mincnt=1, cmap='viridis')
        ax.plot(x,y, 'k')
        ax.set_xlabel('$tIC_{1}$', fontsize=26)
        ax.set_ylabel('$tIC_{2}$', fontsize=26)
        ticks = [ax.xaxis.get_minor_ticks(),ax.yaxis.get_minor_ticks(),
                 ax.xaxis.get_major_ticks(),ax.yaxis.get_major_ticks()]
        for k in range(0,len(ticks)):
            if k == 0 or 1:
                for tick in ticks[k]:
                    tick.label.set_fontsize(22)
            else:
                for tick in ticks[k]:
                    tick.label.set_fontsize(22)
        minorLocator = MultipleLocator(100)
        # for the minor ticks, use no labels; default NullFormatter
        ax.xaxis.set_minor_locator(minorLocator)
        ax.tick_params(which='minor', width=1.75, length=3.5, labelsize=22)
        fig.tight_layout()
        fig.savefig('%s/%s/tICA_%s_frame_%d.png'%(output,name[i],str(name[i]),j))

# }}}
