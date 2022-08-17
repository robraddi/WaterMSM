#!/bin/bash

#########################################################
# **Additional resources:**
# MPH paper: https://mpharrigan.com/research/wetmsm.pdf
# https://github.com/mpharrigan/wetmsm
#########################################################

# Parameters and Directories:{{{
top='sys_1YCR.pdb'
solute='AtomIndices.dat'
trrs='TRAJECTORIES/TRRS/Traj_RUN*_CLONE*.trr'
solvent='solvent.txt'
assign='TRANSFORMED/SolventAssignments/'
n_components=10
lag_time=5
n_shells=4
shell_width=0.3 # Å

#}}}

# 1) Indices:{{{
## Obatining indices for all solvent atoms:
msmb AtomIndices --atoms Atoms\
  --pdb $top\
  --out solvent.txt\
  --water\
#}}}

# 2) SolventShellsFeaturizer:{{{
## Featurizing Solvation shells around solute atoms:
msmb SolventShellsFeaturizer --trjs $trrs\
 --top $top\
 --solute_indices $solute\
 --solvent_indices $solvent\
 --n_shells 4 --shell_width 0.3\
 --transformed TRANSFORMED\
 --out waterfeats
#}}}

mv waterfeats.pkl ./TRANSFORMED/
cd TRANSFORMED

# 3) tICA on SolventFeaurues:{{{
## Perform tICA on solvent feeatures:
msmb tICA -i ./\
  -t ./\
  --n_components $n_components\
  --lag_time $lag_time\
  -o wetmsm_tICA_model
#}}}

# 4) SolventShellsAssigner:{{{
## Implement SolventShellsAssigner to assign solvent atoms to shells.
msmb SolventShellsAssigner --trjs "$trrs"\
 --top $top\
 --solute_indices $solute\
 --solvent_indices $solvent\
 --n_shells $n_shells\
 --shell_width $shell_width\
 --transformed .\
 --out SolventAssignments
#}}}

# 5) SolventApplyComponents:{{{
## Apply tICA to SolventAssignments
msmb SolventApplyComponents --trjs "$trrs"\
 --top $top\
 --solvent_indices $solvent\
 --solute_indices $solute\
 --assignments "$assign"\
 --component 'wetmsm_tICA_model.pkl'\
 --out watercomponent
#}}}

# 6) SolventWriteVMD:{{{
## Writes a file that VMD can parse
msmb SolventWriteVMD --trj "$trrs"\
 --top $top\
 --dataset './watercomponent'\
 -o watervmd

#}}}


