# WaterMSM

**A Markov State Model of Solvent Features Reveals Water Dynamics in Protein-Protein Binding**

**Authors**
* Robert M. Raddi
    - Department of Chemistry, Temple University
* Vincent A. Voelz
    - Department of Chemistry, Temple University
    
Building a Markov State Model (MSM) using strictly solvent features to investigate the role water molecules have in the binding reaction of p53 TAD to the hydrophobic pocket of MDM2.


### Below is a description of what this repository contains:

#### NOTE: This code was written in 2016-2017 using Python 2

- [`scripts/analysis/`](scripts/analysis/): 
  - [`wetMSM.sh`](scripts/analysis/wetMSM.sh): **bash script for generating solvent features**
  - [`build_MSM.py`](scripts/analysis/build_MSM.py): **Python script for using solvent shell features (generated from `wetMSM.sh`) to project onto tICA space and build an MSM**
  - [`solvent_shells.py`](scripts/analysis/solvent_shells.py): **Python script for plotting the 1st and 2nd components and determining the slowest solvent motions**
   
- [`scripts/movies/`](scripts/movies/): **Python scripts for generating images for movie creation**
  - [`movie_making.py`](scripts/movies/movie_making.py): **Python script for concatenating subplots to form a movie**
  - [`tica_movie_subplot.py`](scripts/movies/tica_movie_subplot.py): **Python script to generate images of the tICA subspace along a trajectory**


- [`scripts/miscellaneous/`](scripts/miscellaneous/): **an assortment of scripts mainly used with Chimera such as molmap, calculating water counts and hbond counts**


#### To use most of the scripts inside [`scripts/miscellaneous/`](scripts/miscellaneous/) you must have [UCSF CHIMERA](https://www.cgl.ucsf.edu/chimera/) downloaded and sourced in your `~/.bashrc`

```bash
# First, you must specify the path of Chimera (for MacOS)
$ echo 'chimera=/Applications/Chimera.app/Contents/MacOS/chimera' > ~/.bashrc

# Call on chimera from the command line using the following:
$ chimera --nogui --nostatus --script "{script name}.py"
```











