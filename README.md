# README for the **Coupled-mode Shallow Water model (CSW)**. 
Copyright (c) 2016, Samuel M Kelly (smkelly@d.umn.edu)

This model simulates internal wave generation and propagation using vertical-mode rather than depth coordinates. Publications and presentations that use this code should cite:

**S. M. Kelly, P. F. J. Lermusiaux, T. Duda, and P. J. Haley, Jr. (2016) A Coupled-mode Shallow Water model for tidal analysis: Internal-tide reflection and refraction by the Gulf Stream, J. Phys. Oceanogr., 3661-3679.** 

**A. C. Savage, A. F. Waterhouse, and S. M. Kelly (2020) Internal tide nonstationarity and wave-mesoscale interactions in the Tasman Sea, J. Phys. Oceanogr., 2931-2951.** 

**S. M. Kelly, A. F. Waterhouse, and A. C. Savage (2021) Global Dynamics of the Stationary M2 Mode‐1 Internal Tide, Geophys. Res. Lett., 10.1029/2020GL091692.** 

**L. N. Thomas, S.M. Kelly, T. Klenz, W.R. Young, L. Rainville, H. Simmons, V. Hormann, and I. Stokes (2024) Why Near-Inertial Waves Are Less Affected by Vorticity in the Northeast Pacific Than in the North Atlantic, Oceanogr., 10.5670/oceanog.2024.301.** 

# 1 Creating the model inputs
## 1.1 Create a grid
Step 1: The global grid file uses Smith and Sandwell bathymetry, WOA stratification with ARGO mixed layer depths. It is created on the native bathymetry grid. using /matlab/grid_pre_processing/make_grid.m This script takes a long time (e.g., 3 days). 

The present inputs are:  
(1) Smith and Sandwell bathymetry (The function looks for 'topo_24.1.nc', but see http://topex.ucsd.edu/marine_topo/ for newer data)  
(2) WOA23 temperature and salinity (The function looks for 'woa23_decav91C0_t00_04.nc', see https://www.ncei.noaa.gov/thredds-ocean/catalog/woa23/DATA/temperature/netcdf/decav91C0/0.25/catalog.html)  
(3) ARGO mixed layer depths (The function looks for 'Argo_mixedlayers_monthlyclim_04142022.nc', see http://mixedlayer.ucsd.edu/)
(4) The matlab seawater toolbox 

The user can set:
>Nm=8; (number of vertical modes, can be less during simulations)  
>Nm0=128; (Number of structure functions to use when solving the eigenvalue problem)  
>dz=1; (vertical resolution in eigenvalue problem)  
>H_min=16; (Minimum depth to compute modes)  
>H_max=6000; (Maximum depth to compute modes)  

Step 2: The global grid file is smoothed and interpolated onto the CSW grid using /matlab/make_grid.m This script takes about 1/2 an hour, but requires around 100 GB of RAM. The CSW grid file must be created each time the model resolution is changed. 

The user can set:
>dx=1/50 (Degrees)  
>latlims=[-80 66]; (Domain limits)
>Ns=11; (Number of points in the 2D Guassian filter that smooths the global grid prior to interpolation onto the new resolution)

## 1.2 Create the tidal forcing
You have to create a new tidal forcing file when you (i) change the tidal data set, (ii) add tidal constituents, or (ii) change the tidal filters (smoothness, max values, etc.)

Open Matlab and run ./matlab/make_tides.m

The present inputs are:  
(1) the netCDF grid file created above  
(2) TPXO8 tides (http://volkov.oce.orst.edu/tides/tpxo8_atlas.html) 

The user can set:

>fid.grid='../../../17-6_global_grids/10th_deg_grid.nc'; (the grid file)  
>Nc=1;       (Number of tidal constituents)   
>Ns=3;       (Number of grid points to smooth tidal velocities)  
>H_min=16;	 (Set minimum depth for tides)  
>thresh=1;   (Maximum tidal velocity m/s) 


## 1.3 Create the wind forcing
You have to create a new wind forcing file when you (i) change the wind time period or (ii) change the wind filters (high pass, equatorial gap, etc.)

Open Matlab and run ./matlab/make_tides.m

The present inputs are:  
(1) the netCDF grid file created above  
(2) TPXO8 tides (http://volkov.oce.orst.edu/tides/tpxo8_atlas.html) 

The user can set:

>fid.grid='../../../17-6_global_grids/10th_deg_grid.nc'; (the grid file)  
>Nc=1;       (Number of tidal constituents)   
>Ns=3;       (Number of grid points to smooth tidal velocities)  
>H_min=16;	 (Set minimum depth for tides)  
>thresh=1;   (Maximum tidal velocity m/s) 


# 2 Compiling and running the CSW model:
## 2.1 Edit the compile-time parameters
All of the CSW's parameters are set in the header file: ./src/CSW.h 

Therefore, everytime something is altered, it CSW must be recompiled. This is accomplished by viewing/altering ./src/makefile and following the instructions there. The typical command (launched in the ./src/ directory) is simply

>make clean
>make 

The executable is ./src/cswexec

## 2.2 Running the CSW model
Create a "run" directory and enter it

Submit the job to a scheduler using the command

>qsub ../csw/run_csw.pbs

Make sure the pbs file processor request matches the processors assumed in CSW.h


# 3 Minnesota Supercomputing Institute 
Connect to MSI via:
ssh -Y user@mesabi.msi.umn.edu

Some useful stuff for a .bashrc:

>alias matlab='matlab -nosplash -nodesktop'  
>alias ijob='qsub -I -X -l walltime=4:00:00,nodes=1:ppn=4,mem=64gb'  
>alias isub='isub -n nodes=1:ppn=8 -m 64GB -w 2:00:00'  

The current module loading commands (in .bashrc) are:

>module load ncview
>module load netcdf/4.3.3.1-intel2015update3-serial
>module unload ompi
>module load netcdf/4.6.1-intel2018-parallel
>module load matlab

Then there's some useful pbs commands:

>showq -w user=  
>qstat -f  
>checkjob -v  
>acctinfo  
>groupquota  
>qdel  
>qsub  
>pdsh -j <id>


# 4 Editing the source code using GIT
The source code is contained in a git repository (moved from bitbucket to github 4/2025), which can be used in the following ways

## 4.1 Getting a fresh copy of the source code
This is known as "cloning":

>git clone https://username@github.com/username/csw


## 4.2 Saving changes
First, see if there are any new or changed files:  
>git status  and/or  git diff

Saving changes is a 3-step process. First you add the files with changes to the "index":  
>git add <filename>  
>git add * (git reset HEAD <filename> removes the file from the index)

Then you commit the changes in the "index" to the "head":  
>git commit -m "Commit message"

Lastly, you send the head to the remote repository:  
>git push origin master


## 4.3 Code branches
Start a branch: git checkout -b project_x

Check available branches:  
>git branch

Switch branches:  
>git checkout project_z

delete a branch:  
>git branch -d project_x

Add/commit to a branch: same as usual 

Send the branch head to the remote repository:   
>git push origin <branch>


## 4.4 Reset a file you messed up
>git checkout -- <file>


## 4.5 Update your local index to most recent commit
Say you fixed some bugs somewhere else and committed the update. Now you're working on another area of the code using an older working copy. In order to commit your changes to the newest version of the code, you have to update the index  
>git reset


## 4.6 Merging a branch
When the branch code is ready to be implemented in the master, you could merge, but it is better to rebase (provided nobody is working on the master).   
>git rebase -i master











