# Multi-resolution approximation model, v1.0: Distributed Parallel Implementation

Note: As of February 23, 2021 the Shallow-Tree MRA is depreicated. Users should instead use the [Deep-Tree MRA](https://github.com/lewisblake/DeepTreeMRA).

March 11, 2019

Software authors:
Lewis Blake, Colorado School of Mines (lblake@mines.edu).
Dorit Hammerling, Colorado School of Mines (hammerling@mines.edu).

An associated technical report can be found at https://opensky.ucar.edu/islandora/object/technotes%3A579 .

This code is based on the MRA model described in
"A multi-resolution approximation for massive spatial datasets" by Matthias Katzfuss, 2017 in 
the Journal of American Statistical Association (DOI: 10.1080/01621459.2015.1123632).
Also at arXiv: https://arxiv.org/abs/1507.04789 .
References to this manuscript are found throughout the codebase.

This MRA implementation differs from both the Full-Layer Parallel MRA and the Serial MRA implementation in a few key ways.
It is designed which High-Performance Computing (HPC) applications in mind.
Here, the user specifies the first few levels for the creationg of the prior to be computed in serial. 
After serial computations, workers are assigned specific sections of the data with which to perform parallel computations.
Once parallel computations are complete, the remainder of the posterior is computed in serial.
In doing so, we reduce the memory build up on any given node at any given time. 

Designed and implemented with MATLAB 2018a (Version 9.4)
Previous versions may not be supported.
Required toolboxes:
- Statistics and Machine Learning Toolbox
- Optimization Toolbox
- Parallel Computing Toolbox

## GETTING STARTED:
This MATLAB codebase allows users to apply the multi-resolution approximation model to 2D spatial data sets.

The user_input.m script is where much of the user input can be modified (see USER INPUT and ADDITIONAL USER INPUT below).
If desired, users should modify model parameters within user_input.m.
The main.m script runs the model. 
Within the Matlab Editor Tab, selecting the 'Run' button from main.m will execute the code.

The repository is structured as follows: 
The LICENSE.txt, README.txt, and 'Distributed MRA' folder are contained in the 'MRA-Distributed' repository.
The user_input.m script, the main.m script, the find_num_levels_suggested.m function are contained within the 'Distributed MRA' folder.
Within the 'Distributed MRA' folder, there are four other folders: 'Data', 'Plots', 'Results', and 'subroutines'.
The 'Data' folder contains example data sets. 
The 'Results' folder is the default folder for results to be saved. Initially empty. 
The 'Plots' folder is the default folder for spatial prediction plots to be saved. Initially empty.
The 'subroutines' folder contains the functions used to execute the model which are called from main.m.

## EXAMPLE DATA: 
Two datasets are included in this distribution: satelliteData.mat and simulatedData.mat. 
These files are contained within the 'Data' folder.
Both of these data sets were originally used in Heaton, M.J., Datta, A., Finley, A.O. et al. JABES (2018). https://doi.org/10.1007/s13253-018-00348-w

## PARALLELIZATION:
A significant benefit of the MRA is that it lends itself to execution in parallel. 
For this reason, the portions of the creation of the prior, portions of the posterior inference, and spatial prediction in this codebase were designed to run in parallel using spmd. 
Within user_input.m, users can specify the number of workers in the parallel pool by setting NUM_WORKERS. 
SPMD blocks throughout the code will execute in parallel using NUM_WORKERS workers.
Within the MATLAB Cluster Profile Manager, the user can specify the desired cluster settings.
These settings include the number of nodes, nmber of cores, number of MPI processes, wall times, and many others.
For further reference please see https://www.mathworks.com/help/distcomp/discover-clusters-and-use-cluster-profiles.html .

## PRELIMINARIES:

### find_num_levels_suggested.m 

This is a stand-alone function and thus is not part of the subroutines.
This function estimates the NUM_LEVELS_M for a given dataset size (NUM_DATA_POINTS_n), number of knots (NUM_KNOTS_r), and number of partitions (NUM_PARTITIONS_J).
Note that NUM_LEVELS_M is a positive integer.



## USER INPUT:

### user_input.m 

In user_input.m, the areas requiring user input are as follows:

dataSource: | 'satellite' | 'simulated' |
    - These are the dataSource's for the data provided. Default is 'satellite'.
    - In order to use a different data set, see the section of load_data.m below and feed the case string for the data used to dataSource in user_input.m.

calculationType: | 'prediction' | 'optimize' | 'likelihood' | 'build_structure' |
Default is 'likelihood'.
calculationType can be set to any of the following calculation modes:
	- prediction: Uses given values for the parameters (theta and varEps) and just conducts spatial prediction. Parameters can be changed in load_data.m	
	- optimize: Optimizes over the range, variance and measurement error. The range and variance parameters are stored as a vector: theta. The measurement error is stored as a double: varEps.	
	- likelihood: Calculates the log-likelihood.
	- build_structure: Builds the multi-resolution structure. Reports summary statistics and produces a histogram of the number of observations assigned to regions at the finest resolution.

#### User Input relevant for any calculationType:

NUM_LEVELS_M: Total number of levels in the hierarchical domain-partitioning. By default set to 9.

NUM_PARTITIONS_J: Number of partitions for each region at each level. Only implemented for J = 2 or J = 4. By default set to 2.

NUM_KNOTS_r: Number of knots per partition. By default set to 64.

offsetPercentage: Offset percentage from partition boundaries. Must be between 0 and 1.
This quantity determines the buffer between the boundaries of a region where knots can be placed.
offsetPercentage is also used at the coarsest resolution to extend the maximal x and y domain boundaries as to include data points that may be exactly on the boundary within a region.
The domain boundaries define a rectangular region determined by the minimal and maximal x and y coordinate locations.
Preferably set offsetPercentage to be a small number (e.g. 0.01).

NUM_WORKERS: Number of workers in the parallel pool. Must be set to be a power of J <= nRegions(nLevelsInSerial). See PARALLELIZATION above. Default is 4.

nLevelsInSerial: Number of levels to compute the prior and posterior in serial. Default is *DEFAULT*.

verbose: Boolean variable indicating whether to produce progress indicators. Default is true.

resultsFilePath: Optional file path to save results for each calculationType. 
Set to be a string (e.g. resultsFilesPath = '/Users/JerryGarcia/Desktop/';). 
By default results are saved in the 'Results' folder.

#### User inputs relevant if calculationType = 'prediction' or 'build_structure'

displayPlots: Boolean variable indicating whether to display plots if predicting or building the multi-resolution structure.
savePlots: Boolean variable indicating whether to save plots if predicting or building the multi-resolution structure.
(Note: If not executing the 'prediction' or 'build_structure' calculationType, these booleans are not relevant.)

nXGrid: Number of prediction grid points in x-direction. By default set to 200.
nYGrid: Number of prediction gridpoints in y-direction. By default set to 200.
(Note: These parameters define a nXGrid x nYGrid prediction grid of spatial prediction locations if predicting.
The prediction grid is only defined within rectangular region given by the domain boundaries discussed above.)

plotsFilePath: Optional file path to save prediction plots if plotting.
Set to be a string (e.g. plotsFilesPath = '/Users/JerryGarcia/Pictures/';).
By default plots are saved in the 'Plots' folder.

#### User inputs relevant if calculationType = 'optimize'

lowerBound: Vector of lower-bound values required for the optimization search. Default is [0,0,0].

upperBound: Vector of upper-bound values required for the optimization search. Default is [10,1,5].

initialEstimate: Vector of inital estimates of parameteres required for the optimization search. Default is [5, 0.3, 0.1].


## ADDITIONAL USER INPUT

### load_data.m

In load_data.m the user can specify the type of data being used and the file paths. 
The file paths are presently relative for using the data provided. 
Data in other locations can be loaded using absolute files paths. 
In order to use a different data set, a new case within the switch clause must be added with the case given as a string, a file path to the data set with the load() function, and appropriate values for theta and varEps. 
If these values are not known, they can be given lower and upper bounds and can then be estimated using the 'optimize' mode. 
An example of what a case for a new data set may be is as follows.

e.g., Within the switch clause, specify:
case 'myData'
load('/Users/JerryGarcia/Documents/Data/myData.mat')
theta = [2, 1]; varEps = 0.01;

Data being used must have three columns 'lat', 'lon', and 'obs' denoting latitude, longitude, and the observations or be coerced from their native format into variables with those names.

The user can also change the values of theta and varEps in load_data.m.
Values can determined by the 'optimize' mode. For the 'satellite' and 'simulated' data provided, those values as determined by the 'optimize' mode are set as the default values.

### evaluate_covariance.m 

evaluate_covariance is a general covariance function. By default, it is set as an exponential and can be changed here.


## OUTPUT:

Model output is dependent on the calculationType (computational mode) performed. 

1) For the 'prediction' mode, the output is a .mat file with the MRA results stored within the Results folder. This .mat file contains the prediction locations, prediction mean, and the prediction variance.
If either boolean variables 'displayPlots' or 'savePlots' are set to true, three plots are also produced corresponding to the observations, predicted values, and the prediction variance with the create_plots() function. 
Saving these plots can be accomplished by setting savePlots to true in user_input.m. 

2) For the 'optimize' mode, optimized values for theta and varEps are stored in a .mat file stored witin the Results folder. 

3) The 'likelihood' mode returns the log-likelihood stored in a .mat file within the Results folder.
If verbose is set to true, the log-likelihood will print to the Command Window as well.

4) For the 'build_structure' mode, summary statistics of the distribution of observations to regions at the finest resolution are reported within '/Results/structureSummaryStats.txt'. A histogram is also produced within the Plots folder.

## NOTE: 
If computing on a remote server and file pathing is an issue, comment out the call addpath('subroutines') in user_input.m and copy files in subroutines into the same folder as main.m. 
This may also be necessary for data sets as well.

