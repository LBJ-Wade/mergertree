# Monte Carlo Merger Tree Generation Code in Code
## Purpose:
    This code generates merger trees for a set of halos listed in a catalog file using the algorithm proposed in Parkinson, Cole, & Helly (2007 https://arxiv.org/abs/0708.1382). 
    
## Compile:
    To compile the code, you will need a C compiler and the GNU Scientific Library (GSL) installed in your computer. 
You can compile the code by typing

	Make

in the directory where the source files are located. 

## Files input files:
    To run the code, you need to prepare two text files
	1) a halo catalog file

	2) a snapshot list file

    There is a small example in the directory ./example to show you how they look like. 
In detail, the halo catalog file follows this format:
	1) The first line is an integer, which is the number of snapshots for the final tree saved. Note that the full tree generated when the code is running contains finer temporal resolution, but the final output only takes a set of discrete snapshots. The scale factors at which the snapshots are taken are listed in another file. 
	
    2) The second line provides the name of the snapshot scale factor list file. 
	
    3) The third line is an arbitary string. It is not used in the code, but you need to put in a string to hold the line for the code to read the file correctly. 
	
    4) The fourth line is a string for the name of the directory where the output merger trees are stored. 
	
    5) Number of halos will be generated (Ntrees). 
	
    6) From the sixth line, a list of halo mass starts. You are supposed to put the final mass of the halo in the first collumn, and put '0' in the second and the third column. The halo mass is in unit of Msun/h. If you are generating Ntrees with different masses, you need to have Ntrees lines down the list. If you are generating Ntrees with a same final mass, you just need to write one line with the mass you want. 

    For the snapshot list file, you simple write the scale factor from the earliest time to the latest time, one on each line. The number of lines in this file should match the number in the second entry of the halo catalog file. 

## Run the code:
    You are now ready to run the code, if you have all the files and directories correctly created. To run the code, just type:

	<dir/>mergertree <halocatalogfile> <mode> <random_number_seed>

    In the commend line, <halocatalogfile> is the file name of the halo catalog, <mode> is an integer either 0 or 1, and <random_number_seed> is another integer allowing you to choose a seed for the random number sequence in this run. 

    When <mode>==0, the code assumes that you are generating Ntrees with each tree's final mass provided in an entry in the halo catalog file. 
    
    When <model>==1, the code assumes that you are generating Ntrees with the same final mass that is provided in the first entry of the halo mass list. The code stops reading the list, and keeps generating Ntrees. 

##Other important tweaks you may need:

	1). You may what to change the mass resolution for your trees. You can change it in init.c. Find the variable "Mass_res", which is the mass resolution in units of Msun/h. 
	
    2). "Redshift" in init.c set the redshift where you the final descendent halo is. 
	
    3). "Redshift_max" is the maximum redshift for the code to generate tables for inteplation. It needs to be greater than the redshift of the earliest snapshot. 
	
    4). Mass_min and Mass_max set the lower and higher mass range for generating tables related to halo mass. Mass_mass needs to be smaller than the mass resolution, and Mass_max needs to be greater than the largest halo in the halo mass list. 
	
    5). If the code complaints about interpolation accuracy when it is running, you may need to increase Num_Bin_LnM, Num_Bin_Z, or Num_Bin_Ju. 
	
    6). The parameters Model.g0, Model.g1, and Model.g2 are set to be the values proposed in Parkinson et al. (2006) to best fit Millennium simulation merger trees. It is possible some variations on these values can fit simulations with other cosmology and halo identifiers better. 

Finally, the code was originally written for my own interests. It is not well documented in any way. And, I haven't been using it for quite a few years now. If you find any problem when using it, please do let me know. I am more than happy to help!
