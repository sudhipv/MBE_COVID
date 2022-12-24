# README #

# SEIRD Compartmental Model for COVID-19 Evolution

This repository contains codes used for generating the results in the paper titled "Scalable computational algorithms for geo-spatial COVID-19 spread in high performance computing".
Note that the code has only been tested on DIGITAL RESEARCH ALLIANCE OF CANADA (DRAC) machines  (https://docs.alliancecan.ca/wiki/Technical_documentation).

![south_ON](https://user-images.githubusercontent.com/121322281/209414633-39038460-9c30-49f7-b6ef-e7183f71a912.jpeg)
### The Southern Ontario domain divided in to 200 subdomains.

![inf_50](https://user-images.githubusercontent.com/121322281/209414526-61b6cdd4-d1a6-4135-be34-31f4157b5890.png)
### Infected density of Southern Ontario at 50 days from Sep 1st 2020.

# Packages Needed

1. Gmsh - For pre-processing :  https://gmsh.info/#Documentation 
	
	Gmsh is used for generating the corresponding mesh files necessary for simulation and scalability studies. Currently a mesh file is added to "ontario/mesh/south/southON.msh" which is for the domain of Southern Onatrio having 
	209429 vertices for coarse grid and 834586 vertices in fine grid. Similarly one mesh for square domain can be found at "MMS/2d/mesh/square.msh" with 606 vertices. If the user wishes to try out different desnity of mesh sizes, 
	the above package can be installed in local machine for preprocessing and then the generated .msh files can be transferred to cluster. 		

2. FreeFEM - For model and domain decomposition-based solver : https://doc.freefem.org/introduction/index.html

	This package is necessray to be installed in cluster for compilation and simulation of the model.
	
3. Paraview - For post-processing :  https://www.paraview.org/about/

	This is a post-processing package which can be used for viewing results generated from the FreeFEM solver. User can install this package in local machine and copy files from cluster (such as .vtu and .pvd files).

# This code contains folders as listed below

1. doc - Documentation for the installation of FreeFEM on the DRAC machines.
2. GIS - Contains data on the geography of Southern Ontario concerning the relevant geometry for finite element mesh generation and data on the population of different Public Health Units.
3. MMS - Contains the code necessary to reproduce the model verification results from Appendix B using the method of manufactured solutions.
4. ontario - Contains the code necessary to reproduce the simulation of COVID-19 propagation in Southern Ontario from Section 5.2.
5. PC_covid - Contains the code necessary for the comparison of the preconditioners from Section 5.1.
6. ODE - ODE model for comparing with PDE model in square domain.


# Steps for compilation and Running the code

1. Install Gmsh in your cluster. : Follow instructions from the page : https://gitlab.onelab.info/gmsh/gmsh/
2. Install FreeFEM in your respective cluster. Instructions for installing it in DRAC machines are provided inside "doc/freefem_install_intel2020.txt".
3. Install paraview in your local machine for viewing the generated results.
4. Compile and run the code in cluster following instructions in "doc/instructions_parallel.txt". 
5. Copy results from cluster to local machine and view in paraview.


# Equations from the manuscript:

### 1. Section 5.1

Code : PC_covid/covid_scalar_twolevel.edp 

The comparison of various preconditioners can be run by modifying the code above. 
Different preconditioners can be called by altering few lines of code as specified inside the code : Lines : 216 - 297.

The parameters from Table 1 can be seen on lines 84-97.
The initial condition in Eq. (5.1) can be seen on line 134
The weak form in Eq. (A.1) - Eq. (A.5) in the paper can be seen on lines 376-408.


### 2. Section 5.2

Code : ontario/covid_scalar2L_southON.edp

Eq. (5.3) can be seen on lines 95-96.
The parameters from Table 1 can be seen on lines: 90-109.
Initial conditions in Eq. (5.2) can be seen on lines 138-378. Note that, generating results for Eastern, Central and Western parts as in Fig. (10) can be done by providing the 
initial conditions to only the respective parts of the domain and setting other parts with zero.

The weak form in Eq. (A.1) - Eq. (A.5) can be seen on lines 531-560.


### 3. Section B.2 

Code : MMS/1d/seird_MMS.py


Note that the above code is written in FEniCS which can be installed follwoing instructions from https://fenicsproject.org/download/.

The parameters from Table 3 can be seen on lines 73-83.
Eq. (B.8) - Eq. (B.12) can be seen on lines 107-115.
The weak form in Eq. (A.1) - Eq. (A.5) can be seen on lines 200-219.

### 4. Section B.3.

Code : MMS/2d/2d_Neumann.edp

The parameters from Table 1 are found in lines 36-49
Eq. (B.13) - Eq. (B.17) can be seen on lines 67-71
Eq. (B.18) - Eq. (B-21) can be seen on lines 74-80
The weak form in Eq. (A.1) - Eq. (A.5) can be seen on lines 163-212.

### 5. Section B.4 

Code : ODE/seird_ODE.py


The PDE results for comparison of ODE and PDE can be generated by the same code as in "ontario/covid_scalar2L_southON.edp". 

Initial condition can be seen in line 28-33.
Eq. (B.22) - Eq. (B.26) can be seen in lines 57-61.







