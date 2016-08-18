### Project information
##### Project title: Porting the Essential Dynamics/Molecular Dynamics method for large-scale nucleic acid simulations to ARCHER 
##### Student name: Zhuowei(James) Si
##### EPCC supervisor name(s): Elena Breitmoser (EPCC), Iain Bethune (EPCC)
##### External supervisor name(s) and organisation: Charlie Laughton (The University of Nottingham)

### Brief Introduction
Essential Dynamics/Molecular Dynamics (EDMD) is a method developed by the Laughton group in Nottingham for efficient simulation of biomolecules. Rather than computing forces using a classical force-field on an atomistic level as per conventional Molecular Dynamics, EDMD treats sections of the biomolecule such as DNA tetrads as a single entity (coarse-graining), computes forces on each body independently based on a simplified force-field, and then integrates Newton's equations of motion in a reduced space spanned by the Principal Components of the motion (i.e. capturing the Essential Dynamics rather than the full atomistic Molecular Dynamics).

Two implementations of the EDMD method exist: A serial Python code with poor performance and a parallel Fortran 95/MPI code which lacks of the flexibility. The aim of the project is to evaluate the two existing codes, then implement a new version which combines the best of both worlds: a flexible parallelisation scheme to achieve good load balance on ARCHER, and the improved integrator and user interface from the Python code.

### File Introduction
1. ./src : The folder contains all the source code inside;
2. ./test: The folder with two testing input file: the coordinate file and the tetrad parameter file;
3. ./data: The folder where the simulation results will be stored.
4. "config.txt" : The configuration file contains parameters to be initialised before the ED/MD simulation starts;
5. "Makefile"   : The Makefile to compile the code on ARCHER;
6. "README.md"  : The simple instruction of this code.
7. "edmddna.pbs": The script to submit the job to the back end of ARCHER

### Compile and Run
1. To compile the code on ARCHER, just type the command at the code directory: make.

2. To run the code on the back end of ARCHER, the code needs to be submitted: qsub edmddna.pbs

### Reference
1. [The QCP rotation calculation method](http://theobald.brandeis.edu/qcp/) in src/qcprot/. Developed by <br>  
 Douglas L. Theobald (2005), "Rapid calculation of RMSD using a quaternion-based characteristic polynomial.", Acta Crystallographica A 61(4):478-480. <br>  
 Pu Liu, Dmitris K. Agrafiotis, and Douglas L. Theobald (2009), "Fast determination of the optimal rotational matrix for macromolecular superpositions.", in press, Journal of Computational Chemistry 

2. The random number generator in `calculate_Random_Forces(Tetrad* tetrad)` in src/edmd.cpp <br>  
  Adapted from the following Fortran 77 code:
  ALGORITHM 712, COLLECTED ALGORITHMS FROM ACM. THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE, VOL. 18, NO. 4, DECEMBER, 1992, PP. 434-435. <br>  
  The function returns a normally distributed pseudo-random number with zero mean and unit variance.
The algorithm uses the ratio of uniforms method of A.J. Kinderman and J.F. Monahan augmented with quadratic bounding curves.
