# Project information
## Project title: Porting the Essential Dynamics/Molecular Dynamics method for large-scale nucleic acid simulations to ARCHER 
Student name: Zhuowei(James) Si
EPCC supervisor name(s): Elena Breitmoser, Iain Bethune
External supervisor name(s) and organisation: Charlie Laughton (The University of Nottingham)

Superviser: Elena Breitmoser (EPCC), Iain Bethune (EPCC), Charlie Laughton (The University of Nottingham)

Essential Dynamics/Molecular Dynamics (EDMD) is a method developed by the Laughton group in Nottingham for efficient simulation of biomolecules.  Rather than computing forces using a classical force-field on an atomistic level as per conventional Molecular Dynamics, EDMD treats sections of the biomolecule such as DNA tetrads as a single entity (coarse-graining), computes forces on each body independently based on a simplified force-field, and then integrates Newton's equations of motion in a reduced space spanned by the Principal Components of the motion (i.e. capturing the Essential Dynamics rather than the full atomistic Molecular Dynamics).