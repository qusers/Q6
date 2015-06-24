#Debugging hints for EVB simulations using Q.
**Created: Nov 21, 2014**    
**Author: Beat Amrein**    

    [with a long term goal to offer checklist with energy-offsets and common solutions]

Work carefully trough this checklist to make sure that not you, but your system is the problem.

Energy wrong by >50KCAL  
Are you mapping to the right water run?  
Your profile is ever increasing and has **no TS** [dG* > 50kcal]:  
(you are not running the same reaction in water and protein).
     
    
###Debugging:  

0. Are you mapping to the right water run?
1. Check Movie (compare with your chemical intuition).
2. Check Fep Files. 
3. Compare FEP/Movie from Water- with Protein-run
4. Locate your problem, which is very likely in your FEP files.   
   [Morse, Dihedral, SoftRepulsion, Q- vs PDB numbering, etc.]  

Fix:  
Repair/Adjust your FEP File.

Energy wrong by >5KCAL  
Your **standard deviation is large** [ >5 kcal ]:  
(your replicates are not the same, maybe a flipping water)  

1. Compare the movies of the trajectories with the highest deviation.  
  You likely observe a flipping HOH, GLU, HIS etc.  
  Fix:  
  Use a distance or angle restraint to keep the bad guy under control.  
  Or use a sequence restraint and enable apply it also on the hydrogen of the water...  
  
2. You did not equilibrate your system long enough and/or your frames are too short.  
  Fix:  
  Try to run longer.  
  
3. If you do work on a proton tranfer, you might want to apply a sequence restraint  
  on your transferred protons.  
  Fix:  
  Use a sequence restraint.  
  
**Water-Run specific Problems:**  
Check the movies:  

    Q1) Your water runs are very inconsistent and you have too many degrees of freedom?  
      A1 Apply a sequence restraint of 1.0 on your reaction.  
    Q2) A solvent water molecule interferes with your reaction?  
      A2 Apply a higher sequence restraint, or use additional distance restraints, or run in gas-phase.  

**Protein-Run specific Problems:**  
Check the movies:  

    Q1) You can not obtain the literature energies in protein, but your water run is cool?  
      A: Check the ionization of close-by residues.  
      A: Check hydrogen bonds and group contributions.  
      A: Look at temperature, distances etc.  
      A: This is a very common problem and requires you to play around.   
  
**Shake Failure:**  
The previous run has hot atoms. If you are equilibrating the system:  
try to run the previous step longer (eg double the time).   

**NaN:**  
Your system exploded.   
- Look at the movie (and/or the logfile) to see what happened:   
  If you can't see a thing (immediate explosion), decrease the stepsize (e.g. to 0.01) and the trajectory interval (e.g. to 1).   
- Also the pdb file you used to generate the topology with vmd and try to locate strange bonds that formed, when you e.g., docked some substrate.  
- If this does not help and you use a fep file, you have an error in your fep file.  
