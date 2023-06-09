<!-- #region -->
# A Tutorial and Basic Scripts for <p>FASTDock: A Pipeline for Allosteric Drug Discovery

## Furyal Ahmed and Charles L. Brooks III<p><p>Biophysics Program, University of Michigan, Ann Arbor, MI 48109<p>Department of Chemistry, University of Michigan, Ann Arbor, MI 48109

## Abstract

### Allostery is involved in innumerable biological processes and plays a fundamental role in human disease. Thus, exploration of allosteric modulation is crucial for research on biological mechanisms and in the development of novel therapeutics. Development of small molecule allosteric effectors can be used as tools to probe biological mechanisms of interest. One of the main limitations in targeting allosteric sites, however, is the difficulty in uncovering them for specific receptors. Furthermore, upon discovery of novel allosteric modulation, early lead generation is made more difficult as compared to orthosteric sites, because there is likely no information about the types of molecules that can bind at the site. In the work described here, we present a novel drug discovery pipeline, FASTDock, which allows one to uncover ligandable sites, as well as small-molecules that target the given site, without requiring pre-existing knowledge of ligands that can bind in the targeted site. By using a hierarchical screening strategy, this method has the potential to enable high throughput screens of an exceptionally large database of targeted ligand space.


## Scripts delineating the FASTDock workflow

0. *`script/parameters.py`*

      **File containing all modifiable parameters input into FASTDock. Modify
      this according to your system before running remaining scripts.**<p>

1. *`script/dock_probes.py`*

      **Generates a grid for a given protein structure to use for docking. Protein
      pdb must be in CHARMM format with all ligands, cofactors, and waters stripped. 
      Docks a set of 18 probes to the generated grid structure. 
      _Requires a CUDA enabled GPU._** 
      
      **to run:** `python script/dock_probes.py`<p>

2. *`script/cluster_and_score.py`*

      **Clusters each chemotype based on a 2A radius using MMTSB. Clusters different 
      chemotypes together to identify putative binding sites and scores them.**

      **to run:** `python script/cluster_and_score.py`<p>

3. *`script/MACCS_screen.py`*

      **Performs fingerprint screening of a compound library given a FASTDock 
      cluster for reference. Reference file must be in mol2 format. Database must be 
      given as sdf.** 

      **to run:** `python script/MACCS_screen.py`<p>

## Running the tutorial *`FastDock_Tutorial.ipynb`* 

### The tutorial recapitulates the same run as the basic set of scripts noted above. You can go through the steps in the tutorial and observe the function of each component as well as visualize the results.
           
## For information on installing pyCHARMM and CHARMM in a conda enviroment and installing MMTSB ToolSet see Install_Tools folder 
<!-- #endregion -->
