[![Snakemake](https://img.shields.io/badge/snakemake-≥6.0.2-brightgreen.svg)](https://snakemake.github.io)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Author: Sherine Awad 

VCF ToolBox  
=======================================================================================

This is a vcf toolbox for multiple things/commands to do with vcf. Like index, bgzip, merge, split, quick check for relatedness, Runs of Homozygozity (ROH), etc. 

Change the config.yaml file appropriately according to your data. 
Update samples.tsv to include your samples. You can edit config file to change this name.  


Then run: snakemake -jnumber_of_cores, for example for 5 cores use:

    snakemake -j5 

and for a dry run use: 

    snakemake -j1 -n 


and to print the commands in a dry run use:

    snakemake -j1 -n -p 

For the sake reproducibility, use conda to pull same versions of tools. Snakemake and conda have to be installed in your system:

    snakemake --cores --use-conda

### TODO 

Add everyday quick things we need to do with a VCF occasionally.


