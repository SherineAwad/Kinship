[![Snakemake](https://img.shields.io/badge/snakemake-≥6.0.2-brightgreen.svg)](https://snakemake.github.io)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

Author: Sherine Awad 

Kinship and Relatedness 
=======================================================================================

Various ways to check kinship, relatedness, IBS, and IBS. 

Change the config.yaml file appropriately according to your data. 

Update SAMPLES in config file to include your samples. The pipeline uses vctfools, King, and AKT. More tools will be added.   

A sample name is included in the config file without the suffix 'r_1.fq.gz' or 'r_2.fq.gz' if paired end and without 'fq.gz' if single-end reads. 


### How to run the pipeline 
The pipeline requires snakemake and conda installed. 

Then you use: 

    snakemake -jnumber_of_core

For example for 5 cores use:

    snakemake -j5 

For a dry run use: 

    snakemake -j1 -n 


To print the commands in a dry run use:

    snakemake -j1 -n -p 

For the sake reproducibility, use conda to pull same versions of tools. Snakemake and conda have to be installed in your system:

    snakemake --cores --use-conda

### TODO 

Add more tools as we go 


