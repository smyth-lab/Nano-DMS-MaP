# Nano-DMS-MaP

Code to perform Nanopore-DMS-MaP (mutational profiling) RNA structural probing analysis.

Note: The code accompanying the Nature Methods Nano-DMS-MaP paper is located at this branch: https://github.com/smyth-lab/Nano-DMS-MaP/tree/NMethods-10.1038/s41592-023-01862-7 

In this repository we provide both the code to process samples, as well as for visualization.

Briefly, the pipeline has the following steps:

    Basecalling (and demultiplexing) of fast5 files from native barcoded SQK-NBD112 or SQK-NBD114 (Q20+ kits) runs with guppy.
    Read-to-isoform assignment with Isoquant, followed by sorting.
    Per isoform alignment with LAST.
    RNA-Framework analysis 
      rf-count with optimized settings on LAST-aligned BAM files. 
      rf-norm to quantify reactivity rates (we recommend having a non-modified control to perform Siegfried et al. normalization)
      rf-correlate and rf-combine to evaluate the reproducability of reactivity (>0.9 expected) and calculate the mean reactivity of repplicates.
    
    Option A: Evaluating the agreement with published secondary structures by calculating Receiver-Operator Curves
    Option B: Performing de novo RNA structure prediction using EternaFold (or rf-fold)

Additional workflows in this repository include:

    Generating a custom gtf file for isoquant to take into account primer binding sites
    Mutation profile analysis of BAM files using perbase
    Analyzing the frequency of mutations per read from RNA-Framework mm files
    Generating .varna files of known or predicted structures colored by reactivities


# Installation instructions

Prerequisites: 
- Linux system
- Ana-/Miniconda: https://docs.conda.io/projects/miniconda/en/latest/ 
- Java (for Varna)

  
1. Clone the github repository into the current location. This will automatically clone additional tools available on GitHub, i.e. [RNAFramework](https://rnaframework-docs.readthedocs.io/en/latest/) and [EternaFold](https://eternafold.eternagame.org/) (into the tools folder)
```
git clone https://github.com/smyth-lab/Nano-DMS-MaP
```
2. Move into the Nano-DMS-MaP folder
```
cd Nano-DMS-MaP
```
3. Create a conda environment that contains the majority of tools required for analysis. This may take a while ...
```
conda env create -f environment.yml
```
4. Activate the conda environment
```
conda activate Nano-DMS-MaP
```
5. Compile the EternaFold binary
```
cd tools/EternaFold/src
make multi
cd ../../../
```
7. Download the most current version of [Varna](https://varna.lisn.upsaclay.fr/) (needs to be located in the working directory to function).
```
wget https://varna.lisn.upsaclay.fr/bin/VARNAv3-93.jar
```
8. Start the jupyter lab server. If running on a cluster connected via ssh, add the no-browser and port options. 
```
jupyter lab (--no-browser --port=XXXX)
```
If launching jupyter lab on a cluster, run the following command on your local computer to create an SSH tunnel. 
```
(ssh -N -L XXXX:localhost:XXXX user@server)
```
10. Start analysis in jupypter lab (either opened automatically when launched on local machine, or open by copying the URL shown after the `jupyter lab` command. 
