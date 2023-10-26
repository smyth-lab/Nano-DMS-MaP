# Nano-DMS-MaP

Code to analyze and visualize RNA structural probing by Nanopore DMS mutational profiling data. 

Note: The code accompanying the Nature Methods Nano-DMS-MaP paper is located at this branch: https://github.com/smyth-lab/Nano-DMS-MaP/tree/NMethods-10.1038/s41592-023-01862-7 

## Core Workflow 

_(numbers in square brackets indicate in which jupyter notebook the process is documented)_

1. [1] Basecalling (and demultiplexing) of fast5 files from native ligation barcoded SQK-NBD112 or SQK-NBD114 (Q20+ kits) runs of specific amplicons.
2. [1] Read-to-isoform assignment with [Isoquant](https://github.com/ablab/IsoQuant), followed by sorting.
3. [1] Per isoform alignment with [LAST](https://gitlab.com/mcfrith/last).
4. [1] [RNA Framework](https://rnaframework-docs.readthedocs.io/en/latest/) analysis 
    - rf-count with optimized settings on LAST-aligned BAM files. 
    - rf-norm to quantify reactivity rates (we recommend having a non-modified control to evaluate DMS-induced mutation rates and for normalization)
    - rf-correlate and rf-combine to evaluate the reproducability of reactivity (>0.9 expected) and calculate the mean reactivity of replicates.
5. [4] Plotting DMS reactivities per sample and isoforms and generating varna files of known or predicted structures colored by DMS reactivities. 
6. [4] Performing de novo RNA structure prediction using [EternaFold](https://github.com/eternagame/EternaFold) (or rf-fold)
    
## Optional (QC) steps
- [0] Generating a custom gtf file for Isoquant to take into account primer binding sites
- [2] Visualization of Isoquant's isoform assignment
- [3] Mutation profile analysis of BAM alignment files using [perbase](https://github.com/sstadick/perbase)
- [5] Subsampling with [sambamba](https://github.com/biod/sambamba) to evaluate whether resequencing is worthwhile

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
7. Download the most current version of [VARNA](https://varna.lisn.upsaclay.fr/) (needs to be located in the working directory to function).
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


# References

Tool|Citation|License
-----|-------|-------
LAST|Michiaki Hamada, Yukiteru Ono, Kiyoshi Asai, Martin C Frith, Training alignment parameters for arbitrary sequencers with LAST-TRAIN, Bioinformatics, Volume 33, Issue 6, March 2017, Pages 926–928, https://doi.org/10.1093/bioinformatics/btw742|GNU GPLv3
Isoquant|rjibelski, A.D., Mikheenko, A., Joglekar, A. et al. Accurate isoform discovery with IsoQuant using long reads. Nat Biotechnol 41, 915–918 (2023). https://doi.org/10.1038/s41587-022-01565-y|GNU GPLv2
RNA Framework|Incarnato D, Morandi E, Simon LM, Oliviero S. RNA Framework: an all-in-one toolkit for the analysis of RNA structures and post-transcriptional modifications. Nucleic Acids Res. 2018 Sep 19;46(16):e97. doi: 10.1093/nar/gky486. PMID: 29893890; PMCID: PMC6144828.|GNU GPLv3
EternaFold|Wayment-Steele HK, Kladwang W, Strom AI, Lee J, Treuille A, Becka A; Eterna Participants; Das R. RNA secondary structure packages evaluated and improved by high-throughput experiments. Nat Methods. 2022 Oct;19(10):1234-1242. doi: 10.1038/s41592-022-01605-0|BSD-3-Clause
Perbase|Seth Stadick, Perbase, v0.9.0 released on 20.07.2023 https://github.com/sstadick/perbase|MIT license
Sambamba|Artem Tarasov, Albert J. Vilella, Edwin Cuppen, Isaac J. Nijman, Pjotr Prins, Sambamba: fast processing of NGS alignment formats, Bioinformatics, Volume 31, Issue 12, June 2015, Pages 2032–2034, https://doi.org/10.1093/bioinformatics/btv098|GNU GPLv2+
VARNA|Kévin Darty, Alain Denise, Yann Ponty, VARNA: Interactive drawing and editing of the RNA secondary structure, Bioinformatics, Volume 25, Issue 15, August 2009, Pages 1974–1975, https://doi.org/10.1093/bioinformatics/btp250|GNU GPLv3.0
