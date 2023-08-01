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


