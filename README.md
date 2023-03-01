# Nano-DMS-MaP

Code to perform Nanopore-DMS-MaP RNA structural probing analysis. 

In this repository we provide both the code to process samples, as well as for visualization. 

Briefly, the pipeline has the following steps: 

1. Basecalling (and demultiplexing) of fast5 files from native barcoded SQK-NBD112 or SQK-NBD114 (Q20+ kits) runs with guppy. 
2. Read-to-isoform assignment with Isoquant, followed by sorting. 
3. Per isoform alignment with LAST. 
4. RNA-Framework analysis
  4.a rf-count with optimized settings on LAST-aligned BAM files. 
  4.b rf-norm to quantify reactivity rates (we recommend having a non-modified control to perform Siegfried et al. normalization)
  4. c rf-correlate and rf-combine to evaluate the reproducability of reactivity (>0.9 expected) and calculate the mean reactivity of repplicates. 
5. Option A: Evaluating the agreement with published secondary structures by calculating Receiver-Operator Curves
5. Option B: Performing de novo RNA structure prediction using EternaFold (or rf-fold)

Additional workflows in this repository include:
- Generating a custom gtf file for isoquant to take into account primer binding sites
- Mutation profile analysis of BAM files using perbase
- Analyzing the frequency of mutations per read from RNA-Framework mm files
- Generating .varna files of known or predicted structures colored by reactivities
