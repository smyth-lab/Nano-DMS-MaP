#! /home/pbohn/miniconda3/envs/nanodms/bin/python

import numpy as np
from Bio import SeqIO
import gzip
import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser(
    prog = "split fastq file into one for each isoform",
    description = "Reads in a fastq file and isoquant read_assignment file and sorts all reads uniquely mapped to an isoform into a new fastq file")
parser.add_argument("-fi", "--fastq_infile")
parser.add_argument("-ri", "--read_assignment_file")
parser.add_argument("-rdf", "--read_assignment_df")
parser.add_argument("-rli", "--LAST_read_assignment_file")
parser.add_argument("-o", "--outfile_prefix")
parser.add_argument("-m", "--read_id_mapping_file")
parser.add_argument("-u", "--output_unsorted", action="store_true")



def sort_fastq_by_isoform_LAST(fastq_infile, LAST_read_assignment_file, fastq_outdir, output_unsorted=False):
    
    read_df = pd.read_csv(LAST_read_assignment_file, sep="\t", names=["transcripts", "read_id"], compression="gzip")
    
        
    read_appearances = read_df["read_id"].value_counts().reset_index()
    read_appearances.columns = ["read_id", "count"]
    total_num_reads = len(read_appearances)
    unique_appearances = read_appearances[read_appearances["count"] == 1]
    unique_reads = len(unique_appearances)
    
    print("Detected in total", total_num_reads, "aligned reads.")
    print("Of these,", unique_reads, "reads were assigned uniquely by LAST.")
    
    read_df = read_df[read_df["read_id"].isin(unique_appearances["read_id"].values)]
    
    isoforms = sorted(np.unique(read_df["transcripts"].values))
    
    print("Detected in total", len(isoforms), "RNAs with at least one unique assignment.")
    
    read_dict = read_df[["read_id", "transcripts"]].set_index("read_id").to_dict()["transcripts"]
    
    print("Reading in files")
    
    sorted_reads = {}
    
    for isoform in isoforms:
        sorted_reads[isoform] = []
    
    os.makedirs(fastq_outdir, exist_ok=True)
    
    num_sorted = 0
    num_not_sorted = 0
    unsorted_reads = []

    if ".gz" in fastq_infile:
        
        with gzip.open(fastq_infile, "rt") as infile:
            for read in SeqIO.parse(infile, "fastq"):
                try:
                    sorted_reads[read_dict[read.id]].append(read)
                    num_sorted += 1
                except:
                    num_not_sorted += 1
                    if output_unsorted:
                        unsorted_reads.append(read)

    isoform_counts = []
    for isoform in isoforms:
        isoform_reads = sorted_reads[isoform]
        isoform_count = len(isoform_reads)
        isoform_counts.append((isoform, isoform_count))
        with gzip.open(f"{fastq_outdir}/{isoform}.fastq.gz", "wt") as outfile:
            SeqIO.write(isoform_reads, outfile, "fastq")
        if output_unsorted:
            with gzip.open(f"{fastq_outdir}/unsorted.fastq.gz", "wt") as outfile:
                SeqIO.write(unsorted_reads, outfile, "fastq")
        print("Sorted", isoform_count, "reads for isoform", isoform)
        
    print("Sorted reads in total:", num_sorted)
    print("Not sorted reads in total:", num_not_sorted)


def sort_fastq_by_isoform_isoquant(fastq_infile, read_assignment_file, fastq_outdir, output_unsorted=False, df = False):
    if df:
        read_df = pd.read_pickle(read_assignment_file)
    else:
        read_df = pd.read_csv(read_assignment_file, sep="\t", header=2)
    
    total_num_reads = len(read_df["#read_id"].drop_duplicates())
    
    read_df.drop(read_df.query('assignment_type not in ["unique", "unique_minor_difference"]').index, inplace=True)
    
    print("Detected in total", total_num_reads, "aligned reads.")
    print("Of these,", len(read_df), "reads were assigned uniquely by isoquant.")
    
    isoforms = sorted(np.unique(read_df["isoform_id"].values))

    print("Detected in total", len(isoforms), "isoforms with at least one unique assignment.")
    
    read_dict = read_df[["#read_id", "isoform_id"]].set_index("#read_id").to_dict()["isoform_id"]
    
    print("Reading in files")
    
    sorted_reads = {}
    
    for isoform in isoforms:
        sorted_reads[isoform] = []
    
    os.makedirs(fastq_outdir, exist_ok=True)
    
    num_sorted = 0
    num_not_sorted = 0
    unsorted_reads = []

    if ".gz" in fastq_infile:

        with gzip.open(fastq_infile, "rt") as infile:
            for read in SeqIO.parse(infile, "fastq"):
                try:
                    sorted_reads[read_dict[read.id]].append(read)
                    num_sorted += 1
                except:
                    num_not_sorted += 1
                    if output_unsorted:
                        unsorted_reads.append(read)
    isoform_counts = []
    for isoform in isoforms:
        isoform_reads = sorted_reads[isoform]
        isoform_count = len(isoform_reads)
        isoform_counts.append((isoform, isoform_count))
        with gzip.open(f"{fastq_outdir}/{isoform}.fastq.gz", "wt") as outfile:
            SeqIO.write(isoform_reads, outfile, "fastq")
    if output_unsorted:
        with gzip.open(f"{fastq_outdir}/unsorted.fastq.gz", "wt") as outfile:
            SeqIO.write(unsorted_reads, outfile, "fastq")

        print("Sorted", isoform_count, "reads for isoform", isoform)

    print("Sorted reads in total:", num_sorted)
    print("Not sorted reads in total:", num_not_sorted)
    with open(f"{fastq_outdir}/isoform_counts.csv", "w") as csv_file:
        for isoform, count in isoform_counts:
            csv_file.write(f"{isoform}; {count}\n")

if __name__ == "__main__":
    args = parser.parse_args()

    print("Reading in reads from file", args.fastq_infile)
    if args.output_unsorted:
        print("Also outputting unsorted reads")

    if args.read_assignment_file:
        print("Detected read_assignment_file option, performing fastq per isoform sorting")
        print("Reading read assignment file from",args.read_assignment_file)
        print("Outfile prefix:", args.outfile_prefix)
        sort_fastq_by_isoform_isoquant(args.fastq_infile, args.read_assignment_file, args.outfile_prefix, args.output_unsorted)

    elif args.read_assignment_df:
        print("Detected assignment_df option, performing fastq sorting")
        print("Reading in assignment_df from", args.read_assignment_df)
        print("Outfile prefix:", args.outfile_prefix)
        sort_fastq_by_isoform_isoquant(args.fastq_infile, args.read_assignment_df, args.outfile_prefix, args.output_unsorted, df = True)
    
    elif args.LAST_read_assignment_file:
        print("Detected read_to_RNA file option, performing fastq sorting from LAST assingment file")
        sort_fastq_by_isoform_LAST(args.fastq_infile, args.LAST_read_assignment_file, args.outfile_prefix, args.output_unsorted)
