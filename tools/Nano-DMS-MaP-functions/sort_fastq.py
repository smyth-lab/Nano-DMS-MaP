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
parser.add_argument("-fo", "--fastq_out_prefix")
parser.add_argument("-ir", "--isoquant-read_assignment_file")
parser.add_argument("-rdf", "--read_assignment_df")
parser.add_argument("-tsv", "--read_assignment_file")
parser.add_argument("-o", "--outfile_prefix")
parser.add_argument("-m", "--read_id_mapping_file")
parser.add_argument("-u", "--output_unsorted", action="store_true")

def sort_fastq_by_tsv(fastq_infile, read_assignment_file, fastq_outprefix, output_unsorted=False):
    
    read_df = pd.read_csv(read_assignment_file, sep="\t", usecols=[0,1], names=["assignment", "read_id"], keep_default_na=False, na_values=['_'])
    
    read_appearances = read_df["read_id"].value_counts().reset_index()
    read_appearances.columns = ["read_id", "count"]
    total_num_reads = len(read_appearances)
    unique_appearances = read_appearances[read_appearances["count"] == 1]
    unique_reads = len(unique_appearances)
    
    print("Detected in total", total_num_reads, "aligned reads.")
    print("Of these,", unique_reads, "reads were assigned uniquely.")
    
    read_df = read_df[read_df["read_id"].isin(unique_appearances["read_id"].values)]
    
    assignments = sorted(np.unique(read_df["assignment"].values))
    
    print("Detected in total", len(assignments), "assignments with at least one unique read.")
    
    read_dict = read_df[["read_id", "assignment"]].set_index("read_id").to_dict()["assignment"]
    
    print("Reading in files")
    
    sorted_reads = {}
    
    for assignment in assignments:
        sorted_reads[assignment] = []

    num_sorted = 0
    num_not_sorted = 0
    unsorted_reads = []

    if ".gz" in fastq_infile:
        
        with gzip.open(fastq_infile, "rt") as infile:
            for read in SeqIO.parse(infile, "fastq"):
                if " " in read.id:
                    read_id = read.id.split(" ")[0]
                else:
                    read_id = read.id
                try:
                    sorted_reads[read_dict[read_id]].append(read)
                    num_sorted += 1
                except:
                    num_not_sorted += 1
                    if output_unsorted:
                        unsorted_reads.append(read)
    
    
    outdir = "/".join(fastq_outprefix.split("/")[:-1])
    os.makedirs(outdir, exist_ok=True)
    assignment_counts = []
    for assignment in assignments:
        assignment_reads = sorted_reads[assignment]
        assignment_count = len(assignment_reads)
        assignment_counts.append((assignment, assignment_count))
        with gzip.open(f"{fastq_outprefix}{assignment}.fastq.gz", "wt") as outfile:
            SeqIO.write(assignment_reads, outfile, "fastq")
        if output_unsorted:
            with gzip.open(f"{fastq_outprefix}unsorted.fastq.gz", "wt") as outfile:
                SeqIO.write(unsorted_reads, outfile, "fastq")
        print("Sorted", assignment_count, "reads for assignment", assignment)
        
    print("Sorted reads in total:", num_sorted)
    print("Not sorted reads in total:", num_not_sorted)


def sort_fastq_by_isoform_isoquant(fastq_infile, isoquant_read_assignment_file, fastq_outdir, output_unsorted=False, df = False):
    if df:
        read_df = pd.read_pickle(isoquant_read_assignment_file)
    else:
        read_df = pd.read_csv(isoquant_read_assignment_file, sep="\t", header=2)
    
    total_num_reads = len(read_df)
    
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

def sort_fastq_isoform_by_4sU(fastq_dir, fastq_out_prefix, read_id_mapping_file):
    import pandas as pd

    read_id_mappings_df = pd.read_pickle(read_id_mapping_file)
    isoforms = np.unique(read_id_mappings_df["isoform"].values)

    for isoform in isoforms:
        fastq_file = f"{fastq_dir}/{isoform}.fastq.gz"

        print("Sorting isoform", isoform)

        read_id_mappings = read_id_mappings_df[read_id_mappings_df["isoform"] == isoform]

        read_dict = {}
        sorted_reads = {}

        for _, row in read_id_mappings.iterrows():
            sU_status = row["4sU_status"]
            read_ids = row["read_id"]
            sorted_reads[sU_status] = []

            print(f"Found", len(read_ids), "read ids with 4sU status", sU_status)
            for read_id in read_ids:
                read_dict[read_id] = sU_status

        num_sorted = 0
        num_not_sorted = 0
        if ".gz" in fastq_file:

            with gzip.open(fastq_file, "rt") as infile:
                for read in SeqIO.parse(infile, "fastq"):
                    try:
                        sorted_reads[read_dict[read.id]].append(read)
                        num_sorted += 1
                    except:
                        num_not_sorted += 1

        for sU_status in list(sorted_reads.keys()):
            sU_status_reads = sorted_reads[sU_status]
            sU_status_count = len(sU_status_reads)
            os.makedirs(f"{fastq_out_prefix}/4sU_{sU_status}/", exist_ok=True)
            with gzip.open(f"{fastq_out_prefix}/4sU_{sU_status}/{isoform}.fastq.gz", "wt") as outfile:
                SeqIO.write(sU_status_reads, outfile, "fastq")

            print("Sorted", sU_status_count, "reads for 4sU status", sU_status)

        print("Sorted reads in total:", num_sorted)
        print("Not sorted reads in total:", num_not_sorted)

if __name__ == "__main__":
    args = parser.parse_args()

    print("Reading in reads from file", args.fastq_infile)
    if args.output_unsorted:
        print("Also outputting unsorted reads")

    if args.isoquant_read_assignment_file:
        print("Detected isoquant_read_assignment_file option, performing fastq per isoform sorting according to Isoquant assignments.")
        print("Reading read assignment file from",args.isoquant_read_assignment_file)
        print("Outfile prefix:", args.outfile_prefix)
        sort_fastq_by_isoform_isoquant(args.fastq_infile, args.isoquant_read_assignment_file, args.fastq_out_prefix, args.output_unsorted)

    elif args.read_assignment_df:
        print("Detected assignment_df option, performing fastq sorting")
        print("Reading in assignment_df from", args.read_assignment_df)
        print("Outfile prefix:", args.outfile_prefix)
        sort_fastq_by_isoform_isoquant(args.fastq_infile, args.read_assignment_df, args.fastq_out_prefix, args.output_unsorted, df = True)
    elif args.read_id_mapping_file:
        print("Detected read_id_mappings_file option, resorting per isoform reads according to 4sU status")
        sort_fastq_isoform_by_4sU(args.fastq_infile, args.fastq_out_prefix, args.read_id_mapping_file)
    elif args.read_assignment_file:
        print("Detected tsv read assignment file, performing fastq sorting according to assignment.")
        sort_fastq_by_tsv(args.fastq_infile, args.read_assignment_file, args.outfile_prefix, args.output_unsorted)
