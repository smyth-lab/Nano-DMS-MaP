import numpy as np
from Bio import SeqIO
import gzip
import argparse
import subprocess

parser = argparse.ArgumentParser(
    prog = "calculate per read stats",
    description = "Calculates median or mean qscores per read in (opt. gzipped) fastq files, correctly")
parser.add_argument("-i", "--fastq_infile")
parser.add_argument("-o", "--outfile_prefix")
parser.add_argument("-s", "--statistic")
parser.add_argument("--histogram", action="store_true")#
parser.add_argument("-rl", "--max_readlength", default="150")

def ascii_to_phred(ascii_qscores):
    return np.array([ord(x) -33 for x in list(ascii_qscores)], dtype=np.byte)

def phred_to_prob(phred_scores):
    phred_scores = phred_scores / -10
    probs = np.power(10,phred_scores)
    
    return probs

def prob_to_phred(prob):
    phred = -10 * np.log10(prob)
    return phred

#takes in list or np array
def mean_qscore(phred_scores):
    phred_scores = np.array(phred_scores) / -10
    probs = np.power(10,phred_scores)
    mean_prob = np.mean(probs)
    
    return prob_to_phred(mean_prob)

def calc_per_read_stats(fastq_infile, outfile_prefix, statistic, histogram):
    
    if statistic == "mean":
        read_mean_fastq = []
        if ".gz" in fastq_infile:
            with gzip.open(fastq_infile, "rt") as handle:
                for lineno, line in enumerate(handle):
                    if (1+lineno) % 4 == 0:
                        read_mean_fastq.append(mean_qscore(ascii_to_phred(list(line.strip()))))
        else:
            with open(fastq_infile, "r") as handle:
                for lineno, line in enumerate(handle):
                    if (1+lineno) % 4 == 0:
                        read_mean_fastq.append(mean_qscore(ascii_to_phred(list(line.strip()))))
        
        read_mean_fastq = np.array(read_mean_fastq)
        
        if histogram:
            hist_data = np.histogram(read_mean_fastq, range=[0,40], bins=40)
            hist_data = np.array([hist_data[0], hist_data[1][1:]], dtype=int)
            np.savetxt(f"{outfile_prefix}_read_mean_histogram.csv", hist_data, fmt='%1i', header=f"min:0; max:40")
        
        np.savetxt(f"{outfile_prefix}_read_mean.csv", read_mean_fastq, fmt='%1.3f')
        
        
    elif statistic == "median":
        read_median_fastq = []
        if ".gz" in fastq_infile:
            with gzip.open(fastq_infile, "rt") as handle:
                for lineno, line in enumerate(handle):
                    if (1+lineno) % 4 == 0:
                        read_median_fastq.append(np.median(ascii_to_phred(list(line.strip()))))
        else:
            with open(fastq_infile, "r") as handle:
                for lineno, line in enumerate(handle):
                    if (1+lineno) % 4 == 0:
                        read_median_fastq.append(np.median(ascii_to_phred(list(line.strip()))))
        
        read_median_fastq = np.array(read_median_fastq)
        
        
        if histogram:
            hist_data = np.histogram(read_median_fastq, range=[0,40], bins=40)
            hist_data = np.array([hist_data[0], hist_data[1][1:]], dtype=int)
            np.savetxt(f"{outfile_prefix}_read_median_histogram.csv", hist_data, fmt='%1i', header=f"min:0; max:40")
        
        np.savetxt(f"{outfile_prefix}_read_median.csv", read_median_fastq, fmt='%1i')
        
        
    else:
        raise("Only median or mean per read statistics are available here!")
        

def calc_read_length(fastq_infile, outfile_prefix, histogram):
    
    read_lengths = []
    if ".gz" in fastq_infile:
        with gzip.open(fastq_infile, "rt") as handle:
            for lineno, line in enumerate(handle):
                if (1+lineno) % 4 == 0:
                    read_lengths.append(len(line.strip()))
    else:
        with open(fastq_infile, "r") as handle:
            for lineno, line in enumerate(handle):
                if (3+lineno) % 4 == 0:
                    read_lengths.append(len(line.strip()))

    read_lengths = np.array(read_lengths)

    if histogram:
        hist_data = np.histogram(read_lengths, range=[0,np.max(read_lengths)], bins=np.max(read_lengths))
        
        hist_data = np.array([hist_data[0], hist_data[1][1:]], dtype=int)
        np.savetxt(f"{outfile_prefix}_read_length_histogram.csv", hist_data, delimiter=";", fmt='%1i', header=f"min:0; max:{np.max(read_lengths)}")

    np.savetxt(f"{outfile_prefix}_read_length.csv", read_lengths, fmt='%1i')


import subprocess 

def generate_quality_matrix(fastq_infile, outfile, max_read_length):

    counter = 0
    ignored = 0
    if ".gz" in fastq_infile:
        num_reads = int(int(subprocess.check_output(f"zcat {fastq_infile} | wc -l", shell=True).split()[0])/4)
        
        read_qscore_matrix = np.full((num_reads, int(max_read_length)), fill_value=-1, dtype=int)
        print("Size of matrix:", read_qscore_matrix.shape)
        with gzip.open(fastq_infile, "rt") as handle:
            for lineno, line in enumerate(handle):
                if (1+lineno) % 4 == 0:
                    qscores = ascii_to_phred(list(line.strip()))
                    if len(qscores) > max_read_length:
                        print("Read is longer than max_read_length, ignoring!")
                        ignored += 1
                    else:
                        read_qscore_matrix[counter, :len(qscores)] = qscores.astype(np.byte)
                        counter += 1
    else:
        num_reads = int(int(subprocess.check_output(f"wc -l {fastq_infile}", shell=True).split()[0])/4)
        
        read_qscore_matrix = np.full((num_reads, int(max_read_length)), fill_value=-1, dtype=np.byte)
        print("Size of matrix:", read_qscore_matrix.shape)
        with open(fastq_infile, "r") as handle:
            for lineno, line in enumerate(handle):
                if (1+lineno) % 4 == 0:
                    qscores = ascii_to_phred(list(line.strip()))
                    if len(qscores) > max_read_length:
                        print("Read is longer than max_read_length, ignoring!")
                        ignored += 1
                    else:
                        read_qscore_matrix[counter, :len(qscores)] = qscores.astype(np.byte)
                        counter += 1
    
    read_qscore_matrix = read_qscore_matrix[:counter, :]
    per_pos_hist = np.apply_along_axis(lambda a: np.histogram(a, range=[0,40], bins=40)[0], 0, read_qscore_matrix)
    
    print(f"Ignored in total {ignored} reads that were longer than max_read_length {max_read_length}.")
    
    np.savetxt(f"{outfile}_read_quality_per_position_histogram.csv", per_pos_hist, fmt="%i")
    
    
if __name__ == "__main__":
    args = parser.parse_args()

    print("Reading in reads from file", args.fastq_infile)
    print("Outfile prefix:", args.outfile_prefix)
    
    print("Statistic:", args.statistic)
    if args.histogram:
        print("Also generating a histogram")
    
    if (args.statistic == "median") | (args.statistic == "mean"):
        calc_per_read_stats(args.fastq_infile, args.outfile_prefix, args.statistic, args.histogram)
        
    if (args.statistic == "length"):
        calc_read_length(args.fastq_infile, args.outfile_prefix, args.histogram)
    
    if (args.statistic == "position_hist"):
        generate_quality_matrix(args.fastq_infile, args.outfile_prefix, int(args.max_readlength))
        
        
