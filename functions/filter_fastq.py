import numpy as np
from Bio import SeqIO
import gzip
import argparse


parser = argparse.ArgumentParser(
    prog = "filter_mean_quality",
    description = "Filters all reads in gzipped fastq files for their mean qscore, correctly")
parser.add_argument("-i", "--infile")
parser.add_argument("-o", "--outfile")
parser.add_argument("-i1", "--infile_r1")
parser.add_argument("-o1", "--outfile_r1")
parser.add_argument("-i2", "--infile_r2")
parser.add_argument("-o2", "--outfile_r2")
parser.add_argument("-m", "--min_mean_qscore")
parser.add_argument("-n", "--number_reads", default="all")
parser.add_argument("-p", "--paired", action="store_true")
parser.add_argument("-nt3", "--no_trim3",action="store_true")
parser.add_argument("-fqo", "--read_mean_fastq_outdir")
parser.add_argument("-fqi", "--read_mean_fastq_prefix")


def ascii_to_phred(ascii_qscores):
    return np.array([ord(x) -33 for x in list(ascii_qscores)], dtype=int)

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

def filter_file(infile, outfile, min_mean_qscore, num_reads = "all", trim_3end = True):
    records = []
    count = 0
    with gzip.open(infile, "rt") as handle:
        
        for read in SeqIO.parse(handle, "fastq"):
            if num_reads != "all":
                if count >= int(num_reads):
                    break
            
            phred_scores = list(read.letter_annotations["phred_quality"])
            #print(phred_scores)
            mean_qscore = mean_qscore(phred_scores)
            
            if mean_qscore > min_mean_qscore:
                not_trimmed = True
                while not_trimmed:
                    if phred_scores[-1] > min_mean_qscore:
                        not_trimmed = False
                    else:
                        phred_scores = phred_scores[:-1]
                        read = read[:-1]
                        
                records.append(read)
                count += 1
    print(f"{len(records)} reads passed the mean qscore filter of {min_mean_qscore}")
    with gzip.open(outfile, "wt") as outhandle:
        SeqIO.write(records, outhandle, "fastq")
        
def filter_paired_files(file_read1, file_read2, outfile1, outfile2, min_mean_qscore, num_reads, read_mean_fastq_prefix, read_mean_fastq_outdir, trim_3end = True):
    records = []
    count = 0
    
    if not read_mean_fastq_prefix:
        read1_mean_fastq = []
        with gzip.open(file_read1, "rt") as handle:
            for lineno, line in enumerate(handle):
                if (1+lineno) % 4 == 0:
                    read1_mean_fastq.append(mean_qscore(ascii_to_phred(list(line.strip()))))

        read2_mean_fastq = []
        with gzip.open(file_read2, "rt") as handle:
            for lineno, line in enumerate(handle):
                if (1+lineno) % 4 == 0:
                    read2_mean_fastq.append(mean_qscore(ascii_to_phred(list(line.strip()))))


        read1_mean_fastq = np.array(read1_mean_fastq)
        read2_mean_fastq = np.array(read2_mean_fastq)

        if read_mean_fastq_outdir:
            np.savetxt(f"{read_mean_fastq_outdir}_mean_qscores_r1.csv", read1_mean_fastq)
            np.savetxt(f"{read_mean_fastq_outdir}_mean_qscores_r2.csv", read2_mean_fastq)
    
    else:
        read1_mean_fastq = np.genfromtxt(f"{read_mean_fastq_prefix}_mean_qscores_r1.csv")
        read2_mean_fastq = np.genfromtxt(f"{read_mean_fastq_prefix}_mean_qscores_r2.csv")
    
    passing_reads_r1 = read1_mean_fastq > min_mean_qscore
    passing_reads_r2 = read2_mean_fastq > min_mean_qscore
    
    combined_mean_fastq = np.logical_and(passing_reads_r1, passing_reads_r2)
    
    print(f"{np.count_nonzero(combined_mean_fastq)} reads of {len(combined_mean_fastq)} passed the mean qscore filter of {min_mean_qscore}")
    
    count = 0
    records = []
    with gzip.open(file_read1, "rt") as handle:
        for filtered, read in zip(combined_mean_fastq, SeqIO.parse(handle, "fastq")):
            if num_reads != "all":
                if count >= int(num_reads):
                    break
            if filtered:
                phred_scores = list(read.letter_annotations["phred_quality"])
                
                not_trimmed = True
                while not_trimmed:
                    if phred_scores[-1] > min_mean_qscore:
                        not_trimmed = False
                    else:
                        phred_scores = phred_scores[:-1]
                        read = read[:-1]
                        
                records.append(read)
                count += 1
                
    with gzip.open(outfile1, "wt") as outhandle:
        SeqIO.write(records, outhandle, "fastq")
        
    records = []
    with gzip.open(file_read2, "rt") as handle:
        for filtered, read in zip(combined_mean_fastq, SeqIO.parse(handle, "fastq")):
            if num_reads != "all":
                if count >= int(num_reads):
                    break
            if filtered:
                phred_scores = list(read.letter_annotations["phred_quality"])
                
                not_trimmed = True
                while not_trimmed:
                    if phred_scores[-1] > min_mean_qscore:
                        not_trimmed = False
                    else:
                        phred_scores = phred_scores[:-1]
                        read = read[:-1]
                        
                records.append(read)
                count += 1
                
    with gzip.open(outfile2, "wt") as outhandle:
        SeqIO.write(records, outhandle, "fastq")
        
        

        
        
if __name__ == "__main__":
    args = parser.parse_args()
    if args.paired: 
        print("Paired flag detected.")
        for io_path in ["infile_r1", "infile_r2", "outfile_r1", "outfile_r2"]:
            if io_path not in args:
                raise("Error: specified paired reads, but did not provide pairing-specific input/output files (-i1, -i2, -o1, -o2)")
        print("Reading in read 1 from file", args.infile_r1)
        print("Reading in read 2 from file", args.infile_r2)
        print("Outfile read 1:", args.outfile_r1)
        print("Outfile read 2:", args.outfile_r2)
        
        if args.read_mean_fastq_prefix:
            print("Loading mean fastq per read data from", args.read_mean_fastq_prefix)
            
        if args.read_mean_fastq_outdir:
            print("Storing mean fastq per read data at", args.read_mean_fastq_outdir)
    else:
        for io_path in ["infile", "outfile"]:
            if io_path not in args:
                raise("Error: did not provide pairing-specific input/output files (-i, -o)")
        print("Reading in reads from file", args.infile)
        print("Outfile:", args.outfile)
    
    print("Removing reads with mean qscore below", args.min_mean_qscore)
    print("Outputting", args.number_reads, "reads")
    if args.paired:
        filter_paired_files(args.infile_r1, args.infile_r2, args.outfile_r1, args.outfile_r2, int(args.min_mean_qscore), args.number_reads, args.read_mean_fastq_prefix, args.read_mean_fastq_outdir, ~args.no_trim3)
    else:
        filter_file(args.infile, args.outfile, int(args.min_mean_qscore), args.number_reads, ~args.no_trim3)

        
        
