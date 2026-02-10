#import necesarry libraries
from pathlib import Path
from Bio import SeqIO

"""
This script reads a FASTA file containing small subunit (SSU) rRNA sequences, 
filters out sequences that are shorter than a specified minimum length, and writes the curated sequences to a new FASTA file.
The script also prints out which sequences were kept and which were dropped based on their length."""
#define the main function to perform the curation of the FASTA file
def main ():
    #define input and output FASTA file paths, and create the output directory if it doesn't exist
    input_fasta = Path("data/raw/ssu_sequences.fasta")
    output_fasta = Path("data/curated/ssu_curated_sequences.fasta")
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    #define a minimum length threshold for sequences to be kept, and initialize a list to store the kept sequences   
    min_length = 1000
    kept = []
    #iterate over each sequence record in the input FASTA file, clean and filter the sequences based on length, and write the kept sequences to the output FASTA file
    for r in SeqIO.parse(input_fasta, "fasta"):
        seq = str(r.seq).upper().replace("U", "T")
        seq = "".join(c for c in seq if c in "ACGTN-")
        if len(seq) < min_length:
            print (f" [DROP] {r.id} length {len(seq)} < {min_length}")
            continue
        r.seq = type(r.seq)(seq)
        kept.append(r)
        SeqIO.write(kept, output_fasta, "fasta")
        print (f" [KEEP] {r.id} length {len(seq)} >= {min_length}")

if __name__ == "__main__":
    main()
