from pathlib import Path
from Bio import SeqIO

def main ():
    input_fasta = Path("data/raw/ssu_sequences.fasta")
    output_fasta = Path("data/curated/ssu_curated_sequences.fasta")
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
        
    min_length = 1000
    kept = []
    
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