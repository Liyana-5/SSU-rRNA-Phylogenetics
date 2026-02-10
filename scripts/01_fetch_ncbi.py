#import necessary Libraries 
from pathlib import Path
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

#make connection with entrez and set email
Entrez.email = "liyanaps17@gmail.com"

#define constants for sequence length and maximum number of records to fetch
MIN_LEN = 1000
MAX_LEN = 2500
RETMAX = 10

#function to read taxa from a tsv file and return a list of tuples containing group, taxon, and query
def read_taxa(tsv_path: Path):
    """Read taxa from a TSV file.
    The TSV file should have three columns: group, taxon, and query.
    Returns a list of tuples: (group, taxon, query).
    The query is used to search for SSU sequences in NCBI, and should be a valid Entrez search term.
    """
    rows = []
    with tsv_path.open() as f:
        next(f)  # skip header
        for line in f:
            line = line.strip()
            if not line:
                continue
            group, taxon, query = line.split("\t")
            rows.append((group, taxon, query))
    return rows
#function to search for nucleotide IDs in NCBI using a search term and return a list of IDs
def search_ids(term: str) -> list[str]:
    """Search for nucleotide IDs in NCBI using a search term
    Returns a list of nucleotide IDs (strings).
    """
    with Entrez.esearch(db="nucleotide", term=term, retmax=RETMAX) as h:
        res = Entrez.read(h)
    return list(res.get("IdList", []))
#function to fetch the first record for a given nucleotide ID and return it as a SeqRecord object
def fetch_first_record(nuc_id: str) -> SeqRecord | None:
    """
Fetch the first record for a given nucleotide ID.
Returns a SeqRecord object if a record is found, or None if no records are found
    
    """
    with Entrez.efetch(db="nucleotide", id=nuc_id, rettype="fasta", retmode="text") as h:
        recs = list(SeqIO.parse(h, "fasta"))
    return recs[0] if recs else None
#function to fetch a reasonable 18S rRNA record for a given query, constrained by sequence length
def fetch_reasonable_18s_record(query: str) -> tuple[SeqRecord | None, str | None]:
    """
Fetch a reasonable 18S rRNA record for a given query.
The query is used to search for SSU sequences in NCBI, and should be a valid Entrez search term.
The function constrains the search by sequence length to avoid whole genomes/contigs.
Returns a tuple: (SeqRecord or None, nucleotide ID or None).
    """
    # constrain by sequence length to avoid whole genomes/contigs
    term = f"({query}) AND {MIN_LEN}:{MAX_LEN}[SLEN]"
    ids = search_ids(term)
    if not ids:
        return None, None
    #iterate over the fetched IDs, fetch the corresponding records, and return the first one that meets the length criteria
    for nuc_id in ids:
        rec = fetch_first_record(nuc_id)
        if rec is None:
            continue
        seqlen = len(rec.seq)
        if MIN_LEN <= seqlen <= MAX_LEN:
            return rec, nuc_id

    return None, None

#main function to read taxa, fetch SSU sequences, and write them to a fasta file
def main():
    #define paths for input taxa file and output fasta file, and create output directory if it doesn't exist
    taxa_file = Path("data/taxa.tsv")
    out = Path("data/raw/ssu_sequences.fasta")
    out.parent.mkdir(parents=True, exist_ok=True)

    records: list[SeqRecord] = []
    #iterate over taxa, fetch reasonable 18S records, and store them in a list of SeqRecord objects
    for group, taxon, query in read_taxa(taxa_file):
        rec, nuc_id = fetch_reasonable_18s_record(query)
        #print(f"Query '{query}' for taxon '{taxon}' returned nucleotide ID: {nuc_id}")
        if rec is None:
            print(f"[WARN] No suitable SSU hit for {taxon} (len {MIN_LEN}-{MAX_LEN})")
            continue
        #set the record ID to a combination of group and taxon, and store the record in the list
        rec.id = f"{group}|{taxon}"
        rec.name = rec.id
        rec.description = ""
        records.append(rec)

        print(f"[OK] {taxon} <- id {nuc_id} (len={len(rec.seq)})")
    #write the fetched records to a fasta file
    SeqIO.write(records, out, "fasta")
    print(f"\nWrote {len(records)} sequences to {out}")


if __name__ == "__main__":
    main()
"""Fetch SSU sequences from NCBI for given taxa."""