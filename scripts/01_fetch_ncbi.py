from pathlib import Path
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

Entrez.email = "liyanaps17@gmail.com"

MIN_LEN = 1000
MAX_LEN = 2500
RETMAX = 10

def read_taxa(tsv_path: Path):
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

def search_ids(term: str) -> list[str]:
    with Entrez.esearch(db="nucleotide", term=term, retmax=RETMAX) as h:
        res = Entrez.read(h)
    return list(res.get("IdList", []))

def fetch_first_record(nuc_id: str) -> SeqRecord | None:
    with Entrez.efetch(db="nucleotide", id=nuc_id, rettype="fasta", retmode="text") as h:
        recs = list(SeqIO.parse(h, "fasta"))
    return recs[0] if recs else None

def fetch_reasonable_18s_record(query: str) -> tuple[SeqRecord | None, str | None]:
    # constrain by sequence length to avoid whole genomes/contigs
    term = f"({query}) AND {MIN_LEN}:{MAX_LEN}[SLEN]"
    ids = search_ids(term)
    if not ids:
        return None, None

    for nuc_id in ids:
        rec = fetch_first_record(nuc_id)
        if rec is None:
            continue
        seqlen = len(rec.seq)
        if MIN_LEN <= seqlen <= MAX_LEN:
            return rec, nuc_id

    return None, None


def main():
    taxa_file = Path("data/taxa.tsv")
    out = Path("data/raw/ssu_sequences.fasta")
    out.parent.mkdir(parents=True, exist_ok=True)

    records: list[SeqRecord] = []

    for group, taxon, query in read_taxa(taxa_file):
        rec, nuc_id = fetch_reasonable_18s_record(query)

        if rec is None:
            print(f"[WARN] No suitable SSU hit for {taxon} (len {MIN_LEN}-{MAX_LEN})")
            continue

        rec.id = f"{group}|{taxon}"
        rec.name = rec.id
        rec.description = ""
        records.append(rec)

        print(f"[OK] {taxon} <- id {nuc_id} (len={len(rec.seq)})")

    SeqIO.write(records, out, "fasta")
    print(f"\nWrote {len(records)} sequences to {out}")


if __name__ == "__main__":
    main()
"""Fetch SSU sequences from NCBI for given taxa."""