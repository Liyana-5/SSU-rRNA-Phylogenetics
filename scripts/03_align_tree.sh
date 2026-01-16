#!/usr/bin/env bash
# Align sequences and build a phylogenetic tree
set -euo pipefail

# Usage check
if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <input_fasta> <output_dir> [threads]"
  echo "Example: $0 data/curated/ssu_curated.fasta results/run1 4"
  exit 1
fi
# Input parameters

INPUT_FASTA=$1
THREADS=${3:-4} # Default to 4 threads if not provided

ALIGNMENT_DIR="results/alignments"
TREE_DIR="results/trees"

ALIGNMENT_FILE="${ALIGNMENT_DIR}/aligned_sequences.fasta"
TREE_FILE="${TREE_DIR}/phylogenetic_tree.nwk"
MAFFT_LOG="${ALIGNMENT_DIR}/mafft.log"

# Create output directory if it doesn't exist

mkdir -p "${ALIGNMENT_DIR}" "${TREE_DIR}"

# Step 1: Align sequences using MAFFT
echo "[INFO] Running MAFFT..."
mafft --thread "${THREADS}" --auto "${INPUT_FASTA}" > "${ALIGNMENT_FILE}" 2> "${MAFFT_LOG}"

# Step 2: Build a phylogenetic tree using FastTree
echo "[INFO] Running FastTree..."
FastTree -nt "${ALIGNMENT_FILE}" > "${TREE_FILE}"

echo "[DONE] Alignment and phylogenetic tree construction completed."
echo "Aligned sequences: ${ALIGNMENT_FILE}"
echo "MAFFT log:        ${MAFFT_LOG}"
echo "Tree (Newick):    ${TREE_FILE}"

