# SSU rRNA Phylogenetics of Apicomplexan Parasites

phylogenetic analysis exploring how small-subunit (SSU/18S) rRNA captures evolutionary relationships among Apicomplexan parasites.
The project combines sequence retrieval, curation, alignment, tree inference, and visualisation in a reproducible workflow.

---

# Motivation
SSU rRNA is one of the most widely used molecular markers for studying eukaryotic evolution, yet its strong structural constraints raise questions about how much phylogenetic signal it retains—particularly in parasitic lineages subject to extreme evolutionary pressures.

Apicomplexans (Plasmodium, Toxoplasma, Babesia, Cryptosporidium) provide an interesting system to explore this balance between conservation and divergence.

---

# Question

**To what extent can SSU rRNA resolve evolutionary relationships within Apicomplexa when analysed in a robust, outgroup-aware phylogenetic framework?**

---
# Approach

1. Retrieve SSU rRNA sequences from NCBI

2. Curate and standardise sequences to remove partial or low-quality records

3. Align sequences using MAFFT

4. Infer phylogeny with FastTree

5. Visualise and root trees using multiple eukaryotic outgroups

---
The analysis is organised using Snakemake to ensure reproducibility.

## Structure 
```bash
scripts/   # analysis scripts
data/      # raw and curated sequences
results/   # alignments, trees, figures
Snakefile  # workflow definition
```
## How to run

```bash
snakemake -n
snakemake -j 1
```
Recommended to run on a Linux filesystem in windows

---
## Observations

- Apicomplexan taxa form a coherent, well-supported clade
- Closely related Plasmodium species cluster with high confidence
- Support decreases at deeper nodes, consistent with SSU rRNA conservation

Overall, the inferred topology aligns with established evolutionary relationships, illustrating both the strengths and limitations of SSU rRNA for phylogenetic inference.

---

## Notes

This repository focuses on methodological clarity and biological interpretation, rather than exhaustive taxon sampling. The workflow is general and can be applied to other organismal groups or molecular markers.

---

Author: Liyana Saleem




