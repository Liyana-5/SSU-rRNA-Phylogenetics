
# SSU SSU rRNA Phylogenetics of Apicomplexan Parasites
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


The analysis is organised using Snakemake to ensure reproducibility.
---

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
## Limitations and interpretation

While SSU (18S) rRNA provides a conserved and widely used marker for eukaryotic phylogenetics, its strong structural and functional constraints limit the amount of independent evolutionary signal it can retain, particularly over deep evolutionary timescales. In parasitic lineages such as Apicomplexa, lineage-specific rate heterogeneity and compositional bias further increase susceptibility to phylogenetic artefacts, including long-branch attraction.

Accordingly, the inferred trees are best interpreted as capturing broad evolutionary relationships rather than providing definitive resolution of deep branching order. These results highlight both the utility and the limitations of SSU rRNA, reinforcing the importance of careful outgroup selection, sequence curation, and complementary multi-gene or structure-informed approaches for studying ancient ribosomal systems.

---
## Notes

This repository focuses on methodological clarity and biological interpretation within a defined project scope, rather than on exhaustive taxon sampling or explicit constraint-based hypothesis testing. The workflow is modular and transferable, and can be extended to incorporate additional taxa, markers, or more complex evolutionary models.

---

Author : Liyana Saleem

Author: Liyana Saleem




