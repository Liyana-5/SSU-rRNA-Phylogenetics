rule all:
    input:
        "results/figures/phylogenetic_tree_clean.png"


rule fetch_ssu:
    output:
        "data/raw/ssu_sequences.fasta"
    shell:
        "python scripts/01_fetch_ncbi.py"


rule curate_fasta:
    input:
        "data/raw/ssu_sequences.fasta"
    output:
        "data/curated/ssu_curated_sequences.fasta"
    shell:
        "python scripts/02_curate_fasta.py"


rule align_and_tree:
    input:
        "data/curated/ssu_curated_sequences.fasta"
    output:
        tree="results/trees/phylogenetic_tree.nwk"
    shell:
        """
        bash scripts/03_align_tree.sh \
            {input} \
            results \
            4
        """


rule plot_tree:
    input:
        "results/trees/phylogenetic_tree.nwk"
    output:
        "results/figures/phylogenetic_tree_clean.png"
    shell:
        """
        python scripts/04_plot_tree.py \
            -i {input} \
            -o {output} \
            --outgroup Outgroup \
            --show_group_in_label
        """
