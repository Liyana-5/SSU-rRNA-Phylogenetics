#!/usr/bin/env python3
"""
Plot a phylogenetic tree from a Newick file and save as an image.

Assumes terminal labels look like: Group|Taxon_name
Example: Apicomplexa|Plasmodium_falciparum
"""

from __future__ import annotations

from pathlib import Path
from Bio import Phylo
import matplotlib.pyplot as plt
import argparse
from typing import Optional


def split_label(name: str) -> tuple[str, str]:
    """Return (group, taxon) from 'group|taxon'. If not present, group='Unknown'."""
    if not name:
        return "Unknown", ""
    if "|" in name:
        g, t = name.split("|", 1)
        return g, t
    return "Unknown", name


def clean_terminal_label(name: str, show_group: bool = True) -> str:
    group, taxon = split_label(name)
    taxon = taxon.replace("_", " ").strip()
    if show_group and group != "Unknown":
        return f"{taxon} ({group})"
    return taxon


def terminals_in_group(tree, group_name: str):
    out = []
    for t in tree.get_terminals():
        g, _ = split_label(t.name or "")
        if g == group_name:
            out.append(t)
    return out


def root_by_outgroup_clade(tree, outgroup_name: str) -> Optional[str]:
    """Root on MRCA of all outgroup tips if >=2, else on single tip if ==1."""
    og_terms = terminals_in_group(tree, outgroup_name)
    if len(og_terms) >= 2:
        og_mrca = tree.common_ancestor(og_terms)
        tree.root_with_outgroup(og_mrca)
        return "outgroup_clade(MRCA)"
    if len(og_terms) == 1:
        tree.root_with_outgroup(og_terms[0])
        return "single_outgroup_tip"
    return None


def root_by_most_distant_outgroup_tip(tree, outgroup_name: str) -> Optional[str]:
    """
    Fallback rooting: pick the outgroup tip most distant from all non-outgroup tips.
    Uses branch-length distances in the tree.
    """
    og_terms = terminals_in_group(tree, outgroup_name)
    non_og_terms = [t for t in tree.get_terminals()
                    if split_label(t.name or "")[0] != outgroup_name]

    if not og_terms or not non_og_terms:
        return None

    best_tip = None
    best_score = -1.0

    for og in og_terms:
        # For each outgroup tip, compute its minimum distance to any non-outgroup tip
        dmin = min(tree.distance(og, ing) for ing in non_og_terms)
        if dmin > best_score:
            best_score = dmin
            best_tip = og

    if best_tip is None:
        return None

    tree.root_with_outgroup(best_tip)
    return f"most_distant_outgroup_tip:{best_tip.name}"


def main():
    p = argparse.ArgumentParser(description="Plot a Newick tree to a PNG (auto-rooted, readable).")
    p.add_argument("-i", "--input", required=True, help="Input Newick file")
    p.add_argument("-o", "--output", required=True, help="Output image file (e.g. PNG)")
    p.add_argument("--outgroup", default="Outgroup",
                   help="Group name (left of '|') to use as outgroup (default: Outgroup)")
    p.add_argument("--title", default="SSU rRNA phylogeny (MAFFT + FastTree)")
    p.add_argument("--width", type=float, default=12)
    p.add_argument("--height", type=float, default=10)
    p.add_argument("--hide_support", action="store_true",
                   help="Hide internal node support/confidence values to reduce clutter")
    p.add_argument("--show_group_in_label", action="store_true",
                   help="Show group in tip labels, e.g. 'Homo sapiens (Outgroup)'")
    args = p.parse_args()

    tree_path = Path(args.input)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    tree = Phylo.read(str(tree_path), "newick")

    # Rooting strategy: outgroup clade -> most distant outgroup tip -> midpoint
    method = root_by_outgroup_clade(tree, args.outgroup)
    if method is None:
        method = root_by_most_distant_outgroup_tip(tree, args.outgroup)
    if method is None:
        tree.root_at_midpoint()
        method = "midpoint"

    tree.ladderize()

    # Plot
    fig = plt.figure(figsize=(args.width, args.height))
    ax = fig.add_subplot(1, 1, 1)

    Phylo.draw(
        tree,
        do_show=False,
        axes=ax,
        # Only label terminals (tips)
        label_func=lambda c: clean_terminal_label(c.name or "", show_group=args.show_group_in_label)
        if c.is_terminal() else None,
        show_confidence=(not args.hide_support),
    )

    ax.set_title(f"{args.title}\n(rooting: {method})", fontsize=14)

    # Make it cleaner: remove axes/ticks/spines
    ax.set_xlabel("")
    ax.set_ylabel("")
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()

    print(f"[INFO] Rooting method: {method}")
    print(f"[DONE] Saved: {out_path}")


if __name__ == "__main__":
    main()
