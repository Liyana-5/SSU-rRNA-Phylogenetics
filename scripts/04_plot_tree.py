"""
Plot a phylogenetic tree from a Newick file and save as an image.

Assumes terminal labels look like: Group|Taxon_name
Example: Apicomplexa|Plasmodium_falciparum
"""
# import necessary libraries
from __future__ import annotations

from pathlib import Path
from Bio import Phylo
import matplotlib.pyplot as plt
import argparse
from typing import Optional

# Helper functions for parsing labels and rooting the tree
def split_label(name: str) -> tuple[str, str]:
    """Return (group, taxon) from 'group|taxon'. If not present, group='Unknown'."""
    if not name:
        return "Unknown", ""
    if "|" in name:
        g, t = name.split("|", 1)
        return g, t
    return "Unknown", name

# Clean up terminal labels for plotting, optionally showing group in parentheses
def clean_terminal_label(name: str, show_group: bool = True) -> str:
    # Given a terminal label like "Group|Taxon_name", return a cleaned-up label for plotting.
    group, taxon = split_label(name)
    # Replace underscores with spaces in the taxon name for better readability
    taxon = taxon.replace("_", " ").strip()
    # Optionally include the group in parentheses if it's not "Unknown"
    if show_group and group != "Unknown":
        return f"{taxon} ({group})"
    return taxon

# Get terminal nodes (tips) in the tree that belong to a specific group
def terminals_in_group(tree, group_name: str):
    # Return a list of terminal nodes in the tree that belong to the specified group.
    out = []
    # Iterate over all terminal nodes in the tree and check if their group matches the specified group name.
    for t in tree.get_terminals():
        # Extract the group from the terminal label and compare it to the specified group name.
        g, _ = split_label(t.name or "")
        if g == group_name:
            out.append(t)
    return out

# Root the tree using the specified outgroup group name. First try to root on the MRCA of all outgroup tips, then fallback to the most distant outgroup tip if needed.
def root_by_outgroup_clade(tree, outgroup_name: str) -> Optional[str]:
    """Root on MRCA of all outgroup tips if >=2, else on single tip if ==1."""
    # First try to root on the MRCA of all outgroup tips if there are 2 or more, which is more robust. If only 1 outgroup tip is found, root on that single tip. If no outgroup tips are found, return None to indicate failure.
    og_terms = terminals_in_group(tree, outgroup_name)
    # If there are 2 or more outgroup tips, root on their MRCA (most recent common ancestor) to ensure a more stable rooting. If there is exactly 1 outgroup tip, root on that tip. If there are no outgroup tips, return None to indicate that this rooting strategy failed.
    if len(og_terms) >= 2:
        # Find the MRCA of all outgroup tips and root the tree using that MRCA as the outgroup. This is more robust than rooting on a single tip, as it uses the shared ancestry of multiple outgroup taxa to determine the root.
        og_mrca = tree.common_ancestor(og_terms)
        # Root the tree using the MRCA of the outgroup tips as the outgroup. This will place the root on the branch leading to that MRCA, effectively rooting the tree based on the shared ancestry of the outgroup taxa.
        tree.root_with_outgroup(og_mrca)
        return "outgroup_clade(MRCA)"
    # If there is exactly 1 outgroup tip, root on that single tip. This is less robust than using the MRCA of multiple outgroup tips, but can still be used as a fallback if only one outgroup taxon is available.
    if len(og_terms) == 1:
        tree.root_with_outgroup(og_terms[0])
        return "single_outgroup_tip"
    return None

# If outgroup rooting fails, root on the most distant outgroup tip as a fallback. This uses branch-length distances in the tree to find the outgroup tip that is most distant from all non-outgroup tips, and roots on that tip. This is a heuristic method that can be used when no clear outgroup clade is available.
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

# Main function to parse arguments, read the tree, root it, and plot it
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
