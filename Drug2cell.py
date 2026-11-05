"""Score and analyze drug-target activity in an AnnData object."""

from __future__ import annotations

import argparse
import csv
import logging
import pickle
import re
from pathlib import Path
from typing import Any

import drug2cell as d2c
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

LOGGER = logging.getLogger("drug2cell_analysis")


def parse_arguments() -> argparse.Namespace:
    """Parse input, output, and analysis options."""
    parser = argparse.ArgumentParser(
        description="Score drug-target activity and analyze ranked features."
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Path to the input AnnData file.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        type=Path,
        help="Directory for analysis outputs.",
    )
    parser.add_argument(
        "--groupby",
        required=True,
        nargs="+",
        help="One or more AnnData observation columns for grouping.",
    )
    parser.add_argument(
        "--layer",
        default="data",
        help="Expression layer used for scoring; use 'X' to score adata.X.",
    )
    parser.add_argument(
        "--category",
        action="append",
        default=[],
        help="Optional drug category to include; may be supplied more than once.",
    )
    parser.add_argument(
        "--plot-feature",
        help="Optional drug feature name to plot on the UMAP.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=10,
        help="Number of top ranked features to export per group.",
    )
    parser.add_argument(
        "--chembl-dict",
        type=Path,
        help="Optional trusted pickle containing nested drug-target mappings.",
    )
    return parser.parse_args()


def validate_inputs(args: argparse.Namespace) -> None:
    """Validate paths and options before loading analysis data."""
    if not args.input.is_file():
        raise FileNotFoundError(f"Input AnnData file was not found: {args.input}")
    if args.top_n < 1:
        raise ValueError("--top-n must be a positive integer.")
    if args.chembl_dict is not None and not args.chembl_dict.is_file():
        raise FileNotFoundError(
            f"Drug-target mapping file was not found: {args.chembl_dict}"
        )
    args.output_dir.mkdir(parents=True, exist_ok=True)


def safe_filename(value: str) -> str:
    """Convert a label to a filesystem-safe filename component."""
    normalized = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._")
    return normalized or "group"


def export_drug_target_mapping(
    mapping_file: Path,
    output_file: Path,
) -> None:
    """Export a nested drug-target mapping pickle as a tab-separated table.

    Only load pickle files from trusted sources because pickle data can execute
    arbitrary code during deserialization.
    """
    with mapping_file.open("rb") as input_handle:
        mapping: Any = pickle.load(input_handle)

    if not isinstance(mapping, dict):
        raise TypeError("The drug-target mapping must contain a dictionary.")

    with output_file.open("w", newline="", encoding="utf-8") as output_handle:
        writer = csv.writer(output_handle, delimiter="\t")
        writer.writerow(["class", "drug_id", "drug_name", "genes"])
        for class_name, drugs in mapping.items():
            if not isinstance(drugs, dict):
                raise TypeError(
                    f"Drug entries for class {class_name!r} must be a dictionary."
                )
            for drug_label, genes in drugs.items():
                if not isinstance(drug_label, str) or "|" not in drug_label:
                    raise ValueError(
                        "Each drug label must contain an identifier and name "
                        "separated by '|'."
                    )
                drug_id, drug_name = drug_label.split("|", maxsplit=1)
                if isinstance(genes, str):
                    gene_list = [genes]
                else:
                    gene_list = list(genes)
                writer.writerow(
                    [class_name, drug_id, drug_name, ",".join(map(str, gene_list))]
                )


def get_score_object(adata: sc.AnnData) -> sc.AnnData:
    """Return the drug-score AnnData stored by Drug2Cell."""
    score_object = adata.uns.get("drug2cell")
    if not isinstance(score_object, sc.AnnData):
        raise RuntimeError(
            "Drug2Cell did not add an AnnData score object to adata.uns['drug2cell']."
        )
    return score_object


def matrix_to_frame(score_object: sc.AnnData) -> pd.DataFrame:
    """Convert the drug-score matrix to a labeled DataFrame."""
    matrix = score_object.X
    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    return pd.DataFrame(
        np.asarray(matrix),
        index=score_object.obs_names,
        columns=score_object.var_names,
    )


def export_ranked_features(
    score_object: sc.AnnData,
    groupby: str,
    output_dir: Path,
    top_n: int,
    suffix: str = "",
) -> None:
    """Rank score features by group and export full and top-ranked tables."""
    if groupby not in score_object.obs.columns:
        raise KeyError(f"Grouping column {groupby!r} is missing from score metadata.")

    sc.tl.rank_genes_groups(
        score_object,
        groupby=groupby,
        method="wilcoxon",
    )
    ranked = sc.get.rank_genes_groups_df(score_object, group=None)
    if ranked.empty:
        raise RuntimeError(f"No ranked features were produced for {groupby!r}.")

    suffix_part = f"_{safe_filename(suffix)}" if suffix else ""
    group_part = safe_filename(groupby)
    ranked.to_csv(
        output_dir / f"ranked_features_{group_part}{suffix_part}.tsv",
        sep="\t",
        index=False,
    )
    ranked.groupby("group", sort=False).head(top_n).to_csv(
        output_dir / f"top_{top_n}_features_{group_part}{suffix_part}.tsv",
        sep="\t",
        index=False,
    )


def plot_umap(
    score_object: sc.AnnData,
    output_dir: Path,
    feature: str | None,
) -> None:
    """Save a UMAP colored by the configured feature or available metadata."""
    if "X_umap" not in score_object.obsm:
        LOGGER.warning("No UMAP coordinates found; skipping UMAP plot.")
        return

    color: str | None = feature
    if color is not None and color not in score_object.var_names:
        raise KeyError(f"Requested UMAP feature {color!r} is missing.")
    if color is None:
        for candidate in ("subcluster", "cell_type", "group"):
            if candidate in score_object.obs.columns:
                color = candidate
                break
    if color is None:
        LOGGER.warning("No default UMAP color column is available; skipping plot.")
        return

    sc.pl.umap(score_object, color=color, show=False)
    plt.savefig(output_dir / f"umap_{safe_filename(color)}.pdf", bbox_inches="tight")
    plt.close()


def plot_category_scores(
    score_object: sc.AnnData,
    groupby_columns: list[str],
    categories: list[str],
    output_dir: Path,
) -> sc.AnnData | None:
    """Plot and return a score object restricted to selected categories."""
    if not categories:
        return None

    plot_args = d2c.util.prepare_plot_args(
        score_object,
        categories=categories,
    )
    feature_names = plot_args.get("var_names", [])
    if not feature_names:
        raise ValueError("No drug features matched the requested categories.")
    category_scores = score_object[:, feature_names].copy()

    for groupby in groupby_columns:
        if groupby not in category_scores.obs.columns:
            raise KeyError(
                f"Grouping column {groupby!r} is missing from score metadata."
            )
        plot = sc.pl.dotplot(
            category_scores,
            var_names=feature_names,
            groupby=groupby,
            swap_axes=True,
            show=False,
        )
        plot.savefig(
            output_dir / f"category_scores_{safe_filename(groupby)}.pdf",
            bbox_inches="tight",
        )
        plt.close("all")

    return category_scores


def main() -> None:
    """Run drug scoring, visualizations, and differential ranking."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(message)s",
    )
    args = parse_arguments()
    validate_inputs(args)

    if args.chembl_dict is not None:
        export_drug_target_mapping(
            args.chembl_dict,
            args.output_dir / "drug_target_mapping.tsv",
        )

    LOGGER.info("Loading AnnData from %s", args.input)
    adata = sc.read_h5ad(args.input)
    missing_groups = [column for column in args.groupby if column not in adata.obs]
    if missing_groups:
        raise KeyError(
            "Grouping columns are missing from AnnData observations: "
            + ", ".join(missing_groups)
        )

    if args.layer == "X":
        d2c.score(adata, use_raw=False)
    else:
        if args.layer not in adata.layers:
            raise KeyError(
                f"Expression layer {args.layer!r} is missing. "
                "Use --layer X to score adata.X."
            )
        d2c.score(adata, use_raw=False, layer=args.layer)

    score_object = get_score_object(adata)
    score_object.write_h5ad(args.output_dir / "drug_scores.h5ad", compression="gzip")
    matrix_to_frame(score_object).to_csv(
        args.output_dir / "drug_score_matrix.tsv",
        sep="\t",
        index=True,
    )

    plot_umap(score_object, args.output_dir, args.plot_feature)
    for groupby in args.groupby:
        export_ranked_features(
            score_object,
            groupby,
            args.output_dir,
            args.top_n,
        )

    category_scores = plot_category_scores(
        score_object,
        args.groupby,
        args.category,
        args.output_dir,
    )
    if category_scores is not None:
        for groupby in args.groupby:
            export_ranked_features(
                category_scores,
                groupby,
                args.output_dir,
                args.top_n,
                suffix="selected_categories",
            )

    adata.write_h5ad(args.output_dir / "drug2cell_annotated.h5ad", compression="gzip")
    LOGGER.info("Analysis outputs were written to %s", args.output_dir)


if __name__ == "__main__":
    main()
