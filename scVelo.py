"""Reusable RNA velocity analysis workflow for multiple Loom files."""

from pathlib import Path
from typing import Optional, Sequence

import loompy
import matplotlib.pyplot as plt
import scvelo as scv


# Configure input files, output locations, and metadata fields for the dataset.
LOOM_FILES = [
    "<path_to_loom_file_1.loom>",
    "<path_to_loom_file_2.loom>",
]
OUTPUT_DIRECTORY = Path("<output_directory>")
MERGED_LOOM_FILE = OUTPUT_DIRECTORY / "merged.loom"
OUTPUT_H5AD_FILE = OUTPUT_DIRECTORY / "velocity_analysis.h5ad"

LOOM_KEY = "Accession"
CELL_ID_DELIMITER = ":"
CELL_ID_SUFFIX_TO_REMOVE = "x"
CELL_BARCODE_SUFFIX = "-1"
CELL_ID_SEPARATOR = "___"
SAMPLE_ORDER: Sequence[str] = ()

GROUP_KEY = "cluster_label"
VELOCITY_MODE = "stochastic"
TOP_GENE_COUNT = 100
CELL_CYCLE_ENABLED = True
PAGA_ENABLED = True
ADDITIONAL_GENE_LIST: Sequence[str] = ()
STARTING_CELL_INDEX: Optional[int] = None


def validate_configuration() -> None:
    """Validate configured input files and analysis options."""
    if not LOOM_FILES:
        raise ValueError("Configure at least one Loom input file.")

    unresolved_files = [path for path in LOOM_FILES if path.startswith("<")]
    if unresolved_files:
        raise ValueError("Replace the Loom file placeholders with valid input paths.")

    missing_files = [path for path in LOOM_FILES if not Path(path).is_file()]
    if missing_files:
        raise FileNotFoundError(f"Loom input files were not found: {missing_files}")

    if VELOCITY_MODE not in {"steady_state", "stochastic", "dynamical"}:
        raise ValueError(
            "VELOCITY_MODE must be 'steady_state', 'stochastic', or 'dynamical'."
        )

    if TOP_GENE_COUNT < 1:
        raise ValueError("TOP_GENE_COUNT must be a positive integer.")

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)


def merge_loom_files() -> None:
    """Combine the configured Loom files into one dataset."""
    loompy.combine(
        LOOM_FILES,
        str(MERGED_LOOM_FILE),
        key=LOOM_KEY,
    )


def normalize_cell_id(
    raw_cell_id: str,
    delimiter: str = CELL_ID_DELIMITER,
    suffix_to_remove: str = CELL_ID_SUFFIX_TO_REMOVE,
    barcode_suffix: str = CELL_BARCODE_SUFFIX,
    separator: str = CELL_ID_SEPARATOR,
) -> str:
    """Convert a sample-prefixed Loom cell ID into a unique normalized ID."""
    sample_id, delimiter_found, barcode = raw_cell_id.partition(delimiter)
    if not delimiter_found or not sample_id or not barcode:
        raise ValueError(
            f"Cell ID {raw_cell_id!r} does not match the configured format."
        )

    if suffix_to_remove and barcode.endswith(suffix_to_remove):
        barcode = barcode[: -len(suffix_to_remove)]

    if barcode_suffix and not barcode.endswith(barcode_suffix):
        barcode = f"{barcode}{barcode_suffix}"

    return f"{barcode}{separator}{sample_id}"


def load_and_prepare_data():
    """Load the merged Loom dataset and normalize its cell identifiers."""
    adata = scv.read(str(MERGED_LOOM_FILE))
    adata.var_names_make_unique()

    normalized_cell_ids = [
        normalize_cell_id(cell_id) for cell_id in adata.obs_names.astype(str)
    ]
    if len(normalized_cell_ids) != len(set(normalized_cell_ids)):
        raise ValueError("Normalized cell IDs are not unique.")

    adata.obs_names = normalized_cell_ids

    if SAMPLE_ORDER:
        sample_ids = adata.obs_names.to_series().str.rsplit(
            CELL_ID_SEPARATOR, n=1
        ).str[-1]
        sample_rank = {sample_id: rank for rank, sample_id in enumerate(SAMPLE_ORDER)}
        ordered_indices = sorted(
            range(len(adata.obs_names)),
            key=lambda index: (
                sample_rank.get(sample_ids.iloc[index], len(sample_rank)),
                index,
            ),
        )
        adata = adata[ordered_indices].copy()

    return adata


def save_plot(filename: str) -> None:
    """Save the current Matplotlib figure in the configured output directory."""
    output_path = OUTPUT_DIRECTORY / filename
    plt.gcf().savefig(output_path, bbox_inches="tight")
    plt.close("all")


def run_velocity_analysis(adata):
    """Run preprocessing, velocity estimation, and graph-based analysis."""
    if GROUP_KEY not in adata.obs:
        raise KeyError(
            f"The configured grouping column {GROUP_KEY!r} is missing from adata.obs."
        )

    scv.pp.filter_and_normalize(adata)
    scv.pp.moments(adata)

    if "X_umap" not in adata.obsm:
        scv.tl.umap(adata)

    if VELOCITY_MODE == "dynamical":
        scv.tl.recover_dynamics(adata)

    scv.tl.velocity(adata, mode=VELOCITY_MODE)
    scv.tl.velocity_graph(adata)
    scv.tl.velocity_confidence(adata)
    scv.tl.velocity_pseudotime(adata)

    if CELL_CYCLE_ENABLED:
        scv.tl.score_genes_cell_cycle(adata)

    if PAGA_ENABLED:
        scv.tl.paga(adata, groups=GROUP_KEY)

    if VELOCITY_MODE == "dynamical":
        scv.tl.latent_time(adata)

    return adata


def export_velocity_plots(adata) -> None:
    """Create and save standard velocity and cell-state visualizations."""
    scv.pl.velocity_embedding_stream(
        adata,
        basis="umap",
        color=GROUP_KEY,
        min_mass=0.1,
        show=False,
    )
    save_plot("velocity_embedding_stream.pdf")

    scv.pl.velocity_embedding_grid(
        adata,
        basis="umap",
        color=GROUP_KEY,
        min_mass=0.1,
        show=False,
    )
    save_plot("velocity_embedding_grid.pdf")

    scv.pl.scatter(
        adata,
        color=("velocity_length", "velocity_confidence"),
        cmap="coolwarm",
        perc=[5, 95],
        show=False,
    )
    save_plot("velocity_confidence.pdf")

    scv.pl.scatter(
        adata,
        color="velocity_pseudotime",
        color_map="gnuplot",
        basis="umap",
        show=False,
    )
    save_plot("velocity_pseudotime.pdf")

    if CELL_CYCLE_ENABLED:
        scv.pl.scatter(
            adata,
            color_gradients=["S_score", "G2M_score"],
            smooth=True,
            perc=[5, 95],
            show=False,
        )
        save_plot("cell_cycle_scores.pdf")

    if PAGA_ENABLED:
        scv.pl.paga(
            adata,
            basis="umap",
            size=50,
            alpha=0.02,
            min_edge_width=3,
            node_size_scale=1.5,
            show=False,
        )
        save_plot("paga_velocity.pdf")

    if VELOCITY_MODE == "dynamical":
        scv.pl.scatter(
            adata,
            color="latent_time",
            color_map="gnuplot",
            basis="umap",
            show=False,
        )
        save_plot("latent_time_umap.pdf")

    if STARTING_CELL_INDEX is not None:
        if not 0 <= STARTING_CELL_INDEX < adata.n_obs:
            raise IndexError("STARTING_CELL_INDEX is outside the available cells.")
        x, y = scv.utils.get_cell_transitions(
            adata,
            basis="umap",
            starting_cell=STARTING_CELL_INDEX,
        )
        ax = scv.pl.velocity_graph(
            adata,
            c="lightgrey",
            edge_width=0.05,
            show=False,
        )
        scv.pl.scatter(
            adata,
            x=x,
            y=y,
            s=120,
            c="ascending",
            cmap="gnuplot",
            ax=ax,
            show=False,
        )
        save_plot("cell_transition_example.pdf")


def export_gene_results(adata) -> None:
    """Export ranked genes and optional gene-specific plots."""
    if VELOCITY_MODE != "dynamical":
        return

    if "fit_likelihood" not in adata.var:
        raise KeyError("The dynamical model did not produce fit_likelihood values.")

    ranked_likelihood = adata.var["fit_likelihood"].dropna().sort_values(
        ascending=False
    )
    ranked_likelihood.to_csv(
        OUTPUT_DIRECTORY / "ranked_genes_by_likelihood.tsv",
        sep="\t",
        header=True,
    )

    top_genes = ranked_likelihood.index[:TOP_GENE_COUNT].tolist()
    if top_genes:
        scv.pl.heatmap(
            adata,
            var_names=top_genes,
            sortby="latent_time",
            col_color=GROUP_KEY,
            n_convolve=100,
            color_map="viridis",
            row_cluster=False,
            show=False,
        )
        save_plot("top_likelihood_genes_heatmap.pdf")

    scv.tl.rank_dynamical_genes(adata, groupby=GROUP_KEY)
    ranked_by_group = scv.get_df(adata, "rank_dynamical_genes/names")
    ranked_by_group.to_csv(
        OUTPUT_DIRECTORY / "ranked_genes_by_group.tsv",
        sep="\t",
        index=True,
    )

    if ADDITIONAL_GENE_LIST:
        present_genes = [
            gene for gene in ADDITIONAL_GENE_LIST if gene in adata.var_names
        ]
        missing_genes = sorted(set(ADDITIONAL_GENE_LIST) - set(present_genes))
        if missing_genes:
            raise ValueError(f"Configured genes are absent from the data: {missing_genes}")

        scv.pl.heatmap(
            adata,
            var_names=present_genes,
            sortby="latent_time",
            col_color=GROUP_KEY,
            n_convolve=100,
            color_map="viridis",
            row_cluster=False,
            show=False,
        )
        save_plot("configured_gene_heatmap.pdf")

        scv.pl.scatter(
            adata,
            var_names=present_genes,
            frameon=True,
            color=GROUP_KEY,
            show=False,
        )
        save_plot("configured_gene_phase_portraits.pdf")

        scv.pl.scatter(
            adata,
            x="latent_time",
            y=present_genes,
            frameon=True,
            color=GROUP_KEY,
            show=False,
        )
        save_plot("configured_gene_expression_dynamics.pdf")

    group_scores = adata.obs.groupby(GROUP_KEY, observed=False)[
        ["velocity_length", "velocity_confidence"]
    ].mean()
    group_scores.to_csv(OUTPUT_DIRECTORY / "velocity_scores_by_group.tsv", sep="\t")


def main() -> None:
    """Run the complete Loom-to-velocity analysis workflow."""
    validate_configuration()
    merge_loom_files()

    scv.settings.set_figure_params("scvelo")
    adata = load_and_prepare_data()
    adata = run_velocity_analysis(adata)

    export_velocity_plots(adata)
    export_gene_results(adata)
    adata.write(str(OUTPUT_H5AD_FILE), compression="gzip")


if __name__ == "__main__":
    main()
