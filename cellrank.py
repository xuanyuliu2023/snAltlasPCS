"""Reusable RNA velocity and lineage analysis with Scanpy, scVelo, and CellRank."""

from pathlib import Path
from typing import Optional, Sequence

import cellrank as cr
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import scvelo as scv


# Configure input files, metadata keys, model options, and output locations.
INPUT_H5AD = Path("<input_annotated_data>.h5ad")
OUTPUT_DIRECTORY = Path("<output_directory>")

GROUP_KEY = "cluster_label"
MACROSTATE_KEY = "macrostates_fwd"
LAYOUT_KEY = "umap"
EXPRESSION_LAYER = "normalized_expression"
PLOT_EXPRESSION_LAYER: Optional[str] = None

RUN_VELOCITY_PREPROCESSING = True
RUN_DYNAMICAL_MODEL = True
RUN_CELLRANK_ANALYSIS = True
RUN_DRIVER_ANALYSIS = False
RUN_TREND_ANALYSIS = False

MIN_SHARED_COUNTS = 20
N_TOP_GENES = 2000
N_PCS = 30
N_NEIGHBORS = 30
VELOCITY_JOBS = 1
CONNECTIVITY_WEIGHT = 0.2
VELOCITY_WEIGHT = 0.8
STATE_COUNTS = (5, 12)
MIN_CELLS_PER_STATE = 50

LINEAGE_NAME = "<lineage_name>"
DRIVER_CLUSTERS: Sequence[str] = ()
INITIAL_STATES: Sequence[str] = ()
TERMINAL_STATES: Sequence[str] = ()
GENES_TO_PLOT: Sequence[str] = ()
TOP_DRIVER_GENE_COUNT = 2000
TOP_TF_COUNT = 50

GENE_ANNOTATION_FILE: Optional[Path] = None
TRANSCRIPTION_FACTOR_FILE: Optional[Path] = None

GAM_KNOTS = 6
GAM_SMOOTHING_PENALTY = 10.0
TREND_CLUSTER_RESOLUTION = 0.4
RANDOM_STATE = 0


def validate_configuration() -> None:
    """Validate paths and options before starting the analysis."""
    if str(INPUT_H5AD).startswith("<"):
        raise ValueError("Set INPUT_H5AD to an existing annotated data file.")
    if not INPUT_H5AD.is_file():
        raise FileNotFoundError(f"Annotated data file not found: {INPUT_H5AD}")
    if not 0 <= CONNECTIVITY_WEIGHT <= 1:
        raise ValueError("CONNECTIVITY_WEIGHT must be between zero and one.")
    if not 0 <= VELOCITY_WEIGHT <= 1:
        raise ValueError("VELOCITY_WEIGHT must be between zero and one.")
    if abs(CONNECTIVITY_WEIGHT + VELOCITY_WEIGHT - 1.0) > 1e-8:
        raise ValueError("Kernel weights must sum to one.")
    if TOP_DRIVER_GENE_COUNT < 1 or TOP_TF_COUNT < 1:
        raise ValueError("Driver gene limits must be positive integers.")
    if RUN_DRIVER_ANALYSIS and LINEAGE_NAME.startswith("<"):
        raise ValueError("Set LINEAGE_NAME before enabling driver analysis.")
    if RUN_DRIVER_ANALYSIS and not DRIVER_CLUSTERS:
        raise ValueError("Configure DRIVER_CLUSTERS before enabling driver analysis.")
    if RUN_TREND_ANALYSIS:
        if LINEAGE_NAME.startswith("<"):
            raise ValueError("Set LINEAGE_NAME before enabling trend analysis.")
        if not DRIVER_CLUSTERS:
            raise ValueError("Configure DRIVER_CLUSTERS before enabling trend analysis.")
    if GENE_ANNOTATION_FILE is not None and not GENE_ANNOTATION_FILE.is_file():
        raise FileNotFoundError(f"Gene annotation file not found: {GENE_ANNOTATION_FILE}")
    if TRANSCRIPTION_FACTOR_FILE is not None and not TRANSCRIPTION_FACTOR_FILE.is_file():
        raise FileNotFoundError(
            f"Transcription factor file not found: {TRANSCRIPTION_FACTOR_FILE}"
        )
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)


def configure_plotting() -> None:
    """Set consistent plotting defaults and vector-friendly font output."""
    scv.settings.verbosity = 3
    scv.settings.set_figure_params("scvelo")
    cr.settings.verbosity = 2
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    mpl.rcParams["savefig.dpi"] = 300


def save_current_figure(file_name: str) -> None:
    """Save and close the current Matplotlib figure."""
    plt.gcf().savefig(OUTPUT_DIRECTORY / file_name, bbox_inches="tight")
    plt.close("all")


def load_data():
    """Read the annotated data and validate required metadata."""
    adata = sc.read_h5ad(INPUT_H5AD)
    if GROUP_KEY not in adata.obs:
        raise KeyError(f"Required observation column {GROUP_KEY!r} is missing.")
    return adata


def preprocess_and_estimate_velocity(adata):
    """Optionally preprocess the data and estimate RNA velocity."""
    if not RUN_VELOCITY_PREPROCESSING:
        return adata

    scv.pp.filter_and_normalize(
        adata,
        min_shared_counts=MIN_SHARED_COUNTS,
        n_top_genes=N_TOP_GENES,
        subset_highly_variable=False,
    )
    sc.tl.pca(adata)
    sc.pp.neighbors(
        adata,
        n_pcs=N_PCS,
        n_neighbors=N_NEIGHBORS,
        random_state=RANDOM_STATE,
    )
    scv.pp.moments(adata, n_pcs=None, n_neighbors=None)

    if RUN_DYNAMICAL_MODEL:
        scv.tl.recover_dynamics(adata, n_jobs=VELOCITY_JOBS)
        scv.tl.velocity(adata, mode="dynamical")
    else:
        scv.tl.velocity(adata, mode="stochastic")

    return adata


def build_transition_kernel(adata):
    """Combine expression-connectivity and RNA-velocity transition kernels."""
    connectivity_kernel = cr.kernels.ConnectivityKernel(adata)
    connectivity_kernel.compute_transition_matrix()

    velocity_kernel = cr.kernels.VelocityKernel(adata)
    velocity_kernel.compute_transition_matrix()

    combined_kernel = (
        CONNECTIVITY_WEIGHT * connectivity_kernel
        + VELOCITY_WEIGHT * velocity_kernel
    )
    combined_kernel.write_to_adata()
    return combined_kernel


def fit_lineage_estimator(transition_kernel):
    """Fit macrostates and infer initial, terminal, and fate states."""
    estimator = cr.estimators.GPCCA(transition_kernel)
    estimator.fit(
        cluster_key=GROUP_KEY,
        n_states=list(STATE_COUNTS),
        n_cells=MIN_CELLS_PER_STATE,
    )
    estimator.predict_terminal_states(allow_overlap=True)
    estimator.predict_initial_states(allow_overlap=True)

    if INITIAL_STATES:
        estimator.set_initial_states(states=list(INITIAL_STATES), allow_overlap=True)
    if TERMINAL_STATES:
        estimator.set_terminal_states(states=list(TERMINAL_STATES))

    estimator.compute_fate_probabilities()
    return estimator


def export_core_plots(adata, transition_kernel, estimator) -> None:
    """Create standard transition, state, and fate probability plots."""
    transition_kernel.plot_projection()
    save_current_figure("transition_kernel_projection.pdf")

    estimator.plot_macrostates(
        which="all",
        discrete=True,
        legend_loc="right",
        basis=LAYOUT_KEY,
    )
    save_current_figure("macrostates.pdf")

    estimator.plot_macrostates(
        which="terminal",
        discrete=False,
        basis=LAYOUT_KEY,
    )
    save_current_figure("terminal_states.pdf")

    estimator.plot_macrostates(
        which="initial",
        discrete=False,
        basis=LAYOUT_KEY,
    )
    save_current_figure("initial_states.pdf")

    estimator.plot_fate_probabilities(same_plot=False, basis=LAYOUT_KEY)
    save_current_figure("fate_probabilities.pdf")

    cr.pl.circular_projection(
        adata,
        keys=[GROUP_KEY],
        legend_loc="right",
    )
    save_current_figure("circular_projection.pdf")


def load_gene_annotation() -> dict[str, str]:
    """Load optional gene-to-biotype annotations from a two-column TSV file."""
    if GENE_ANNOTATION_FILE is None:
        return {}

    annotation = pd.read_csv(
        GENE_ANNOTATION_FILE,
        sep="\t",
        header=None,
        dtype=str,
    )
    if annotation.shape[1] < 2:
        raise ValueError("Gene annotation input must contain at least two columns.")
    return dict(zip(annotation.iloc[:, 0], annotation.iloc[:, 1]))


def load_transcription_factors() -> set[str]:
    """Load optional transcription factor identifiers from a TSV file."""
    if TRANSCRIPTION_FACTOR_FILE is None:
        return set()

    factors = pd.read_csv(
        TRANSCRIPTION_FACTOR_FILE,
        sep="\t",
        header=None,
        dtype=str,
    )
    if factors.shape[1] < 1:
        raise ValueError("Transcription factor input must contain at least one column.")
    return set(factors.iloc[:, 0])


def compute_and_export_drivers(adata, estimator):
    """Compute lineage driver genes and export ranked annotations."""
    drivers = estimator.compute_lineage_drivers(
        lineages=LINEAGE_NAME,
        clusters=list(DRIVER_CLUSTERS),
        cluster_key=MACROSTATE_KEY,
        layer=EXPRESSION_LAYER,
        use_raw=False,
    )

    fate_probability_key = f"fate_probability_{LINEAGE_NAME}"
    adata.obs[fate_probability_key] = (
        estimator.fate_probabilities[LINEAGE_NAME].X.flatten()
    )

    gene_types = load_gene_annotation()
    transcription_factors = load_transcription_factors()
    drivers["gene_type"] = [
        gene_types.get(gene, "unknown") for gene in drivers.index
    ]
    drivers["is_transcription_factor"] = [
        gene in transcription_factors for gene in drivers.index
    ]

    q_value_column = f"{LINEAGE_NAME}_qval"
    correlation_column = f"{LINEAGE_NAME}_corr"
    if q_value_column not in drivers.columns:
        raise KeyError(f"Expected driver statistics column {q_value_column!r} is missing.")
    if correlation_column not in drivers.columns:
        raise KeyError(
            f"Expected driver statistics column {correlation_column!r} is missing."
        )

    drivers = drivers.loc[drivers[q_value_column] < 0.05].sort_values(
        by=q_value_column,
        ascending=True,
    )
    drivers.to_csv(
        OUTPUT_DIRECTORY / "lineage_driver_genes.tsv",
        sep="\t",
        index=True,
    )

    protein_coding = drivers.loc[drivers["gene_type"] == "protein_coding"]
    top_protein_coding = protein_coding.sort_values(
        by=correlation_column,
        ascending=False,
    ).head(TOP_DRIVER_GENE_COUNT)
    top_protein_coding.to_csv(
        OUTPUT_DIRECTORY / "top_protein_coding_driver_genes.tsv",
        sep="\t",
        index=True,
    )

    top_factors = drivers.loc[
        (drivers["gene_type"] == "protein_coding")
        & drivers["is_transcription_factor"]
    ].sort_values(by=correlation_column, ascending=False).head(TOP_TF_COUNT)
    top_factors.to_csv(
        OUTPUT_DIRECTORY / "top_transcription_factor_drivers.tsv",
        sep="\t",
        index=True,
    )
    return drivers


def plot_driver_genes(adata, drivers) -> None:
    """Plot expression of the top lineage drivers in the configured embedding."""
    selected_genes = list(drivers.index[: min(8, len(drivers))])
    if not selected_genes:
        return

    sc.pl.embedding(
        adata,
        basis=LAYOUT_KEY,
        color=selected_genes,
        ncols=3,
        frameon=False,
        layer=PLOT_EXPRESSION_LAYER,
        use_raw=False,
        show=False,
    )
    save_current_figure("top_driver_gene_expression.pdf")


def run_gene_trend_analysis(adata, drivers) -> None:
    """Fit a GAM and export gene trend visualizations and cluster summaries."""
    model = cr.models.GAMR(
        adata,
        n_knots=GAM_KNOTS,
        smoothing_penalty=GAM_SMOOTHING_PENALTY,
    )
    selected_genes = list(GENES_TO_PLOT) or list(drivers.index[:10])
    if not selected_genes:
        raise ValueError("Configure genes or run driver analysis before trend plotting.")

    cr.pl.gene_trends(
        adata,
        model=model,
        data_key=EXPRESSION_LAYER,
        genes=selected_genes,
        same_plot=True,
        ncols=5,
        time_key="palantir_pseudotime",
        hide_cells=True,
        weight_threshold=(1e-3, 1e-3),
    )
    save_current_figure("gene_trends.pdf")

    cr.pl.heatmap(
        adata,
        model=model,
        lineages=LINEAGE_NAME,
        cluster_key=MACROSTATE_KEY,
        show_fate_probabilities=True,
        data_key=EXPRESSION_LAYER,
        genes=list(drivers.index),
        time_key="palantir_pseudotime",
        show_all_genes=False,
        weight_threshold=(1e-3, 1e-3),
        n_jobs=1,
    )
    save_current_figure("lineage_driver_heatmap.pdf")

    trend_key = f"lineage_{LINEAGE_NAME}_trend"
    adata.uns.pop(trend_key, None)
    cr.pl.cluster_trends(
        adata,
        model=model,
        lineage=LINEAGE_NAME,
        data_key=EXPRESSION_LAYER,
        genes=list(drivers.index),
        time_key="palantir_pseudotime",
        weight_threshold=(1e-3, 1e-3),
        n_jobs=1,
        random_state=RANDOM_STATE,
        clustering_kwargs={
            "resolution": TREND_CLUSTER_RESOLUTION,
            "random_state": RANDOM_STATE,
        },
        neighbors_kwargs={"random_state": RANDOM_STATE},
    )

    if trend_key not in adata.uns:
        raise KeyError(f"CellRank did not create the expected trend result {trend_key!r}.")
    adata.uns[trend_key].obs.to_csv(
        OUTPUT_DIRECTORY / "lineage_cluster_trends.tsv",
        sep="\t",
        index=True,
    )


def main() -> None:
    """Run the configured CellRank workflow and save the annotated dataset."""
    validate_configuration()
    configure_plotting()
    adata = load_data()
    adata = preprocess_and_estimate_velocity(adata)

    if RUN_CELLRANK_ANALYSIS:
        transition_kernel = build_transition_kernel(adata)
        estimator = fit_lineage_estimator(transition_kernel)
        export_core_plots(adata, transition_kernel, estimator)
    else:
        transition_kernel = None
        estimator = None

    drivers = None
    if RUN_DRIVER_ANALYSIS:
        if estimator is None:
            raise RuntimeError("Driver analysis requires CellRank state estimation.")
        drivers = compute_and_export_drivers(adata, estimator)
        plot_driver_genes(adata, drivers)

    if RUN_TREND_ANALYSIS:
        if estimator is None:
            raise RuntimeError("Trend analysis requires CellRank state estimation.")
        if drivers is None:
            drivers = estimator.compute_lineage_drivers(
                lineages=LINEAGE_NAME,
                clusters=list(DRIVER_CLUSTERS),
                cluster_key=MACROSTATE_KEY,
                layer=EXPRESSION_LAYER,
                use_raw=False,
            )
        run_gene_trend_analysis(adata, drivers)

    if transition_kernel is not None:
        transition_kernel.write_to_adata()
    output_file = OUTPUT_DIRECTORY / "cellrank_analysis.h5ad"
    adata.write_h5ad(output_file, compression="gzip")


if __name__ == "__main__":
    main()
