"""Per-sample gene-expression import and Scrublet doublet removal."""

from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc
import scrublet as scr


CELLRANGER_DIR = Path("data/cellranger")
SAMPLE_METADATA = Path("config/sample_metadata.csv")
OUTPUT_H5AD = Path("results/gex_doublet_removed.h5ad")
QC_SUMMARY = Path("results/gex_doublet_removal_summary.csv")
EXPECTED_DOUBLET_RATE = 0.10


samples = pd.read_csv(SAMPLE_METADATA)
adatas = []
summary = []

for sample_id in samples["sample"]:
    matrix_file = (
        CELLRANGER_DIR
        / sample_id
        / "sample_filtered_feature_bc_matrix.h5"
    )

    adata = sc.read_10x_h5(matrix_file)
    adata.var_names_make_unique()
    adata.obs["sample"] = sample_id
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    scrublet = scr.Scrublet(
        adata.X,
        expected_doublet_rate=EXPECTED_DOUBLET_RATE,
    )
    doublet_score, predicted_doublet = scrublet.scrub_doublets(verbose=False)
    adata.obs["doublet_score"] = doublet_score
    adata.obs["predicted_doublet"] = predicted_doublet

    n_input = adata.n_obs
    adata = adata[~adata.obs["predicted_doublet"]].copy()
    adata.obs_names = sample_id + "_" + adata.obs_names
    adata.layers["counts"] = adata.X.copy()

    summary.append(
        {
            "sample": sample_id,
            "input_cells": n_input,
            "predicted_doublets": n_input - adata.n_obs,
            "retained_cells": adata.n_obs,
        }
    )
    adatas.append(adata)


gex = ad.concat(adatas, join="outer", merge="same", index_unique=None)
OUTPUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
gex.write_h5ad(OUTPUT_H5AD, compression="gzip")
pd.DataFrame(summary).to_csv(QC_SUMMARY, index=False)
