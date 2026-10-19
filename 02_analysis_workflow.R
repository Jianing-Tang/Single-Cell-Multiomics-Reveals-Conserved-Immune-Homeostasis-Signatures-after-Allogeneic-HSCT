#!/usr/bin/env Rscript

# Core scRNA-seq and scTCR-seq analysis workflow.
# Run the sections sequentially and complete both annotation checkpoints.

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(harmony)
  library(zellkonverter)
  library(scRepertoire)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(stringdist)
  library(ggplot2)
})

set.seed(42)
dir.create("results", showWarnings = FALSE)

sample_metadata <- read.csv(
  "config/sample_metadata.csv",
  stringsAsFactors = FALSE
)


# =============================================================================
# Part 1. PBMC atlas
# =============================================================================

sce <- readH5AD("results/gex_doublet_removed.h5ad")
mainExpName(sce) <- "RNA"
pbmc <- as.Seurat(sce, counts = "counts", data = NULL)
DefaultAssay(pbmc) <- "RNA"

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(
  pbmc,
  subset = nFeature_RNA >= 500 &
    nFeature_RNA <= 4000 &
    nCount_RNA < 20000 &
    percent.mt < 15
)

pbmc_metadata <- pbmc[[]] |>
  rownames_to_column("cell_barcode") |>
  left_join(sample_metadata, by = "sample") |>
  column_to_rownames("cell_barcode")
pbmc@meta.data <- pbmc_metadata[colnames(pbmc), , drop = FALSE]

pbmc <- NormalizeData(pbmc, scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 3000)
pbmc <- ScaleData(pbmc, features = VariableFeatures(pbmc))
pbmc <- RunPCA(pbmc, npcs = 50, seed.use = 42)
pbmc <- RunHarmony(pbmc, group.by.vars = "sample")
pbmc <- FindNeighbors(pbmc, reduction = "harmony", dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.8)
pbmc <- RunUMAP(
  pbmc,
  reduction = "harmony",
  dims = 1:30,
  n.neighbors = 30,
  min.dist = 0.3,
  metric = "cosine",
  reduction.name = "umap_allcell",
  seed.use = 42
)

pbmc_markers <- FindAllMarkers(
  pbmc,
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25
)
write.csv(pbmc_markers, "results/pbmc_cluster_markers.csv", row.names = FALSE)

pbmc_marker_panel <- c(
  "CD3D", "CD3E", "CD79A", "MS4A1", "CD14", "LYZ", "FCGR3A",
  "GNLY", "NKG7", "CSF3R", "JCHAIN", "MZB1", "TPSAB1", "PPBP"
)

pdf("results/pbmc_marker_dotplot.pdf", width = 11, height = 6)
print(DotPlot(pbmc, features = pbmc_marker_panel) + RotatedAxis())
dev.off()

# Annotation checkpoint 1. Cell-level labels were curated after inspection of
# canonical markers and the
# exported marker table. The annotation file contains: cell_barcode,
# cellanno_broad.
pbmc_annotations <- read.csv(
  "config/pbmc_cell_annotations.csv",
  stringsAsFactors = FALSE
)
pbmc$cellanno_broad <- pbmc_annotations$cellanno_broad[
  match(colnames(pbmc), pbmc_annotations$cell_barcode)
]

pdf("results/pbmc_umap.pdf", width = 8, height = 6)
print(
  DimPlot(
    pbmc,
    reduction = "umap_allcell",
    group.by = "cellanno_broad",
    label = TRUE,
    repel = TRUE
  ) + NoLegend()
)
dev.off()


# =============================================================================
# Part 2. T-cell landscape
# =============================================================================

Idents(pbmc) <- "cellanno_broad"
tcell <- subset(pbmc, idents = "T")
DefaultAssay(tcell) <- "RNA"
tcell[["RNA"]] <- split(tcell[["RNA"]], f = tcell$sample)

tcell <- NormalizeData(tcell, scale.factor = 10000)
tcell <- FindVariableFeatures(tcell, selection.method = "vst", nfeatures = 2000)
tcell <- ScaleData(tcell, features = VariableFeatures(tcell))
tcell <- RunPCA(tcell, npcs = 50, seed.use = 42)
tcell <- IntegrateLayers(
  object = tcell,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE
)
tcell[["RNA"]] <- JoinLayers(tcell[["RNA"]])
tcell <- FindNeighbors(tcell, reduction = "integrated.cca", dims = 1:30)
tcell <- FindClusters(tcell, resolution = 0.8)
tcell <- RunUMAP(
  tcell,
  reduction = "integrated.cca",
  dims = 1:30,
  n.neighbors = 30,
  min.dist = 0.3,
  metric = "cosine",
  reduction.name = "umap_cca_T",
  seed.use = 42
)

tcell_markers <- FindAllMarkers(
  tcell,
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25
)
write.csv(tcell_markers, "results/tcell_cluster_markers.csv", row.names = FALSE)

tcell_marker_panel <- c(
  "CCR7", "LTB", "MAL", "IL7R", "GZMK", "CCL5", "NKG7", "GNLY",
  "GZMB", "PRF1", "ZNF683", "HLA-DRA", "FOXP3", "IL32", "TRDC"
)

pdf("results/tcell_marker_dotplot.pdf", width = 12, height = 7)
print(DotPlot(tcell, features = tcell_marker_panel) + RotatedAxis())
dev.off()

# Annotation checkpoint 2. The annotation file contains: cell_barcode,
# T_cell_anno_V1.
tcell_annotations <- read.csv(
  "config/tcell_cell_annotations.csv",
  stringsAsFactors = FALSE
)
tcell$T_cell_anno_V1 <- tcell_annotations$T_cell_anno_V1[
  match(colnames(tcell), tcell_annotations$cell_barcode)
]

pbmc$cellanno_T_detailed <- pbmc$cellanno_broad
pbmc$cellanno_T_detailed[match(colnames(tcell), colnames(pbmc))] <-
  tcell$T_cell_anno_V1

pdf("results/tcell_umap.pdf", width = 8, height = 6)
print(
  DimPlot(
    tcell,
    reduction = "umap_cca_T",
    group.by = "T_cell_anno_V1",
    label = TRUE,
    repel = TRUE
  ) + NoLegend()
)
dev.off()

saveRDS(pbmc, "results/pbmc_harmony.rds")
saveRDS(tcell, "results/tcell_cca.rds")


# =============================================================================
# Part 3. scTCR-seq integration and Consistent-clonotype identification
# =============================================================================

vdj_samples <- sample_metadata |>
  filter(has_vdj) |>
  pull(sample)

contig_list <- lapply(
  file.path(
    "data/cellranger",
    vdj_samples,
    "filtered_contig_annotations.csv"
  ),
  read.csv
)
names(contig_list) <- vdj_samples

combined_tcr <- combineTCR(
  contig_list,
  samples = vdj_samples,
  removeNA = FALSE,
  removeMulti = FALSE,
  filterMulti = FALSE
)

tcell <- combineExpression(
  combined_tcr,
  tcell,
  cloneCall = "gene",
  group.by = "sample",
  proportion = FALSE
)

conventional_abt <- subset(
  tcell,
  subset = !T_cell_anno_V1 %in% c("MAIT", "NKT", "gdT")
)

clone_metadata <- conventional_abt[[]] |>
  rownames_to_column("cell_barcode") |>
  mutate(
    TRA_gene = str_split_fixed(CTgene, "_", 2)[, 1],
    TRB_gene = str_split_fixed(CTgene, "_", 2)[, 2],
    TRA_gene = na_if(TRA_gene, "NA"),
    TRB_gene = na_if(TRB_gene, "NA"),
    alpha_strict = str_split_fixed(CTstrict, "_", 2)[, 1],
    beta_strict = str_split_fixed(CTstrict, "_", 2)[, 2],
    valid_CTstrict =
      str_count(CTstrict, fixed("_")) == 1 &
      str_detect(alpha_strict, fixed(";")) &
      str_detect(beta_strict, fixed(";")) &
      !is.na(TRA_gene) & TRA_gene != "" &
      !is.na(TRB_gene) & TRB_gene != ""
  ) |>
  filter(valid_CTstrict)

consistent_clonotypes <- clone_metadata |>
  distinct(pair, CTstrict, origin) |>
  count(pair, CTstrict, name = "n_origins") |>
  filter(n_origins == 2) |>
  select(pair, CTstrict)

clone_metadata <- clone_metadata |>
  left_join(
    consistent_clonotypes |> mutate(consistent = TRUE),
    by = c("pair", "CTstrict")
  ) |>
  mutate(
    consistent = replace_na(consistent, FALSE),
    clone_status = case_when(
      origin == "Donor" & consistent ~ "Donor_consistent",
      origin == "Donor" ~ "Donor_unique",
      origin == "Recipient" & consistent ~ "Recipient_consistent",
      origin == "Recipient" ~ "Recipient_unique"
    )
  )

clone_frequency <- clone_metadata |>
  count(sample, CTstrict, name = "clonal_frequency") |>
  mutate(
    clone_size = cut(
      clonal_frequency,
      breaks = c(0, 1, 5, 20, 100, Inf),
      labels = c("Single", "Small", "Medium", "Large", "Hyperexpanded")
    )
  )

clone_metadata <- clone_metadata |>
  left_join(clone_frequency, by = c("sample", "CTstrict"))

clone_annotations <- clone_metadata |>
  select(cell_barcode, clone_status, clonal_frequency, clone_size) |>
  column_to_rownames("cell_barcode")

conventional_abt <- AddMetaData(conventional_abt, clone_annotations)

write.table(
  consistent_clonotypes,
  "results/consistent_CTstrict_clonotypes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
saveRDS(conventional_abt, "results/conventional_abT_CTstrict.rds")


# =============================================================================
# Part 4. Predicted antigen-reactivity annotation with VDJdb
# =============================================================================

clonotypes <- clone_metadata |>
  filter(!is.na(CTaa)) |>
  mutate(
    alpha_aa = str_remove(CTaa, "_.*$"),
    beta_aa = str_remove(CTaa, "^.*_")
  ) |>
  distinct(pair, CTstrict, CTaa, TRA_gene, TRB_gene, alpha_aa, beta_aa)

trb_query <- clonotypes |>
  transmute(
    pair,
    CTstrict,
    TRB_CDR3 = beta_aa,
    TRB_call = TRB_gene
  ) |>
  separate_longer_delim(c(TRB_CDR3, TRB_call), delim = ";") |>
  mutate(
    TRBV = str_remove(str_extract(TRB_call, "TRBV[^.]+"), "\\*.*$"),
    TRBJ = str_remove(str_extract(TRB_call, "TRBJ[^.]+"), "\\*.*$")
  ) |>
  filter(
    !is.na(TRB_CDR3), TRB_CDR3 != "",
    !is.na(TRBV), TRBV != "",
    !is.na(TRBJ), TRBJ != ""
  ) |>
  distinct()

vdjdb <- read.delim(
  "data/reference/vdjdb_2026-06-03.txt",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = ""
) |>
  filter(
    species == "HomoSapiens",
    gene == "TRB",
    vdjdb.score >= 1,
    !is.na(cdr3),
    !is.na(v.segm),
    !is.na(j.segm)
  ) |>
  transmute(
    TRB_CDR3 = cdr3,
    TRBV = str_remove(v.segm, "\\*.*$"),
    TRBJ = str_remove(j.segm, "\\*.*$"),
    antigen_species = antigen.species,
    antigen_epitope = antigen.epitope,
    antigen_gene = antigen.gene,
    reported_MHC = mhc.a,
    vdjdb_score = vdjdb.score
  )

exact_hits <- trb_query |>
  inner_join(
    vdjdb,
    by = c("TRB_CDR3", "TRBV", "TRBJ"),
    relationship = "many-to-many"
  ) |>
  mutate(match_rule = "Exact CDR3beta + V/J")

hamming1_hits <- trb_query |>
  inner_join(
    vdjdb,
    by = c("TRBV", "TRBJ"),
    suffix = c("_query", "_vdjdb"),
    relationship = "many-to-many"
  ) |>
  filter(nchar(TRB_CDR3_query) == nchar(TRB_CDR3_vdjdb)) |>
  mutate(
    amino_acid_substitutions = stringdist(
      TRB_CDR3_query,
      TRB_CDR3_vdjdb,
      method = "hamming"
    )
  ) |>
  filter(amino_acid_substitutions <= 1) |>
  mutate(match_rule = "CDR3beta Hamming distance <= 1 + V/J")

collapse_hits <- function(hits, prefix) {
  hits |>
    group_by(pair, CTstrict) |>
    summarise(
      predicted = TRUE,
      antigen_species = paste(sort(unique(antigen_species)), collapse = ";"),
      antigen_epitopes = paste(sort(unique(antigen_epitope)), collapse = ";"),
      antigen_genes = paste(sort(unique(antigen_gene)), collapse = ";"),
      reported_MHC_not_filtered = paste(
        sort(unique(reported_MHC)),
        collapse = ";"
      ),
      max_vdjdb_score = max(vdjdb_score),
      .groups = "drop"
    ) |>
    rename_with(~ paste0(prefix, "_", .x), -c(pair, CTstrict))
}

exact_annotation <- collapse_hits(exact_hits, "VDJdb_exact")
hamming1_annotation <- collapse_hits(hamming1_hits, "VDJdb_hamming1")

vdjdb_metadata <- conventional_abt[[]] |>
  rownames_to_column("cell_barcode") |>
  left_join(exact_annotation, by = c("pair", "CTstrict")) |>
  left_join(hamming1_annotation, by = c("pair", "CTstrict")) |>
  mutate(
    VDJdb_exact_predicted = replace_na(VDJdb_exact_predicted, FALSE),
    VDJdb_hamming1_predicted = replace_na(
      VDJdb_hamming1_predicted,
      FALSE
    )
  ) |>
  column_to_rownames("cell_barcode")

conventional_abt@meta.data <- vdjdb_metadata[
  colnames(conventional_abt),
  ,
  drop = FALSE
]

write.table(
  exact_hits,
  "results/VDJdb_exact_TRB_hits.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  hamming1_hits,
  "results/VDJdb_hamming1_TRB_hits.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
saveRDS(
  conventional_abt,
  "results/conventional_abT_CTstrict_VDJdb.rds"
)
