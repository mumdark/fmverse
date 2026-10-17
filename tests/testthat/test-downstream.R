# Numerical-kernel regressions; public S4 APIs are tested in test-fmobj.R.
testthat::skip_if_not_installed("Seurat")
testthat::skip_if_not_installed("SeuratObject")
testthat::skip_if_not_installed("ggplot2")

make_downstream_result <- function(seed = 42L, n_cells = 48L) {
  set.seed(seed)
  n_genes <- 16L
  labels <- rep(c("A", "B"), each = n_cells / 2L)
  batch <- rep(rep(c("b1", "b2"), each = n_cells / 4L), 2L)
  counts <- matrix(stats::rpois(n_genes * n_cells, lambda = 2), nrow = n_genes,
                   dimnames = list(paste0("gene", seq_len(n_genes)), paste0("cell", seq_len(n_cells))))
  counts[1:3, labels == "B"] <- counts[1:3, labels == "B"] + 4L
  object <- Seurat::CreateSeuratObject(counts = counts, meta.data = data.frame(
    cell_type = labels, batch = batch, treatment = ifelse(labels == "B", "drug_x", "control"),
    dose = ifelse(labels == "B", 1, 0), row.names = colnames(counts)
  ))
  embedding <- cbind(
    ifelse(labels == "B", 2, -2) + stats::rnorm(n_cells, sd = 0.25),
    stats::rnorm(n_cells), matrix(stats::rnorm(n_cells * 4L), ncol = 4L)
  )
  rownames(embedding) <- colnames(counts)
  colnames(embedding) <- paste0("scGPT_", seq_len(ncol(embedding)))
  object[["scGPT"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding, key = "scGPT_", assay = "RNA"
  )
  vocab <- c("<pad>" = 0L, "<cls>" = 1L, gene1 = 2L, gene2 = 3L, gene3 = 4L,
             gene4 = 5L, gene5 = 6L, gene6 = 7L)
  gene_ids <- matrix(rep(c(1L, 2L, 3L, 4L, 5L), n_cells), nrow = n_cells, byrow = TRUE)
  values <- matrix(stats::runif(n_cells * 5L, 1, 50), nrow = n_cells)
  pad_mask <- matrix(FALSE, nrow = n_cells, ncol = 5L)
  structure(list(
    model_name = "scGPT", embeddings = embedding, seurat = object,
    model = NULL, metadata = list(checkpoint_name = "whole-human", device = "cpu"),
    inputs = list(gene_ids = gene_ids, values = values, pad_mask = pad_mask,
                  cell_names = rownames(embedding), vocab = vocab, assay = "RNA"),
    analyses = list(), plots = list(), history = list()
  ), class = c("fmverse_result", "list"))
}

test_that("legacy Seurat assays are upgraded to the v5 representation", {
  counts <- matrix(
    1:12, nrow = 3L,
    dimnames = list(paste0("gene", 1:3), paste0("cell", 1:4))
  )
  object <- Seurat::CreateSeuratObject(counts = counts)
  object@assays$RNA <- SeuratObject::CreateAssayObject(counts = counts)
  expect_true(inherits(object[["RNA"]], "Assay"))
  converted <- fmverse_to_seurat_v5(object)
  expect_true(inherits(converted[["RNA"]], "Assay5"))
  expect_true("counts" %in% SeuratObject::Layers(converted[["RNA"]]))
})

test_that("downstream functions keep the result class and expose editable ggplots", {
  testthat::skip_if_not_installed("uwot")
  object <- make_downstream_result()
  clustered <- .fm_cluster_legacy(object, method = "louvain", dims = 1:5, neighbors = 8L)
  expect_s3_class(clustered, "fmverse_result")
  expect_identical(rownames(clustered$embeddings), rownames(object$embeddings))
  expect_equal(clustered$embeddings, object$embeddings)
  expect_true("fm_cluster" %in% colnames(clustered$seurat[[]]))
  expect_s3_class(clustered$plots$cluster_umap, "ggplot")
  expect_s3_class(clustered$plots$cluster_umap + ggplot2::theme_bw(), "ggplot")

  reference <- make_downstream_result(seed = 4L)
  reference$metadata$checkpoint_name <- "scGPT-whole-human"
  reference$seurat$cell_type <- rep(c("T", "B"), each = 24L)
  annotated <- .fm_annotate_legacy(object, reference, label = "cell_type", k = 3L, truth = "cell_type")
  expect_s3_class(annotated, "fmverse_result")
  expect_true("fm_annotation" %in% colnames(annotated$seurat[[]]))
  expect_s3_class(annotated$plots$annotation_umap, "ggplot")
  expect_s3_class(annotated$plots$annotation_confusion, "ggplot")
})

test_that("Leiden clustering and GO enrichment return the common result object", {
  testthat::skip_if_not_installed("leidenbase")
  object <- make_downstream_result()
  clustered <- .fm_cluster_legacy(object, method = "leiden", dims = 1:5,
                           neighbors = 8L, resolution = 0.6)
  expect_s3_class(clustered, "fmverse_result")
  expect_s3_class(clustered$plots$cluster_umap, "ggplot")

  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("org.Hs.eg.db")
  marker_genes <- c("CD3D", "CD3E", "CD2", "TRAC", "TRBC1", "TRBC2", "CD4",
                    "IL7R", "LCK", "CD247", "CD28", "CCR7", "LTB", "MAL", "TCF7")
  markers <- data.frame(gene = marker_genes, p_val = rep(1e-8, length(marker_genes)),
                        p_val_adj = rep(1e-6, length(marker_genes)),
                        cluster = "T_cells", stringsAsFactors = FALSE)
  enriched <- .fm_enrich_legacy(clustered, markers = markers, ontology = "BP", pvalue_cutoff = 0.1)
  expect_s3_class(enriched, "fmverse_result")
  expect_true(is.data.frame(enriched$analyses$enrichment$table))
  expect_s3_class(enriched$plots$enrichment_dotplot, "ggplot")
  expect_s3_class(enriched$plots$enrichment_barplot, "ggplot")
  expect_lte(nrow(enriched$plots$enrichment_dotplot$data), 30L)
})

test_that("markers and drug-response result tables are retained", {
  object <- make_downstream_result()
  object$seurat$fm_cluster <- factor(rep(c("0", "1"), each = 24L))
  marked <- .fm_markers_legacy(object, group_by = "fm_cluster", min_pct = 0,
                       logfc_threshold = 0, top_n = 6L)
  expect_s3_class(marked, "fmverse_result")
  expect_gt(nrow(marked$analyses$markers$table), 0L)
  expect_s3_class(marked$plots$markers_dotplot, "ggplot")
  expect_s3_class(marked$plots$markers_volcano, "ggplot")
  plot_file <- tempfile(fileext = ".png")
  ggplot2::ggsave(plot_file, plot = marked$plots$markers_volcano,
                  width = 4, height = 3, dpi = 72)
  expect_true(file.exists(plot_file) && file.info(plot_file)$size > 0)
  unlink(plot_file)

  drug <- .fm_drug_response_legacy(object, treatment = "treatment", drug = "drug_x",
                           control_label = "control", dose = "dose", seed = 3L)
  expect_s3_class(drug, "fmverse_result")
  expect_true(all(c("auc", "auprc") %in% names(drug$analyses$drug_response$metrics)))
  expect_s3_class(drug$plots$drug_response_roc, "ggplot")
  expect_s3_class(drug$plots$drug_response_dose, "ggplot")
})

testthat::skip_if_not_installed("torch")

test_downstream_fake_model <- torch::nn_module(
  classname = "test_downstream_fake_model",
  initialize = function() {},
  forward = function(gene_ids, values, pad_mask = NULL, return_expression = FALSE,
                     return_attention = FALSE, ...) {
    result <- list()
    if (return_expression) result$expression_prediction <- values +
      torch::torch_mean(values, dim = 2, keepdim = TRUE)
    if (return_attention) {
      n <- as.integer(gene_ids$size()[[1L]])
      l <- as.integer(gene_ids$size()[[2L]])
      result$attention <- torch::torch_ones(c(n, l, l)) / l
    }
    result
  }
)

test_that("decoder task interfaces produce ggplot objects from saved tokens", {
  object <- make_downstream_result()
  object$model <- NULL
  object$metadata$checkpoint <- NULL
  expect_error(fmverse:::fmverse_ensure_model(object), "no live model")

  object <- make_downstream_result()
  object$model <- test_downstream_fake_model()
  reconstructed <- .fm_reconstruct_legacy(object, max_cells = 8L)
  expect_s3_class(reconstructed$plots$reconstruction_scatter, "ggplot")
  expect_gt(nrow(reconstructed$analyses$reconstruct$comparison), 0L)

  perturbed <- .fm_perturb_legacy(object, genes = "gene1", max_cells = 8L, top_n = 4L)
  expect_s3_class(perturbed$plots$perturbation_heatmap, "ggplot")
  expect_true(all(c("perturbation", "gene", "predicted_effect") %in%
                    colnames(perturbed$analyses$perturbation$effects)))
  expect_true(any(abs(perturbed$analyses$perturbation$effects$predicted_effect) > 0))
  expect_length(perturbed$analyses$perturbation$cells_by_perturbation$gene1, 8L)

  network <- .fm_grn_legacy(object, max_cells = 4L, max_length = 5L, top_n = 10L)
  expect_s3_class(network$plots$grn_network, "ggplot")
  expect_true(all(c("source", "target", "weight") %in%
                    colnames(network$analyses$grn$edges)))
})

test_that("Harmony and Slingshot add reductions, metadata, and plots", {
  object <- make_downstream_result()
  testthat::skip_if_not_installed("harmony")
  testthat::skip_if_not_installed("uwot")
  integrated <- .fm_integrate_legacy(object, batch = "batch", dims = 1:5,
                             group_by = "cell_type", seed = 4L)
  expect_equal(integrated$embeddings, object$embeddings)
  expect_true(all(c("batch_same_neighbor_fraction_before",
                    "batch_same_neighbor_fraction_after") %in%
                    names(integrated$analyses$integration$metrics)))
  expect_true("fm_harmony" %in% SeuratObject::Reductions(integrated$seurat))
  expect_s3_class(integrated$plots$integration_batch_before, "ggplot")
  expect_s3_class(integrated$plots$integration_batch_after, "ggplot")

  testthat::skip_if_not_installed("slingshot")
  trajectory <- .fm_trajectory_legacy(object, group_by = "cell_type", root = "A",
                              dims = 1:5, seed = 4L)
  expect_true("fm_pseudotime_1" %in% colnames(trajectory$seurat[[]]))
  expect_s3_class(trajectory$plots$trajectory_pseudotime, "ggplot")
  expect_s3_class(trajectory$plots$trajectory_curves, "ggplot")
  expect_identical(trajectory$plots$trajectory_curves$labels$x, "UMAP_1")
  expect_identical(trajectory$plots$trajectory_curves$labels$y, "UMAP_2")
})
