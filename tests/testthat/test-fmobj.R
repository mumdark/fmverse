make_fmobj_fixture <- function(model = "scGPT", seed = 4L) {
  set.seed(seed)
  genes <- c("CD3D", "CD3E", "CD2", "TRAC", "TRBC1", "CD4", "IL7R", "LCK",
    "CD247", "CD28", "CCR7", "LTB", "MAL", "TCF7", "MS4A1", "CD79A")
  labels <- rep(c("A", "B"), each = 24L)
  counts <- matrix(rpois(16L * 48L, 3), 16L,
    dimnames = list(genes, paste0("cell", seq_len(48L))))
  counts[1:6, labels == "B"] <- counts[1:6, labels == "B"] + 12L
  s <- Seurat::CreateSeuratObject(counts = counts, meta.data = data.frame(
    cell_type = labels, batch = rep(c("b1", "b2"), 24L),
    treatment = ifelse(labels == "B", "drug", "control"), dose = as.integer(labels == "B"),
    row.names = colnames(counts)))
  x <- cbind(ifelse(labels == "B", 3, -3) + rnorm(48L, sd = .1), matrix(rnorm(48L * 5L), 48L))
  dimnames(x) <- list(colnames(s), paste0(model, "_", seq_len(ncol(x))))
  s[[model]] <- SeuratObject::CreateDimReducObject(x, key = paste0(model, "_"), assay = "RNA")
  ids <- matrix(rep(2:17, each = 48L), 48L)
  structure(list(model_name = model, seurat = s, embeddings = x, model = NULL,
    metadata = list(checkpoint_name = if (model == "scGPT") "whole-human" else "gf-12L-95M-i4096"),
    audit = list(matched_key_count = 1L), inputs = list(gene_ids = ids,
      values = matrix(runif(length(ids), 1, 5), 48L), pad_mask = matrix(FALSE, 48L, 16L),
      attention_mask = matrix(1L, 48L, 16L), cell_names = rownames(x), assay = "RNA",
      vocab = c("<pad>" = 0L, stats::setNames(1L, if (model == "scGPT") "<cls>" else "<mask>"),
        stats::setNames(2:17, genes))),
    analyses = list(), plots = list(), history = list()), class = c("fmverse_result", "list"))
}

test_that("S4 fmobj migrates a legacy run and stores embedding only in Seurat", {
  old <- make_fmobj_fixture()
  obj <- as_fmobj(old, run_id = "s")
  expect_s4_class(obj, "fmobj")
  expect_setequal(methods::slotNames(obj), c("seurat", "runs", "active_run", "history", "version"))
  expect_true(methods::validObject(obj))
  expect_identical(fm_embeddings(obj), old$embeddings)
  expect_false("embeddings" %in% names(fm_run(obj)))
  expect_false("seurat" %in% names(fm_run(obj)))
  expect_false("model" %in% names(fm_run(obj)$metadata))
  expect_output(show(obj), "1 model runs")
  expect_error(fm_embeddings(obj, run_id = "unknown"), "Unknown run_id")
  expect_identical(dim(fm_seurat(obj)), dim(old$seurat))
  broken <- obj; broken@runs$s$cells <- rev(broken@runs$s$cells)
  expect_error(methods::validObject(broken), "row order")
})

test_that("empty containers upgrade v4 assays and preserve input metadata", {
  s <- make_fmobj_fixture()$seurat
  counts <- SeuratObject::LayerData(s, layer = "counts")
  s@assays$RNA <- SeuratObject::CreateAssayObject(counts = counts)
  s$fm_batch <- s$batch
  empty <- fmobj(s)
  expect_s4_class(empty, "fmobj")
  expect_s4_class(fm_seurat(empty)[["RNA"]], "Assay5")
  expect_equal(nrow(fm_runs(empty)), 0L)
  expect_equal(SeuratObject::LayerData(fm_seurat(empty), layer = "counts"), counts)
  old <- make_fmobj_fixture(); old$seurat$fm_batch <- old$seurat$batch
  obj <- as_fmobj(old)
  expect_identical(fm_seurat(obj)$fm_batch, old$seurat$fm_batch)
  old$seurat@assays$RNA <- SeuratObject::CreateAssayObject(counts = counts)
  upgraded <- as_fmobj(old)
  expect_s4_class(fm_seurat(upgraded)[["RNA"]], "Assay5")
  expect_error(fmverse(obj, run_id = NA_character_), "model_dir.*missing")
})

test_that("cell-subset runs preserve the full shared dataset and other embeddings", {
  old <- make_fmobj_fixture()
  obj <- as_fmobj(old, "full")
  partial <- old
  ix <- seq_len(16L)
  partial$embeddings <- partial$embeddings[ix, , drop = FALSE]
  for (field in c("gene_ids", "values", "pad_mask", "attention_mask")) {
    partial$inputs[[field]] <- partial$inputs[[field]][ix, , drop = FALSE]
  }
  partial$inputs$cell_names <- rownames(partial$embeddings)
  partial$seurat[["scGPT"]] <- SeuratObject::CreateDimReducObject(partial$embeddings, key = "scGPT_", assay = "RNA")
  obj <- .fm_add_legacy(obj, partial, "subset")
  obj <- fm_cluster(obj, run_id = "subset", method = "louvain", dims = 1:5, neighbors = 4L)
  expect_equal(ncol(fm_seurat(obj)), 48L)
  expect_equal(nrow(fm_embeddings(obj, "subset")), 16L)
  expect_identical(fm_embeddings(obj, "full"), old$embeddings)
  column <- fm_analysis(obj, "subset", "cluster")$cluster_column
  expect_true(all(is.na(fm_seurat(obj)[[]][setdiff(colnames(old$seurat), partial$inputs$cell_names), column])))
})

test_that("multiple runs, overwrite cleanup, automatic IDs and save/load are isolated", {
  a <- make_fmobj_fixture(); b <- make_fmobj_fixture("Geneformer", seed = 9L)
  obj <- .fm_add_legacy(as_fmobj(a, "s"), b, "g")
  expect_equal(fm_runs(obj)$run_id, c("s", "g"))
  expect_identical(fm_embeddings(obj, "s"), a$embeddings)
  expect_identical(fm_embeddings(obj, "g"), b$embeddings)
  expect_identical(obj@active_run, "g")
  expect_error(fmverse(obj, model_dir = tempdir(), run_id = "s"),
    "missing `fmverse.json`")
  expect_error(.fm_add_legacy(obj, a, "s"), "already exists")
  obj <- fm_cluster(obj, run_id = "s", method = "louvain", dims = 1:5, neighbors = 8L)
  stale_names <- fm_run(obj, "s")$artifacts$metadata
  replaced <- .fm_add_legacy(obj, a, "s", overwrite = TRUE)
  expect_length(fm_run(replaced, "s")$analyses, 0L)
  expect_false(any(stale_names %in% colnames(fm_seurat(replaced)[[]])))
  expect_identical(fm_embeddings(replaced, "g"), b$embeddings)
  auto <- .fm_add_legacy(as_fmobj(a), a)
  expect_equal(length(unique(fm_runs(auto)$run_id)), 2L)
  p <- tempfile(fileext = ".rds"); saveRDS(obj, p); restored <- readRDS(p)
  expect_identical(fm_embeddings(restored, "g"), b$embeddings)
  expect_s3_class(fm_plot(restored, "s", "cluster", "umap"), "ggplot")
  expect_error(fm_reconstruct(restored, run_id = "s", max_cells = 1L), "no live model")
})

test_that("run and active-run selection agree and neither pollutes another run", {
  a <- make_fmobj_fixture(); b <- make_fmobj_fixture("Geneformer", 10L)
  obj <- .fm_add_legacy(as_fmobj(a, "s"), b, "g")
  counts <- SeuratObject::LayerData(fm_seurat(obj), layer = "counts")
  explicit <- fm_cluster(obj, run_id = "s", method = "louvain", dims = 1:5, neighbors = 8L, seed = 8L)
  active <- fm_use(obj, "s") |> fm_cluster(method = "louvain", dims = 1:5, neighbors = 8L, seed = 8L)
  expect_equal(fm_analysis(explicit, "s", "cluster")$table, fm_analysis(active, "s", "cluster")$table)
  expect_identical(explicit@active_run, "g")
  both <- fm_cluster(explicit, run_id = "g", method = "louvain", dims = 1:5, neighbors = 8L)
  expect_false(identical(fm_analysis(both, "s", "cluster")$cluster_column,
    fm_analysis(both, "g", "cluster")$cluster_column))
  expect_equal(SeuratObject::LayerData(fm_seurat(both), layer = "counts"), counts)
  expect_identical(fm_embeddings(both, "s"), a$embeddings)
  expect_identical(fm_embeddings(both, "g"), b$embeddings)
  expect_error(fm_cluster(both, run_id = "g", reduction = "scGPT"), "does not belong")
  expect_error(fm_cluster(both, run_id = "s", embedding = "integrated"), "Analysis not found")
  old <- .fm_cluster_legacy(a, method = "louvain", dims = 1:5, neighbors = 8L)
  migrated <- as_fmobj(old, "old")
  expect_s3_class(fm_plot(migrated, task = "cluster", plot = "umap") + ggplot2::theme_bw(), "ggplot")
})

test_that("marker dependencies are invalidated, including downstream enrichment", {
  obj <- as_fmobj(make_fmobj_fixture(), "s") |>
    fm_cluster(method = "louvain", dims = 1:5, neighbors = 8L) |>
    fm_markers(min_pct = 0, logfc_threshold = 0)
  expect_gt(nrow(fm_analysis(obj, task = "markers")$table), 0L)
  expect_s3_class(fm_plot(obj, task = "markers", plot = "dotplot"), "ggplot")
  if (requireNamespace("clusterProfiler", quietly = TRUE) && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    obj <- fm_enrich(obj, pvalue_cutoff = 1)
    expect_s3_class(fm_plot(obj, task = "enrichment", plot = "dotplot"), "ggplot")
  }
  obj <- fm_cluster(obj, method = "louvain", dims = 1:5, neighbors = 8L, resolution = 1)
  expect_error(fm_analysis(obj, task = "markers"), "stale")
  expect_true(fm_analysis(obj, task = "markers", allow_stale = TRUE)$stale)
  expect_error(fm_enrich(obj), "stale")
  if (!is.null(fm_run(obj)$analyses$enrichment)) expect_error(fm_plot(obj, task = "enrichment"), "stale")
})

test_that("annotation, Harmony, trajectory and drug response use selected run", {
  obj <- as_fmobj(make_fmobj_fixture(), "q")
  ref <- as_fmobj(make_fmobj_fixture(seed = 8L), "reference")
  obj <- fm_annotate(obj, ref, label = "cell_type", truth = "cell_type", reference_run_id = "reference", k = 3L)
  expect_s3_class(fm_plot(obj, task = "annotation", plot = "confusion"), "ggplot")
  expect_identical(fm_analysis(obj, task = "annotation")$parameters$reference_run_id, "reference")
  wrong <- ref; wrong@runs$reference$metadata$checkpoint_name <- "blood"
  expect_error(fm_annotate(obj, wrong, label = "cell_type"), "same scGPT checkpoint")
  if (requireNamespace("harmony", quietly = TRUE)) {
    obj <- fm_integrate(obj, batch = "batch", dims = 1:5, group_by = "cell_type", reduction_name = "custom_integrated")
    expect_equal(nrow(fm_embeddings(obj, embedding = "integrated")), 48L)
    expect_null(fm_analysis(obj, task = "integration")$embeddings)
    obj <- fm_cluster(obj, reduction = "custom_integrated", method = "louvain", dims = 1:5, neighbors = 8L)
    expect_equal(fm_analysis(obj, task = "cluster")$depends_on, "integration")
    obj <- fm_integrate(obj, batch = "batch", dims = 1:5)
    expect_false("custom_integrated" %in% names(fm_run(obj)$artifacts$reductions))
    expect_error(fm_integrate(obj, batch = "batch", reduction_name = "scGPT"), "original embedding is protected")
    expect_error(fm_analysis(obj, task = "cluster"), "stale")
  }
  if (requireNamespace("slingshot", quietly = TRUE)) {
    obj <- fm_trajectory(obj, group_by = "cell_type", root = "A", dims = 1:5)
    expect_equal(nrow(fm_analysis(obj, task = "trajectory")$pseudotime), 48L)
    expect_s3_class(fm_plot(obj, task = "trajectory", plot = "curves"), "ggplot")
  }
  obj <- fm_drug_response(obj, treatment = "treatment", drug = "drug", dose = "dose")
  expect_s3_class(fm_plot(obj, task = "drug_response", plot = "roc"), "ggplot")
  gf <- make_fmobj_fixture("Geneformer", 10L)
  both <- .fm_add_legacy(obj, gf, "g")
  both <- fm_drug_response(both, run_id = "g", treatment = "treatment", drug = "drug")
  alone <- fm_drug_response(as_fmobj(gf, "g"), treatment = "treatment", drug = "drug")
  expect_equal(fm_analysis(both, "g", "drug_response")$predictions,
    fm_analysis(alone, "g", "drug_response")$predictions)
})

test_that("decoder APIs rehydrate only the requested run model", {
  old <- make_fmobj_fixture()
  fake <- torch::nn_module(initialize = function() {}, forward = function(gene_ids, values,
    pad_mask = NULL, return_expression = FALSE, return_attention = FALSE, ...) {
      out <- list()
      if (return_expression) out$expression_prediction <- values + 1
      if (return_attention) out$attention <- torch::torch_ones(c(gene_ids$size(1), gene_ids$size(2), gene_ids$size(2)))
      out
  })()
  testthat::local_mocked_bindings(fmverse_ensure_model = function(result) fake, .package = "fmverse")
  obj <- as_fmobj(old, "s")
  obj <- fm_reconstruct(obj, max_cells = 2L)
  obj <- fm_perturb(obj, genes = "CD3D", max_cells = 2L, top_n = 4L)
  obj <- fm_grn(obj, max_cells = 2L, max_length = 8L, top_n = 4L)
  expect_s3_class(fm_plot(obj, task = "reconstruct", plot = "scatter"), "ggplot")
  expect_s3_class(fm_plot(obj, task = "perturbation", plot = "heatmap"), "ggplot")
  expect_s3_class(fm_plot(obj, task = "grn", plot = "network"), "ggplot")
  expect_false(any(vapply(fm_run(obj), inherits, logical(1), "nn_module")))
})
