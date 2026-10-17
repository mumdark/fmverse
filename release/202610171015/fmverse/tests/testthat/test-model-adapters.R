test_that("scBERT and scFoundation expose their supported official checkpoints", {
  sb <- fmverse_scbert_checkpoints()
  sf <- fmverse_scfoundation_checkpoints()
  expect_equal(sb$checkpoint, "panglao-human")
  expect_equal(sb$genes, 16906L)
  expect_equal(sf$checkpoint, "01B-resolution")
  expect_equal(sf$genes, 19264L)
  expect_error(fmverse:::fmverse_resolve_scbert_checkpoint("unknown"),
    "Unknown scBERT checkpoint")
  expect_error(fmverse:::fmverse_resolve_scfoundation_checkpoint("unknown"),
    "Unknown scFoundation checkpoint")
})

test_that("scBERT state audit and hidden-state forward work in R torch", {
  model <- fmverse:::fmverse_scbert_model(
    num_tokens = 7L, max_seq_len = 8L, dim = 8L,
    depth = 1L, heads = 2L, dim_head = 4L
  )
  path <- tempfile(fileext = ".safetensors")
  safetensors::safe_save_file(model$state_dict(), path)
  loaded <- fmverse_load_scbert_weights(model, path)
  expect_equal(loaded$audit$matched_key_count, loaded$audit$model_key_count)
  expect_equal(loaded$audit$shape_mismatch_count, 0L)

  # 0 is a valid Python token ID and maps to the first row of the R embedding.
  ids <- torch::torch_tensor(matrix(c(0L, 2L, 6L), nrow = 1L),
    dtype = torch::torch_long())
  hidden <- loaded$model(ids, return_encodings = TRUE)
  expect_identical(as.integer(hidden$size()), c(1L, 3L, 8L))
})

test_that("scGPT keeps the batch dimension for singleton forward passes", {
  model <- fmverse_scgpt_model(ntoken = 32L, pad_id = 1L, d_model = 8L,
    nhead = 2L, d_hid = 16L, nlayers = 2L, dropout = 0)
  model$eval()
  gene_ids <- torch::torch_tensor(matrix(c(2L, 3L, 4L, 5L), nrow = 1L),
    dtype = torch::torch_long())
  values <- torch::torch_tensor(matrix(c(0, 1, 2, 3), nrow = 1L),
    dtype = torch::torch_float())
  pad_mask <- torch::torch_tensor(matrix(FALSE, nrow = 1L, ncol = 4L),
    dtype = torch::torch_bool())
  output <- torch::with_no_grad(model(gene_ids, values, pad_mask))
  embedding <- torch::as_array(output$cell_embedding_l2)
  expect_identical(dim(embedding), c(1L, 8L))
  expect_true(all(is.finite(embedding)))
})

test_that("scBERT tokenization retains fixed Panglao positions for absent genes", {
  genes <- paste0("g", seq_len(16906L))
  gene_file <- tempfile(fileext = ".tsv")
  writeLines(c("gene_name", genes), gene_file)
  on.exit(unlink(gene_file), add = TRUE)
  counts <- Matrix::Matrix(matrix(c(3L, 0L, 0L, 5L), ncol = 2L,
    dimnames = list(c("g100", "g10"), c("cell1", "cell2"))), sparse = TRUE)
  object <- Seurat::CreateSeuratObject(counts = counts)
  inputs <- fmverse:::fmverse_scbert_tokenize(object, gene_index_file = gene_file)
  expect_identical(sort(unique(Matrix::summary(inputs$gene_ids)$j)), c(10L, 100L))
})

test_that("scFoundation state audit, special-token embeddings, and forward work", {
  config <- list(
    seq_len = 12L, n_class = 7L, pad_token_id = 5L, mask_token_id = 4L,
    bin_num = 3L, bin_alpha = 1,
    encoder = list(module_type = "performer", hidden_dim = 8L, depth = 1L,
      heads = 2L, dim_head = 4L),
    decoder = list(module_type = "transformer", hidden_dim = 8L, depth = 1L,
      heads = 2L)
  )
  model <- fmverse_scfoundation_build_model(config)
  path <- tempfile(fileext = ".safetensors")
  safetensors::safe_save_file(model$state_dict(), path)
  loaded <- fmverse_load_scfoundation_weights(model, path)
  expect_equal(loaded$audit$matched_key_count, loaded$audit$model_key_count)

  values <- torch::torch_tensor(matrix(c(2, 4, 5), nrow = 1L),
    dtype = torch::torch_float())
  positions <- torch::torch_tensor(matrix(c(1L, 10L, 11L), nrow = 1L),
    dtype = torch::torch_long())
  hidden <- loaded$model(values = values, gene_positions = positions)
  expect_identical(as.integer(hidden$size()), c(1L, 3L, 8L))
})

test_that("scFoundation count preprocessing preserves fixed gene order and total tokens", {
  counts <- Matrix::Matrix(matrix(c(4L, 2L, 0L, 5L), nrow = 2L), sparse = TRUE,
    dimnames = list(c("A1BG", "A1CF"), c("cell1", "cell2")))
  object <- Seurat::CreateSeuratObject(counts = counts)
  inputs <- fmverse:::fmverse_scfoundation_tokenize(object, tgthighres = "a5")
  expect_identical(length(inputs$gene_order), 19264L)
  expect_identical(inputs$gene_order[1:2], c("A1BG", "A1CF"))
  expect_identical(inputs$cell_names, c("cell1", "cell2"))
  expect_identical(inputs$cells[[1L]]$gene_positions,
    c(0L, 1L, 19264L, 19265L))
  expect_identical(inputs$cells[[2L]]$gene_positions,
    c(1L, 19264L, 19265L))
  expect_equal(length(inputs$cells[[1L]]$values), 4L)
})

test_that("scFoundation dispatch creates a shared S4 fmobj run", {
  config <- list(
    seq_len = 19266L, n_class = 7L, pad_token_id = 5L, mask_token_id = 4L,
    bin_num = 3L, bin_alpha = 1,
    encoder = list(module_type = "performer", hidden_dim = 8L, depth = 1L,
      heads = 2L, dim_head = 4L),
    decoder = list(module_type = "transformer", hidden_dim = 8L, depth = 1L,
      heads = 2L)
  )
  checkpoint_dir <- tempfile("scfoundation-checkpoint-")
  dir.create(checkpoint_dir)
  on.exit(unlink(checkpoint_dir, recursive = TRUE), add = TRUE)
  model <- fmverse_scfoundation_build_model(config)
  safetensors::safe_save_file(model$state_dict(),
    file.path(checkpoint_dir, "scfoundation_cell.safetensors"))
  jsonlite::write_json(config,
    file.path(checkpoint_dir, "scfoundation_cell_config.json"),
    auto_unbox = TRUE, pretty = TRUE)
  gene_index <- system.file("extdata", "scfoundation_gene_index_19264.tsv",
    package = "fmverse")
  file.copy(gene_index, file.path(checkpoint_dir, "scfoundation_gene_index_19264.tsv"))
  jsonlite::write_json(list(schema_version = 1L, model = "scFoundation",
    checkpoint = "01B-resolution", files = list(
      weights = "scfoundation_cell.safetensors",
      config = "scfoundation_cell_config.json",
      gene_index = "scfoundation_gene_index_19264.tsv")),
    file.path(checkpoint_dir, "fmverse.json"), auto_unbox = TRUE, pretty = TRUE)
  counts <- Matrix::Matrix(matrix(c(4L, 2L, 0L, 5L), nrow = 2L), sparse = TRUE,
    dimnames = list(c("A1BG", "A1CF"), c("cell1", "cell2")))
  seurat <- Seurat::CreateSeuratObject(counts = counts)

  result <- fmverse(seurat, model_dir = checkpoint_dir, verbose = FALSE)
  expect_s4_class(result, "fmobj")
  expect_equal(fm_run(result)$model, "scFoundation")
  expect_equal(fm_run(result)$metadata$batch_size, 1L)
  expect_identical(rownames(fm_embeddings(result)), c("cell1", "cell2"))
  expect_equal(ncol(fm_embeddings(result)), 32L)
})

test_that("both model reductions flow through the common fmobj downstream interface", {
  set.seed(12L)
  n_cells <- 48L
  genes <- paste0("gene", seq_len(24L))
  labels <- rep(c("A", "B"), each = n_cells / 2L)
  counts <- matrix(rpois(length(genes) * n_cells, 2), nrow = length(genes),
    dimnames = list(genes, paste0("cell", seq_len(n_cells))))
  counts[1:5, labels == "B"] <- counts[1:5, labels == "B"] + 8L
  seurat <- Seurat::CreateSeuratObject(counts = Matrix::Matrix(counts, sparse = TRUE))

  make_run <- function(model_name, seed) {
    set.seed(seed)
    embedding <- matrix(rnorm(n_cells * 8L), nrow = n_cells,
      dimnames = list(colnames(seurat), paste0(model_name, "_", seq_len(8L))))
    embedding[, 1L] <- ifelse(labels == "B", 4, -4) + rnorm(n_cells, sd = 0.05)
    local <- seurat
    local[[model_name]] <- SeuratObject::CreateDimReducObject(
      embeddings = embedding, key = paste0(model_name, "_"), assay = "RNA"
    )
    legacy <- structure(list(model_name = model_name, embeddings = embedding,
      seurat = local, model = NULL, metadata = list(checkpoint_name = "test"),
      audit = list(), inputs = list(cell_names = rownames(embedding), assay = "RNA"),
      analyses = list(), plots = list(), history = list()),
      class = c("fmverse_result", "list"))
    legacy
  }

  obj <- as_fmobj(make_run("scBERT", 21L), run_id = "sb")
  obj <- .fm_add_legacy(obj, make_run("scFoundation", 22L), run_id = "sf")
  for (id in c("sb", "sf")) {
    obj <- fm_cluster(obj, run_id = id, method = "louvain", dims = 1:5,
      neighbors = 8L, seed = 3L)
    expect_s3_class(fm_plot(obj, run_id = id, task = "cluster"), "ggplot")
    obj <- suppressWarnings(fm_markers(obj, run_id = id, top_n = 8L,
      min_pct = 0, logfc_threshold = 0))
    expect_true(nrow(fm_analysis(obj, run_id = id, task = "markers")$table) > 0L)
    expect_s3_class(fm_plot(obj, run_id = id, task = "markers", plot = "dotplot"), "ggplot")
  }
  expect_identical(fm_runs(obj)$run_id, c("sb", "sf"))
})

