test_that("UCE catalogs the official checkpoints and bundled species", {
  checkpoints <- fmverse_uce_checkpoints()
  expect_equal(checkpoints$checkpoint, c("4layer", "33layer"))
  expect_equal(checkpoints$layers, c(4L, 33L))
  expect_equal(nrow(fmverse_uce_species()), 8L)
  expect_true(all(c("object", "model_dir", "...", "run_id", "overwrite") %in%
    names(formals(fmverse))))
  expect_false(any(c("model", "checkpoint") %in% names(formals(fmverse))))
  expect_true(all(c("human", "mouse", "frog", "zebrafish") %in%
    fmverse_uce_species()$species))
})

test_that("small UCE module returns normalized CLS embeddings and layer states", {
  model <- fmverse_uce_model(
    nlayers = 4L, token_dim = 8L, d_model = 16L, nhead = 4L,
    d_hid = 32L, output_dim = 16L, dropout = 0,
    max_len = 20L, n_tokens = 30L
  )
  state <- model$state_dict()
  expect_true("pos_encoder.pe" %in% names(state))
  expect_true("pe_embedding.weight" %in% names(state))
  model$eval()
  token_ids <- torch::torch_tensor(
    matrix(c(3L, 4L, 2L, 3L, 5L, 0L), nrow = 3L),
    dtype = torch::torch_long()
  )
  attention_mask <- torch::torch_tensor(
    matrix(c(1, 1, 1, 1, 1, 0), nrow = 2L), dtype = torch::torch_float()
  )
  result <- model(token_ids, attention_mask, return_layers = TRUE)
  expect_equal(as.integer(result$gene_output$size()), c(3L, 2L, 16L))
  expect_equal(as.integer(result$cell_embedding$size()), c(2L, 16L))
  expect_length(result$layers, 4L)
  expected_tokens <- model$pe_embedding(token_ids + 1L)
  expected_tokens <- expected_tokens / torch::torch_clamp(
    torch::torch_sqrt(torch::torch_sum(expected_tokens^2, dim = 3L, keepdim = TRUE)),
    min = 1e-12
  )
  expect_equal(torch::as_array(result$token_embedding), torch::as_array(expected_tokens),
    tolerance = 1e-6)
  norms <- torch::torch_sqrt(torch::torch_sum(result$cell_embedding^2, dim = 2L))
  expect_equal(as.numeric(torch::as_array(norms)), c(1, 1), tolerance = 1e-5)
})

test_that("UCE loader strictly accepts a CPU safetensors state dictionary", {
  make_small_model <- function() fmverse_uce_model(
    nlayers = 4L, token_dim = 8L, d_model = 8L, nhead = 2L,
    d_hid = 16L, output_dim = 8L, dropout = 0,
    max_len = 8L, n_tokens = 12L
  )
  source <- make_small_model()
  checkpoint <- tempfile(fileext = ".torch")
  file.create(checkpoint)
  safe_checkpoint <- sub("\\.torch$", ".cpu.safetensors", checkpoint)
  safetensors::safe_save_file(source$state_dict(), safe_checkpoint)
  target <- make_small_model()
  loaded <- fmverse_load_uce_weights(target, checkpoint)
  expect_equal(loaded$audit$weight_format, "safetensors")
  expect_identical(normalizePath(loaded$audit$weight_source), normalizePath(safe_checkpoint))
  expect_equal(loaded$audit$checkpoint_tensors, loaded$audit$model_tensors)
  expect_length(loaded$audit$missing_keys, 0L)
  expect_length(loaded$audit$unexpected_keys, 0L)
  expect_equal(
    torch::as_array(loaded$model$pe_embedding$weight),
    torch::as_array(source$pe_embedding$weight), tolerance = 1e-6
  )
})

test_that("UCE tokenization retains cell order and the official special tokens", {
  temp <- tempfile("fmverse-uce-")
  dir.create(temp)
  index_path <- file.path(temp, "species_gene_index.tsv")
  chrom_path <- file.path(temp, "species_chrom.csv")
  utils::write.table(
    data.frame(species = "human", gene_symbol = paste0("G", 1:4),
      pe_row_idx = 0:3), index_path, sep = "\t", quote = FALSE,
    row.names = FALSE
  )
  utils::write.csv(
    data.frame(gene_symbol = paste0("G", 1:4), chromosome = c("1", "1", "2", "2"),
      start = c(10, 20, 30, 40), species = "human"),
    chrom_path, row.names = FALSE
  )
  counts <- Matrix::Matrix(
    matrix(c(2, 1, 0, 3, 1, 0, 4, 1, 2, 2, 1, 1, 3, 2, 2, 1), nrow = 4,
      dimnames = list(paste0("G", 1:4), paste0("cell_", letters[1:4]))),
    sparse = TRUE
  )
  pbmc <- Seurat::CreateSeuratObject(counts = counts)
  assets <- list(gene_index = index_path, chromosome = chrom_path)
  tokens <- fmverse_uce_tokenize(
    pbmc, species = "human", sample_size = 4L, pad_length = 20L,
    seed = 42L, assets = assets
  )
  expect_identical(tokens$cell_names, paste0("cell_", letters[1:4]))
  expect_equal(unname(tokens$token_ids[, 1L]), rep(3L, 4L))
  expect_true(all(tokens$attention_mask[, 1L] == 1L))
  expect_equal(tokens$token_ids, fmverse_uce_tokenize(
    pbmc, species = "human", sample_size = 4L, pad_length = 20L,
    seed = 42L, assets = assets
  )$token_ids)
  split_pbmc <- pbmc
  split_pbmc[["RNA"]] <- split(pbmc[["RNA"]], f = factor(c("batch_a", "batch_a", "batch_b", "batch_b")))
  expect_false("counts" %in% SeuratObject::Layers(split_pbmc[["RNA"]]))
  split_tokens <- fmverse_uce_tokenize(
    split_pbmc, species = "human", sample_size = 4L, pad_length = 20L,
    seed = 42L, assets = assets
  )
  expect_equal(split_tokens$token_ids, tokens$token_ids)
})
