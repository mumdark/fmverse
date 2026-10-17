# Numerical-kernel regressions; public S4 APIs are tested in test-fmobj.R.
testthat::skip_if_not_installed("Seurat")
testthat::skip_if_not_installed("SeuratObject")
testthat::skip_if_not_installed("torch")
testthat::skip_if_not_installed("ggplot2")

test_that("current Geneformer checkpoints expose versioned input specifications", {
  registry <- fmverse_geneformer_checkpoints()
  testthat::expect_identical(registry$checkpoint, c(
    "Geneformer-V1-10M", "Geneformer-V2-104M", "Geneformer-V2-316M",
    "Geneformer-V2-104M_CLcancer"
  ))
  testthat::expect_identical(registry$dictionary_release,
    c("gc30M", "gc104M", "gc104M", "gc104M"))
  testthat::expect_identical(registry$input_size, c(2048L, 4096L, 4096L, 4096L))
  testthat::expect_identical(registry$special_tokens, c(FALSE, TRUE, TRUE, TRUE))
  testthat::expect_identical(
    fmverse:::fmverse_normalize_geneformer_checkpoint_name("V2-104M_CLcancer"),
    "Geneformer-V2-104M_CLcancer"
  )
})

make_geneformer_downstream_result <- function(seed = 17L, n_cells = 24L) {
  set.seed(seed)
  config <- list(
    vocab_size = 32L, hidden_size = 8L, num_hidden_layers = 2L,
    num_attention_heads = 2L, intermediate_size = 16L,
    max_position_embeddings = 32L, type_vocab_size = 2L, pad_token_id = 0L,
    hidden_dropout_prob = 0, attention_probs_dropout_prob = 0,
    layer_norm_eps = 1e-12
  )
  n_genes <- 16L
  cell_type <- rep(c("A", "B"), each = n_cells / 2L)
  batch <- rep(rep(c("b1", "b2"), each = n_cells / 4L), 2L)
  counts <- matrix(stats::rpois(n_genes * n_cells, lambda = 2), nrow = n_genes,
    dimnames = list(paste0("GENE", seq_len(n_genes)), paste0("cell", seq_len(n_cells))))
  counts[1:4, cell_type == "B"] <- counts[1:4, cell_type == "B"] + 7L
  object <- Seurat::CreateSeuratObject(counts = counts, meta.data = data.frame(
    cell_type = cell_type, batch = batch,
    treatment = ifelse(cell_type == "B", "drug_x", "control"),
    dose = ifelse(cell_type == "B", 1, 0), row.names = colnames(counts)
  ))
  vocab <- c("<pad>" = 0L, "<mask>" = 1L, stats::setNames(2:17, rownames(counts)),
             stats::setNames(18:31, paste0("OTHER", 1:14)))
  gene_ids <- matrix(rep(2:17, each = n_cells), nrow = n_cells)
  attention_mask <- matrix(TRUE, nrow = n_cells, ncol = n_genes)
  model <- fmverse_geneformer_build_model(config)
  model$eval()
  input_tensor <- torch::torch_tensor(gene_ids, dtype = torch::torch_long())
  mask_tensor <- torch::torch_tensor(attention_mask, dtype = torch::torch_long())
  output <- torch::with_no_grad(model(input_ids = input_tensor,
    attention_mask = mask_tensor, return_layers = TRUE))
  embedding <- torch::as_array(output$cell_embedding$cpu())
  embedding <- matrix(embedding, nrow = n_cells, ncol = config$hidden_size)
  rownames(embedding) <- colnames(counts)
  colnames(embedding) <- paste0("Geneformer_", seq_len(ncol(embedding)))
  object[["Geneformer"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding, key = "Geneformer_", assay = "RNA"
  )
  token_labels <- stats::setNames(names(vocab), as.character(unname(vocab)))
  symbol_to_token <- stats::setNames(unname(vocab[rownames(counts)]), rownames(counts))
  structure(list(
    model_name = "Geneformer", embeddings = embedding, seurat = object,
    model = model, audit = list(checkpoint_key_count = 0L),
    inputs = list(
      gene_ids = gene_ids, values = matrix(0, nrow = n_cells, ncol = n_genes),
      pad_mask = !attention_mask, attention_mask = attention_mask,
      cell_names = rownames(embedding), vocab = vocab,
      gene_symbols = token_labels, symbol_to_token = symbol_to_token,
      assay = "RNA", token_counts = rep.int(n_genes, n_cells)
    ),
    analyses = list(), plots = list(), history = list(),
    metadata = list(checkpoint_name = "gf-12L-95M-i4096", device = "cpu", config = config)
  ), class = c("fmverse_result", "list"))
}

test_that("Geneformer R model exposes BERT layers and a tied-vocabulary MLM head", {
  result <- make_geneformer_downstream_result(n_cells = 8L)
  ids <- torch::torch_tensor(result$inputs$gene_ids[1:2, 1:5, drop = FALSE],
                             dtype = torch::torch_long())
  mask <- torch::torch_ones_like(ids)
  output <- torch::with_no_grad(result$model(
    input_ids = ids, attention_mask = mask, return_layers = TRUE, return_logits = TRUE
  ))
  testthat::expect_length(output$hidden_states, 3L)
  testthat::expect_equal(as.integer(output$cell_embedding$size()), c(2L, 8L))
  testthat::expect_equal(as.integer(output$logits$size()), c(2L, 5L, 32L))
  testthat::expect_true(all(is.finite(torch::as_array(output$cell_embedding))))
})

test_that("Geneformer vocabulary keeps special-token names for downstream decoding", {
  dictionary <- data.frame(
    gene_id = c("<pad>", "<mask>", "ENSG00000000003"),
    token_id = 0:2,
    stringsAsFactors = FALSE
  )
  vocab <- fmverse:::fmverse_geneformer_vocab(dictionary)
  testthat::expect_identical(unname(vocab[["<mask>"]]), 1L)
  testthat::expect_identical(unname(vocab[["<pad>"]]), 0L)
})

test_that("Geneformer V2 pooling excludes CLS and EOS boundary tokens", {
  config <- list(
    vocab_size = 12L, hidden_size = 8L, num_hidden_layers = 2L,
    num_attention_heads = 2L, intermediate_size = 16L,
    max_position_embeddings = 8L, type_vocab_size = 2L, pad_token_id = 0L,
    hidden_dropout_prob = 0, attention_probs_dropout_prob = 0,
    layer_norm_eps = 1e-12, cls_token_id = 2L, eos_token_id = 3L
  )
  model <- fmverse_geneformer_build_model(config)
  model$eval()
  ids <- torch::torch_tensor(matrix(c(2L, 5L, 7L, 3L), nrow = 1L),
                             dtype = torch::torch_long())
  mask <- torch::torch_ones_like(ids)
  output <- torch::with_no_grad(model(input_ids = ids, attention_mask = mask,
                                       return_layers = TRUE))
  expected <- torch::torch_mean(output$hidden_states[[2L]][, 2:3, ], dim = 2L)
  testthat::expect_equal(
    as.numeric(torch::as_array(output$cell_embedding)),
    as.numeric(torch::as_array(expected)), tolerance = 1e-6
  )
})

test_that("Geneformer result chains through common downstream functions", {
  testthat::skip_if_not_installed("uwot")
  result <- make_geneformer_downstream_result()
  original <- result$embeddings
  clustered <- .fm_cluster_legacy(result, method = "louvain", dims = 1:6, neighbors = 6L)
  testthat::expect_s3_class(clustered, "fmverse_result")
  testthat::expect_equal(clustered$embeddings, original)
  testthat::expect_true("fm_cluster" %in% colnames(clustered$seurat[[]]))
  testthat::expect_s3_class(clustered$plots$cluster_umap, "ggplot")
  testthat::expect_s3_class(clustered$plots$cluster_umap + ggplot2::theme_bw(), "ggplot")

  reference <- make_geneformer_downstream_result(seed = 22L)
  annotated <- .fm_annotate_legacy(clustered, reference, label = "cell_type", truth = "cell_type",
                           reduction = "Geneformer", k = 3L)
  testthat::expect_true("fm_annotation" %in% colnames(annotated$seurat[[]]))
  testthat::expect_s3_class(annotated$plots$annotation_umap, "ggplot")
  testthat::expect_s3_class(annotated$plots$annotation_confusion, "ggplot")

  marked <- .fm_markers_legacy(clustered, group_by = "cell_type", min_pct = 0,
                       logfc_threshold = 0, top_n = 8L)
  testthat::expect_gt(nrow(marked$analyses$markers$table), 0L)
  testthat::expect_s3_class(marked$plots$markers_dotplot, "ggplot")
  testthat::expect_s3_class(marked$plots$markers_volcano, "ggplot")

  reconstructed <- .fm_reconstruct_legacy(clustered, cells = clustered$inputs$cell_names[1:2],
                                  max_cells = 2L, seed = 2L)
  testthat::expect_gt(nrow(reconstructed$analyses$reconstruct$comparison), 0L)
  testthat::expect_s3_class(reconstructed$plots$reconstruction_scatter, "ggplot")

  perturbed <- .fm_perturb_legacy(clustered, genes = "GENE1",
                          cells = clustered$inputs$cell_names[1:2], max_cells = 2L, top_n = 5L)
  testthat::expect_gt(nrow(perturbed$analyses$perturbation$effects), 0L)
  testthat::expect_s3_class(perturbed$plots$perturbation_heatmap, "ggplot")

  network <- .fm_grn_legacy(clustered, cells = clustered$inputs$cell_names[1:2],
                    max_cells = 2L, max_length = 8L, top_n = 12L)
  testthat::expect_true(all(c("source", "target", "weight") %in%
                              colnames(network$analyses$grn$edges)))
  testthat::expect_s3_class(network$plots$grn_network, "ggplot")

  drug <- .fm_drug_response_legacy(clustered, treatment = "treatment", drug = "drug_x",
                           control_label = "control", dose = "dose", seed = 3L)
  testthat::expect_true(all(c("auc", "auprc") %in%
                              names(drug$analyses$drug_response$metrics)))
  testthat::expect_s3_class(drug$plots$drug_response_roc, "ggplot")
  testthat::expect_s3_class(drug$plots$drug_response_dose, "ggplot")
})
