# R implementation of the selected scGPT cell-embedding input path.
fmverse_scgpt_tokenize <- function(object, vocab_file, assay = NULL,
                                   cells = NULL, max_length = 1200L,
                                   n_bins = 51L,
                                   preprocessing = c("official", "deterministic"),
                                   seed = 20261013L) {
  preprocessing <- match.arg(preprocessing)
  if (!file.exists(vocab_file)) stop("Vocabulary file does not exist: ", vocab_file)
  if (max_length < 2L || n_bins < 2L) stop("Invalid max_length or n_bins.")
  if (preprocessing == "official") {
    if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || !is.finite(seed) ||
        seed < 0 || seed > .Machine$integer.max || seed != as.integer(seed)) {
      stop("seed must be one non-negative integer no larger than .Machine$integer.max.")
    }
    old_rng_kind <- RNGkind()
    had_rng_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_rng_seed) old_rng_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (!identical(RNGkind(), old_rng_kind)) {
        do.call(RNGkind, as.list(old_rng_kind))
      }
      if (had_rng_seed) {
        assign(".Random.seed", old_rng_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
  }
  object <- fmverse_to_seurat_v5(object)
  if (is.null(assay)) assay <- SeuratObject::DefaultAssay(object)
  if (!assay %in% SeuratObject::Assays(object)) stop("Assay not found: ", assay)
  assay_object <- object[[assay]]
  if (!"counts" %in% SeuratObject::Layers(assay_object)) {
    stop("Assay ", assay, " has no raw counts layer.")
  }
  counts <- SeuratObject::LayerData(object, assay = assay, layer = "counts")
  counts <- methods::as(counts, "dgCMatrix")
  feature_metadata <- assay_object[[]]
  gene_names <- if ("scGPT_gene_name" %in% colnames(feature_metadata)) {
    as.character(feature_metadata$scGPT_gene_name)
  } else {
    rownames(counts)
  }
  if (length(gene_names) != nrow(counts) || anyNA(gene_names)) {
    stop("Feature-level `scGPT_gene_name` metadata must match the count rows.")
  }
  vocab <- jsonlite::fromJSON(vocab_file, simplifyVector = TRUE)
  if (is.list(vocab)) vocab <- unlist(vocab, use.names = TRUE)
  if (is.null(names(vocab))) stop("vocab.json did not decode to a named vocabulary.")
  vocab <- stats::setNames(as.integer(vocab), names(vocab))
  if (!all(c("<pad>", "<cls>") %in% names(vocab))) {
    stop("Checkpoint vocabulary must include <pad> and <cls>.")
  }
  pad_id <- unname(vocab[["<pad>"]])
  cls_id <- unname(vocab[["<cls>"]])
  if (is.null(cells)) {
    cell_names <- colnames(counts)
  } else {
    cell_names <- as.character(cells)
    if (anyDuplicated(cell_names)) stop("cells contains duplicate barcodes.")
    if (any(!cell_names %in% colnames(counts))) stop("Some requested cells are absent.")
    counts <- counts[, cell_names, drop = FALSE]
  }
  feature_ids <- unname(vocab[gene_names])
  matched <- !is.na(feature_ids)
  if (!any(matched)) stop("No feature names matched the scGPT vocabulary.")
  counts <- counts[which(matched), , drop = FALSE]
  feature_ids <- as.integer(feature_ids[matched])

  n_cells <- ncol(counts)
  max_length <- as.integer(max_length)
  gene_matrix <- matrix(as.integer(pad_id), nrow = n_cells, ncol = max_length)
  value_matrix <- matrix(-2, nrow = n_cells, ncol = max_length)
  token_counts <- integer(n_cells)
  max_genes <- max_length - 1L
  pointers <- counts@p
  for (j in seq_len(n_cells)) {
    if (preprocessing == "official") {
      cell_seed_hex <- digest::digest(
        paste0(as.integer(seed), "|", cell_names[[j]]),
        algo = "crc32", serialize = FALSE
      )
      cell_seed <- strtoi(substr(cell_seed_hex, 1L, 7L), base = 16L)
      set.seed(cell_seed)
    }
    first <- pointers[[j]] + 1L
    last <- pointers[[j + 1L]]
    gene_matrix[j, 1L] <- as.integer(cls_id)
    if (first > last) {
      token_counts[[j]] <- 1L
      next
    }
    offsets <- seq.int(first, last)
    row_index <- counts@i[offsets] + 1L
    expression <- as.numeric(counts@x[offsets])
    keep <- expression != 0
    row_index <- row_index[keep]
    expression <- expression[keep]
    if (preprocessing == "official") {
      # scGPT bins every expressed gene before its collator truncates/samples.
      bins <- stats::quantile(
        expression, probs = seq(0, 1, length.out = as.integer(n_bins) - 1L),
        names = FALSE, type = 7
      )
      left <- findInterval(expression, bins)
      right <- findInterval(expression, bins, left.open = TRUE)
      binned_expression <- ceiling(
        stats::runif(length(expression)) * (right - left) + left
      )
      if (length(expression) > max_genes) {
        # Equivalent in distribution to torch.randperm(n)[:max_genes]:
        # sample without replacement and preserve the sampled random order.
        selected <- sample.int(length(expression), max_genes)
        row_index <- row_index[selected]
        binned_expression <- binned_expression[selected]
      }
    } else {
      if (length(expression) > max_genes) {
        token_ids <- feature_ids[row_index]
        hash <- vapply(token_ids, function(token_id) {
          digest::digest(
            paste0(cell_names[[j]], "|", token_id),
            algo = "crc32", serialize = FALSE
          )
        }, character(1))
        selected <- order(hash)[seq_len(max_genes)]
        selected <- selected[order(row_index[selected])]
        row_index <- row_index[selected]
        expression <- expression[selected]
      }
      bins <- stats::quantile(
        expression, probs = seq(0, 1, length.out = as.integer(n_bins) - 1L),
        names = FALSE, type = 7
      )
      left <- findInterval(expression, bins)
      right <- findInterval(expression, bins, left.open = TRUE)
      binned_expression <- ceiling((left + right) / 2)
    }
    n <- length(row_index)
    if (n > 0L) {
      gene_matrix[j, seq_len(n) + 1L] <- feature_ids[row_index]
      value_matrix[j, seq_len(n) + 1L] <- binned_expression
    }
    token_counts[[j]] <- n + 1L
  }
  list(
    gene_ids = torch::torch_tensor(gene_matrix, dtype = torch::torch_long()),
    values = torch::torch_tensor(value_matrix, dtype = torch::torch_float()),
    pad_mask = torch::torch_tensor(gene_matrix == pad_id, dtype = torch::torch_bool()),
    cell_names = cell_names,
    assay = assay,
    matched_features = sum(matched),
    total_features = nrow(object[[assay]]),
    token_counts = token_counts,
    preprocessing = list(
      value_source = "raw counts layer",
      method = preprocessing,
      binning = if (preprocessing == "official") {
        "scGPT quantiles; random digitize at repeated boundaries (R runif)"
      } else {
        "scGPT quantiles; deterministic midpoint for tied values"
      },
      truncation = if (preprocessing == "official") {
        "bin first, then sample and permute max_length - 1 matched nonzero features"
      } else {
        "stable CRC32 sample of max_length - 1 matched nonzero features"
      },
      seed = if (preprocessing == "official") as.integer(seed) else NULL,
      seed_scheme = if (preprocessing == "official") {
        "R RNG seed from CRC32(seed | cell barcode)"
      } else {
        NULL
      },
      max_length = max_length,
      n_bins = as.integer(n_bins)
    )
  )
}

#' Embed a Seurat object with the official scGPT checkpoint.
#'
#' @param object Seurat v4 or v5 object.
#' @param model Upstream adapter: `scGPT`, `Geneformer`, `scBERT`,
#'   `scFoundation`, or `UCE`.
#' @param batch_size Cells per forward pass.
#' @param checkpoint Official checkpoint name (for example `whole-human`,
#'   `Geneformer-V2-104M`, `panglao-human`, or `01B-resolution`) or a local
#'   checkpoint directory / weight-file path.
#' @param checkpoint_dir Optional root directory containing the registered
#'   `scGPT-<checkpoint>` folders. Also read from `options(fmverse.checkpoint_dir)`
#'   or the `FMVERSE_CHECKPOINT_DIR` environment variable.
#' @param config_file Internal path to a model configuration file from the
#'   selected folder manifest.
#' @param assets_dir Internal path to UCE's self-contained assets directory.
#' @param vocab_file Optional vocab.json path; inferred beside checkpoint.
#' @param max_length Maximum input token sequence length; inferred from the
#'   selected checkpoint configuration when omitted.
#' @param preprocessing Input rules: `official` follows scGPT's stochastic
#'   binning and sampling order; `deterministic` uses midpoint ties and CRC32 truncation.
#' @param seed Random seed used by the official-style R preprocessing.
#' @param gene_index_file Ordered gene index for scBERT or scFoundation.
#' @param gene_id Optional feature-metadata column containing Ensembl IDs for
#'   Geneformer input. Ignored by scGPT.
#' @param species Species key for UCE (for example `human` or `mouse`).
#' @param gene_symbol Optional feature metadata column containing gene symbols
#'   for UCE. Otherwise common symbol fields or feature names are used.
#' @param sample_size Number of expression-weighted gene samples per UCE cell.
#' @param pad_length Maximum token sequence length for UCE.
#' @param pre_normalized Whether scFoundation input is raw counts (`F`) or
#'   already normalized (`T`).
#' @param tgthighres scFoundation total-count resolution token.
#' @param pool_type scFoundation output pooling mode, `all` or `max`.
#' @param cells Optional cell barcodes to embed.
#' @param assay Counts assay.
#' @param device Torch device, CPU by default for reproducible parity.
#' @param verbose Report inference progress.
#' @return A structured fmverse_result with embeddings, Seurat object, model, and audit.
#' @export
.fmverse_infer_legacy <- function(object, model = "scGPT", batch_size = 16L,
                    checkpoint = "whole-human", vocab_file = NULL, max_length = NULL,
                    cells = NULL, assay = NULL, device = "cpu",
                    verbose = interactive(),
                    preprocessing = c("official", "deterministic"),
                    seed = 20261013L, checkpoint_dir = NULL, gene_id = NULL,
                    gene_index_file = NULL, species = "human", gene_symbol = NULL,
                    config_file = NULL, assets_dir = NULL,
                    sample_size = 1024L, pad_length = 1536L, pre_normalized = "F",
                    tgthighres = "a5", pool_type = "all") {
  model_key <- tolower(gsub("[ _-]", "", as.character(model[[1L]])))
  if (identical(model_key, "scbert")) {
    if (missing(checkpoint) || is.null(checkpoint) || identical(checkpoint, "whole-human")) {
      checkpoint <- "panglao-human"
    }
    if (missing(batch_size)) batch_size <- 1L
    return(.fmverse_scbert_legacy(
      object = object, checkpoint = checkpoint, checkpoint_dir = checkpoint_dir,
      gene_index_file = gene_index_file, batch_size = batch_size, cells = cells,
      assay = assay, device = device, verbose = verbose
    ))
  }
  if (identical(model_key, "scfoundation")) {
    if (missing(checkpoint) || is.null(checkpoint) || identical(checkpoint, "whole-human")) {
      checkpoint <- "01B-resolution"
    }
    if (missing(batch_size)) batch_size <- 1L
    return(.fmverse_scfoundation_legacy(
      object = object, checkpoint = checkpoint, checkpoint_dir = checkpoint_dir,
      gene_index_file = gene_index_file, config_file = config_file,
      batch_size = batch_size, cells = cells,
      assay = assay, device = device, verbose = verbose,
      pre_normalized = pre_normalized, tgthighres = tgthighres,
      pool_type = pool_type
    ))
  }
  if (identical(model_key, "uce")) {
    if (missing(checkpoint) || is.null(checkpoint) || identical(checkpoint, "whole-human")) {
      checkpoint <- "4layer"
    }
    if (missing(batch_size)) batch_size <- 1L
    return(.fmverse_uce_legacy(
      object = object, checkpoint = checkpoint, species = species,
      checkpoint_dir = checkpoint_dir, batch_size = batch_size,
      assets_dir = assets_dir,
      sample_size = sample_size, pad_length = pad_length, seed = seed,
      gene_symbol = gene_symbol, cells = cells, assay = assay,
      device = device, verbose = verbose
    ))
  }
  if (identical(model_key, "geneformer")) {
    if (missing(checkpoint)) checkpoint <- "gf-12L-95M-i4096"
    return(.fmverse_geneformer_legacy(
      object = object, checkpoint = checkpoint, checkpoint_dir = checkpoint_dir,
      batch_size = batch_size, max_length = max_length, cells = cells,
      assay = assay, gene_id = gene_id, device = device, verbose = verbose
    ))
  }
  if (!identical(model_key, "scgpt")) {
    stop("Supported models are `scGPT`, `Geneformer`, `scBERT`, `scFoundation`, and `UCE`.")
  }
  if (missing(checkpoint) || is.null(checkpoint) || !length(checkpoint)) {
    stop("Provide checkpoint as best_model.pt or its directory.")
  }
  if (length(batch_size) != 1L || is.na(batch_size) || batch_size < 1L ||
      batch_size != as.integer(batch_size)) {
    stop("batch_size must be a positive integer.")
  }
  preprocessing <- match.arg(preprocessing)
  resolved_checkpoint <- fmverse_resolve_scgpt_checkpoint(
    checkpoint, checkpoint_dir = checkpoint_dir
  )
  checkpoint_path <- resolved_checkpoint$path
  checkpoint_file <- if (dir.exists(checkpoint_path)) {
    file.path(checkpoint_path, "best_model.pt")
  } else {
    checkpoint_path
  }
  if (is.null(vocab_file)) vocab_file <- file.path(dirname(checkpoint_file), "vocab.json")
  args_file <- if (is.null(config_file)) file.path(dirname(checkpoint_file), "args.json") else
    config_file
  if (!file.exists(args_file)) stop("Expected args.json beside the checkpoint.")
  config <- jsonlite::fromJSON(args_file, simplifyVector = TRUE)
  if (is.null(max_length)) {
    max_length <- config$max_seq_len
    if (is.null(max_length) || !length(max_length) || is.na(max_length)) {
      max_length <- 1200L
    }
  }
  if (length(max_length) != 1L || is.na(max_length) || max_length < 2L ||
      max_length != as.integer(max_length)) {
    stop("max_length must be an integer of at least 2.")
  }
  vocab <- jsonlite::fromJSON(vocab_file, simplifyVector = TRUE)
  if (is.list(vocab)) vocab <- unlist(vocab, use.names = TRUE)
  vocab <- stats::setNames(as.integer(vocab), names(vocab))
  object <- fmverse_to_seurat_v5(object)
  tokenized <- fmverse_scgpt_tokenize(
    object, vocab_file = vocab_file, assay = assay, cells = cells,
    max_length = max_length, n_bins = config$n_bins,
    preprocessing = preprocessing, seed = seed
  )
  model_module <- fmverse_scgpt_model(
    ntoken = length(vocab), pad_id = unname(vocab[["<pad>"]]),
    d_model = config$embsize, nhead = config$nheads,
    d_hid = config$d_hid, nlayers = config$nlayers, dropout = config$dropout
  )
  loaded <- fmverse_load_scgpt_weights(
    model_module, checkpoint_file, nlayers_cls = config$n_layers_cls
  )
  model_module <- loaded$model
  torch_device <- torch::torch_device(device)
  model_module$to(device = torch_device)
  model_module$eval()
  n_cells <- length(tokenized$cell_names)
  embeddings <- matrix(NA_real_, nrow = n_cells, ncol = config$embsize)
  starts <- seq.int(1L, n_cells, by = as.integer(batch_size))
  progress_step <- max(1L, ceiling(length(starts) / 20L))
  for (batch_number in seq_along(starts)) {
    first <- starts[[batch_number]]
    index <- first:min(first + as.integer(batch_size) - 1L, n_cells)
    sequence_end <- max(tokenized$token_counts[index])
    sequence_index <- seq_len(sequence_end)
    output <- torch::with_no_grad({
      model_module(
        tokenized$gene_ids[index, sequence_index, drop = FALSE]$to(device = torch_device),
        tokenized$values[index, sequence_index, drop = FALSE]$to(device = torch_device),
        tokenized$pad_mask[index, sequence_index, drop = FALSE]$to(device = torch_device),
        return_layers = FALSE
      )
    })
    batch_embeddings <- torch::as_array(output$cell_embedding_l2$cpu())
    embeddings[index, ] <- matrix(batch_embeddings, nrow = length(index), ncol = config$embsize)
    if (isTRUE(verbose) && (batch_number %% progress_step == 0L || batch_number == length(starts))) {
      message(sprintf("Embedded %d/%d cells", max(index), n_cells))
    }
    rm(output, batch_embeddings)
    # Torch tensors own native memory; collect their R wrappers periodically so
    # long CPU runs over thousands of small batches do not retain stale buffers.
    if (batch_number %% 10L == 0L) invisible(gc(verbose = FALSE))
  }
  rownames(embeddings) <- tokenized$cell_names
  colnames(embeddings) <- paste0("scGPT_", seq_len(ncol(embeddings)))
  object[["scGPT"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embeddings, key = "scGPT_", assay = tokenized$assay
  )
  structure(
    list(
      model_name = "scGPT",
      embeddings = embeddings,
      seurat = object,
      model = model_module,
      audit = loaded$audit,
      inputs = list(
        gene_ids = as.matrix(torch::as_array(tokenized$gene_ids$cpu())),
        values = as.matrix(torch::as_array(tokenized$values$cpu())),
        pad_mask = as.matrix(torch::as_array(tokenized$pad_mask$cpu())),
        cell_names = tokenized$cell_names,
        vocab = vocab,
        assay = tokenized$assay
      ),
      analyses = list(),
      plots = list(),
      history = list(),
      metadata = list(
        checkpoint_name = resolved_checkpoint$name,
        checkpoint = normalizePath(checkpoint_file, winslash = "/", mustWork = TRUE),
        vocabulary = normalizePath(vocab_file, winslash = "/", mustWork = TRUE),
        device = device,
        batch_size = as.integer(batch_size),
        max_length = as.integer(max_length),
        matched_features = tokenized$matched_features,
        total_features = tokenized$total_features,
        preprocessing_mode = preprocessing,
        seed = if (preprocessing == "official") as.integer(seed) else NULL,
        n_cells = n_cells,
        token_counts = tokenized$token_counts,
        preprocessing = tokenized$preprocessing
      )
    ),
    class = c("fmverse_result", "list")
  )
}
