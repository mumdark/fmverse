# Geneformer rank-value tokenization and inference adapter.

fmverse_geneformer_checkpoint_catalog <- function() {
  data.frame(
    checkpoint = c(
      "Geneformer-V1-10M", "Geneformer-V2-104M", "Geneformer-V2-316M",
      "Geneformer-V2-104M_CLcancer"
    ),
    family = c("V1", "V2", "V2", "V2 continual-learning"),
    parameters = c("10M", "104M", "316M", "104M"),
    training_data = c(
      "Genecorpus-30M: ~30M human single-cell transcriptomes",
      "Genecorpus-104M: ~104M non-cancer human single-cell transcriptomes",
      "Genecorpus-104M: ~104M non-cancer human single-cell transcriptomes",
      "V2-104M continued on ~14M cancer-domain transcriptomes"
    ),
    input_size = c(2048L, 4096L, 4096L, 4096L),
    dictionary_release = c("gc30M", "gc104M", "gc104M", "gc104M"),
    special_tokens = c(FALSE, TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
}

#' List the currently released official Geneformer base checkpoints.
#'
#' @return A data frame describing the four official checkpoints, training
#'   corpus, input size, and Geneformer tokenization version.
#' @export
fmverse_geneformer_checkpoints <- function() {
  fmverse_geneformer_checkpoint_catalog()
}

fmverse_geneformer_checkpoint_aliases <- c(
  "v1-10m" = "Geneformer-V1-10M",
  "geneformer-v1-10m" = "Geneformer-V1-10M",
  "gf-6l-30m-i2048" = "Geneformer-V1-10M",
  "v2-104m" = "Geneformer-V2-104M",
  "geneformer-v2-104m" = "Geneformer-V2-104M",
  "v2-316m" = "Geneformer-V2-316M",
  "geneformer-v2-316m" = "Geneformer-V2-316M",
  "v2-104m-clcancer" = "Geneformer-V2-104M_CLcancer",
  "v2-104m_clcancer" = "Geneformer-V2-104M_CLcancer",
  "geneformer-v2-104m-clcancer" = "Geneformer-V2-104M_CLcancer",
  "geneformer-v2-104m_clcancer" = "Geneformer-V2-104M_CLcancer"
)

fmverse_normalize_geneformer_checkpoint_name <- function(checkpoint) {
  value <- trimws(as.character(checkpoint[[1L]]))
  aliases <- fmverse_geneformer_checkpoint_aliases
  alias <- unname(aliases[tolower(value)])
  if (length(alias) && !is.na(alias)) return(alias)
  canonical <- fmverse_geneformer_checkpoint_catalog()$checkpoint
  canonical_hit <- canonical[tolower(canonical) == tolower(value)]
  if (length(canonical_hit)) return(canonical_hit[[1L]])
  value
}

fmverse_resolve_geneformer_checkpoint <- function(checkpoint,
                                                   checkpoint_dir = NULL) {
  if (is.null(checkpoint) || !length(checkpoint)) {
    checkpoint <- "gf-12L-95M-i4096"
  }
  checkpoint <- path.expand(as.character(checkpoint[[1L]]))
  normalized_name <- fmverse_normalize_geneformer_checkpoint_name(checkpoint)
  is_path <- file.exists(checkpoint) || dir.exists(checkpoint)
  if (file.exists(checkpoint) && !dir.exists(checkpoint)) {
    manifest_root <- !is.null(checkpoint_dir) &&
      file.exists(file.path(checkpoint_dir, "fmverse.json"))
    model_dir <- if (manifest_root) checkpoint_dir else dirname(checkpoint)
    checkpoint_file <- checkpoint
  } else if (dir.exists(checkpoint)) {
    model_dir <- checkpoint
    checkpoint_file <- file.path(model_dir, "model.safetensors")
  } else {
    roots <- c(
      checkpoint_dir,
      getOption("fmverse.geneformer_checkpoint_dir"),
      Sys.getenv("FMVERSE_GENEFORMER_CHECKPOINT_DIR", unset = ""),
      file.path(getwd(), "00Geneformer", "checkpoint"),
      file.path(dirname(getwd()), "00Geneformer", "checkpoint")
    )
    roots <- unique(roots[!is.na(roots) & nzchar(roots)])
    model_dir <- NULL
    for (root in roots) {
      candidates <- c(
        file.path(root, normalized_name), file.path(root, checkpoint),
        file.path(root, "ctheodoris-Geneformer", normalized_name),
        file.path(root, "ctheodoris-Geneformer", checkpoint)
      )
      hit <- candidates[dir.exists(candidates)]
      if (length(hit)) {
        model_dir <- hit[[1L]]
        break
      }
    }
    if (is.null(model_dir)) {
      stop("Geneformer checkpoint `", checkpoint,
           "` was not found. Set checkpoint_dir or FMVERSE_GENEFORMER_CHECKPOINT_DIR.",
           call. = FALSE)
    }
    checkpoint_file <- file.path(model_dir, "model.safetensors")
  }
  if (!file.exists(checkpoint_file)) stop("Missing model.safetensors: ", checkpoint_file)
  manifest_file <- file.path(model_dir, "fmverse.json")
  folder_info <- NULL
  if (file.exists(manifest_file)) {
    folder_info <- .fm_model_manifest(model_dir)
    if (!identical(folder_info$model, "Geneformer")) {
      stop("The selected model folder is not marked as Geneformer.", call. = FALSE)
    }
    model_name <- folder_info$checkpoint
  }
  config_file <- if (is.null(folder_info)) file.path(model_dir, "config.json") else
    folder_info$files$config
  if (!file.exists(config_file)) stop("Missing Geneformer config.json: ", config_file)
  config <- jsonlite::fromJSON(config_file, simplifyVector = TRUE)
  if (is.null(folder_info)) model_name <- if (is_path) basename(model_dir) else normalized_name
  catalog <- fmverse_geneformer_checkpoint_catalog()
  catalog_row <- match(tolower(model_name), tolower(catalog$checkpoint))
  dictionary_release <- if (!is.na(catalog_row)) {
    catalog$dictionary_release[[catalog_row]]
  } else if (grepl("30m|v1", model_name, ignore.case = TRUE) ||
             config$vocab_size > 25000L) {
    "gc30M"
  } else if (grepl("95m", model_name, ignore.case = TRUE)) {
    "gc95M"
  } else {
    "gc104M"
  }
  special_tokens <- identical(dictionary_release, "gc104M")
  dictionary_stem <- switch(dictionary_release,
    gc30M = "gc30M", gc95M = "gc95M", gc104M = "gc104M"
  )
  roots <- c(
    checkpoint_dir,
    getOption("fmverse.geneformer_checkpoint_dir"),
    Sys.getenv("FMVERSE_GENEFORMER_CHECKPOINT_DIR", unset = ""),
    file.path(getwd(), "00Geneformer", "checkpoint"),
    file.path(dirname(getwd()), "00Geneformer", "checkpoint"),
    dirname(model_dir), dirname(dirname(model_dir))
  )
  dictionary_dir <- NULL
  if (!is.null(folder_info)) {
    dictionary_dir <- dirname(folder_info$files$token_dictionary)
  } else {
    for (root in unique(roots[!is.na(roots) & nzchar(roots)])) {
      candidates <- if (identical(dictionary_release, "gc30M")) {
        c(file.path(root, "geneformer", "gene_dictionaries_30m"),
          file.path(root, "gene_dictionaries_30m"))
      } else {
        c(file.path(root, "geneformer"), root)
      }
      token_file <- paste0("token_dictionary_", dictionary_stem, ".tsv")
      hit <- candidates[dir.exists(candidates) &
        file.exists(file.path(candidates, token_file))]
      if (length(hit)) {
        dictionary_dir <- hit[[1L]]
        break
      }
    }
  }
  if (is.null(dictionary_dir)) {
    stop("Geneformer token/rank dictionaries were not found beside the checkpoint.", call. = FALSE)
  }
  list(
    name = model_name,
    model_dir = normalizePath(model_dir, winslash = "/", mustWork = TRUE),
    checkpoint = normalizePath(checkpoint_file, winslash = "/", mustWork = TRUE),
    config = config,
    config_file = normalizePath(config_file, winslash = "/", mustWork = TRUE),
    dictionary_dir = normalizePath(dictionary_dir, winslash = "/", mustWork = TRUE),
    dictionary_release = dictionary_release,
    special_tokens = special_tokens,
    token_dictionary = if (is.null(folder_info)) file.path(dictionary_dir,
      paste0("token_dictionary_", dictionary_stem, ".tsv")) else
      folder_info$files$token_dictionary,
    gene_median_dictionary = if (is.null(folder_info)) file.path(dictionary_dir,
      paste0("gene_median_dictionary_", dictionary_stem, ".tsv")) else
      folder_info$files$gene_median_dictionary,
    gene_name_id_dictionary = if (is.null(folder_info)) file.path(dictionary_dir,
      paste0("gene_name_id_dict_", dictionary_stem, ".tsv")) else
      folder_info$files$gene_name_id_dictionary,
    ensembl_mapping_dictionary = if (is.null(folder_info)) file.path(dictionary_dir,
      paste0("ensembl_mapping_dict_", dictionary_stem, ".tsv")) else
      folder_info$files$ensembl_mapping_dictionary
  )
}

fmverse_geneformer_read_dictionaries <- function(paths, model_vocab_size) {
  token <- utils::read.delim(paths$token_dictionary, check.names = FALSE,
                             stringsAsFactors = FALSE)
  median <- utils::read.delim(paths$gene_median_dictionary, check.names = FALSE,
                              stringsAsFactors = FALSE)
  symbol <- utils::read.delim(paths$gene_name_id_dictionary, check.names = FALSE,
                              stringsAsFactors = FALSE)
  required <- list(token = c("gene_id", "token_id"),
                   median = c("gene_id", "median"),
                   symbol = c("gene_name", "gene_id"))
  for (name in names(required)) {
    if (!all(required[[name]] %in% colnames(get(name)))) {
      stop("Geneformer ", name, " dictionary is missing columns: ",
           paste(setdiff(required[[name]], colnames(get(name))), collapse = ", "))
    }
  }
  token$token_id <- as.integer(token$token_id)
  median$median <- as.numeric(median$median)
  token <- token[!is.na(token$token_id) & token$token_id < model_vocab_size, , drop = FALSE]
  if (!all(c("<pad>", "<mask>") %in% token$gene_id)) {
    stop("Geneformer token dictionary must include <pad> and <mask>.")
  }
  if (isTRUE(paths$special_tokens) && !all(c("<cls>", "<eos>") %in% token$gene_id)) {
    stop("This Geneformer checkpoint requires <cls> and <eos> tokens.")
  }
  token <- token[order(token$token_id), , drop = FALSE]
  if (nrow(token) != model_vocab_size || anyDuplicated(token$token_id) ||
      !identical(token$token_id, seq.int(0L, model_vocab_size - 1L))) {
    stop("Geneformer token dictionary IDs do not match config.json vocab_size.")
  }
  list(token = token, median = median, symbol = symbol, config = paths$config)
}

fmverse_geneformer_vocab <- function(token_dictionary) {
  stats::setNames(
    as.integer(token_dictionary$token_id),
    as.character(token_dictionary$gene_id)
  )
}

fmverse_geneformer_map_features <- function(object, assay, counts, symbol_dictionary,
                                            gene_id = NULL) {
  feature_metadata <- object[[assay]][[]]
  original <- rownames(counts)
  if (!all(original %in% rownames(feature_metadata))) {
    stop("The selected assay feature metadata is not aligned with the raw count rows.")
  }
  feature_metadata <- feature_metadata[original, , drop = FALSE]
  if (is.null(gene_id)) {
    auto_id <- intersect(c("geneformer_ensembl_id", "ensembl_id", "gene_id"),
                         colnames(feature_metadata))
    if (length(auto_id)) gene_id <- auto_id[[1L]]
  }
  if (!is.null(gene_id)) {
    if (length(gene_id) != 1L || !gene_id %in% colnames(feature_metadata)) {
      stop("`gene_id` must name a feature metadata column in the selected assay.")
    }
    ensembl <- as.character(feature_metadata[[gene_id]])
  } else {
    gene_names <- if ("geneformer_gene_name" %in% colnames(feature_metadata)) {
      as.character(feature_metadata$geneformer_gene_name)
    } else if ("scGPT_gene_name" %in% colnames(feature_metadata)) {
      as.character(feature_metadata$scGPT_gene_name)
    } else {
      original
    }
    ensembl <- gene_names
    needs_mapping <- !grepl("^ENSG[0-9]+(?:\\.[0-9]+)?$", ensembl)
    if (any(needs_mapping)) {
      symbol <- symbol_dictionary$gene_name
      ids <- symbol_dictionary$gene_id
      names(ids) <- symbol
      symbols <- original[needs_mapping]
      mapped <- unname(ids[match(symbols, names(ids))])
      # Gene symbols in 10x feature tables are case sensitive in official
      # Geneformer mapping; this fallback handles common uppercase aliases.
      absent <- is.na(mapped)
      if (any(absent)) {
        symbol_upper <- toupper(symbol)
        ids_upper <- ids[!duplicated(symbol_upper)]
        names(ids_upper) <- symbol_upper[!duplicated(symbol_upper)]
        mapped[absent] <- unname(ids_upper[match(toupper(symbols[absent]), names(ids_upper))])
      }
      ensembl[needs_mapping] <- mapped
    }
  }
  ensembl <- sub("\\.[0-9]+$", "", ensembl)
  ensembl[!grepl("^ENSG[0-9]+$", ensembl)] <- NA_character_
  ensembl
}

fmverse_geneformer_tokenize <- function(object, dictionaries, assay = NULL, cells = NULL,
                                         max_length = 4096L, gene_id = NULL) {
  object <- fmverse_to_seurat_v5(object)
  if (is.null(assay)) assay <- SeuratObject::DefaultAssay(object)
  if (!assay %in% SeuratObject::Assays(object)) stop("Assay not found: ", assay)
  if (!"counts" %in% SeuratObject::Layers(object[[assay]])) {
    stop("Assay ", assay, " has no raw counts layer.")
  }
  counts <- methods::as(SeuratObject::LayerData(object, assay = assay, layer = "counts"),
                        "dgCMatrix")
  if (is.null(cells)) {
    cell_names <- colnames(counts)
  } else {
    cell_names <- as.character(cells)
    if (!length(cell_names) || anyDuplicated(cell_names) ||
        any(!cell_names %in% colnames(counts))) {
      stop("`cells` must contain unique barcodes present in the Seurat object.")
    }
    counts <- counts[, cell_names, drop = FALSE]
  }
  raw_cell_totals <- Matrix::colSums(counts)
  ensembl <- fmverse_geneformer_map_features(
    object, assay, counts, dictionaries$symbol, gene_id = gene_id
  )
  valid <- !is.na(ensembl) & ensembl %in% dictionaries$token$gene_id &
    ensembl %in% dictionaries$median$gene_id
  if (!any(valid)) {
    stop("No features map to Geneformer's selected token/rank dictionaries.")
  }
  mapped_features <- which(valid)
  unique_ids <- unique(ensembl[mapped_features])
  feature_to_gene <- match(ensembl[mapped_features], unique_ids)
  collapse <- Matrix::sparseMatrix(
    i = feature_to_gene, j = seq_along(mapped_features), x = 1,
    dims = c(length(unique_ids), length(mapped_features))
  )
  gene_counts <- collapse %*% counts[mapped_features, , drop = FALSE]
  rownames(gene_counts) <- unique_ids
  colnames(gene_counts) <- cell_names
  token_index <- match(unique_ids, dictionaries$token$gene_id)
  token_ids <- dictionaries$token$token_id[token_index]
  median_index <- match(unique_ids, dictionaries$median$gene_id)
  medians <- dictionaries$median$median[median_index]
  if (any(!is.finite(medians) | medians <= 0)) {
    stop("Geneformer median dictionary contains a non-positive or non-finite value.")
  }
  max_length <- as.integer(max_length)
  special_tokens <- isTRUE(dictionaries$special_tokens)
  token_capacity <- max_length - if (special_tokens) 2L else 0L
  if (length(max_length) != 1L || is.na(max_length) || max_length < 1L ||
      max_length > dictionaries$config$max_position_embeddings || token_capacity < 1L) {
    stop("`max_length` must be between 1 and the checkpoint's max_position_embeddings.")
  }
  pointers <- gene_counts@p
  cell_tokens <- vector("list", length(cell_names))
  n_truncated <- 0L
  token_counts <- integer(length(cell_names))
  for (j in seq_along(cell_names)) {
    first <- pointers[[j]] + 1L
    last <- pointers[[j + 1L]]
    if (first > last || !is.finite(raw_cell_totals[[j]]) || raw_cell_totals[[j]] <= 0) {
      stop("Cell `", cell_names[[j]], "` has no counts for Geneformer tokenization.")
    }
    offsets <- seq.int(first, last)
    values <- gene_counts@x[offsets]
    nonzero <- values > 0
    offsets <- offsets[nonzero]
    values <- values[nonzero]
    rows <- gene_counts@i[offsets] + 1L
    normalized <- values / raw_cell_totals[[j]] * 10000 / medians[rows]
    order_index <- order(-normalized, rows, method = "radix")
    if (length(order_index) > token_capacity) {
      order_index <- order_index[seq_len(token_capacity)]
      n_truncated <- n_truncated + 1L
    }
    token_counts[[j]] <- length(order_index)
    ranked_ids <- as.integer(token_ids[rows[order_index]])
    if (special_tokens) {
      cls_id <- dictionaries$token$token_id[match("<cls>", dictionaries$token$gene_id)]
      eos_id <- dictionaries$token$token_id[match("<eos>", dictionaries$token$gene_id)]
      if (is.na(cls_id) || is.na(eos_id)) {
        stop("The selected Geneformer dictionary is missing required <cls>/<eos> tokens.")
      }
      cell_tokens[[j]] <- c(as.integer(cls_id), ranked_ids, as.integer(eos_id))
    } else {
      cell_tokens[[j]] <- ranked_ids
    }
    token_counts[[j]] <- length(cell_tokens[[j]])
  }
  width <- max(token_counts)
  gene_matrix <- matrix(as.integer(dictionaries$token$token_id[
    match("<pad>", dictionaries$token$gene_id)
  ]), nrow = length(cell_names), ncol = width)
  value_matrix <- matrix(0, nrow = length(cell_names), ncol = width)
  for (j in seq_along(cell_tokens)) {
    n <- length(cell_tokens[[j]])
    gene_matrix[j, seq_len(n)] <- cell_tokens[[j]]
    value_matrix[j, seq_len(n)] <- seq_len(n)
    if (special_tokens) {
      value_matrix[j, 1L] <- 0
      value_matrix[j, n] <- 0
    }
  }
  pad_id <- dictionaries$token$token_id[match("<pad>", dictionaries$token$gene_id)]
  attention_mask <- gene_matrix != pad_id
  vocab <- fmverse_geneformer_vocab(dictionaries$token)
  symbol_by_ensembl <- dictionaries$symbol$gene_name[
    match(unique_ids, dictionaries$symbol$gene_id)
  ]
  names(symbol_by_ensembl) <- unique_ids
  token_to_symbol <- setNames(dictionaries$token$gene_id,
                              as.character(dictionaries$token$token_id))
  official_symbol_by_gene <- dictionaries$symbol$gene_name[
    match(dictionaries$token$gene_id, dictionaries$symbol$gene_id)
  ]
  use_symbols <- !is.na(official_symbol_by_gene)
  token_to_symbol[as.character(dictionaries$token$token_id[use_symbols])] <-
    official_symbol_by_gene[use_symbols]
  input_gene_labels <- stats::setNames(rownames(counts), ensembl)
  input_gene_labels <- input_gene_labels[!is.na(names(input_gene_labels))]
  input_gene_labels <- input_gene_labels[!duplicated(names(input_gene_labels))]
  input_token <- token_to_symbol
  input_token[as.character(token_ids)] <- unname(input_gene_labels[unique_ids])
  token_to_symbol[as.character(token_ids)] <- ifelse(
    is.na(input_token[as.character(token_ids)]) | !nzchar(input_token[as.character(token_ids)]),
    token_to_symbol[as.character(token_ids)], input_token[as.character(token_ids)]
  )
  symbol_to_token <- dictionaries$token$token_id[
    match(dictionaries$symbol$gene_id, dictionaries$token$gene_id)
  ]
  names(symbol_to_token) <- dictionaries$symbol$gene_name
  valid_symbol <- !is.na(names(symbol_to_token)) & !is.na(symbol_to_token) &
    symbol_to_token < dictionaries$config$vocab_size
  symbol_to_token <- symbol_to_token[valid_symbol]
  symbol_to_token <- symbol_to_token[!duplicated(names(symbol_to_token))]
  list(
    gene_ids = gene_matrix, values = value_matrix,
    pad_mask = !attention_mask, attention_mask = attention_mask,
    cell_names = cell_names, vocab = vocab, gene_symbols = token_to_symbol,
    symbol_to_token = symbol_to_token,
    assay = assay, token_counts = token_counts,
    matched_features = length(mapped_features), total_features = nrow(counts),
    mapped_ensembl_ids = unique_ids, truncated_cells = n_truncated,
    preprocessing = list(
      method = "Geneformer TranscriptomeTokenizer rank-value encoding",
      dictionary_release = dictionaries$dictionary_release,
      normalization = paste0("raw_counts / total_cell_counts * 10000 / ",
                             dictionaries$dictionary_release, " gene median"),
      ordering = "descending normalized expression; deterministic feature-order tie break for equal ranks",
      special_tokens = special_tokens, padding = "right with <pad>", max_length = max_length
    )
  )
}

#' Embed a Seurat object with an official Geneformer checkpoint.
#' @param object Seurat v4/v5 object.
#' @param checkpoint Official model name or local model directory/file.
#' @param checkpoint_dir Root containing model and gene dictionaries.
#' @param batch_size Cells per forward pass.
#' @param max_length Maximum rank-ordered gene tokens.
#' @param cells Optional cell barcodes to embed.
#' @param assay Counts assay.
#' @param gene_id Optional assay feature-metadata field containing Ensembl IDs.
#' @param device Torch device.
#' @param verbose Report progress.
#' @return An `fmobj` containing a Geneformer run, tokenized inputs, and embedding.
#' @export
.fmverse_geneformer_legacy <- function(object, checkpoint = "gf-12L-95M-i4096",
                               checkpoint_dir = NULL, batch_size = 8L,
                               max_length = NULL, cells = NULL, assay = NULL,
                               gene_id = NULL, device = "cpu", verbose = interactive()) {
  if (length(batch_size) != 1L || is.na(batch_size) || batch_size < 1L ||
      batch_size != as.integer(batch_size)) {
    stop("`batch_size` must be a positive integer.")
  }
  resolved <- fmverse_resolve_geneformer_checkpoint(checkpoint, checkpoint_dir)
  config <- resolved$config
  if (is.null(max_length)) max_length <- config$max_position_embeddings
  dictionary_paths <- resolved
  dictionary_paths$config <- config
  dictionaries <- fmverse_geneformer_read_dictionaries(dictionary_paths, config$vocab_size)
  dictionaries$dictionary_release <- resolved$dictionary_release
  dictionaries$special_tokens <- resolved$special_tokens
  cls_id <- dictionaries$token$token_id[match("<cls>", dictionaries$token$gene_id)]
  eos_id <- dictionaries$token$token_id[match("<eos>", dictionaries$token$gene_id)]
  config$cls_token_id <- if (length(cls_id) && !is.na(cls_id)) as.integer(cls_id) else -1L
  config$eos_token_id <- if (length(eos_id) && !is.na(eos_id)) as.integer(eos_id) else -1L
  tokenized <- fmverse_geneformer_tokenize(
    object, dictionaries, assay = assay, cells = cells,
    max_length = max_length, gene_id = gene_id
  )
  model <- fmverse_geneformer_build_model(config)
  loaded <- fmverse_load_geneformer_weights(model, resolved$checkpoint)
  model <- loaded$model
  torch_device <- torch::torch_device(device)
  model$to(device = torch_device)
  model$eval()
  n_cells <- length(tokenized$cell_names)
  embeddings <- matrix(NA_real_, nrow = n_cells, ncol = config$hidden_size)
  starts <- seq.int(1L, n_cells, by = as.integer(batch_size))
  progress_step <- max(1L, ceiling(length(starts) / 20L))
  for (batch_number in seq_along(starts)) {
    first <- starts[[batch_number]]
    index <- first:min(first + as.integer(batch_size) - 1L, n_cells)
    sequence_end <- max(tokenized$token_counts[index])
    slice <- seq_len(sequence_end)
    input_ids <- torch::torch_tensor(
      tokenized$gene_ids[index, slice, drop = FALSE], dtype = torch::torch_long(),
      device = torch_device
    )
    attention_mask <- torch::torch_tensor(
      tokenized$attention_mask[index, slice, drop = FALSE], dtype = torch::torch_long(),
      device = torch_device
    )
    output <- torch::with_no_grad(model(
      input_ids = input_ids, attention_mask = attention_mask,
      embedding_layer = -1L
    ))
    values <- torch::as_array(output$cell_embedding$cpu())
    embeddings[index, ] <- matrix(values, nrow = length(index), ncol = config$hidden_size)
    if (isTRUE(verbose) && (batch_number %% progress_step == 0L || batch_number == length(starts))) {
      message(sprintf("Geneformer embedded %d/%d cells", max(index), n_cells))
    }
    rm(output, input_ids, attention_mask, values)
    if (batch_number %% 10L == 0L) invisible(gc(verbose = FALSE))
  }
  rownames(embeddings) <- tokenized$cell_names
  colnames(embeddings) <- paste0("Geneformer_", seq_len(ncol(embeddings)))
  object <- fmverse_to_seurat_v5(object)
  if (!is.null(cells)) object <- subset(object, cells = tokenized$cell_names)
  object[["Geneformer"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embeddings, key = "Geneformer_", assay = tokenized$assay
  )
  structure(
    list(
      model_name = "Geneformer", embeddings = embeddings, seurat = object,
      model = model, audit = loaded$audit,
      inputs = list(
        gene_ids = tokenized$gene_ids, values = tokenized$values,
        pad_mask = tokenized$pad_mask, attention_mask = tokenized$attention_mask,
        cell_names = tokenized$cell_names, vocab = tokenized$vocab,
        gene_symbols = tokenized$gene_symbols,
        symbol_to_token = tokenized$symbol_to_token, assay = tokenized$assay,
        token_counts = tokenized$token_counts
      ),
      analyses = list(), plots = list(), history = list(),
      metadata = list(
        checkpoint_name = resolved$name, checkpoint = resolved$checkpoint,
        config = config, config_file = resolved$config_file,
        dictionary_dir = resolved$dictionary_dir,
        dictionary_release = resolved$dictionary_release, device = device,
        batch_size = as.integer(batch_size), max_length = as.integer(max_length),
        embedding_layer = -1L, pooling = "mean over non-padding token states",
        matched_features = tokenized$matched_features,
        total_features = tokenized$total_features,
        truncated_cells = tokenized$truncated_cells,
        n_cells = n_cells, token_counts = tokenized$token_counts,
        preprocessing = tokenized$preprocessing
      )
    ),
    class = c("fmverse_result", "list")
  )
}
