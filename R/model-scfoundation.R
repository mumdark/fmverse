# R-native xTrimoGene/scFoundation cell checkpoint architecture.

fmverse_scfoundation_auto_bin <- torch::nn_module(
  classname = "fmverse_scfoundation_auto_bin",
  initialize = function(dim, bin_num, bin_alpha, mask_token_id, pad_token_id) {
    self$dim <- as.integer(dim)
    self$bin_num <- as.integer(bin_num)
    self$bin_alpha <- as.numeric(bin_alpha)
    self$mlp <- torch::nn_linear(1L, bin_num)
    self$mlp2 <- torch::nn_linear(bin_num, bin_num)
    self$LeakyReLU <- torch::nn_leaky_relu(negative_slope = 0.1)
    self$emb <- torch::nn_embedding(bin_num, dim)
    self$emb_mask <- torch::nn_embedding(1L, dim)
    self$emb_pad <- torch::nn_embedding(1L, dim)
    self$mask_token_id <- as.numeric(mask_token_id)
    self$pad_token_id <- as.numeric(pad_token_id)
  },
  forward = function(x) {
    values <- x$squeeze(3L)
    mask_flag <- (values == self$mask_token_id)$unsqueeze(3L)$to(
      dtype = torch::torch_float()
    )
    pad_flag <- (values == self$pad_token_id)$unsqueeze(3L)$to(
      dtype = torch::torch_float()
    )
    x <- self$mlp(x)
    x <- self$LeakyReLU(x)
    x_cross <- self$mlp2(x)
    weight <- (self$bin_alpha * x + x_cross)$softmax(dim = 3L)
    # R torch uses one-based Embedding indices; this retrieves Python rows 0..n-1.
    bin_index <- torch::torch_arange(1L, self$bin_num,
      dtype = torch::torch_long(), device = x$device)
    embedded <- torch::torch_matmul(weight, self$emb(bin_index))
    first <- torch::torch_ones(1L, dtype = torch::torch_long(), device = x$device)
    mask_embedding <- self$emb_mask(first)$view(c(1L, 1L, self$dim))
    pad_embedding <- self$emb_pad(first)$view(c(1L, 1L, self$dim))
    embedded <- embedded * (1 - mask_flag - pad_flag) +
      mask_embedding * mask_flag + pad_embedding * pad_flag
    embedded
  }
)

fmverse_scfoundation_cell_module <- torch::nn_module(
  classname = "fmverse_scfoundation_cell_module",
  initialize = function(config) {
    seq_len <- as.integer(config$seq_len)
    encoder_config <- config$encoder
    decoder_config <- config$decoder
    self$seq_len <- seq_len
    self$num_tokens <- as.integer(config$n_class)
    self$pad_token_id <- as.integer(config$pad_token_id)
    self$mask_token_id <- as.integer(config$mask_token_id)
    self$token_emb <- fmverse_scfoundation_auto_bin(
      encoder_config$hidden_dim, config$bin_num, config$bin_alpha,
      config$mask_token_id, config$pad_token_id
    )
    self$pos_emb <- torch::nn_embedding(seq_len + 1L, encoder_config$hidden_dim)
    self$encoder <- fmverse_scfoundation_module(
      encoder_config, seq_len, encoder = TRUE
    )
    self$decoder <- fmverse_scfoundation_module(
      decoder_config, seq_len, encoder = FALSE
    )
    self$decoder_embed <- torch::nn_linear(
      encoder_config$hidden_dim, decoder_config$hidden_dim, bias = TRUE
    )
    self$norm <- torch::nn_layer_norm(decoder_config$hidden_dim)
    self$to_final <- torch::nn_linear(decoder_config$hidden_dim, 1L)
  },
  forward = function(values, gene_positions) {
    x <- self$token_emb(values$unsqueeze(3L))
    x <- x + self$pos_emb(gene_positions + 1L)
    self$encoder(x)
  }
)

fmverse_scfoundation_module <- function(config, seq_len, encoder = TRUE) {
  module_type <- tolower(as.character(config$module_type))
  if (module_type == "performer") {
    fmverse_scfoundation_performer(
      max_seq_len = seq_len, dim = as.integer(config$hidden_dim),
      depth = as.integer(config$depth), heads = as.integer(config$heads),
      dim_head = as.integer(config$dim_head %||% 64L),
      nb_features = if (is.null(config$nb_features)) NULL else as.integer(config$nb_features),
      ff_mult = as.integer(config$ff_mult %||% 4L),
      ff_dropout = as.numeric(config$ff_dropout %||% 0),
      attn_dropout = as.numeric(config$attn_dropout %||% 0)
    )
  } else if (module_type == "transformer") {
    fmverse_scfoundation_transformer(
      max_seq_len = seq_len, dim = as.integer(config$hidden_dim),
      depth = as.integer(config$depth), heads = as.integer(config$heads),
      dropout = as.numeric(config$dropout %||% 0.1)
    )
  } else {
    stop("Unsupported scFoundation module type: ", module_type, call. = FALSE)
  }
}

#' Construct the R-native scFoundation cell model from its exported config.
#' @param config JSON-decoded scFoundation `cell` configuration.
#' @return An R torch module containing the encoder and checkpoint decoder.
#' @export
fmverse_scfoundation_build_model <- function(config) {
  required <- c("seq_len", "n_class", "pad_token_id", "mask_token_id", "bin_num",
                "bin_alpha", "encoder", "decoder")
  absent <- setdiff(required, names(config))
  if (length(absent)) stop("scFoundation config is missing: ", paste(absent, collapse = ", "))
  fmverse_scfoundation_cell_module(config)
}

fmverse_scfoundation_keymap <- function(checkpoint_names, model_state) {
  target <- names(model_state)
  mapped <- stats::setNames(rep(NA_character_, length(checkpoint_names)), checkpoint_names)
  for (key in checkpoint_names) {
    candidates <- unique(c(key, sub("^module\\.", "", key)))
    # torch for R's TransformerEncoder includes an explicit `layers` component;
    # the official xTrimoGene Python wrapper stores its ModuleList directly.
    candidates <- c(candidates, sub("(transformer_encoder)[.]([0-9]+[.])",
      "\\1.layers.\\2", candidates))
    hit <- intersect(unique(candidates), target)
    if (length(hit) == 1L) mapped[[key]] <- hit
  }
  mapped
}

#' Strictly load an exported scFoundation `cell` tensor state dictionary.
#' @param model Model from `fmverse_scfoundation_build_model()`.
#' @param path Flat Python state dict or safetensors produced from the official
#'   nested MMF checkpoint.
#' @return The loaded module and tensor audit.
#' @export
fmverse_load_scfoundation_weights <- function(model, path) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop("scFoundation weight file does not exist: ", path, call. = FALSE)
  }
  weights <- if (grepl("[.]safetensors$", path, ignore.case = TRUE)) {
    safetensors::safe_load_file(path, framework = "torch")
  } else {
    torch::load_state_dict(path)
  }
  if (!is.list(weights) || is.null(names(weights))) {
    stop("scFoundation weights must be a flat named tensor state dictionary.", call. = FALSE)
  }
  # PyTorch stores this Performer bookkeeping buffer as a rank-0 scalar,
  # while torch for R represents it as a one-element tensor. It is not used
  # by the fixed-projection evaluation path, but needs a shape-compatible view
  # for strict state-dict loading.
  scalar_buffer <- "decoder.performer.calls_since_last_redraw"
  if (scalar_buffer %in% names(weights) &&
      length(weights[[scalar_buffer]]$size()) == 0L) {
    weights[[scalar_buffer]] <- weights[[scalar_buffer]]$view(c(1L))
  }
  state <- model$state_dict()
  keymap <- fmverse_scfoundation_keymap(names(weights), state)
  unmapped <- names(keymap)[is.na(keymap)]
  targets <- unname(keymap[!is.na(keymap)])
  missing <- setdiff(names(state), targets)
  duplicates <- unique(targets[duplicated(targets)])
  if (length(unmapped) || length(missing) || length(duplicates)) {
    stop("scFoundation state-dict audit failed. Unmapped checkpoint keys: ",
      paste(unmapped, collapse = ", "), "; missing R keys: ",
      paste(missing, collapse = ", "), "; duplicate targets: ",
      paste(duplicates, collapse = ", "), call. = FALSE)
  }
  shape <- function(x) paste(as.integer(x$size()), collapse = "x")
  mismatch <- names(weights)[vapply(names(weights), function(key) {
    !identical(shape(weights[[key]]), shape(state[[keymap[[key]]]]))
  }, logical(1))]
  if (length(mismatch)) stop("scFoundation tensor shapes differ for: ",
    paste(mismatch, collapse = ", "), call. = FALSE)
  mapped <- stats::setNames(vector("list", length(state)), names(state))
  for (key in names(weights)) mapped[[keymap[[key]]]] <- weights[[key]]
  model$load_state_dict(mapped, strict = TRUE)
  list(model = model, audit = list(
    checkpoint = normalizePath(path, winslash = "/", mustWork = TRUE),
    checkpoint_key_count = length(weights), model_key_count = length(state),
    matched_key_count = length(weights), missing_model_keys = character(),
    unmapped_checkpoint_keys = character(), shape_mismatch_count = 0L,
    shape_adaptations = if (scalar_buffer %in% names(weights))
      paste0(scalar_buffer, ": rank-0 scalar to one-element R tensor") else character(),
    key_mapping = data.frame(checkpoint_key = names(weights), model_key = unname(keymap),
      shape = vapply(weights, shape, character(1)), stringsAsFactors = FALSE)
  ))
}

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x
