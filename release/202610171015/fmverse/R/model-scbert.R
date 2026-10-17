# R-native scBERT PerformerLM and strict checkpoint adapter.

fmverse_scbert_model <- torch::nn_module(
  classname = "fmverse_scbert_model",
  initialize = function(num_tokens = 7L, max_seq_len = 16907L, dim = 200L,
                        depth = 6L, heads = 10L, dim_head = 64L,
                        nb_features = NULL) {
    self$max_seq_len <- as.integer(max_seq_len)
    self$token_emb <- torch::nn_embedding(num_tokens, dim)
    self$pos_emb <- torch::nn_module(
      classname = "fmverse_scbert_gene_position_embedding",
      initialize = function(sequence_length, embedding_dim) {
        self$emb <- torch::nn_embedding(sequence_length, embedding_dim)
      },
      forward = function(input_ids) {
        # R torch's Embedding API uses one-based indices; row 1 matches Python ID 0.
        pos <- torch::torch_arange(1L, input_ids$size(2L),
          dtype = torch::torch_long(), device = input_ids$device)
        self$emb(pos)
      }
    )(max_seq_len, dim)
    self$dropout <- torch::nn_dropout(0)
    self$performer <- fmverse_performer_core(
      dim, depth, heads, dim_head, nb_features, ff_mult = 4L,
      attn_dropout = 0, ff_dropout = 0, qkv_bias = FALSE
    )
    self$norm <- torch::nn_layer_norm(dim)
    self$to_out <- torch::nn_linear(dim, num_tokens)
  },
  forward = function(input_ids, return_encodings = FALSE) {
    if (input_ids$size(2L) > self$max_seq_len) stop("scBERT input exceeds max_seq_len.")
    hidden <- self$token_emb(input_ids + 1L) + self$pos_emb(input_ids)
    hidden <- self$performer(self$dropout(hidden))
    hidden <- self$norm(hidden)
    if (return_encodings) hidden else self$to_out(hidden)
  }
)

#' Construct the R-native scBERT Performer model
#' @return The fixed 200-dimensional, six-layer scBERT model module.
#' @export
fmverse_scbert_build_model <- function() fmverse_scbert_model()

fmverse_scbert_keymap <- function(checkpoint_names, model_state) {
  target <- names(model_state)
  mapped <- stats::setNames(rep(NA_character_, length(checkpoint_names)), checkpoint_names)
  for (key in checkpoint_names) {
    candidates <- unique(c(key, sub("^module\\.", "", key)))
    hit <- intersect(candidates, target)
    if (length(hit) == 1L) mapped[[key]] <- hit
  }
  mapped
}

#' Build and strictly load the official scBERT Performer checkpoint.
#' @param model R-native scBERT module.
#' @param path Python-exported flat state dict or safetensors file.
#' @return Loaded module and a complete tensor-name/shape audit.
#' @export
fmverse_load_scbert_weights <- function(model, path) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop("scBERT weight file does not exist: ", path, call. = FALSE)
  }
  weights <- if (grepl("[.]safetensors$", path, ignore.case = TRUE)) {
    safetensors::safe_load_file(path, framework = "torch")
  } else {
    torch::load_state_dict(path)
  }
  if (!is.list(weights) || is.null(names(weights))) {
    stop("scBERT checkpoint must be a flat named tensor state dictionary.", call. = FALSE)
  }
  # PyTorch serializes this inference bookkeeping buffer as a rank-0 scalar;
  # torch for R represents it as a one-element tensor. It does not participate
  # in the fixed-projection evaluation path, but reshape it for strict loading.
  scalar_buffer <- "performer.calls_since_last_redraw"
  if (scalar_buffer %in% names(weights) &&
      length(weights[[scalar_buffer]]$size()) == 0L) {
    weights[[scalar_buffer]] <- weights[[scalar_buffer]]$view(c(1L))
  }
  state <- model$state_dict()
  keymap <- fmverse_scbert_keymap(names(weights), state)
  unmapped <- names(keymap)[is.na(keymap)]
  targets <- unname(keymap[!is.na(keymap)])
  missing <- setdiff(names(state), targets)
  duplicate <- unique(targets[duplicated(targets)])
  if (length(unmapped) || length(missing) || length(duplicate)) {
    stop("scBERT state-dict audit failed. Unmapped checkpoint keys: ",
      paste(unmapped, collapse = ", "), "; missing R keys: ",
      paste(missing, collapse = ", "), "; duplicate targets: ",
      paste(duplicate, collapse = ", "), call. = FALSE)
  }
  shape <- function(x) paste(as.integer(x$size()), collapse = "x")
  mismatch <- names(weights)[vapply(names(weights), function(key) {
    !identical(shape(weights[[key]]), shape(state[[keymap[[key]]]]))
  }, logical(1))]
  if (length(mismatch)) {
    desc <- vapply(mismatch, function(key) paste0(key, " (checkpoint ",
      shape(weights[[key]]), ", R ", shape(state[[keymap[[key]]]]), ")"), character(1))
    stop("scBERT tensor shape mismatch: ", paste(desc, collapse = "; "), call. = FALSE)
  }
  mapped <- stats::setNames(vector("list", length(state)), names(state))
  for (key in names(weights)) mapped[[keymap[[key]]]] <- weights[[key]]
  model$load_state_dict(mapped, strict = TRUE)
  list(model = model, audit = list(
    checkpoint = normalizePath(path, winslash = "/", mustWork = TRUE),
    checkpoint_key_count = length(weights), model_key_count = length(state),
    matched_key_count = length(weights), missing_model_keys = character(),
    unmapped_checkpoint_keys = character(), shape_mismatch_count = 0L,
    shape_adaptations = if (scalar_buffer %in% names(weights))
      "performer.calls_since_last_redraw: rank-0 scalar to one-element R tensor" else character(),
    key_mapping = data.frame(checkpoint_key = names(weights),
      model_key = unname(keymap), shape = vapply(weights, shape, character(1)),
      stringsAsFactors = FALSE)
  ))
}
