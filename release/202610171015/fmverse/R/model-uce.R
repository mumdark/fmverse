# R-native UCE architecture. Parameter and buffer names mirror the official
# snap-stanford/UCE TransformerModel so its .torch state dictionaries load
# without conversion.

fmverse_uce_positional_encoding <- torch::nn_module(
  classname = "fmverse_uce_positional_encoding",
  initialize = function(d_model = 1280L, dropout = 0.05, max_len = 1536L) {
    self$dropout <- torch::nn_dropout(dropout)
    # The checkpoint contains the official sinusoidal position table. A
    # persistent buffer of the right shape is enough here; the checkpoint
    # loader replaces it before inference.
    self$pe <- torch::nn_buffer(
      torch::torch_zeros(c(as.integer(max_len), 1L, as.integer(d_model)))
    )
  },
  forward = function(x) {
    # The selected rows already have shape (sequence, 1, dimension), which
    # broadcasts over the cell/batch axis in the official (sequence,batch,dim)
    # layout.
    x <- x + self$pe[seq_len(x$size(1L)), , ]
    self$dropout(x)
  }
)

fmverse_uce_full_block <- function(in_features, out_features, dropout = 0.05) {
  torch::nn_sequential(
    torch::nn_linear(in_features, out_features, bias = TRUE),
    torch::nn_layer_norm(out_features),
    torch::nn_gelu(),
    torch::nn_dropout(p = dropout)
  )
}

fmverse_uce_model <- function(nlayers = 4L, token_dim = 5120L,
                              d_model = 1280L, nhead = 20L,
                              d_hid = 5120L, output_dim = 1280L,
                              dropout = 0.05, max_len = 1536L,
                              n_tokens = 145469L) {
  if (length(nlayers) != 1L || !nlayers %in% c(4L, 33L)) {
    stop("UCE supports the official 4- and 33-layer architectures.", call. = FALSE)
  }
  fmverse_uce_module(
    nlayers = as.integer(nlayers), token_dim = as.integer(token_dim),
    d_model = as.integer(d_model), nhead = as.integer(nhead),
    d_hid = as.integer(d_hid), output_dim = as.integer(output_dim),
    dropout = dropout, max_len = as.integer(max_len),
    n_tokens = as.integer(n_tokens)
  )
}

fmverse_uce_module <- torch::nn_module(
  classname = "fmverse_uce_module",
  initialize = function(nlayers, token_dim, d_model, nhead, d_hid,
                        output_dim, dropout, max_len, n_tokens) {
    self$pos_encoder <- fmverse_uce_positional_encoding(
      d_model = d_model, dropout = dropout, max_len = max_len
    )
    self$d_model <- d_model
    self$encoder <- torch::nn_sequential(
      torch::nn_linear(token_dim, d_model),
      torch::nn_gelu(),
      torch::nn_layer_norm(d_model)
    )
    layer <- torch::nn_transformer_encoder_layer(
      d_model = d_model, nhead = nhead, dim_feedforward = d_hid,
      # The upstream TransformerEncoderLayer call leaves activation at its
      # PyTorch default, ReLU. GELU is used only in the separate MLP blocks.
      dropout = dropout, activation = "relu", batch_first = FALSE,
      norm_first = FALSE, layer_norm_eps = 1e-5
    )
    self$transformer_encoder <- torch::nn_transformer_encoder(layer, nlayers)
    self$decoder <- torch::nn_sequential(
      fmverse_uce_full_block(d_model, 1024L, dropout),
      fmverse_uce_full_block(1024L, output_dim, dropout),
      fmverse_uce_full_block(output_dim, output_dim, dropout),
      torch::nn_linear(output_dim, output_dim)
    )
    self$binary_decoder <- torch::nn_sequential(
      fmverse_uce_full_block(output_dim + 1280L, 2048L, dropout),
      fmverse_uce_full_block(2048L, 512L, dropout),
      fmverse_uce_full_block(512L, 128L, dropout),
      torch::nn_linear(128L, 1L)
    )
    self$gene_embedding_layer <- torch::nn_sequential(
      torch::nn_linear(token_dim, d_model),
      torch::nn_gelu(),
      torch::nn_layer_norm(d_model)
    )
    # UCE stores its frozen ESM2/protein and special-token rows as a regular
    # embedding in the released checkpoint.
    self$pe_embedding <- torch::nn_embedding(n_tokens, token_dim)
    self$nlayers <- as.integer(nlayers)
    self$output_dim <- as.integer(output_dim)
  },
  forward = function(token_ids, attention_mask, return_layers = FALSE) {
    # Official Python indices are zero-based; R torch embedding indices are
    # one-based. The stored rows themselves remain in their original order.
    x <- self$pe_embedding(token_ids + 1L)
    raw_token_embedding <- x
    # The official inference runner L2-normalizes each 5,120-dimensional
    # protein/special-token vector before the learned encoder. Keep this at
    # the model boundary so every R entry point follows the same operation.
    token_norm <- torch::torch_sqrt(torch::torch_sum(x * x, dim = 3L, keepdim = TRUE))
    x <- x / torch::torch_clamp(token_norm, min = 1e-12)
    token_embedding <- x
    x <- self$encoder(x) * sqrt(self$d_model)
    encoder_output <- x
    x <- self$pos_encoder(x)
    transformer_input <- x
    padding_mask <- 1 - attention_mask
    if (return_layers) layer_outputs <- vector("list", self$nlayers)
    for (i in seq_len(self$nlayers)) {
      x <- self$transformer_encoder$layers[[i]](
        x, src_key_padding_mask = padding_mask
      )
      if (return_layers) layer_outputs[[i]] <- x
    }
    gene_output <- self$decoder(x)
    # UCE defines the cell vector as the decoded CLS token, L2-normalized.
    cell_embedding <- gene_output[1L, , ]
    cell_embedding <- cell_embedding / torch::torch_clamp(
      torch::torch_sqrt(torch::torch_sum(cell_embedding * cell_embedding, dim = 2L,
                                         keepdim = TRUE)), min = 1e-12
    )
    result <- list(gene_output = gene_output, cell_embedding = cell_embedding)
    if (return_layers) {
      names(layer_outputs) <- sprintf("layer_%02d", 0:(self$nlayers - 1L))
      result$raw_token_embedding <- raw_token_embedding
      result$token_embedding <- token_embedding
      result$encoder_output <- encoder_output
      result$transformer_input <- transformer_input
      result$layers <- layer_outputs
    }
    result
  }
)

#' Construct the R-native UCE model architecture
#'
#' @param nlayers Number of Transformer layers in an official UCE checkpoint.
#' @return A torch \code{nn_module} compatible with the official UCE state dict.
#' @export
fmverse_uce_build_model <- function(nlayers = 4L) {
  fmverse_uce_model(nlayers = nlayers)
}

#' Load and strictly audit an official UCE checkpoint
#'
#' @param model An \code{fmverse_uce_build_model()} instance.
#' @param path Path to an official \code{.torch} model file or CPU safetensors copy.
#' @return The loaded model and a tensor-name/shape audit.
#' @export
fmverse_load_uce_weights <- function(model, path) {
  if (!file.exists(path)) stop("UCE checkpoint does not exist: ", path, call. = FALSE)
  safe_path <- if (grepl("\\.safetensors$", path, ignore.case = TRUE)) path else {
    sub("\\.torch$", ".cpu.safetensors", path, ignore.case = TRUE)
  }
  if (file.exists(safe_path)) {
    weights <- safetensors::safe_load_file(safe_path, framework = "torch")
    weight_source <- safe_path
  } else {
    weights <- tryCatch(torch::load_state_dict(path), error = function(error) {
      stop("R torch could not read the UCE .torch checkpoint. Convert the official state ",
        "dictionary to an adjacent .safetensors file (see 00UCE/prepare_assets.py). ",
        "Original error: ", conditionMessage(error), call. = FALSE)
    })
    weight_source <- path
  }
  state <- model$state_dict()
  checkpoint_keys <- names(weights)
  missing <- setdiff(names(state), checkpoint_keys)
  unexpected <- setdiff(checkpoint_keys, names(state))
  shape_of <- function(x) paste(as.integer(x$size()), collapse = "x")
  common <- intersect(names(state), checkpoint_keys)
  bad_shape <- common[vapply(common, function(key) {
    !identical(as.integer(state[[key]]$size()), as.integer(weights[[key]]$size()))
  }, logical(1))]
  if (length(missing) || length(unexpected) || length(bad_shape)) {
    detail <- paste0(
      "UCE checkpoint audit failed; missing keys: ", paste(missing, collapse = ", "),
      "; unexpected keys: ", paste(unexpected, collapse = ", "),
      "; shape mismatches: ", paste(bad_shape, collapse = ", ")
    )
    stop(detail, call. = FALSE)
  }
  model$load_state_dict(weights, strict = TRUE)
  list(
    model = model,
    audit = list(
      checkpoint = normalizePath(path, winslash = "/", mustWork = TRUE),
      weight_source = normalizePath(weight_source, winslash = "/", mustWork = TRUE),
      weight_format = if (grepl("\\.safetensors$", weight_source, ignore.case = TRUE)) {
        "safetensors"
      } else "pytorch state_dict",
      checkpoint_tensors = length(weights), model_tensors = length(state),
      missing_keys = character(), unexpected_keys = character(),
      shape_mismatches = character(),
      tensor_table = data.frame(
        name = common, shape = vapply(common, function(key) shape_of(weights[[key]]), ""),
        stringsAsFactors = FALSE
      )
    )
  )
}
