# R-native architecture matching the pinned scGPT whole-human checkpoint.
# The saved classifier head is not part of this pretraining checkpoint.

fmverse_scgpt_gene_encoder <- torch::nn_module(
  classname = "fmverse_scgpt_gene_encoder",
  initialize = function(ntoken, d_model, pad_id) {
    # scGPT's checkpoint IDs are zero-based; torch for R exposes embedding
    # indices as one-based values, while the checkpoint rows remain unchanged.
    self$embedding <- torch::nn_embedding(
      ntoken, d_model, padding_idx = as.integer(pad_id) + 1L
    )
    self$enc_norm <- torch::nn_layer_norm(d_model)
  },
  forward = function(gene_ids) {
    self$enc_norm(self$embedding(gene_ids + 1L))
  }
)

fmverse_scgpt_value_encoder <- torch::nn_module(
  classname = "fmverse_scgpt_value_encoder",
  initialize = function(d_model, dropout = 0.2, max_value = 512) {
    self$linear1 <- torch::nn_linear(1, d_model)
    self$activation <- torch::nn_relu()
    self$linear2 <- torch::nn_linear(d_model, d_model)
    self$norm <- torch::nn_layer_norm(d_model)
    self$dropout <- torch::nn_dropout(dropout)
    self$max_value <- max_value
  },
  forward = function(values) {
    values <- torch::torch_clamp(values$unsqueeze(3), max = self$max_value)
    values <- self$activation(self$linear1(values))
    values <- self$linear2(values)
    self$dropout(self$norm(values))
  }
)

fmverse_scgpt_expr_decoder <- torch::nn_module(
  classname = "fmverse_scgpt_expr_decoder",
  initialize = function(d_model) {
    self$fc <- torch::nn_sequential(
      torch::nn_linear(d_model, d_model), torch::nn_leaky_relu(),
      torch::nn_linear(d_model, d_model), torch::nn_leaky_relu(),
      torch::nn_linear(d_model, 1)
    )
  },
  forward = function(x) self$fc(x)
)

fmverse_scgpt_cls_decoder <- torch::nn_module(
  classname = "fmverse_scgpt_cls_decoder",
  initialize = function(d_model, n_cls, nlayers = 3L) {
    layers <- list()
    if (nlayers > 1L) {
      for (i in seq_len(nlayers - 1L)) {
        layers[[length(layers) + 1L]] <- torch::nn_linear(d_model, d_model)
        layers[[length(layers) + 1L]] <- torch::nn_relu()
        layers[[length(layers) + 1L]] <- torch::nn_layer_norm(d_model)
      }
    }
    self$`_decoder` <- do.call(torch::nn_sequential, layers)
    self$out_layer <- torch::nn_linear(d_model, n_cls)
  },
  forward = function(x) self$out_layer(self$`_decoder`(x))
)

fmverse_scgpt_mvc_decoder <- torch::nn_module(
  classname = "fmverse_scgpt_mvc_decoder",
  initialize = function(d_model) {
    self$gene2query <- torch::nn_linear(d_model, d_model)
    self$query_activation <- torch::nn_sigmoid()
    self$W <- torch::nn_linear(d_model, d_model, bias = FALSE)
  },
  forward = function(cell_embedding, gene_embedding) {
    query <- self$query_activation(self$gene2query(gene_embedding))
    torch::torch_bmm(self$W(query), cell_embedding$unsqueeze(3))$squeeze(3)
  }
)

#' Construct the R-native scGPT whole-human model
#'
#' @param ntoken Vocabulary size, including the checkpoint's special tokens.
#' @param pad_id Vocabulary ID for `<pad>`.
#' @param d_model Embedding dimension.
#' @param nhead Number of attention heads.
#' @param d_hid Feed-forward hidden dimension.
#' @param nlayers Number of transformer blocks.
#' @param dropout Dropout used by the checkpoint architecture.
#' @return An R torch `nn_module` with checkpoint-compatible components.
#' @export
fmverse_scgpt_model <- function(ntoken = 60697L, pad_id = 60694L,
                                d_model = 512L, nhead = 8L,
                                d_hid = 512L, nlayers = 12L,
                                dropout = 0.2) {
  fmverse_scgpt_module(
    ntoken = ntoken, pad_id = pad_id, d_model = d_model, nhead = nhead,
    d_hid = d_hid, nlayers = nlayers, dropout = dropout
  )
}

fmverse_scgpt_module <- torch::nn_module(
  classname = "fmverse_scgpt_module",
  initialize = function(ntoken, pad_id, d_model, nhead, d_hid, nlayers, dropout) {
    self$encoder <- fmverse_scgpt_gene_encoder(ntoken, d_model, pad_id)
    self$value_encoder <- fmverse_scgpt_value_encoder(d_model, dropout)
    layer <- torch::nn_transformer_encoder_layer(
      d_model = d_model, nhead = nhead, dim_feedforward = d_hid,
      dropout = dropout, activation = "gelu", layer_norm_eps = 1e-5,
      batch_first = TRUE, norm_first = FALSE
    )
    self$transformer_encoder <- torch::nn_transformer_encoder(layer, nlayers)
    self$decoder <- fmverse_scgpt_expr_decoder(d_model)
    self$mvc_decoder <- fmverse_scgpt_mvc_decoder(d_model)
    # The upstream model does not use this table for the selected inference
    # path, but the official checkpoint includes it and it must be retained.
    self$flag_encoder <- torch::nn_embedding(2L, d_model)
    self$pad_id <- pad_id
    self$d_model <- d_model
    self$nlayers <- nlayers
  },
  forward = function(gene_ids, values, pad_mask = NULL, return_layers = FALSE,
                     return_expression = FALSE, return_attention = FALSE) {
    if (is.null(pad_mask)) {
      pad_mask <- gene_ids == self$pad_id
    }
    gene_embedding <- self$encoder(gene_ids)
    value_embedding <- self$value_encoder(values)
    x <- gene_embedding + value_embedding
    if (return_layers) layer_outputs <- vector("list", self$nlayers)
    if (return_attention) attention_sum <- NULL
    for (i in seq_len(self$nlayers)) {
      if (return_attention) {
        # R torch exposes the same per-head weights as PyTorch's
        # MultiheadAttention. The regular encoder-layer call below remains the
        # source of hidden states, so the parity-critical path is unchanged.
        attention_output <- self$transformer_encoder$layers[[i]]$self_attn(
          x, x, x, key_padding_mask = pad_mask,
          need_weights = TRUE, avg_weights = FALSE
        )
        layer_attention <- attention_output[[2L]]
        if (is.null(attention_sum)) {
          attention_sum <- torch::torch_mean(layer_attention, dim = 2)
        } else {
          attention_sum <- attention_sum + torch::torch_mean(layer_attention, dim = 2)
        }
      }
      x <- self$transformer_encoder$layers[[i]](
        x, src_key_padding_mask = pad_mask
      )
      if (return_layers) layer_outputs[[i]] <- x
    }
    # Preserve the batch axis when a forward pass contains one cell.
    cell_embedding <- x[, 1, , drop = FALSE]$squeeze(2)
    norm <- torch::torch_sqrt(
      torch::torch_sum(cell_embedding * cell_embedding, dim = 2, keepdim = TRUE)
    )
    cell_embedding_l2 <- cell_embedding / torch::torch_clamp(norm, min = 1e-12)
    result <- list(
      cell_embedding = cell_embedding,
      cell_embedding_l2 = cell_embedding_l2
    )
    if (return_layers) {
      names(layer_outputs) <- sprintf("layer_%02d", 0:(self$nlayers - 1L))
      result <- c(list(
        gene_embedding = gene_embedding,
        value_embedding = value_embedding,
        transformer_input = gene_embedding + value_embedding
      ), result, layer_outputs)
    }
    if (return_expression) {
      result$expression_prediction <- self$decoder(x)$squeeze(3)
    }
    if (return_attention) {
      result$attention <- attention_sum / self$nlayers
    }
    result
  }
)

#' Map checkpoint keys to the R module's parameter names
#'
#' Both R torch and PyTorch name transformer layers from zero. This adapter
#' maps only the upstream fused QKV names to LibTorch's in-projection names.
#' @param checkpoint_names Character vector of official checkpoint keys.
#' @param model_state Named list returned by `model$state_dict()`.
#' @return A named character vector from checkpoint keys to R keys.
#' @export
fmverse_scgpt_keymap <- function(checkpoint_names, model_state) {
  target <- names(model_state)
  mapped <- setNames(rep(NA_character_, length(checkpoint_names)), checkpoint_names)
  for (key in checkpoint_names) {
    candidates <- key
    candidates <- c(
      candidates,
      sub("self_attn\\.Wqkv\\.weight$", "self_attn.in_proj_weight", key),
      sub("self_attn\\.Wqkv\\.bias$", "self_attn.in_proj_bias", key)
    )
    candidates <- unique(candidates)
    hit <- intersect(candidates, target)
    if (length(hit) == 1L) mapped[[key]] <- hit
  }
  mapped
}

#' Load and audit the official scGPT checkpoint in R torch
#'
#' @param model An `fmverse_scgpt_model()` instance.
#' @param path Path to `best_model.pt`.
#' @param nlayers_cls Optional number of layers in an extra cell-type
#'   classification decoder; inferred from checkpoint keys when omitted.
#' @return A list with the loaded model and a strict key/shape audit.
#' @export
fmverse_load_scgpt_weights <- function(model, path, nlayers_cls = NULL) {
  if (!file.exists(path)) stop("Checkpoint file does not exist: ", path)
  load_path <- fmverse_prepare_scgpt_checkpoint(path)
  checkpoint <- torch::load_state_dict(load_path)

  # The continual-pretraining checkpoint includes the supervised cell-type
  # decoder. Keep that task head in the R module and load it strictly, even
  # though the current fmverse embedding path only uses the shared encoder.
  if (any(startsWith(names(checkpoint), "cls_decoder.")) &&
      !any(startsWith(names(model$state_dict()), "cls_decoder."))) {
    output_key <- "cls_decoder.out_layer.weight"
    if (is.null(checkpoint[[output_key]])) {
      stop("Checkpoint has a cls_decoder but no out_layer.weight.", call. = FALSE)
    }
    n_cls <- as.integer(checkpoint[[output_key]]$size()[[1L]])
    if (is.null(nlayers_cls) || !is.finite(nlayers_cls) || nlayers_cls < 1L) {
      decoder_indices <- names(checkpoint)[grepl(
        "^cls_decoder\\._decoder\\.[0-9]+\\.weight$", names(checkpoint)
      )]
      index <- suppressWarnings(as.integer(sub(
        "^cls_decoder\\._decoder\\.([0-9]+)\\.weight$", "\\1", decoder_indices
      )))
      n_hidden <- sum(index %% 3L == 0L, na.rm = TRUE)
      nlayers_cls <- n_hidden + 1L
    }
    model$cls_decoder <- fmverse_scgpt_cls_decoder(
      d_model = as.integer(checkpoint[[output_key]]$size()[[2L]]),
      n_cls = n_cls,
      nlayers = as.integer(nlayers_cls)
    )
  }

  model_state <- model$state_dict()
  keymap <- fmverse_scgpt_keymap(names(checkpoint), model_state)
  unmapped <- names(keymap)[is.na(keymap)]
  if (length(unmapped)) {
    stop("Checkpoint keys could not be mapped: ", paste(unmapped, collapse = ", "))
  }
  if (anyDuplicated(unname(keymap))) {
    stop("Checkpoint key mapping is not one-to-one.")
  }
  shape_of <- function(x) paste(as.integer(x$size()), collapse = "x")
  shape_errors <- names(keymap)[vapply(names(keymap), function(key) {
    shape_of(checkpoint[[key]]) != shape_of(model_state[[keymap[[key]]]])
  }, logical(1))]
  if (length(shape_errors)) {
    stop("Checkpoint/model shapes differ for: ", paste(shape_errors, collapse = ", "))
  }
  mapped <- checkpoint
  names(mapped) <- unname(keymap)
  model$load_state_dict(mapped)
  list(
    model = model,
    audit = list(
      checkpoint_tensors = length(checkpoint),
      loaded_tensors = length(mapped),
      original_checkpoint = normalizePath(path, winslash = "/", mustWork = TRUE),
      load_checkpoint = normalizePath(load_path, winslash = "/", mustWork = TRUE),
      cpu_rewritten_tensors = if (identical(path, load_path)) 0L else length(checkpoint),
      checkpoint_scalar_parameters = sum(vapply(checkpoint, function(x) {
        prod(as.integer(x$size()))
      }, numeric(1))),
      checkpoint_only_keys = setdiff(names(checkpoint), names(keymap)[!is.na(keymap)]),
      model_only_keys = setdiff(names(model_state), unname(keymap)),
      shape_mismatches = character(),
      keymap = data.frame(
        checkpoint = names(keymap), model = unname(keymap),
        stringsAsFactors = FALSE
      )
    )
  )
}
