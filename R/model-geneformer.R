# R-native BERT architecture and strict checkpoint adapter for Geneformer.

fmverse_geneformer_embeddings <- torch::nn_module(
  classname = "fmverse_geneformer_embeddings",
  initialize = function(vocab_size, hidden_size, max_position_embeddings,
                        type_vocab_size, pad_token_id, hidden_dropout_prob,
                        layer_norm_eps) {
    self$word_embeddings <- torch::nn_embedding(vocab_size, hidden_size,
                                                padding_idx = as.integer(pad_token_id) + 1L)
    self$position_embeddings <- torch::nn_embedding(max_position_embeddings, hidden_size)
    self$token_type_embeddings <- torch::nn_embedding(type_vocab_size, hidden_size)
    self$LayerNorm <- torch::nn_layer_norm(hidden_size, eps = layer_norm_eps)
    self$dropout <- torch::nn_dropout(hidden_dropout_prob)
  },
  forward = function(input_ids) {
    batch_size <- input_ids$size(1L)
    seq_length <- input_ids$size(2L)
    positions <- torch::torch_arange(
      0L, seq_length - 1L, dtype = torch::torch_long(), device = input_ids$device
    )$unsqueeze(1L)$expand(c(batch_size, seq_length))
    token_types <- torch::torch_zeros_like(input_ids)
    hidden <- self$word_embeddings(input_ids + 1L) +
      self$position_embeddings(positions + 1L) +
      self$token_type_embeddings(token_types + 1L)
    self$dropout(self$LayerNorm(hidden))
  }
)

fmverse_geneformer_self_attention <- torch::nn_module(
  classname = "fmverse_geneformer_self_attention",
  initialize = function(hidden_size, num_attention_heads, attention_probs_dropout_prob) {
    self$num_attention_heads <- as.integer(num_attention_heads)
    self$attention_head_size <- as.integer(hidden_size / num_attention_heads)
    self$all_head_size <- as.integer(hidden_size)
    self$query <- torch::nn_linear(hidden_size, hidden_size)
    self$key <- torch::nn_linear(hidden_size, hidden_size)
    self$value <- torch::nn_linear(hidden_size, hidden_size)
    self$dropout <- torch::nn_dropout(attention_probs_dropout_prob)
  },
  forward = function(hidden, attention_mask) {
    batch_size <- hidden$size(1L)
    seq_length <- hidden$size(2L)
    split_heads <- function(x) {
      x$view(c(batch_size, seq_length, self$num_attention_heads,
               self$attention_head_size))$permute(c(1L, 3L, 2L, 4L))
    }
    query <- split_heads(self$query(hidden))
    key <- split_heads(self$key(hidden))
    value <- split_heads(self$value(hidden))
    scores <- torch::torch_matmul(query, key$transpose(3L, 4L)) /
      sqrt(self$attention_head_size)
    masked <- attention_mask$eq(0L)$unsqueeze(2L)$unsqueeze(3L)
    scores <- scores$masked_fill(masked, -10000)
    probabilities <- scores$softmax(dim = 4L)
    probabilities <- self$dropout(probabilities)
    context <- torch::torch_matmul(probabilities, value)$permute(c(1L, 3L, 2L, 4L))
    context <- context$contiguous()$view(c(batch_size, seq_length, self$all_head_size))
    list(context = context, attention = probabilities)
  }
)

fmverse_geneformer_attention_output <- torch::nn_module(
  classname = "fmverse_geneformer_attention_output",
  initialize = function(hidden_size, hidden_dropout_prob, layer_norm_eps) {
    self$dense <- torch::nn_linear(hidden_size, hidden_size)
    self$dropout <- torch::nn_dropout(hidden_dropout_prob)
    self$LayerNorm <- torch::nn_layer_norm(hidden_size, eps = layer_norm_eps)
  },
  forward = function(context, residual) self$LayerNorm(self$dropout(self$dense(context)) + residual)
)

fmverse_geneformer_attention <- torch::nn_module(
  classname = "fmverse_geneformer_attention",
  initialize = function(hidden_size, num_attention_heads, attention_probs_dropout_prob,
                        hidden_dropout_prob, layer_norm_eps) {
    self$self_attention <- fmverse_geneformer_self_attention(
      hidden_size, num_attention_heads, attention_probs_dropout_prob
    )
    self$output <- fmverse_geneformer_attention_output(
      hidden_size, hidden_dropout_prob, layer_norm_eps
    )
  },
  forward = function(hidden, attention_mask) {
    output <- self$self_attention(hidden, attention_mask)
    list(hidden = self$output(output$context, hidden), attention = output$attention)
  }
)

fmverse_geneformer_intermediate <- torch::nn_module(
  classname = "fmverse_geneformer_intermediate",
  initialize = function(hidden_size, intermediate_size) {
    self$dense <- torch::nn_linear(hidden_size, intermediate_size)
    self$activation <- torch::nn_relu()
  },
  forward = function(hidden) self$activation(self$dense(hidden))
)

fmverse_geneformer_layer_output <- torch::nn_module(
  classname = "fmverse_geneformer_layer_output",
  initialize = function(intermediate_size, hidden_size, hidden_dropout_prob, layer_norm_eps) {
    self$dense <- torch::nn_linear(intermediate_size, hidden_size)
    self$dropout <- torch::nn_dropout(hidden_dropout_prob)
    self$LayerNorm <- torch::nn_layer_norm(hidden_size, eps = layer_norm_eps)
  },
  forward = function(intermediate, residual) {
    self$LayerNorm(self$dropout(self$dense(intermediate)) + residual)
  }
)

fmverse_geneformer_layer <- torch::nn_module(
  classname = "fmverse_geneformer_layer",
  initialize = function(hidden_size, num_attention_heads, intermediate_size,
                        attention_probs_dropout_prob, hidden_dropout_prob, layer_norm_eps) {
    self$attention <- fmverse_geneformer_attention(
      hidden_size, num_attention_heads, attention_probs_dropout_prob,
      hidden_dropout_prob, layer_norm_eps
    )
    self$intermediate <- fmverse_geneformer_intermediate(hidden_size, intermediate_size)
    self$output <- fmverse_geneformer_layer_output(
      intermediate_size, hidden_size, hidden_dropout_prob, layer_norm_eps
    )
  },
  forward = function(hidden, attention_mask) {
    attention <- self$attention(hidden, attention_mask)
    output <- self$output(self$intermediate(attention$hidden), attention$hidden)
    list(hidden = output, attention = attention$attention)
  }
)

fmverse_geneformer_encoder <- torch::nn_module(
  classname = "fmverse_geneformer_encoder",
  initialize = function(num_hidden_layers, hidden_size, num_attention_heads,
                        intermediate_size, attention_probs_dropout_prob,
                        hidden_dropout_prob, layer_norm_eps) {
    layers <- lapply(seq_len(num_hidden_layers), function(i) fmverse_geneformer_layer(
      hidden_size, num_attention_heads, intermediate_size,
      attention_probs_dropout_prob, hidden_dropout_prob, layer_norm_eps
    ))
    self$layer <- torch::nn_module_list(layers)
    self$num_hidden_layers <- as.integer(num_hidden_layers)
  },
  forward = function(hidden, attention_mask, return_layers = FALSE,
                     return_attention = FALSE, embedding_layer = 0L) {
    if (return_layers) hidden_states <- vector("list", self$num_hidden_layers + 1L)
    if (return_layers) hidden_states[[1L]] <- hidden
    selected_hidden <- if (embedding_layer == 0L) hidden else NULL
    if (return_attention) attention_sum <- NULL
    for (i in seq_len(self$num_hidden_layers)) {
      output <- self$layer[[i]](hidden, attention_mask)
      hidden <- output$hidden
      if (return_layers) hidden_states[[i + 1L]] <- hidden
      if (embedding_layer == i) selected_hidden <- hidden
      if (return_attention) {
        layer_attention <- torch::torch_mean(output$attention, dim = 2L)
        attention_sum <- if (is.null(attention_sum)) layer_attention else attention_sum + layer_attention
      }
    }
    out <- list(last_hidden_state = hidden, embedding_hidden_state = selected_hidden)
    if (return_layers) {
      names(hidden_states) <- c("embeddings", sprintf("layer_%02d", seq_len(self$num_hidden_layers)))
      out$hidden_states <- hidden_states
    }
    if (return_attention) out$attention <- attention_sum / self$num_hidden_layers
    out
  }
)

fmverse_geneformer_bert <- torch::nn_module(
  classname = "fmverse_geneformer_bert",
  initialize = function(config) {
    self$embeddings <- fmverse_geneformer_embeddings(
      config$vocab_size, config$hidden_size, config$max_position_embeddings,
      config$type_vocab_size, config$pad_token_id, config$hidden_dropout_prob,
      config$layer_norm_eps
    )
    self$encoder <- fmverse_geneformer_encoder(
      config$num_hidden_layers, config$hidden_size, config$num_attention_heads,
      config$intermediate_size, config$attention_probs_dropout_prob,
      config$hidden_dropout_prob, config$layer_norm_eps
    )
  },
  forward = function(input_ids, attention_mask, return_layers = FALSE,
                     return_attention = FALSE, embedding_layer = 0L) {
    hidden <- self$embeddings(input_ids)
    self$encoder(hidden, attention_mask, return_layers, return_attention,
                 embedding_layer = embedding_layer)
  }
)

fmverse_geneformer_predictions <- torch::nn_module(
  classname = "fmverse_geneformer_predictions",
  initialize = function(hidden_size, vocab_size, layer_norm_eps) {
    self$transform <- torch::nn_module_list(list(
      torch::nn_linear(hidden_size, hidden_size),
      torch::nn_relu(),
      torch::nn_layer_norm(hidden_size, eps = layer_norm_eps)
    ))
    self$bias <- torch::nn_parameter(torch::torch_zeros(vocab_size))
  },
  forward = function(hidden, embedding_weight) {
    hidden <- self$transform[[1L]](hidden)
    hidden <- self$transform[[2L]](hidden)
    hidden <- self$transform[[3L]](hidden)
    torch::nnf_linear(hidden, embedding_weight, self$bias)
  }
)

fmverse_geneformer_mlm_head <- torch::nn_module(
  classname = "fmverse_geneformer_mlm_head",
  initialize = function(hidden_size, vocab_size, layer_norm_eps) {
    self$predictions <- fmverse_geneformer_predictions(
      hidden_size, vocab_size, layer_norm_eps
    )
  },
  forward = function(hidden, embedding_weight) self$predictions(hidden, embedding_weight)
)

fmverse_geneformer_model <- function(config) {
  model <- torch::nn_module(
    classname = "fmverse_geneformer_model",
    initialize = function(config) {
      self$bert <- fmverse_geneformer_bert(config)
      self$cls <- fmverse_geneformer_mlm_head(
        config$hidden_size, config$vocab_size, config$layer_norm_eps
      )
      self$hidden_size <- as.integer(config$hidden_size)
      self$num_hidden_layers <- as.integer(config$num_hidden_layers)
      self$cls_token_id <- if (is.null(config$cls_token_id)) -1L else as.integer(config$cls_token_id)
      self$eos_token_id <- if (is.null(config$eos_token_id)) -1L else as.integer(config$eos_token_id)
    },
    forward = function(input_ids, attention_mask, return_layers = FALSE,
                       return_attention = FALSE, return_logits = FALSE,
                       embedding_layer = -1L) {
      # Geneformer's default EmbExtractor uses emb_layer=-1: the second-to-last
      # transformer block, then averages non-padding token embeddings.
      layer_index <- if (embedding_layer < 0L) self$num_hidden_layers + embedding_layer else embedding_layer
      if (layer_index < 0L || layer_index > self$num_hidden_layers) {
        stop("embedding_layer is outside the available hidden-state range.")
      }
      bert <- self$bert(input_ids, attention_mask, return_layers, return_attention,
                        embedding_layer = as.integer(layer_index))
      selected_hidden <- if (isTRUE(return_layers)) {
        bert$hidden_states[[as.integer(layer_index) + 1L]]
      } else {
        bert$embedding_hidden_state
      }
      if (is.null(selected_hidden)) stop("Requested Geneformer hidden state was not produced.")
      pooling_mask <- attention_mask$to(dtype = torch::torch_float())
      if (self$cls_token_id >= 0L) {
        pooling_mask <- pooling_mask * input_ids$ne(as.integer(self$cls_token_id))$to(
          dtype = torch::torch_float()
        )
      }
      if (self$eos_token_id >= 0L) {
        pooling_mask <- pooling_mask * input_ids$ne(as.integer(self$eos_token_id))$to(
          dtype = torch::torch_float()
        )
      }
      valid <- pooling_mask$unsqueeze(3L)
      pooled <- torch::torch_sum(selected_hidden * valid, dim = 2L) /
        torch::torch_clamp(torch::torch_sum(valid, dim = 2L), min = 1)
      output <- list(
        last_hidden_state = bert$last_hidden_state,
        hidden_states = bert$hidden_states,
        cell_embedding = pooled
      )
      if (return_attention) output$attention <- bert$attention
      if (return_logits) {
        output$logits <- self$cls(
          selected_hidden, self$bert$embeddings$word_embeddings$weight
        )
      }
      output
    }
  )
  model(config)
}

#' Construct an R-native Geneformer BERT model from checkpoint config.
#' @param config A decoded HuggingFace config.json list.
#' @return An R torch module with BERT encoder and masked-language-model head.
#' @export
fmverse_geneformer_build_model <- function(config) {
  fmverse_geneformer_model(config)
}

# Map HuggingFace's `attention.self` and prediction-head names to the legal R
# module fields. R torch module lists also expose zero-based state-dict names.
#' @export
fmverse_geneformer_keymap <- function(checkpoint_names, model_state) {
  target <- names(model_state)
  mapped <- stats::setNames(rep(NA_character_, length(checkpoint_names)), checkpoint_names)
  for (key in checkpoint_names) {
    candidates <- key
    candidates <- c(candidates, sub("attention\\.self\\.", "attention.self_attention.", key))
    candidates <- c(candidates,
      sub("cls\\.predictions\\.transform\\.dense\\.",
          "cls.predictions.transform.0.", key),
      sub("cls\\.predictions\\.transform\\.LayerNorm\\.",
          "cls.predictions.transform.2.", key),
      sub("cls\\.predictions\\.decoder\\.bias$", "cls.predictions.bias", key)
    )
    hit <- intersect(unique(candidates), target)
    if (length(hit) == 1L) mapped[[key]] <- hit
  }
  mapped
}

#' Load and audit a pinned Geneformer safetensors checkpoint.
#' @param model An R-native Geneformer model.
#' @param path Checkpoint `model.safetensors` path.
#' @return Loaded model and exact key/shape audit.
#' @export
fmverse_load_geneformer_weights <- function(model, path) {
  if (!file.exists(path)) stop("Geneformer checkpoint does not exist: ", path)
  weights <- safetensors::safe_load_file(path, framework = "torch")
  state <- model$state_dict()
  # Recent Hugging Face exports may persist this fixed BERT positional-index
  # buffer even though it is not a learned parameter and torch for R does not
  # register it in the module state. Validate it explicitly, then account for
  # it in the audit while loading all learned tensors strictly.
  buffer_keys <- intersect("bert.embeddings.position_ids", names(weights))
  if (length(buffer_keys)) {
    positions <- torch::as_array(weights[[buffer_keys[[1L]]]]$cpu())
    expected_positions <- matrix(
      seq.int(0L, as.integer(model$bert$embeddings$position_embeddings$weight$size(1L)) - 1L),
      nrow = 1L
    )
    if (!identical(dim(positions), dim(expected_positions)) ||
        !identical(as.integer(positions), as.integer(expected_positions))) {
      stop("Geneformer position_ids buffer does not match max_position_embeddings.")
    }
  }
  parameter_keys <- setdiff(names(weights), buffer_keys)
  keymap <- fmverse_geneformer_keymap(parameter_keys, state)
  unmapped <- names(keymap)[is.na(keymap)]
  mapped_targets <- unname(keymap[!is.na(keymap)])
  missing <- setdiff(names(state), mapped_targets)
  duplicates <- unique(mapped_targets[duplicated(mapped_targets)])
  if (length(unmapped) || length(missing) || length(duplicates)) {
    stop(
      "Geneformer state-dict key audit failed. Unmapped checkpoint keys: ",
      paste(unmapped, collapse = ", "), "; missing R keys: ",
      paste(missing, collapse = ", "), "; duplicate targets: ",
      paste(duplicates, collapse = ", ")
    )
  }
  checkpoint_dims <- lapply(weights, function(x) as.integer(unlist(x$size())))
  model_dims <- lapply(state, function(x) as.integer(unlist(x$size())))
  shape_ok <- vapply(parameter_keys, function(key) {
    identical(checkpoint_dims[[key]], model_dims[[keymap[[key]]]])
  }, logical(1))
  if (any(!shape_ok)) {
    bad <- names(weights)[!shape_ok]
    stop("Geneformer tensor shapes differ for: ", paste(bad, collapse = ", "))
  }
  mapped_state <- stats::setNames(vector("list", length(state)), names(state))
  for (key in parameter_keys) mapped_state[[keymap[[key]]]] <- weights[[key]]
  model$load_state_dict(mapped_state, strict = TRUE)
  audit_keys <- c(names(keymap), buffer_keys)
  audit_models <- c(unname(keymap), if (length(buffer_keys)) {
    "<buffer:bert.embeddings.position_ids>"
  })
  list(
    model = model,
    audit = list(
      checkpoint = normalizePath(path, winslash = "/", mustWork = TRUE),
      checkpoint_key_count = length(weights), model_key_count = length(state),
      parameter_key_count = length(parameter_keys), buffer_key_count = length(buffer_keys),
      matched_key_count = length(parameter_keys), shape_mismatch_count = 0L,
      missing_model_keys = character(), unmapped_checkpoint_keys = character(),
      key_mapping = data.frame(checkpoint_key = audit_keys,
                               model_key = audit_models,
                               shape = vapply(audit_keys, function(key) {
                                 paste(checkpoint_dims[[key]], collapse = "x")
                               }, character(1)),
                               stringsAsFactors = FALSE)
    )
  )
}
