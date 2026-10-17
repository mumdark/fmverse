# Shared R torch primitives for the Performer blocks used by scBERT and
# scFoundation. Parameter names mirror the official PyTorch module hierarchy.

fmverse_performer_pre_layer_norm <- torch::nn_module(
  classname = "fmverse_performer_pre_layer_norm",
  initialize = function(dim, fn) {
    self$norm <- torch::nn_layer_norm(dim)
    self$fn <- fn
  },
  forward = function(x) self$fn(self$norm(x))
)

fmverse_performer_fast_attention <- torch::nn_module(
  classname = "fmverse_performer_fast_attention",
  initialize = function(dim_head, nb_features = NULL) {
    if (is.null(nb_features)) nb_features <- as.integer(dim_head * log(dim_head))
    self$dim_heads <- as.integer(dim_head)
    self$nb_features <- as.integer(nb_features)
    # PyTorch saves this as a buffer. A frozen parameter has the same state-dict
    # key and storage shape, and lets R torch restore the original projection.
    self$projection_matrix <- torch::nn_parameter(
      torch::torch_zeros(c(self$nb_features, self$dim_heads)), requires_grad = FALSE
    )
  },
  forward = function(q, k, v) {
    data_normalizer <- self$dim_heads ^ (-0.25)
    ratio <- self$nb_features ^ (-0.5)
    project <- function(x, is_query) {
      projected <- torch::torch_matmul(x * data_normalizer,
                                       self$projection_matrix$transpose(1L, 2L))
      diag <- torch::torch_sum(x * x, dim = 4L, keepdim = TRUE) *
        (data_normalizer ^ 2 / 2)
      if (is_query) {
        maximum <- torch::torch_max(projected, dim = 4L, keepdim = TRUE)[[1L]]
      } else {
        maximum <- torch::torch_max(projected)
      }
      ratio * (torch::torch_exp(projected - diag - maximum) + 1e-4)
    }
    q_features <- project(q, TRUE)
    k_features <- project(k, FALSE)
    k_sum <- torch::torch_sum(k_features, dim = 3L)
    denominator <- torch::torch_einsum("...nd,...d->...n", list(q_features, k_sum))
    context <- torch::torch_einsum("...nd,...ne->...de", list(k_features, v))
    torch::torch_einsum("...de,...nd,...n->...ne", list(context, q_features,
      1 / denominator))
  }
)

fmverse_performer_self_attention <- torch::nn_module(
  classname = "fmverse_performer_self_attention",
  initialize = function(dim, heads, dim_head = 64L, nb_features = NULL,
                        qkv_bias = TRUE, dropout = 0) {
    inner_dim <- as.integer(heads * dim_head)
    self$fast_attention <- fmverse_performer_fast_attention(dim_head, nb_features)
    self$heads <- as.integer(heads)
    self$to_q <- torch::nn_linear(dim, inner_dim, bias = qkv_bias)
    self$to_k <- torch::nn_linear(dim, inner_dim, bias = qkv_bias)
    self$to_v <- torch::nn_linear(dim, inner_dim, bias = qkv_bias)
    self$to_out <- torch::nn_linear(inner_dim, dim)
    self$dropout <- torch::nn_dropout(dropout)
  },
  forward = function(x) {
    batch <- x$size(1L); seq <- x$size(2L)
    split_heads <- function(layer) {
      layer(x)$view(c(batch, seq, self$heads, -1L))$permute(c(1L, 3L, 2L, 4L))
    }
    attended <- self$fast_attention(split_heads(self$to_q), split_heads(self$to_k),
                                    split_heads(self$to_v))
    attended <- attended$permute(c(1L, 3L, 2L, 4L))$contiguous()$view(
      c(batch, seq, -1L)
    )
    self$dropout(self$to_out(attended))
  }
)

fmverse_performer_feed_forward <- torch::nn_module(
  classname = "fmverse_performer_feed_forward",
  initialize = function(dim, mult = 4L, dropout = 0) {
    self$w1 <- torch::nn_linear(dim, dim * mult)
    self$act <- torch::nn_gelu()
    self$dropout <- torch::nn_dropout(dropout)
    self$w2 <- torch::nn_linear(dim * mult, dim)
  },
  forward = function(x) self$w2(self$dropout(self$act(self$w1(x))))
)

fmverse_performer_chunk <- torch::nn_module(
  classname = "fmverse_performer_chunk",
  initialize = function(fn) self$fn <- fn,
  forward = function(x) self$fn(x)
)

fmverse_performer_block <- torch::nn_module(
  classname = "fmverse_performer_block",
  initialize = function(dim, heads, dim_head = 64L, nb_features = NULL,
                        ff_mult = 4L, attn_dropout = 0, ff_dropout = 0,
                        qkv_bias = TRUE) {
    attention <- fmverse_performer_self_attention(
      dim, heads, dim_head, nb_features, qkv_bias = qkv_bias, dropout = attn_dropout
    )
    feed_forward <- fmverse_performer_chunk(
      fmverse_performer_feed_forward(dim, ff_mult, ff_dropout)
    )
    self$`0` <- fmverse_performer_pre_layer_norm(dim, attention)
    self$`1` <- fmverse_performer_pre_layer_norm(dim, feed_forward)
  },
  forward = function(x) {
    x <- x + self$`0`(x)
    x + self$`1`(x)
  }
)

fmverse_performer_sequence <- torch::nn_module(
  classname = "fmverse_performer_sequence",
  initialize = function(dim, depth, heads, dim_head = 64L, nb_features = NULL,
                        ff_mult = 4L, attn_dropout = 0, ff_dropout = 0,
                        qkv_bias = TRUE) {
    self$layers <- torch::nn_module_list(lapply(seq_len(depth), function(i) {
      fmverse_performer_block(dim, heads, dim_head, nb_features, ff_mult,
                              attn_dropout, ff_dropout, qkv_bias)
    }))
  },
  forward = function(x) {
    for (i in seq_len(length(self$layers))) x <- self$layers[[i]](x)
    x
  }
)

fmverse_performer_core <- torch::nn_module(
  classname = "fmverse_performer_core",
  initialize = function(dim, depth, heads, dim_head = 64L, nb_features = NULL,
                        ff_mult = 4L, attn_dropout = 0, ff_dropout = 0,
                        qkv_bias = TRUE) {
    self$net <- fmverse_performer_sequence(
      dim, depth, heads, dim_head, nb_features, ff_mult, attn_dropout, ff_dropout,
      qkv_bias
    )
    # Registered in upstream Performer and present in the serialized state dict.
    self$calls_since_last_redraw <- torch::nn_parameter(
      torch::torch_tensor(0L), requires_grad = FALSE
    )
  },
  forward = function(x) self$net(x)
)

fmverse_scfoundation_transformer <- torch::nn_module(
  classname = "fmverse_scfoundation_transformer",
  initialize = function(max_seq_len, dim, depth, heads, dropout = 0.1) {
    self$max_seq_len <- as.integer(max_seq_len)
    layers <- lapply(seq_len(depth), function(i) torch::nn_transformer_encoder_layer(
      d_model = dim, nhead = heads, dim_feedforward = dim * 4L,
      dropout = dropout, activation = "relu", batch_first = TRUE,
      norm_first = FALSE
    ))
    self$transformer_encoder <- torch::nn_module_list(layers)
    self$norm <- torch::nn_layer_norm(dim)
  },
  forward = function(x, padding_mask = NULL) {
    if (x$size(2L) > self$max_seq_len) stop("Input exceeds transformer max_seq_len.")
    for (i in seq_len(length(self$transformer_encoder))) {
      x <- self$transformer_encoder[[i]](x, src_key_padding_mask = padding_mask)
    }
    self$norm(x)
  }
)

fmverse_scfoundation_performer <- torch::nn_module(
  classname = "fmverse_scfoundation_performer",
  initialize = function(max_seq_len, dim, depth, heads, dim_head = 64L,
                        nb_features = NULL, ff_mult = 4L, ff_dropout = 0,
                        attn_dropout = 0) {
    self$max_seq_len <- as.integer(max_seq_len)
    self$performer <- fmverse_performer_core(
      dim, depth, heads, dim_head, nb_features, ff_mult, attn_dropout, ff_dropout
    )
    self$norm <- torch::nn_layer_norm(dim)
  },
  forward = function(x, padding_mask = NULL) {
    if (x$size(2L) > self$max_seq_len) stop("Input exceeds Performer max_seq_len.")
    # The official xTrimoGene Performer routes `mask`, while its wrapper passes
    # `padding_mask`; reproducing that path keeps the cell-wise reference behavior.
    self$norm(self$performer(x))
  }
)
