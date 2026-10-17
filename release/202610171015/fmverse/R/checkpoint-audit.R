#' Inventory a Python PyTorch state dictionary
#'
#' Reads a plain Python state-dictionary file with torch for R and returns the
#' parameter names, shapes, and dtypes without modifying the tensors.
#'
#' @param path Path to a Python `.pt` or `.pth` state-dictionary file.
#' @param normalize_cuda Prepare an R-cached CPU-tagged copy when the checkpoint
#'   stores CUDA locations.
#' @return A data frame with one row per tensor.
#' @export
fmverse_checkpoint_inventory <- function(path, normalize_cuda = TRUE) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop("`path` must identify an existing checkpoint file.", call. = FALSE)
  }

  load_path <- if (isTRUE(normalize_cuda)) {
    fmverse_prepare_scgpt_checkpoint(path)
  } else {
    path
  }
  weights <- torch::load_state_dict(load_path)
  if (!is.list(weights) || is.null(names(weights))) {
    stop("The checkpoint did not load as a named state dictionary.", call. = FALSE)
  }

  shape_of <- function(x) {
    if (!inherits(x, "torch_tensor")) {
      stop("Checkpoint entry is not a torch tensor.", call. = FALSE)
    }
    paste(as.integer(x$size()), collapse = "x")
  }

  dtype_of <- function(x) {
    sub("^torch_", "", as.character(x$dtype))
  }
  numel_of <- function(x) prod(as.numeric(x$size()))

  data.frame(
    name = names(weights),
    shape = vapply(weights, shape_of, character(1)),
    dtype = vapply(weights, dtype_of, character(1)),
    numel = vapply(weights, numel_of, numeric(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' Compare state-dictionary keys and tensor shapes
#'
#' @param reference Named list of reference tensors, usually the checkpoint.
#' @param candidate Named list of R model tensors, usually `model$state_dict()`.
#' @return A list containing missing keys, unexpected keys, shape mismatches,
#'   and a logical `ok` value.
#' @export
fmverse_compare_state_dict <- function(reference, candidate) {
  if (!is.list(reference) || is.null(names(reference)) ||
      !is.list(candidate) || is.null(names(candidate))) {
    stop("Both inputs must be named lists of tensors.", call. = FALSE)
  }

  missing <- setdiff(names(reference), names(candidate))
  unexpected <- setdiff(names(candidate), names(reference))
  shared <- intersect(names(reference), names(candidate))
  shape_of <- function(x) paste(as.integer(x$size()), collapse = "x")
  mismatch <- shared[vapply(shared, function(key) {
    !identical(shape_of(reference[[key]]), shape_of(candidate[[key]]))
  }, logical(1))]

  shape_mismatches <- data.frame(
    name = mismatch,
    reference = vapply(mismatch, function(key) shape_of(reference[[key]]), character(1)),
    candidate = vapply(mismatch, function(key) shape_of(candidate[[key]]), character(1)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  list(
    ok = length(missing) == 0L && length(unexpected) == 0L && nrow(shape_mismatches) == 0L,
    missing = missing,
    unexpected = unexpected,
    shape_mismatches = shape_mismatches
  )
}
