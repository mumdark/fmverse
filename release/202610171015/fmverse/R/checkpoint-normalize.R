# This checkpoint was serialized with CUDA device tags. Rewrite only those
# pickle location strings in an R-managed cache; tensor payloads stay intact.
fmverse_prepare_scgpt_checkpoint <- function(path, cache_dir = NULL) {
  if (length(path) != 1L || is.na(path) || !file.exists(path)) {
    stop("`path` must identify an existing checkpoint file.", call. = FALSE)
  }
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  if (is.null(cache_dir)) {
    cache_dir <- file.path(tools::R_user_dir("fmverse", "cache"), "checkpoints")
  }
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  source_hash <- digest::digest(file = path, algo = "sha256")
  output <- file.path(cache_dir, paste0("scgpt-", source_hash, "-cpu.pt"))
  if (file.exists(output)) return(output)

  work <- tempfile("fmverse-scgpt-cpu-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE, force = TRUE), add = TRUE)
  utils::unzip(path, exdir = work)
  pickle_files <- list.files(
    work, pattern = "^data[.]pkl$", recursive = TRUE, full.names = TRUE
  )
  if (length(pickle_files) != 1L) {
    stop("Expected one PyTorch data.pkl entry in the checkpoint ZIP.", call. = FALSE)
  }

  pickle_path <- pickle_files[[1L]]
  pickle <- readBin(pickle_path, what = "raw", n = file.info(pickle_path)$size)
  # BINUNICODE (0x58), little-endian length 6, and the string "cuda:0".
  needle <- as.raw(c(0x58, 0x06, 0x00, 0x00, 0x00, as.integer(charToRaw("cuda:0"))))
  replacement <- as.raw(c(0x58, 0x03, 0x00, 0x00, 0x00, as.integer(charToRaw("cpu"))))
  starts <- which(pickle == needle[[1L]])
  starts <- starts[starts <= length(pickle) - length(needle) + 1L]
  hits <- starts[vapply(starts, function(i) {
    identical(pickle[i:(i + length(needle) - 1L)], needle)
  }, logical(1))]
  cuda_tag <- charToRaw("cuda:")
  tag_starts <- which(pickle == cuda_tag[[1L]])
  tag_starts <- tag_starts[tag_starts <= length(pickle) - length(cuda_tag) + 1L]
  tag_hits <- tag_starts[vapply(tag_starts, function(i) {
    identical(pickle[i:(i + length(cuda_tag) - 1L)], cuda_tag)
  }, logical(1))]
  if (length(tag_hits) != length(hits)) {
    stop("Checkpoint contains unsupported or mixed CUDA location tags.", call. = FALSE)
  }
  if (!length(hits)) return(path)

  pieces <- vector("list", 2L * length(hits) + 1L)
  cursor <- 1L
  piece <- 1L
  for (start in hits) {
    if (start > cursor) {
      pieces[[piece]] <- pickle[cursor:(start - 1L)]
      piece <- piece + 1L
    }
    pieces[[piece]] <- replacement
    piece <- piece + 1L
    cursor <- start + length(needle)
  }
  if (cursor <= length(pickle)) pieces[[piece]] <- pickle[cursor:length(pickle)]
  writeBin(do.call(c, pieces[seq_len(piece)]), pickle_path)

  archive_root <- basename(dirname(pickle_path))
  output_tmp <- tempfile("fmverse-checkpoint-", tmpdir = cache_dir, fileext = ".pt")
  on.exit(unlink(output_tmp, force = TRUE), add = TRUE)
  zip::zipr(
    zipfile = output_tmp,
    files = archive_root,
    root = work,
    recurse = TRUE,
    compression_level = 1,
    include_directories = FALSE
  )
  if (!file.rename(output_tmp, output) && !file.copy(output_tmp, output, overwrite = TRUE)) {
    stop("Could not place the CPU-normalized checkpoint in the R cache.", call. = FALSE)
  }
  output
}
