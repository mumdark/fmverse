# Model-folder metadata and validation for the folder-based public API.

.fm_model_manifest <- function(model_dir) {
  if (!is.character(model_dir) || length(model_dir) != 1L || is.na(model_dir) ||
      !nzchar(model_dir) || !dir.exists(model_dir)) {
    stop("`model_dir` must name an existing model directory.", call. = FALSE)
  }
  root <- normalizePath(path.expand(model_dir), winslash = "/", mustWork = TRUE)
  manifest_file <- file.path(root, "fmverse.json")
  if (!file.exists(manifest_file)) {
    stop("Model directory is missing `fmverse.json`: ", manifest_file,
         call. = FALSE)
  }
  manifest <- tryCatch(
    jsonlite::fromJSON(manifest_file, simplifyVector = FALSE),
    error = function(e) stop("Cannot read `fmverse.json`: ", conditionMessage(e),
                             call. = FALSE)
  )
  if (!is.list(manifest) || !identical(as.integer(manifest$schema_version), 1L)) {
    stop("`fmverse.json` must declare `schema_version: 1`.", call. = FALSE)
  }
  model <- manifest$model
  checkpoint <- manifest$checkpoint
  if (!is.character(model) || length(model) != 1L || is.na(model) || !nzchar(model)) {
    stop("`fmverse.json` is missing a valid `model` string.", call. = FALSE)
  }
  if (!model %in% c("scGPT", "Geneformer", "scBERT", "scFoundation", "UCE")) {
    stop("Unsupported model in `fmverse.json`: `", model, "`.", call. = FALSE)
  }
  if (!is.character(checkpoint) || length(checkpoint) != 1L || is.na(checkpoint) ||
      !nzchar(checkpoint)) {
    stop("`fmverse.json` is missing a valid `checkpoint` string.", call. = FALSE)
  }

  allowed <- switch(model,
    scGPT = fmverse_scgpt_checkpoint_catalog()$checkpoint,
    Geneformer = c(fmverse_geneformer_checkpoint_catalog()$checkpoint,
                   "gf-12L-95M-i4096", "geneformer-12L-30M"),
    scBERT = "panglao-human",
    scFoundation = "01B-resolution",
    UCE = fmverse_uce_checkpoint_catalog()$checkpoint
  )
  if (!checkpoint %in% allowed) {
    stop("Unsupported `", model, "` checkpoint in `fmverse.json`: `", checkpoint,
         "`. Supported values: ", paste(allowed, collapse = ", "), call. = FALSE)
  }

  files <- manifest$files
  if (!is.list(files)) stop("`fmverse.json` must contain a `files` object.", call. = FALSE)
  required <- switch(model,
    scGPT = c("weights", "config", "vocab"),
    Geneformer = c("weights", "config", "token_dictionary", "gene_median_dictionary",
                   "gene_name_id_dictionary", "ensembl_mapping_dictionary"),
    scBERT = c("weights", "gene_index"),
    scFoundation = c("weights", "config", "gene_index"),
    UCE = c("weights", "assets_dir")
  )
  absent <- setdiff(required, names(files))
  if (length(absent)) {
    stop("`fmverse.json` for ", model, " is missing file entries: ",
         paste(absent, collapse = ", "), call. = FALSE)
  }
  resolved <- list()
  root_lower <- tolower(sub("/$", "", root))
  for (role in required) {
    relative <- files[[role]]
    if (!is.character(relative) || length(relative) != 1L || is.na(relative) ||
        !nzchar(relative)) {
      stop("`fmverse.json` field `files.", role, "` must be a relative path.",
           call. = FALSE)
    }
    relative <- gsub("\\\\", "/", relative)
    if (grepl("^(?:/|[A-Za-z]:/|//)", relative)) {
      stop("`files.", role, "` must be relative to the model directory.", call. = FALSE)
    }
    path <- normalizePath(file.path(root, relative), winslash = "/", mustWork = FALSE)
    path_lower <- tolower(path)
    if (!startsWith(path_lower, paste0(root_lower, "/"))) {
      stop("`files.", role, "` points outside the model directory.", call. = FALSE)
    }
    exists <- if (identical(role, "assets_dir")) dir.exists(path) else
      file.exists(path) && !dir.exists(path)
    if (!exists) stop("Missing model resource `files.", role, "`: ", path,
                      call. = FALSE)
    resolved[[role]] <- normalizePath(path, winslash = "/", mustWork = TRUE)
  }
  if (identical(model, "UCE")) {
    assets <- resolved$assets_dir
    asset_files <- c("species_gene_index.tsv", "species_offsets.tsv", "species_chrom.csv",
                     "all_tokens.safetensors")
    missing_assets <- file.path(assets, asset_files)[!file.exists(file.path(assets, asset_files))]
    species <- fmverse_uce_species_catalog()
    embedding_files <- file.path(assets, "protein_embeddings", species$embedding_file)
    missing_embeddings <- embedding_files[!file.exists(embedding_files)]
    if (length(missing_assets) || length(missing_embeddings)) {
      missing <- c(missing_assets, missing_embeddings)
      stop("UCE model directory is incomplete; missing assets: ",
           paste(missing, collapse = ", "), call. = FALSE)
    }
  }

  defaults <- manifest$defaults
  if (identical(model, "UCE")) {
    default_species <- if (is.list(defaults)) defaults$species else NULL
    if (is.null(default_species)) default_species <- "human"
    if (!is.character(default_species) || length(default_species) != 1L ||
        !default_species %in% fmverse_uce_species_catalog()$species) {
      stop("`defaults.species` in the UCE `fmverse.json` is not supported.", call. = FALSE)
    }
    defaults <- list(species = default_species)
  }
  list(root = root, manifest = manifest, model = model, checkpoint = checkpoint,
       files = resolved, defaults = defaults)
}
