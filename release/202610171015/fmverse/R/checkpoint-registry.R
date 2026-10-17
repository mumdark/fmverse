fmverse_scgpt_checkpoint_catalog <- function() {
  data.frame(
    checkpoint = c(
      "whole-human", "continual-pretrained", "brain", "blood",
      "heart", "lung", "kidney", "pan-cancer"
    ),
    official_name = c(
      "scGPT_human", "scGPT_CP", "scGPT_brain", "scGPT_bc",
      "scGPT_heart", "scGPT_lung", "scGPT_kidney", "scGPT_pancancer"
    ),
    description = c(
      "33 million normal human cells; recommended general-purpose model",
      "Continual pretraining for zero-shot cell embedding",
      "13.2 million brain cells", "10.3 million blood and bone marrow cells",
      "1.8 million heart cells", "2.1 million lung cells",
      "814 thousand kidney cells", "5.7 million cancer cells"
    ),
    folder = paste0("scGPT-", c(
      "whole-human", "continual-pretrained", "brain", "blood",
      "heart", "lung", "kidney", "pan-cancer"
    )),
    drive_folder_id = c(
      "1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y",
      "1_GROJTzXiAV8HB4imruOTk6PEGuNOcgB",
      "1vf1ijfQSk7rGdDGpBntR5bi5g6gNt-Gx",
      "1kkug5C7NjvXIwQGGaGoqXTk_Lb_pDrBU",
      "1GcgXrd7apn6y4Ze_iSCncskX3UsWPY2r",
      "16A1DJ30PT6bodt4bWLa4hpS7gbWZQFBG",
      "1S-1AR65DF120kNFpEbWCvRHPhpkGK3kK",
      "13QzLHilYUd0v3HTwa_9n4G4yEF-hdkqa"
    ),
    official_url = paste0("https://drive.google.com/drive/folders/", c(
      "1oWh_-ZRdhtoGQ2Fw24HP41FgLoomVo-y",
      "1_GROJTzXiAV8HB4imruOTk6PEGuNOcgB",
      "1vf1ijfQSk7rGdDGpBntR5bi5g6gNt-Gx",
      "1kkug5C7NjvXIwQGGaGoqXTk_Lb_pDrBU",
      "1GcgXrd7apn6y4Ze_iSCncskX3UsWPY2r",
      "16A1DJ30PT6bodt4bWLa4hpS7gbWZQFBG",
      "1S-1AR65DF120kNFpEbWCvRHPhpkGK3kK",
      "13QzLHilYUd0v3HTwa_9n4G4yEF-hdkqa"
    )),
    stringsAsFactors = FALSE
  )
}

#' List the supported official scGPT checkpoint names
#'
#' @return A data frame with canonical checkpoint names, descriptions, and
#'   official Google Drive folder identifiers.
#' @export
fmverse_scgpt_checkpoints <- function() {
  fmverse_scgpt_checkpoint_catalog()
}

fmverse_scgpt_checkpoint_aliases <- function() {
  aliases <- c(
    "whole human" = "whole-human", "human" = "whole-human",
    "scgpt_human" = "whole-human", "scgpt-whole-human" = "whole-human",
    "continual pretrained" = "continual-pretrained",
    "continual" = "continual-pretrained", "scgpt_cp" = "continual-pretrained",
    "scgpt-continual-pretrained" = "continual-pretrained",
    "scgpt_brain" = "brain", "scgpt_bc" = "blood", "bc" = "blood",
    "scgpt_blood" = "blood",
    "scgpt_heart" = "heart", "scgpt_lung" = "lung",
    "scgpt_kidney" = "kidney", "pancancer" = "pan-cancer",
    "pan cancer" = "pan-cancer", "scgpt_pancancer" = "pan-cancer"
  )
  aliases
}

fmverse_normalize_scgpt_checkpoint_name <- function(x) {
  x <- trimws(tolower(as.character(x)))
  aliases <- fmverse_scgpt_checkpoint_aliases()
  if (x %in% names(aliases)) return(unname(aliases[[x]]))
  canonical <- fmverse_scgpt_checkpoint_catalog()$checkpoint
  if (x %in% canonical) return(x)
  normalized <- gsub("[_[:space:]]+", "-", x)
  normalized <- gsub("-+", "-", normalized)
  normalized <- sub("^scgpt-", "", normalized)
  if (normalized %in% canonical) return(normalized)
  NA_character_
}

fmverse_scgpt_checkpoint_roots <- function(checkpoint_dir = NULL) {
  roots <- c(
    checkpoint_dir,
    getOption("fmverse.checkpoint_dir"),
    Sys.getenv("FMVERSE_CHECKPOINT_DIR", unset = ""),
    file.path(getwd(), "00scGPT", "checkpoint"),
    file.path(getwd(), "..", "00scGPT", "checkpoint"),
    file.path(getwd(), "checkpoint"),
    file.path(tools::R_user_dir("fmverse", "data"), "checkpoints")
  )
  roots <- roots[!is.na(roots) & nzchar(roots)]
  unique(normalizePath(roots, winslash = "/", mustWork = FALSE))
}

fmverse_resolve_scgpt_checkpoint <- function(checkpoint, checkpoint_dir = NULL) {
  if (length(checkpoint) != 1L || is.na(checkpoint) || !nzchar(checkpoint)) {
    stop("`checkpoint` must be one checkpoint name or a local path.", call. = FALSE)
  }
  checkpoint <- path.expand(as.character(checkpoint))
  if (file.exists(checkpoint) || dir.exists(checkpoint)) {
    return(list(
      path = normalizePath(checkpoint, winslash = "/", mustWork = TRUE),
      name = NA_character_
    ))
  }

  name <- fmverse_normalize_scgpt_checkpoint_name(checkpoint)
  if (is.na(name)) {
    stop(
      "Unknown scGPT checkpoint `", checkpoint,
      "`. Use one of: ", paste(fmverse_scgpt_checkpoint_catalog()$checkpoint, collapse = ", "),
      ", or provide a local checkpoint path.",
      call. = FALSE
    )
  }
  record <- fmverse_scgpt_checkpoint_catalog()[
    fmverse_scgpt_checkpoint_catalog()$checkpoint == name, , drop = FALSE
  ]
  roots <- fmverse_scgpt_checkpoint_roots(checkpoint_dir)
  candidates <- unique(c(
    file.path(roots, record$folder[[1L]]),
    file.path(roots, name)
  ))
  found <- candidates[dir.exists(candidates)]
  if (!length(found)) {
    expected <- file.path(roots, record$folder[[1L]])
    stop(
      "Checkpoint `", name, "` is registered but not present locally. Expected a folder ",
      "containing best_model.pt, args.json, and vocab.json under one of: ",
      paste(expected, collapse = "; "),
      ". Download it from the official folder: https://drive.google.com/drive/folders/",
      record$drive_folder_id[[1L]],
      call. = FALSE
    )
  }
  list(path = normalizePath(found[[1L]], winslash = "/", mustWork = TRUE), name = name)
}

fmverse_scbert_checkpoint_catalog <- function() {
  data.frame(
    checkpoint = "panglao-human",
    official_name = "panglao_pretrain.pth",
    description = "Official 16,906-gene human Panglao pretraining checkpoint",
    genes = 16906L,
    hidden_size = 200L,
    layers = 6L,
    official_url = "https://drive.weixin.qq.com/s?k=AJEAIQdfAAoUxhXE7r",
    stringsAsFactors = FALSE
  )
}

#' List the scBERT checkpoint supported by fmverse
#' @return A data frame describing the official human Panglao checkpoint.
#' @export
fmverse_scbert_checkpoints <- function() fmverse_scbert_checkpoint_catalog()

fmverse_scbert_roots <- function(checkpoint_dir = NULL) {
  roots <- c(checkpoint_dir, getOption("fmverse.scbert_checkpoint_dir"),
    Sys.getenv("FMVERSE_SCBERT_CHECKPOINT_DIR", unset = ""),
    file.path(getwd(), "00scBERT", "checkpoint"),
    file.path(getwd(), "..", "00scBERT", "checkpoint"),
    file.path(tools::R_user_dir("fmverse", "data"), "checkpoints", "scBERT"))
  unique(normalizePath(roots[!is.na(roots) & nzchar(roots)],
    winslash = "/", mustWork = FALSE))
}

fmverse_resolve_scbert_checkpoint <- function(checkpoint, checkpoint_dir = NULL) {
  if (length(checkpoint) != 1L || is.na(checkpoint) || !nzchar(checkpoint)) {
    stop("`checkpoint` must be one scBERT checkpoint name or local path.", call. = FALSE)
  }
  checkpoint <- path.expand(as.character(checkpoint))
  if (file.exists(checkpoint)) {
    return(list(path = normalizePath(checkpoint, winslash = "/", mustWork = TRUE),
      name = "panglao-human"))
  }
  aliases <- c("panglao" = "panglao-human", "human" = "panglao-human",
    "panglao_pretrain.pth" = "panglao-human")
  name <- unname(aliases[tolower(checkpoint)])
  if (!length(name) || is.na(name)) name <- checkpoint
  if (!identical(name, "panglao-human")) {
    stop("Unknown scBERT checkpoint `", checkpoint, "`; supported checkpoint: panglao-human.",
      call. = FALSE)
  }
  roots <- fmverse_scbert_roots(checkpoint_dir)
  candidates <- unlist(lapply(roots, function(root) file.path(root, c(
    "scbert_state.safetensors", "panglao_pretrained.safetensors",
    "scbert_state_dict.pt", "panglao_pretrained_state_dict.pt"))), use.names = FALSE)
  found <- candidates[file.exists(candidates)]
  if (!length(found)) stop("scBERT weights not found. Put the official checkpoint converted to a flat ",
    "state dict at `scbert_state.safetensors` in 00scBERT/checkpoint, or set ",
    "FMVERSE_SCBERT_CHECKPOINT_DIR. Official source: ",
    fmverse_scbert_checkpoint_catalog()$official_url[[1L]], call. = FALSE)
  list(path = normalizePath(found[[1L]], winslash = "/", mustWork = TRUE), name = name)
}

fmverse_scfoundation_checkpoint_catalog <- function() {
  data.frame(
    checkpoint = "01B-resolution",
    version = "ce",
    checkpoint_key = "cell",
    description = "Official human scFoundation general cell-embedding checkpoint",
    genes = 19264L,
    official_url = paste0("https://hopebio2020.sharepoint.com/:f:/s/PublicSharedfiles/",
      "IgBlEJ72TBE5Q76AmgXbgjXiAR69fzcrgzqgUYdSThPLrqk"),
    stringsAsFactors = FALSE
  )
}

#' List scFoundation checkpoints supported by fmverse
#' @return A data frame describing the official general cell checkpoint.
#' @export
fmverse_scfoundation_checkpoints <- function() fmverse_scfoundation_checkpoint_catalog()

fmverse_scfoundation_roots <- function(checkpoint_dir = NULL) {
  roots <- c(checkpoint_dir, getOption("fmverse.scfoundation_checkpoint_dir"),
    Sys.getenv("FMVERSE_SCF_FOUNDATION_CHECKPOINT_DIR", unset = ""),
    file.path(getwd(), "00scFoundation", "checkpoint"),
    file.path(getwd(), "..", "00scFoundation", "checkpoint"),
    file.path(tools::R_user_dir("fmverse", "data"), "checkpoints", "scFoundation"))
  unique(normalizePath(roots[!is.na(roots) & nzchar(roots)],
    winslash = "/", mustWork = FALSE))
}

fmverse_resolve_scfoundation_checkpoint <- function(checkpoint, checkpoint_dir = NULL) {
  if (length(checkpoint) != 1L || is.na(checkpoint) || !nzchar(checkpoint)) {
    stop("`checkpoint` must be one scFoundation checkpoint name or local path.", call. = FALSE)
  }
  checkpoint <- path.expand(as.character(checkpoint))
  if (file.exists(checkpoint)) {
    weight <- normalizePath(checkpoint, winslash = "/", mustWork = TRUE)
    return(list(path = weight, directory = dirname(weight), name = "01B-resolution"))
  }
  name <- tolower(gsub("[_ ]+", "-", checkpoint))
  if (!name %in% c("01b-resolution", "cell", "ce")) {
    stop("Unknown scFoundation checkpoint `", checkpoint,
      "; supported checkpoint: 01B-resolution (cell embedding).", call. = FALSE)
  }
  roots <- fmverse_scfoundation_roots(checkpoint_dir)
  candidates <- unlist(lapply(roots, function(root) c(
    file.path(root, "scfoundation_cell.safetensors"),
    file.path(root, "01B-resolution", "scfoundation_cell.safetensors"),
    file.path(root, "model.safetensors"))), use.names = FALSE)
  found <- candidates[file.exists(candidates)]
  if (!length(found)) stop("scFoundation cell weights not found. Export the official nested `cell` ",
    "weight group to `scfoundation_cell.safetensors` and place it under 00scFoundation/checkpoint ",
      "or set FMVERSE_SCF_FOUNDATION_CHECKPOINT_DIR. Official source: ",
    fmverse_scfoundation_checkpoint_catalog()$official_url[[1L]], call. = FALSE)
  path <- normalizePath(found[[1L]], winslash = "/", mustWork = TRUE)
  list(path = path, directory = dirname(path), name = "01B-resolution")
}

fmverse_uce_checkpoint_catalog <- function() {
  data.frame(
    checkpoint = c("4layer", "33layer"),
    official_name = c("4layer_model.torch", "33l_8ep_1024t_1280.torch"),
    layers = c(4L, 33L),
    description = c(
      "Official four-layer example checkpoint",
      "Official 33-layer paper checkpoint"
    ),
    figshare_file_id = c("42706576", "43423236"),
    official_url = c(
      "https://ndownloader.figshare.com/files/42706576",
      "https://ndownloader.figshare.com/files/43423236"
    ),
    stringsAsFactors = FALSE
  )
}

#' List official UCE checkpoints supported by fmverse
#'
#' @return A data frame with official checkpoint names, layer counts, and URLs.
#' @export
fmverse_uce_checkpoints <- function() fmverse_uce_checkpoint_catalog()

fmverse_uce_species_catalog <- function() {
  data.frame(
    species = c("human", "mouse", "frog", "zebrafish", "mouse_lemur", "pig",
                "macaca_fascicularis", "macaca_mulatta"),
    scientific_name = c(
      "Homo sapiens", "Mus musculus", "Xenopus tropicalis", "Danio rerio",
      "Microcebus murinus", "Sus scrofa", "Macaca fascicularis", "Macaca mulatta"
    ),
    genome_assembly = c(
      "GRCh38", "GRCm39", "Xenopus_tropicalis_v9.1", "GRCz11", "Mmur_3.0",
      "Sscrofa11.1", "Macaca_fascicularis_6.0", "Mmul_10"
    ),
    embedding_file = c(
      "Homo_sapiens.GRCh38.gene_symbol_to_embedding_ESM2.pt",
      "Mus_musculus.GRCm39.gene_symbol_to_embedding_ESM2.pt",
      "Xenopus_tropicalis.Xenopus_tropicalis_v9.1.gene_symbol_to_embedding_ESM2.pt",
      "Danio_rerio.GRCz11.gene_symbol_to_embedding_ESM2.pt",
      "Microcebus_murinus.Mmur_3.0.gene_symbol_to_embedding_ESM2.pt",
      "Sus_scrofa.Sscrofa11.1.gene_symbol_to_embedding_ESM2.pt",
      "Macaca_fascicularis.Macaca_fascicularis_6.0.gene_symbol_to_embedding_ESM2.pt",
      "Macaca_mulatta.Mmul_10.gene_symbol_to_embedding_ESM2.pt"
    ),
    stringsAsFactors = FALSE
  )
}

#' List UCE species bundled with the official checkpoint release
#'
#' @return A data frame with species keys and genome assemblies.
#' @export
fmverse_uce_species <- function() fmverse_uce_species_catalog()

fmverse_uce_roots <- function(checkpoint_dir = NULL) {
  roots <- c(
    checkpoint_dir,
    getOption("fmverse.uce_checkpoint_dir"),
    Sys.getenv("FMVERSE_UCE_CHECKPOINT_DIR", unset = ""),
    file.path(getwd(), "00UCE", "checkpoint"),
    file.path(getwd(), "..", "00UCE", "checkpoint"),
    file.path(tools::R_user_dir("fmverse", "data"), "checkpoints", "UCE")
  )
  unique(normalizePath(roots[!is.na(roots) & nzchar(roots)],
    winslash = "/", mustWork = FALSE))
}

fmverse_resolve_uce_checkpoint <- function(checkpoint = "4layer", checkpoint_dir = NULL) {
  if (length(checkpoint) != 1L || is.na(checkpoint) || !nzchar(checkpoint)) {
    stop("`checkpoint` must be one UCE version or a local checkpoint path.", call. = FALSE)
  }
  checkpoint <- path.expand(as.character(checkpoint))
  catalog <- fmverse_uce_checkpoint_catalog()
  if (file.exists(checkpoint)) {
    path <- normalizePath(checkpoint, winslash = "/", mustWork = TRUE)
    hit <- match(basename(path), catalog$official_name)
    if (is.na(hit)) stop("Cannot infer UCE layer count from checkpoint filename.", call. = FALSE)
    return(list(path = path, directory = dirname(path), name = catalog$checkpoint[[hit]],
      layers = catalog$layers[[hit]]))
  }
  name <- tolower(gsub("[-_ ]", "", checkpoint))
  aliases <- c("4" = "4layer", "4layer" = "4layer", "sample" = "4layer",
    "wholehuman" = "4layer", "33" = "33layer", "33layer" = "33layer", "paper" = "33layer")
  alias_hit <- unname(aliases[name])
  if (!length(alias_hit) || is.na(alias_hit)) {
    stop("Unknown UCE checkpoint `", checkpoint, "`; choose `4layer` or `33layer`.", call. = FALSE)
  }
  name <- alias_hit[[1L]]
  record <- catalog[catalog$checkpoint == name, , drop = FALSE]
  roots <- fmverse_uce_roots(checkpoint_dir)
  version_folder <- paste0("uce-", name)
  candidates <- unique(c(
    unlist(lapply(roots, function(root) file.path(root, version_folder, record$official_name)),
      use.names = FALSE),
    unlist(lapply(roots, function(root) file.path(root, record$official_name)), use.names = FALSE)
  ))
  found <- candidates[file.exists(candidates)]
  if (!length(found)) {
    stop("UCE checkpoint `", name, "` is not present. Expected `", record$official_name[[1L]],
      "` under 00UCE/checkpoint/", version_folder, ". Official download: ",
      record$official_url[[1L]], call. = FALSE)
  }
  path <- normalizePath(found[[1L]], winslash = "/", mustWork = TRUE)
  list(path = path, directory = dirname(path), name = name, layers = record$layers[[1L]])
}

fmverse_resolve_uce_assets <- function(checkpoint_path, checkpoint_dir = NULL,
                                      assets_dir = NULL) {
  roots <- fmverse_uce_roots(checkpoint_dir)
  version_dir <- dirname(normalizePath(checkpoint_path, winslash = "/", mustWork = TRUE))
  candidates <- unique(c(
    assets_dir,
    file.path(version_dir, "assets"),
    unlist(lapply(roots, file.path, "shared"), use.names = FALSE)
  ))
  required <- c("species_gene_index.tsv", "species_offsets.tsv", "species_chrom.csv")
  for (candidate in candidates) {
    if (all(file.exists(file.path(candidate, required)))) {
      return(list(
        directory = normalizePath(candidate, winslash = "/", mustWork = TRUE),
        gene_index = file.path(candidate, required[[1L]]),
        offsets = file.path(candidate, required[[2L]]),
        chromosome = file.path(candidate, required[[3L]]),
        token_file = file.path(candidate, "all_tokens.safetensors"),
        embedding_directory = file.path(candidate, "protein_embeddings")
      ))
    }
  }
  stop("UCE species assets are incomplete. Expected species_gene_index.tsv, ",
    "species_offsets.tsv, species_chrom.csv, and protein_embeddings/ beside the checkpoint.",
    call. = FALSE)
}
