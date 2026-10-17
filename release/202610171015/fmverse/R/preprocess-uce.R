# Official UCE preprocessing and R-native embedding path.

fmverse_uce_tokenize <- function(object, species, gene_symbol = NULL, assay = NULL,
                                 cells = NULL, sample_size = 1024L,
                                 pad_length = 1536L, seed = 1L,
                                 assets = NULL) {
  species <- tolower(as.character(species))
  species_catalog <- fmverse_uce_species_catalog()
  if (length(species) != 1L || is.na(species) || !species %in% species_catalog$species) {
    stop("Unsupported UCE species. Choose one of: ",
      paste(species_catalog$species, collapse = ", "), call. = FALSE)
  }
  if (length(sample_size) != 1L || is.na(sample_size) || sample_size < 1L ||
      sample_size != as.integer(sample_size) || length(pad_length) != 1L ||
      is.na(pad_length) || pad_length < 2L || pad_length != as.integer(pad_length)) {
    stop("sample_size and pad_length must be positive whole numbers.", call. = FALSE)
  }
  if (length(seed) != 1L || is.na(seed) || !is.finite(seed) || seed < 0 ||
      seed > .Machine$integer.max || seed != as.integer(seed)) {
    stop("seed must be one non-negative integer no larger than .Machine$integer.max.",
      call. = FALSE)
  }
  if (is.null(assets)) stop("UCE assets must be resolved before tokenization.", call. = FALSE)
  if (!file.exists(assets$gene_index) || !file.exists(assets$chromosome)) {
    stop("UCE gene index or species chromosome table is missing.", call. = FALSE)
  }

  object <- fmverse_to_seurat_v5(object)
  if (is.null(assay)) assay <- SeuratObject::DefaultAssay(object)
  if (!assay %in% SeuratObject::Assays(object)) stop("Assay not found: ", assay, call. = FALSE)
  count_layers <- grep("^counts($|\\.)", SeuratObject::Layers(object[[assay]]), value = TRUE)
  if (!length(count_layers)) {
    stop("UCE requires raw counts in the Seurat `counts` layer.", call. = FALSE)
  }
  if (length(count_layers) != 1L || !identical(count_layers, "counts")) {
    object[[assay]] <- SeuratObject::JoinLayers(
      object[[assay]], layers = "counts", new = "counts"
    )
  }
  counts <- methods::as(SeuratObject::LayerData(object, assay = assay, layer = "counts"),
    "dgCMatrix")
  if (is.null(cells)) {
    cell_names <- colnames(counts)
  } else {
    cell_names <- as.character(cells)
    if (!length(cell_names) || anyNA(cell_names) || anyDuplicated(cell_names)) {
      stop("cells must be unique, non-missing Seurat cell names.", call. = FALSE)
    }
    if (any(!cell_names %in% colnames(counts))) stop("Some requested cells are absent.", call. = FALSE)
    counts <- counts[, cell_names, drop = FALSE]
  }
  if (any(!is.finite(counts@x)) || any(counts@x < 0)) {
    stop("UCE counts must be finite and non-negative.", call. = FALSE)
  }

  feature_metadata <- object[[assay]][[]]
  feature_metadata <- feature_metadata[match(rownames(counts), rownames(feature_metadata)), , drop = FALSE]
  if (nrow(feature_metadata) != nrow(counts) || anyNA(rownames(feature_metadata))) {
    stop("Could not align Seurat feature metadata to the counts layer.", call. = FALSE)
  }
  if (!is.null(gene_symbol)) {
    if (length(gene_symbol) != 1L || !gene_symbol %in% colnames(feature_metadata)) {
      stop("gene_symbol must name a feature metadata column in the selected assay.", call. = FALSE)
    }
    feature_names <- as.character(feature_metadata[[gene_symbol]])
  } else {
    candidates <- intersect(c("gene_symbol", "gene_symbols", "gene_name", "symbol", "scGPT_gene_name"),
      colnames(feature_metadata))
    feature_names <- if (length(candidates)) as.character(feature_metadata[[candidates[[1L]]]]) else
      rownames(counts)
  }
  if (length(feature_names) != nrow(counts) || anyNA(feature_names)) {
    stop("Gene-symbol metadata must contain one non-missing symbol per count row.", call. = FALSE)
  }
  feature_names <- toupper(trimws(feature_names))
  # Match official Python behavior for repeated symbols: each feature row is
  # retained and points to the first corresponding ESM2/token row.

  gene_index <- utils::read.delim(assets$gene_index, stringsAsFactors = FALSE,
    check.names = FALSE, quote = "", comment.char = "")
  required_index <- c("species", "gene_symbol", "pe_row_idx")
  if (!all(required_index %in% colnames(gene_index))) {
    stop("species_gene_index.tsv must contain columns: ",
      paste(required_index, collapse = ", "), call. = FALSE)
  }
  gene_index <- gene_index[gene_index$species == species, , drop = FALSE]
  gene_index$gene_symbol <- toupper(gene_index$gene_symbol)
  if (!nrow(gene_index) || anyDuplicated(gene_index$gene_symbol)) {
    stop("UCE gene index is empty or has duplicate symbols for species `", species, "`.", call. = FALSE)
  }
  chromosome_table <- utils::read.csv(assets$chromosome, stringsAsFactors = FALSE,
    check.names = FALSE)
  required_chrom <- c("gene_symbol", "chromosome", "start", "species")
  if (!all(required_chrom %in% colnames(chromosome_table))) {
    stop("species_chrom.csv is missing one of: ", paste(required_chrom, collapse = ", "),
      call. = FALSE)
  }
  chromosome_table$gene_symbol <- toupper(chromosome_table$gene_symbol)
  chromosome_table$species_chrom <- paste(chromosome_table$species,
    chromosome_table$chromosome, sep = "_")
  categories <- sort(unique(chromosome_table$species_chrom), method = "radix")
  chromosome_table$chrom_code <- match(chromosome_table$species_chrom, categories) - 1L
  chromosome_table <- chromosome_table[chromosome_table$species == species, , drop = FALSE]
  chromosome_table <- chromosome_table[!duplicated(chromosome_table$gene_symbol), , drop = FALSE]
  coords <- chromosome_table[match(gene_index$gene_symbol, chromosome_table$gene_symbol), , drop = FALSE]
  gene_index$chrom_code <- coords$chrom_code
  gene_index$start <- suppressWarnings(as.numeric(coords$start))
  gene_index$pe_row_idx <- suppressWarnings(as.integer(gene_index$pe_row_idx))
  usable_gene <- !is.na(gene_index$chrom_code) & is.finite(gene_index$start) &
    !is.na(gene_index$pe_row_idx)
  gene_index <- gene_index[usable_gene, , drop = FALSE]
  if (!nrow(gene_index)) stop("No UCE ESM2 gene rows have chromosome coordinates for `",
    species, "`.", call. = FALSE)
  matched <- match(feature_names, gene_index$gene_symbol)
  keep <- !is.na(matched)
  if (!any(keep)) stop("No Seurat feature symbols match UCE's `", species,
    "` ESM2 gene table. Provide gene symbols or set `gene_symbol` to the correct feature field.",
    call. = FALSE)
  matched_rows <- which(keep)
  counts <- counts[matched_rows, , drop = FALSE]
  matched_index <- matched[keep]
  pe_row_idx <- gene_index$pe_row_idx[matched_index]
  chrom_code <- as.integer(gene_index$chrom_code[matched_index])
  genomic_start <- gene_index$start[matched_index]
  if (any(pe_row_idx < 0L | pe_row_idx >= 143574L)) {
    stop("UCE protein-token row indices fall outside the official gene-token range.", call. = FALSE)
  }

  original_rng_kind <- RNGkind()
  had_rng_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_rng_seed) old_rng_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (!identical(RNGkind(), original_rng_kind)) do.call(RNGkind, as.list(original_rng_kind))
    if (had_rng_seed) assign(".Random.seed", old_rng_seed, envir = .GlobalEnv) else
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
  }, add = TRUE)
  set.seed(as.integer(seed))

  per_cell <- vector("list", ncol(counts))
  token_counts <- integer(ncol(counts))
  for (j in seq_len(ncol(counts))) {
    first <- counts@p[[j]] + 1L
    last <- counts@p[[j + 1L]]
    if (first > last) stop("Cell `", cell_names[[j]], "` has no counts in UCE-mapped genes.",
      call. = FALSE)
    offsets <- seq.int(first, last)
    row_index <- counts@i[offsets] + 1L
    expression <- as.numeric(counts@x[offsets])
    active <- expression > 0
    row_index <- row_index[active]
    expression <- expression[active]
    if (!length(expression)) stop("Cell `", cell_names[[j]],
      "` has zero counts in UCE-mapped genes.", call. = FALSE)
    weights <- log1p(expression)
    weights <- weights / sum(weights)
    choice_idx <- sample.int(length(weights), size = as.integer(sample_size),
      replace = TRUE, prob = weights)
    chosen_genes <- row_index[choice_idx]
    chosen_chrom <- chrom_code[chosen_genes]
    chosen_start <- genomic_start[chosen_genes]
    chosen_pe <- pe_row_idx[chosen_genes]
    chrom_order <- order(chosen_chrom, method = "radix")
    chosen_genes <- chosen_genes[chrom_order]
    chosen_chrom <- chosen_chrom[chrom_order]
    chosen_start <- chosen_start[chrom_order]
    chosen_pe <- chosen_pe[chrom_order]
    unique_chrom <- sort(unique(chosen_chrom), method = "radix")
    unique_chrom <- unique_chrom[sample.int(length(unique_chrom), length(unique_chrom))]
    sentence <- 3L # CLS
    for (chrom in unique_chrom) {
      in_chrom <- which(chosen_chrom == chrom)
      in_chrom <- in_chrom[order(chosen_start[in_chrom], method = "radix")]
      sentence <- c(sentence, as.integer(chrom + 143574L), chosen_pe[in_chrom], 2L)
    }
    if (length(sentence) > as.integer(pad_length)) {
      stop("UCE sequence exceeded pad_length for cell `", cell_names[[j]], "`.", call. = FALSE)
    }
    per_cell[[j]] <- as.integer(sentence)
    token_counts[[j]] <- length(sentence)
  }
  width <- max(token_counts)
  token_ids <- matrix(0L, nrow = length(cell_names), ncol = width)
  attention_mask <- matrix(0L, nrow = length(cell_names), ncol = width)
  for (j in seq_along(per_cell)) {
    n <- token_counts[[j]]
    token_ids[j, seq_len(n)] <- per_cell[[j]]
    attention_mask[j, seq_len(n)] <- 1L
  }
  rownames(token_ids) <- rownames(attention_mask) <- cell_names
  list(
    token_ids = token_ids, attention_mask = attention_mask,
    token_counts = token_counts, cell_names = cell_names, assay = assay,
    species = species, gene_names = feature_names[keep],
    matched_features = sum(keep), total_features = nrow(object[[assay]]),
    preprocessing = list(
      value_source = "raw counts layer", method = "expression-weighted gene sampling with replacement",
      sample_size = as.integer(sample_size), pad_length = as.integer(pad_length),
      sequence = "CLS, shuffled chromosome groups; genes ordered by genomic start; chromosome separator; padding",
      seed = as.integer(seed), rng = "R RNG; official Python uses NumPy RNG",
      species = species
    )
  )
}

#' Run UCE inference on a Seurat object
#'
#' @param object Seurat v4/v5 input.
#' @param checkpoint Official `4layer` or `33layer` model.
#' @param species One of the species bundled with UCE.
#' @param batch_size Number of cells evaluated per forward pass.
#' @param sample_size Number of expression-weighted gene samples per cell.
#' @param pad_length Maximum UCE token sequence length.
#' @param seed Seed for R's UCE gene sampling.
#' @param checkpoint_dir Optional root for local UCE checkpoint resources.
#' @param assets_dir Internal explicit path to UCE's self-contained assets.
#' @param gene_symbol Optional assay feature metadata column containing symbols.
#' @param cells Optional Seurat cell barcodes to embed.
#' @param assay Counts assay.
#' @param device Torch device (`cpu` by default).
#' @param verbose Report progress.
#' @return A structured `fmverse_result` with a 1,280-dimensional UCE embedding.
.fmverse_uce_legacy <- function(object, checkpoint = "4layer", species = "human",
                                checkpoint_dir = NULL, batch_size = 1L,
                                assets_dir = NULL,
                                sample_size = 1024L, pad_length = 1536L,
                                seed = 1L, gene_symbol = NULL, cells = NULL,
                                assay = NULL, device = "cpu", verbose = interactive()) {
  if (length(batch_size) != 1L || is.na(batch_size) || batch_size < 1L ||
      batch_size != as.integer(batch_size)) {
    stop("batch_size must be a positive integer.", call. = FALSE)
  }
  resolved <- fmverse_resolve_uce_checkpoint(checkpoint, checkpoint_dir)
  assets <- fmverse_resolve_uce_assets(resolved$path, checkpoint_dir, assets_dir)
  species_record <- fmverse_uce_species_catalog()
  species_record <- species_record[species_record$species == tolower(species), , drop = FALSE]
  if (!nrow(species_record)) stop("Unsupported UCE species: ", species, call. = FALSE)
  tokenized <- fmverse_uce_tokenize(
    object, species = species, gene_symbol = gene_symbol, assay = assay,
    cells = cells, sample_size = sample_size, pad_length = pad_length,
    seed = seed, assets = assets
  )
  model <- fmverse_uce_build_model(resolved$layers)
  loaded <- fmverse_load_uce_weights(model, resolved$path)
  model <- loaded$model
  torch_device <- torch::torch_device(device)
  model$to(device = torch_device)
  model$eval()
  n_cells <- length(tokenized$cell_names)
  embeddings <- matrix(NA_real_, nrow = n_cells, ncol = 1280L,
    dimnames = list(tokenized$cell_names, paste0("UCE_", seq_len(1280L))))
  starts <- seq.int(1L, n_cells, by = as.integer(batch_size))
  progress_step <- max(1L, ceiling(length(starts) / 20L))
  for (batch_number in seq_along(starts)) {
    first <- starts[[batch_number]]
    index <- first:min(first + as.integer(batch_size) - 1L, n_cells)
    sequence_end <- max(tokenized$token_counts[index])
    sequence_index <- seq_len(sequence_end)
    input_ids <- torch::torch_tensor(
      t(tokenized$token_ids[index, sequence_index, drop = FALSE]),
      dtype = torch::torch_long(), device = torch_device
    )
    attention_mask <- torch::torch_tensor(
      tokenized$attention_mask[index, sequence_index, drop = FALSE],
      dtype = torch::torch_float(), device = torch_device
    )
    output <- torch::with_no_grad(model(input_ids, attention_mask))
    values <- torch::as_array(output$cell_embedding$cpu())
    embeddings[index, ] <- matrix(values, nrow = length(index), ncol = 1280L)
    if (isTRUE(verbose) && (batch_number %% progress_step == 0L || batch_number == length(starts))) {
      message(sprintf("UCE embedded %d/%d cells (%s, %s)", max(index), n_cells,
        resolved$name, species))
    }
    rm(output, input_ids, attention_mask, values)
    if (batch_number %% 5L == 0L) invisible(gc(verbose = FALSE))
  }
  object <- fmverse_to_seurat_v5(object)
  if (!is.null(cells)) object <- subset(object, cells = tokenized$cell_names)
  object[["UCE"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embeddings, key = "UCE_", assay = tokenized$assay
  )
  list(
    model_name = "UCE", embeddings = embeddings, seurat = object, model = model,
    audit = loaded$audit,
    inputs = list(
      token_ids = tokenized$token_ids,
      attention_mask = tokenized$attention_mask,
      token_counts = tokenized$token_counts,
      cell_names = tokenized$cell_names,
      gene_names = tokenized$gene_names,
      assay = tokenized$assay, species = tokenized$species
    ),
    analyses = list(), plots = list(), history = list(),
    metadata = list(
      checkpoint_name = resolved$name, checkpoint = resolved$path,
      layers = resolved$layers, species = species_record$species[[1L]],
      scientific_name = species_record$scientific_name[[1L]],
      genome_assembly = species_record$genome_assembly[[1L]],
      device = device, batch_size = as.integer(batch_size),
      embedding_dim = 1280L, n_cells = n_cells,
      matched_features = tokenized$matched_features,
      total_features = tokenized$total_features,
      token_counts = tokenized$token_counts,
      preprocessing = tokenized$preprocessing,
      assets = assets[c("gene_index", "offsets", "chromosome", "directory")]
    )
  ) |>
    structure(class = c("fmverse_result", "list"))
}
