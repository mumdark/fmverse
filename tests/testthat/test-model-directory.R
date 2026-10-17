make_model_directory_fixture <- function(model, checkpoint) {
  root <- tempfile(paste0("fmverse-", model, "-"))
  dir.create(root)
  touch <- function(relative) {
    path <- file.path(root, relative)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines("fixture", path)
    relative
  }
  files <- switch(model,
    scGPT = list(weights = touch("best_model.pt"), config = touch("args.json"),
                 vocab = touch("vocab.json")),
    Geneformer = list(weights = touch("model.safetensors"), config = touch("config.json"),
      token_dictionary = touch("token_dictionary_gc104M.tsv"),
      gene_median_dictionary = touch("gene_median_dictionary_gc104M.tsv"),
      gene_name_id_dictionary = touch("gene_name_id_dict_gc104M.tsv"),
      ensembl_mapping_dictionary = touch("ensembl_mapping_dict_gc104M.tsv")),
    scBERT = list(weights = touch("scbert_state.safetensors"),
                  gene_index = touch("panglao_gene_index.tsv")),
    scFoundation = list(weights = touch("scfoundation_cell.safetensors"),
      config = touch("scfoundation_cell_config.json"),
      gene_index = touch("scfoundation_gene_index_19264.tsv")),
    UCE = {
      assets <- "assets"
      touch(file.path(assets, "species_gene_index.tsv"))
      touch(file.path(assets, "species_offsets.tsv"))
      touch(file.path(assets, "species_chrom.csv"))
      touch(file.path(assets, "all_tokens.safetensors"))
      for (file in fmverse_uce_species()$embedding_file) {
        touch(file.path(assets, "protein_embeddings", file))
      }
      list(weights = touch("4layer_model.torch"), assets_dir = assets)
    }
  )
  manifest <- list(schema_version = 1L, model = model, checkpoint = checkpoint,
                   files = files)
  if (identical(model, "UCE")) manifest$defaults <- list(species = "human")
  jsonlite::write_json(manifest, file.path(root, "fmverse.json"),
                       auto_unbox = TRUE, pretty = TRUE)
  root
}

make_dispatch_result <- function(object, model) {
  object <- fmverse_to_seurat_v5(object)
  cells <- colnames(object)
  embedding <- matrix(seq_len(length(cells) * 4L), nrow = length(cells),
    dimnames = list(cells, paste0(model, "_", seq_len(4L))))
  object[[model]] <- SeuratObject::CreateDimReducObject(embedding,
    key = paste0(gsub("[^A-Za-z0-9]", "", model), "_"),
    assay = SeuratObject::DefaultAssay(object))
  structure(list(model_name = model, embeddings = embedding, seurat = object,
    model = NULL, audit = list(), inputs = list(cell_names = cells),
    metadata = list(checkpoint_name = "fixture"), analyses = list(),
    plots = list(), history = list()), class = c("fmverse_result", "list"))
}

test_that("model manifests validate all registered adapters and relative assets", {
  cases <- list(scGPT = "whole-human", Geneformer = "Geneformer-V2-104M",
    scBERT = "panglao-human", scFoundation = "01B-resolution", UCE = "4layer")
  for (model in names(cases)) {
    root <- make_model_directory_fixture(model, cases[[model]])
    on.exit(unlink(root, recursive = TRUE), add = TRUE)
    info <- fmverse:::.fm_model_manifest(root)
    expect_identical(info$model, model)
    expect_identical(info$checkpoint, cases[[model]])
    expect_true(all(file.exists(unlist(info$files))))
  }
})

test_that("model directory errors identify missing, unsupported, and escaping resources", {
  expect_error(fmverse:::.fm_model_manifest(tempdir()), "missing `fmverse.json`")
  root <- make_model_directory_fixture("scGPT", "whole-human")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  manifest <- jsonlite::fromJSON(file.path(root, "fmverse.json"), simplifyVector = FALSE)
  manifest$checkpoint <- "not-a-checkpoint"
  jsonlite::write_json(manifest, file.path(root, "fmverse.json"), auto_unbox = TRUE)
  expect_error(fmverse:::.fm_model_manifest(root), "Unsupported `scGPT` checkpoint")
  manifest$checkpoint <- "whole-human"
  manifest$files$weights <- "../outside.pt"
  jsonlite::write_json(manifest, file.path(root, "fmverse.json"), auto_unbox = TRUE)
  expect_error(fmverse:::.fm_model_manifest(root), "must be relative|outside the model directory")
})

test_that("adapters honor declared nested resource paths in independent folders", {
  root <- make_model_directory_fixture("Geneformer", "Geneformer-V2-104M")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  manifest_path <- file.path(root, "fmverse.json")
  manifest <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
  for (role in names(manifest$files)) {
    old <- file.path(root, manifest$files[[role]])
    nested <- file.path("declared", paste0(role, ".resource"))
    new <- file.path(root, nested)
    dir.create(dirname(new), recursive = TRUE, showWarnings = FALSE)
    file.copy(old, new)
    unlink(old)
    manifest$files[[role]] <- nested
  }
  jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE, pretty = TRUE)
  jsonlite::write_json(list(vocab_size = 6L), file.path(root, manifest$files$config),
                       auto_unbox = TRUE)
  info <- fmverse:::.fm_model_manifest(root)
  expect_true(all(startsWith(unlist(info$files), normalizePath(root, winslash = "/"))))
  resolved <- fmverse:::fmverse_resolve_geneformer_checkpoint(
    info$files$weights, checkpoint_dir = info$root)
  expect_identical(resolved$config_file, info$files$config)
  expect_identical(resolved$token_dictionary, info$files$token_dictionary)
  expect_identical(resolved$gene_median_dictionary, info$files$gene_median_dictionary)
})

test_that("fmverse infers adapters from model_dir and creates automatic run IDs", {
  cases <- list(scGPT = "whole-human", Geneformer = "Geneformer-V2-104M",
    scBERT = "panglao-human", scFoundation = "01B-resolution", UCE = "4layer")
  calls <- list()
  fake_engine <- function(object, model, checkpoint, batch_size, ...) {
    calls[[model]] <<- list(model = model, checkpoint = checkpoint,
      batch_size = batch_size, options = list(...))
    make_dispatch_result(object, model)
  }
  testthat::local_mocked_bindings(.fmverse_infer_legacy = fake_engine,
                                  .package = "fmverse")
  counts <- Matrix::Matrix(matrix(c(1L, 2L, 3L, 4L, 1L, 2L), nrow = 2L,
    dimnames = list(c("g1", "g2"), paste0("c", seq_len(3L)))), sparse = TRUE)
  seurat <- Seurat::CreateSeuratObject(counts = counts)
  obj <- fmobj(seurat)
  for (model in names(cases)) {
    root <- make_model_directory_fixture(model, cases[[model]])
    on.exit(unlink(root, recursive = TRUE), add = TRUE)
    obj <- fmverse(obj, model_dir = root)
    run <- fm_runs(obj)[nrow(fm_runs(obj)), , drop = FALSE]
    expect_identical(run$model, model)
    expect_identical(run$checkpoint, cases[[model]])
    expect_identical(calls[[model]]$batch_size,
      if (model == "scGPT") 16L else if (model == "Geneformer") 8L else 1L)
    expect_true(file.exists(calls[[model]]$checkpoint))
    if (model == "scGPT") {
      expect_true(file.exists(calls[[model]]$options$vocab_file))
      expect_true(file.exists(calls[[model]]$options$config_file))
    }
    if (model %in% c("scBERT", "scFoundation")) {
      expect_true(file.exists(calls[[model]]$options$gene_index_file))
    }
    if (model == "scFoundation") expect_true(file.exists(calls[[model]]$options$config_file))
    if (model == "UCE") expect_true(dir.exists(calls[[model]]$options$assets_dir))
    if (model == "UCE") expect_identical(calls[[model]]$options$species, "human")
  }
  expect_identical(obj@active_run, "UCE_4layer")
  expect_error(fmverse(obj, model_dir = root, model = "UCE"),
    "Model selection and resource paths are read from `fmverse.json`")
  expect_error(fmverse(obj, model_dir = root, checkpoint = "4layer"),
    "Model selection and resource paths are read from `fmverse.json`")
})

test_that("folder API supports explicit IDs, collision suffixes, and typed shortcuts", {
  root <- make_model_directory_fixture("scBERT", "panglao-human")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  fake_engine <- function(object, model, checkpoint, batch_size, ...) {
    make_dispatch_result(object, model)
  }
  testthat::local_mocked_bindings(.fmverse_infer_legacy = fake_engine,
                                  .package = "fmverse")
  counts <- Matrix::Matrix(matrix(1:6, nrow = 2L,
    dimnames = list(c("g1", "g2"), paste0("c", seq_len(3L)))), sparse = TRUE)
  seurat <- Seurat::CreateSeuratObject(counts = counts)
  first <- fmverse(seurat, model_dir = root, run_id = "named")
  expect_identical(fm_runs(first)$run_id, "named")
  expect_error(fmverse(first, model_dir = root, run_id = "named"),
    "run_id already exists")
  replaced <- fmverse(first, model_dir = root, run_id = "named", overwrite = TRUE)
  expect_identical(fm_runs(replaced)$run_id, "named")
  second <- fmverse(first, model_dir = root)
  expect_identical(tail(fm_runs(second)$run_id, 1L), "scBERT_panglao-human")
  third <- fmverse(second, model_dir = root)
  expect_identical(tail(fm_runs(third)$run_id, 1L), "scBERT_panglao-human_2")
  expect_s4_class(fmverse_scbert(seurat, model_dir = root), "fmobj")
  expect_error(fmverse_geneformer(seurat, model_dir = root),
    "requires a Geneformer model directory")
})
