test_that("all official scGPT model-zoo checkpoints are registered", {
  catalog <- fmverse_scgpt_checkpoints()
  expect_setequal(catalog$checkpoint, c(
    "whole-human", "continual-pretrained", "brain", "blood",
    "heart", "lung", "kidney", "pan-cancer"
  ))
  expect_true(all(nzchar(catalog$drive_folder_id)))
})

test_that("checkpoint aliases resolve consistently", {
  expect_equal(fmverse:::fmverse_normalize_scgpt_checkpoint_name("whole human"), "whole-human")
  expect_equal(fmverse:::fmverse_normalize_scgpt_checkpoint_name("scGPT_brain"), "brain")
  expect_equal(fmverse:::fmverse_normalize_scgpt_checkpoint_name("scGPT_bc"), "blood")
  expect_equal(fmverse:::fmverse_normalize_scgpt_checkpoint_name("pancancer"), "pan-cancer")
})

test_that("registered checkpoint paths can be selected from one checkpoint root", {
  root <- tempfile("fmverse-checkpoints-")
  dir.create(file.path(root, "scGPT-blood"), recursive = TRUE)
  resolved <- fmverse:::fmverse_resolve_scgpt_checkpoint("blood", checkpoint_dir = root)
  expect_equal(resolved$name, "blood")
  expect_equal(basename(resolved$path), "scGPT-blood")
})
