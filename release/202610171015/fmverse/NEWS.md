# fmverse 0.0.0.9005

- Select models through a self-contained folder with a versioned `fmverse.json`;
  infer model and checkpoint without public `model` or `checkpoint` arguments.
- Add manifests to all loadable local checkpoint folders. Bundle Geneformer
  dictionaries and scBERT/scFoundation gene indexes beside their weights.
- Preserve the batch axis in scGPT singleton inference batches and add a
  single-cell forward-pass regression test.
- Preserve model-specific inference defaults, optional `run_id` generation,
  overwrite behavior, S4 runs, and all downstream analysis functions.
- Add manifest/path validation tests and update public examples to `model_dir`.

# fmverse 0.0.0.9003

- Add R-native Performer adapters for scBERT `panglao-human` and scFoundation
  `01B-resolution`, including checkpoint catalogs, strict tensor audits,
  fixed-gene preprocessing and embeddings registered as independent `fmobj` runs.
- Add one-time checkpoint-format conversion and Panglao gene-order export helpers
  under `00scBERT/` and `00scFoundation/`; model weights remain user-provided.
- Reuse the common cluster, annotation, integration, marker, enrichment,
  trajectory and supervised drug-response interfaces for both new reductions.
- Test R-side architecture audits/forward passes and the S4/downstream path;
  official-checkpoint parity remains pending local access to the official weights.

# fmverse 0.0.0.9002

- Add native adapters for Geneformer-V1-10M, Geneformer-V2-104M,
  Geneformer-V2-316M and Geneformer-V2-104M_CLcancer. Select version-matched
  gc30M/gc104M dictionaries, V1/V2 sequence rules, and model config while
  preserving the earlier gf-12L-95M-i4096 workflow.
- Export `fmverse_geneformer_checkpoints()` and document each checkpoint's
  training corpus, token capacity and dictionary release.
- Inference and downstream APIs return the formal S4 `fmobj` class with one shared
  Seurat v5 object, independent model runs, an active run, history and schema version.
- Add multiple checkpoint/model runs with `fmverse(obj, ..., run_id=...)`; use
  `fm_use()`, `fm_runs()`, `fm_run()`, `fm_seurat()`, `fm_embeddings()`,
  `fm_analysis()` and `fm_plot()` for explicit selection and extraction.
- Namespace run-owned metadata, graphs and reductions. Integration retains the
  original embedding. Repeated analyses replace their own outputs and invalidate
  dependent results. Explicit run replacement requires `overwrite=TRUE`.
- Migrate legacy `fmverse_result` objects with `as_fmobj()` and retain the existing
  Seurat v4-to-v5 path. Saved objects read without checkpoint files; only model-dependent
  computation needs the checkpoint. Live model modules are not serialized in runs.
- Preserve scGPT/Geneformer numerical kernels, tokenization, batch inference and
  plotting algorithms. The 64-cell pre/post-refactor check has identical input tensors
  and zero embedding error for both official checkpoints.
