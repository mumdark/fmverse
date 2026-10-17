# fmverse

`fmverse` is an R-native runtime for selected single-cell foundation models. The package supports scGPT, Geneformer's four currently listed core pretrained variants (plus compatibility for `gf-12L-95M-i4096`), scBERT `panglao-human`, scFoundation `01B-resolution`, and both official UCE depths (`4layer` and `33layer`) across eight species. Public inference and analysis functions return an S4 `fmobj` with one shared Seurat object and independent model runs. Geneformer inference uses the same R-native BERT and downstream interface across variants:

1. Read official checkpoint tensors in R and compare every key and shape with the R model.
2. Feed Python and R identical token IDs and masks; compare hidden states at every transformer layer in CPU float32 evaluation mode.
3. Reproduce each version's raw-count tokenization in R and compare embeddings with its Python implementation using the same inputs.

The current Geneformer registry includes `Geneformer-V1-10M`, `Geneformer-V2-104M`, `Geneformer-V2-316M`, and `Geneformer-V2-104M_CLcancer`. V1 selects Genecorpus-30M dictionaries and its 2,048-token, no-boundary-token input format. V2 variants select Genecorpus-104M dictionaries, add CLS/EOS boundary tokens, and use up to 4,096 tokens. Each checkpoint's config builds the correct hidden size, layer count, and vocabulary. The older `gf-12L-95M-i4096` path retains its historical tokenizer behavior. Checkpoint metadata and training-set descriptions are available with `fmverse_geneformer_checkpoints()`.

Select a downloaded checkpoint folder with `model_dir`; its `fmverse.json` declares the model, checkpoint, and required relative file paths. The public inference call does not take `model` or `checkpoint`. For example, select `00scGPT/checkpoint/scGPT-whole-human` or `00Geneformer/checkpoint/Geneformer-V2-104M`. The `fmverse_*_checkpoints()` functions still list supported versions.

`fmverse()` returns an S4 `fmobj`. Its five slots are `seurat`, `runs`, `active_run`, `history`, and `version`. Expression data is shared once; each model run references its own Seurat reduction and stores checkpoint metadata, the unchanged token inputs, weight audit, and downstream analyses. Live torch modules are loaded from their original checkpoints when needed and are not stored in the object. Use the accessors below rather than the old list fields.

The scGPT R model preserves all 163 checkpoint tensors, including the legacy `flag_encoder.weight` table (loaded and audited, though not used by the pinned public inference path). It maps scGPT's fused QKV weights to torch for R's attention projection and converts zero-based scGPT token IDs to the one-based indices required by R torch embeddings. The cell-embedding path stores CLS embeddings and L2-normalizes them as scGPT does.

The Geneformer adapter reads each original `model.safetensors` checkpoint and strictly audits learned tensor names and shapes. It applies the matching Genecorpus-30M, Genecorpus-95M, or Genecorpus-104M token/rank dictionaries, ranks detected genes from raw counts using the official normalization, and returns mean-pooled hidden states from the second-to-last transformer block as that run's reduction. The V1 positional-index buffer is validated separately as a non-learned tensor. The PBMC3K preparation script retains original 10x Ensembl IDs in feature metadata so tokenization does not depend on symbol-to-ID ambiguity. All variants share the same `fmobj` and downstream functions as scGPT.

Two additional R-native encoders are available through the same `fmobj` run interface. `scBERT` supports the human `panglao-human` checkpoint (16,906 ordered Panglao genes, 200-dimensional Performer output); its checkpoint directory now includes the ordered Panglao gene index. `scFoundation` supports the human `01B-resolution` cell-embedding weights (19,264 fixed genes) and exposes its official `pool_type = "all"` or `"max"` output. Both models use the common cluster, annotation, integration, marker, enrichment, trajectory, and supervised drug-response interfaces. Reconstruction, perturbation, and attention-derived GRN methods remain limited to the model adapters that expose those corresponding outputs.

The scBERT `.pth` and scFoundation nested MMF checkpoint are not flat tensor dictionaries. The helper scripts in `../00scBERT/` and `../00scFoundation/` make a one-time format conversion to safetensors and, for scFoundation, export its model config. These conversions do not run inference; all tokenization, weight loading, and inference then run in R torch. Each selectable folder contains `fmverse.json` and the runtime files needed for that checkpoint. See those folders' READMEs for source links and model-specific licenses.

```r
scbert <- fmverse(
  pbmc3k, model_dir = "00scBERT/checkpoint", batch_size = 1
)

scfoundation <- fmverse(
  pbmc3k, model_dir = "00scFoundation/checkpoint", batch_size = 1,
  tgthighres = "a5", pool_type = "all"
)
fm_runs(scfoundation)
scfoundation <- fm_cluster(scfoundation, method = "leiden")
```

The model-specific checkpoint tables are available through `fmverse_scbert_checkpoints()` and `fmverse_scfoundation_checkpoints()`. Each folder is validated before inference; missing labels or assets produce a direct error naming the required file. The model-specific embedding paths have been covered by R-side architecture, strict state-dictionary self-audits, preprocessing, and forward-pass tests.

All eight official model-zoo folders have now been downloaded and audited. Seven checkpoints contain 163 tensors and load through the shared encoder path. The `continual-pretrained` checkpoint adds a 10-tensor cell-type classification decoder; the loader detects its output dimension and layer count, constructs the matching R module, and loads those tensors too. The `blood` model has its own smaller 36,574-token vocabulary, while the other seven use 60,697 tokens. Each run reads that checkpoint's own vocabulary and preprocessing configuration instead of reusing whole-human metadata.

The pinned PyTorch file stores its tensor locations as `cuda:0`. CPU-only torch for R cannot remap those tags during `load_state_dict()`, so `fmverse` makes a CPU-tagged copy in the R user cache by rewriting only the device strings in the checkpoint's pickle metadata; all tensor payload bytes remain unchanged. Inference itself runs in R and does not use Python.

By default, raw-count preprocessing follows the pinned scGPT path: it computes quantile bins for all expressed genes, randomly resolves repeated quantile boundaries, then samples and permutes genes if a sequence exceeds 1,200 tokens. R uses a cell-specific seed derived from the base seed (default `20261013`) and barcode; its prior RNG state is restored afterward, so results do not depend on chunk boundaries. Since R, NumPy, and PyTorch use different RNG streams, the sampled genes and tie outcomes are reproducible in R but are not guaranteed to be the exact draws from a particular Python run. Set `preprocessing = "deterministic"` to use midpoint handling for tied bins and a stable CRC32 subset; this keeps the original fixed-input parity fixture reproducible.

## Environment

Use the existing Windows R 4.4.0 installation and restore the project library with `renv::restore()`. `renv.lock` is generated by `bootstrap.R`. Initial numerical parity runs on CPU float32; the package API also accepts a torch device for later performance checks. The eight official checkpoint folders are cached under `00scGPT/checkpoint` by the sibling PowerShell helper; already-complete files are skipped, and SHA256 values are recorded in `official_scgpt_manifest.json`. The model weights stay outside the R package source.

```r
# Run once from the fmverse project directory to restore dependencies and install it.
source("bootstrap.R")
library(fmverse)
```

For the optional downstream methods, install their backends into the project library:

```r
renv::install(c("harmony", "leidenbase", "slingshot", "clusterProfiler", "org.Hs.eg.db"))
```

After loading the package, pass the downloaded model-version folder:

```r
result <- fmverse(
  pbmc3k,
  model_dir = "00scGPT/checkpoint/scGPT-whole-human",
  batch_size = 16,
  preprocessing = "official",
  seed = 20261013
)
fm_embeddings(result)
fm_seurat(result)[[fm_run(result)$reduction]]

# Select another official model-zoo checkpoint with the same API.
result_blood <- fmverse(
  pbmc3k,
  model_dir = "00scGPT/checkpoint/scGPT-blood",
  batch_size = 8
)

# Load the official Geneformer checkpoint through the same result API.
geneformer_result <- fmverse(
  pbmc3k,
  model_dir = "00Geneformer/checkpoint/Geneformer-V2-104M",
  batch_size = 2
)
fm_embeddings(geneformer_result)

# The same run interface selects the cancer-domain continued-pretraining model.
geneformer_cancer <- fmverse(
  pbmc3k,
  model_dir = "00Geneformer/checkpoint/Geneformer-V2-104M_CLcancer",
  batch_size = 2
)
```

Every supported model folder must contain a `fmverse.json` manifest. All file
paths inside it are relative to that folder and are validated before loading;
model weights and preprocessing dictionaries are not resolved from adjacent
folders. When `max_length` is omitted, the selected checkpoint's config supplies
the sequence length.

For example, the whole-human scGPT folder is marked like this:

```json
{
  "schema_version": 1,
  "model": "scGPT",
  "checkpoint": "whole-human",
  "files": {
    "weights": "best_model.pt",
    "config": "args.json",
    "vocab": "vocab.json"
  }
}
```

Geneformer manifests additionally locate four tokenizer/rank dictionaries;
scBERT and scFoundation locate their ordered gene index (and scFoundation its
config); UCE locates the weights and `assets/` directory and defaults
`defaults.species` to `human`. The runtime reads these fields as paths, so the
folder may be moved or renamed without relying on repository-relative assets.

## Seurat compatibility

`fmverse_to_seurat_v5()` updates legacy objects and casts legacy assays to Seurat v5 `Assay5` objects while retaining the other object components. Model input is read from the raw `counts` layer. For scGPT, tokenization uses the assay feature metadata column `scGPT_gene_name` when present. For Geneformer, the adapter first looks for an Ensembl-ID column (`geneformer_ensembl_id`, `ensembl_id`, or `gene_id`), then falls back to the official symbol-to-Ensembl mapping. The PBMC3K preparation script stores both original 10x symbols and Ensembl IDs.

Padding is trimmed per batch to the longest active sequence in that batch. The fixed-input layer-parity fixture uses 8 cells; full 2,700-cell PBMC3K runs have also been completed for all eight official checkpoints in both R and the pinned Python implementation. The model-zoo embeddings and per-cell comparisons are saved in the sibling `00scGPT/results/model-zoo` directory. The R and Python preprocessing implementations follow the same method but use independent RNG streams, so stochastic tie handling and sequence truncation can select different genes in individual cells.

## Multiple runs and downstream analysis

An `fmobj` owns one Seurat v5 object. Repeated `fmverse()` calls append model runs without duplicating that object or replacing another run's embeddings. `fmobj(pbmc)` creates an empty container. The new run becomes active; an explicit `run_id` on a downstream call does not change the active run.

```r
obj <- fmverse(pbmc, model_dir = "00scGPT/checkpoint/scGPT-whole-human",
               run_id = "scgpt_whole", batch_size = 2)
obj <- fmverse(obj, model_dir = "00Geneformer/checkpoint/Geneformer-V2-104M",
               run_id = "geneformer_95m", batch_size = 2)
fm_runs(obj)

obj <- fm_cluster(obj, run_id = "scgpt_whole", method = "leiden")
obj <- fm_use(obj, run_id = "geneformer_95m")
obj <- obj |>
  fm_cluster(method = "leiden") |>
  fm_markers() |>
  fm_enrich()

fm_embeddings(obj, run_id = "scgpt_whole")
fm_run(obj, run_id = "geneformer_95m")
fm_analysis(obj, run_id = "geneformer_95m", task = "markers")$table
p <- fm_plot(obj, run_id = "geneformer_95m", task = "cluster", plot = "umap")
p + ggplot2::theme_bw()
ggplot2::ggsave("geneformer_clusters.png", p)
```

`run_id` selects a model/checkpoint run. If omitted, the active run is used. Omitted inference IDs are generated from the model and checkpoint with a numeric suffix on collision. An explicitly repeated ID raises an error unless `overwrite = TRUE`; replacement removes only that run's owned results.

Cluster labels, annotations, neighbor graphs, UMAPs, and integration reductions have unique internal names per run. Obtain the actual cluster column from `fm_analysis(obj, run_id, task = "cluster")$cluster_column` and the original reduction name from `fm_run(obj, run_id)$reduction`. `fm_markers()` and `fm_trajectory()` automatically resolve the selected run's cluster column when `group_by` is omitted. Explicit existing metadata columns remain supported.

```r
obj <- fm_annotate(obj, reference = reference_obj, label = "cell_type",
  run_id = "geneformer_95m", reference_run_id = "reference_geneformer")
obj <- fm_integrate(obj, run_id = "geneformer_95m", batch = "batch", method = "harmony")
obj <- fm_cluster(obj, run_id = "geneformer_95m", embedding = "integrated")
obj <- fm_trajectory(obj, run_id = "geneformer_95m", root = "0")
obj <- fm_drug_response(obj, run_id = "geneformer_95m",
  treatment = "treatment", drug = "drug_x", control_label = "control")
```

The original embedding remains available after integration. Embedding-based downstream calls default to `embedding = "original"`; use `"integrated"` explicitly to analyze a current integration result. The compatibility `reduction` argument accepts only a selected run's own original or integration reduction (including its model-name alias). A reduction from another run raises an error.

Each run keeps one current result per task. Repeating a task replaces that result and appends history. Dependent analyses become stale; `fm_analysis()` and `fm_plot()` require recomputation unless `allow_stale = TRUE` is explicitly supplied for inspection. Analysis records contain tables/metrics, parameters, revision, dependencies, and editable ggplots. Package functions never save figures.

The same object supports `fm_reconstruct()`, `fm_perturb()`, and `fm_grn()`, selected by `run_id`. For scGPT, reconstruction and perturbation use the expression decoder. Geneformer reconstruction masks gene tokens and scores the MLM head; perturbation masks a target token and reports contextual hidden-state shifts. Geneformer diagnostics do not represent raw-count reconstruction or causal expression effects. Attention edges are candidate associations. Drug response uses supervised held-out-cell classification of observed conditions, not zero-shot prediction for unseen drugs.

## Save and migrate

```r
saveRDS(obj, "pbmc_fmobj.rds")
obj <- readRDS("pbmc_fmobj.rds")
obj <- as_fmobj(old_result, run_id = "scgpt_whole")
```

Old `fmverse_result` objects are migrated without repeating inference. Downstream functions accept them through this conversion. Saved embedding/table/plot access does not require checkpoint files; model-dependent tasks require their original checkpoint and configuration. Seurat v4 inputs continue through the existing v5 conversion path.

Batch inference, tokenization and numerical algorithms are unchanged. This release supports R pipelines and multiple runs; it does not add HDF5 storage or interrupt/resume execution.

For reproducible examples see `../00fmobj/run_pbmc64.R`. It runs the same 64 PBMC3K cells through both official checkpoints, compares pre/post-refactor tokens and embeddings, saves an S4 multi-run object, and explicitly calls `ggsave()` for each run's plots. Results live under `../00fmobj/results/`. Existing full scGPT and Geneformer runners use the accessor API; their historical result files remain in place.
# Universal Cell Embedding (UCE)

`fmverse` supports both official UCE checkpoint depths through the same S4
`fmobj` run interface:

```r
pbmc <- fmverse_to_seurat_v5(pbmc)
obj <- fmverse(
  pbmc, model_dir = "00UCE/checkpoint/uce-4layer",
  species = "human", batch_size = 1, seed = 20261013
)
fm_embeddings(obj)
```

Use `00UCE/checkpoint/uce-33layer` for the paper model. Available species can be
listed with `fmverse_uce_species()`; the official release currently includes
human, mouse, frog, zebrafish, mouse lemur, pig, cynomolgus macaque, and rhesus
macaque. The default checkpoint is the 4-layer example model. The 4-layer and
33-layer embeddings come from different model depths and should be treated as
separate runs.

UCE requires raw counts and gene symbols. Split Seurat v5 `counts.*` layers are
joined for tokenization. It samples 1,024 genes per cell with
expression-weighted replacement, then uses the selected species' bundled ESM2
gene-token order and chromosome/start-position annotations. The release assets
are organized under `00UCE/checkpoint/`; the main large assets are the model
weights (3.40 GB and 5.69 GB), token matrix (2.98 GB), and multi-species ESM2
table archive (2.74 GB). `00UCE/prepare_assets.py` verifies the official
Figshare checksums and creates the species lookup tables used by R.

The R torch runtime available during development is CPU-only. Keep `batch_size`
small for local inference; this controls the number of cells per forward pass,
not the 1,536-token maximum sequence length.

The selected UCE directory includes its `fmverse.json` and `assets/` directory
generated by `00UCE/prepare_assets.py`; it can be moved independently of the
workspace and passed directly as `model_dir`.
