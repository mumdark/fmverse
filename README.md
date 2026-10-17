# fmverse

`fmverse` runs selected pretrained single-cell foundation models in R and provides a shared interface for embeddings and downstream analyses.

## How it works

Each checkpoint is stored in a model directory with an `fmverse.json` manifest. The manifest identifies the model, checkpoint, and required files. `fmverse()` reads that metadata, loads the checkpoint into a native R `torch` implementation, and stores the run in an S4 `fmobj`. The package does not use Python or `reticulate` for inference. Model weights are distributed separately and are not included in the R package.

## Supported models

| Model | Parameters / checkpoint options | Training data | Description |
| --- | --- | --- | --- |
| **scGPT** | `whole-human`, `continual-pretrained`, `brain`, `blood`, `heart`, `lung`, `kidney`, `pan-cancer` | `whole-human`: 33M normal human cells; tissue/cancer checkpoints: brain 13.2M, blood and bone marrow 10.3M, heart 1.8M, lung 2.1M, kidney 814K, and pan-cancer 5.7M cells. `continual-pretrained` starts from whole-human and adds a supervised cell-type stage using Tabula Sapiens. | Gene-expression token Transformer. The whole-human checkpoint is a general-purpose option; tissue-specific checkpoints match their corresponding context. [Project](https://github.com/bowang-lab/scGPT) · [continual-pretraining tutorial](https://github.com/bowang-lab/scGPT/tree/main/tutorials/zero-shot) |
| **Geneformer** | `Geneformer-V1-10M` (10M parameters); `Geneformer-V2-104M` (104M); `Geneformer-V2-316M` (316M); `Geneformer-V2-104M_CLcancer` (104M) | V1: Genecorpus-30M (~30M human transcriptomes); V2: Genecorpus-104M (~104M non-cancer human transcriptomes); cancer variant: continued training on ~14M cancer transcriptomes. | Human rank-value Transformer. V1 uses up to 2,048 tokens; V2 uses up to 4,096. [Official model card](https://huggingface.co/ctheodoris/Geneformer) |
| **scBERT** | `panglao-human` (6 layers; 200-dimensional hidden representation) | PanglaoDB single-cell collection; the supported checkpoint uses a human 16,906-gene index. | Performer-based masked-gene pretraining for cell representation and annotation. [Paper](https://www.nature.com/articles/s42256-022-00534-z) · [code](https://github.com/TencentAILabHealthcare/scBERT) |
| **scFoundation** | `01B-resolution` (100M parameters) | More than 50M human single-cell transcriptomic profiles. | xTrimoGene encoder-decoder; uses a fixed 19,264-gene index and supports cell-level embeddings. [Official project](https://github.com/biomap-research/scFoundation) |
| **UCE** | `4layer` example; `33layer` flagship (>650M parameters) | Integrated Mega-scale Atlas: over 36M cells from more than 300 datasets, spanning eight species. | Cross-species cell embedding using ESM2 protein embeddings and genomic positions. Supported species: human, mouse, frog, zebrafish, mouse lemur, pig, crab-eating macaque, and rhesus macaque. [Paper](https://www.nature.com/articles/s41586-026-10689-z) · [code](https://github.com/snap-stanford/UCE) |

The two legacy Geneformer folder identifiers `gf-12L-95M-i4096` and `geneformer-12L-30M` remain accepted for compatibility; the four named Geneformer releases above are the core catalog.

## Installation

Requirements: R 4.3.0 or later, Seurat 5, and torch for R. Install the torch runtime once and follow any platform-specific prompts.

```r
install.packages(c("remotes", "torch"))
torch::install_torch()

remotes::install_github("mumdark/fmverse")
```

See the [torch for R installation guide](https://torch.mlverse.org/docs/articles/installation) for platform-specific setup.

## Example

This example assumes `pbmc` is a Seurat v4 or v5 object and that the selected model directory has been downloaded with its required files and `fmverse.json` manifest.

```r
library(fmverse)

obj <- fmverse(
  pbmc,
  model_dir = "path/to/scGPT-whole-human",
  run_id = "scgpt_pbmc",
  batch_size = 2,
  device = "cpu"
)

fm_runs(obj)
embedding <- fm_embeddings(obj, run_id = "scgpt_pbmc")

obj <- fm_cluster(obj, run_id = "scgpt_pbmc", method = "louvain")
p <- fm_plot(obj, run_id = "scgpt_pbmc", task = "cluster", plot = "umap")
p + ggplot2::theme_minimal()
```

`model_dir` must point to a directory containing `fmverse.json` and every resource listed in that file. For UCE, `species` can be supplied to select one of the bundled species assets; the manifest default is human.

## Main object

`fmverse()` returns an S4 `fmobj`. Repeated calls can add runs for different models or checkpoints to the same object.

| Component / accessor | Contents |
| --- | --- |
| `obj@seurat` / `fm_seurat(obj)` | One shared Seurat v5 object containing the expression data, metadata, and reductions for each run. Seurat v4 inputs are converted internally. |
| `obj@runs` / `fm_run(obj, run_id)` | Per-run model and checkpoint metadata, cell order, preprocessing and weight audit information, and downstream analysis records. |
| `obj@active_run` / `fm_use(obj, run_id)` | The run selected by default when a downstream function does not receive `run_id`. |
| `obj@history` / `fm_runs(obj)` | Recorded operations and a summary of available runs. |
| `fm_embeddings(obj, run_id)` | Embedding matrix for the selected run. |
| `fm_analysis(obj, run_id, task)` | Stored analysis tables, parameters, and metrics. |
| `fm_plot(obj, run_id, task, plot)` | A ggplot2 object that can be further customized and saved by the user. |

Analyses and plots are stored by run. Plot functions return ggplot2 objects and do not write files.
