# fmverse

`fmverse` 在 R 中加载并运行预训练单细胞基础模型，为不同模型提供统一的推理和下游分析接口。

## 原理

用户将模型 checkpoint 及其推理资源放在同一目录中。目录内的 `fmverse.json` 声明模型类型、版本和所需文件；`fmverse()` 据此选择对应的 R 实现，通过 torch for R 加载权重并计算细胞 embedding。运行结果保存在 S4 `fmobj` 中。推理在 R 内执行，不依赖 Python 或 reticulate。

模型目录需要包含 `fmverse.json` 和该文件列出的全部资源。模型权重不随 R 包发布。

## 环境与安装

需要 R 4.3.0 或更高版本、Seurat 5 和 torch for R。安装 torch 后按提示重启 R 会话。

```r
install.packages(c("remotes", "torch"))
torch::install_torch()

remotes::install_github("mumdark/fmverse")
```

torch 安装说明见 [torch for R installation guide](https://torch.mlverse.org/docs/articles/installation)。

## 示例

以下示例假设 `pbmc` 是 Seurat v4 或 v5 对象，且模型目录已下载并包含完整资源：

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

obj <- fm_cluster(
  obj,
  run_id = "scgpt_pbmc",
  method = "louvain"
)

p <- fm_plot(
  obj,
  run_id = "scgpt_pbmc",
  task = "cluster",
  plot = "umap"
)
```

`fm_plot()` 返回可继续调整的 ggplot2 对象。

## 主要对象与访问函数

`fmverse()` 返回 S4 `fmobj`。同一个对象可以追加多个模型或 checkpoint 的运行结果。

- `obj@seurat`：唯一的 Seurat v5 对象，保存原始表达数据、细胞元数据和各 run 的 embedding reduction。输入 Seurat v4 对象会在内部转换。
- `obj@runs`：按 `run_id` 保存模型与 checkpoint 信息、细胞顺序、预处理和权重审计信息，以及该 run 的下游分析结果。
- `obj@active_run`：未显式指定 `run_id` 时，下游函数使用的 run。
- `obj@history`：模型推理及下游分析的调用记录。

常用访问函数：

| 函数 | 返回内容 |
| --- | --- |
| `fm_runs(obj)` | 所有 run 的摘要 |
| `fm_run(obj, run_id)` | 指定 run 的记录 |
| `fm_seurat(obj)` | 共享的 Seurat 对象 |
| `fm_embeddings(obj, run_id)` | 指定 run 的 embedding |
| `fm_analysis(obj, run_id, task)` | 指定分析的结果与参数 |
| `fm_plot(obj, run_id, task, plot)` | 可编辑的 ggplot2 图 |
| `fm_use(obj, run_id)` | 设置默认 active run |


