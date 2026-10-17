\name{fmobj}
\alias{fmobj}
\alias{fmobj-class}
\alias{show,fmobj-method}
\alias{as_fmobj}
\alias{fm_runs}
\alias{fm_seurat}
\alias{fm_run}
\alias{fm_use}
\alias{fm_embeddings}
\alias{fm_analysis}
\alias{fm_plot}
\title{S4 container for single-cell foundation-model runs}
\description{Store one shared Seurat object and independently addressable model runs. Public inference and downstream functions return an fmobj.}
\usage{
fmobj(object)
as_fmobj(object, run_id = NULL)
fm_runs(object)
fm_seurat(object)
fm_run(object, run_id = NULL)
fm_use(object, run_id)
fm_embeddings(object, run_id = NULL, embedding = "original")
fm_analysis(object, run_id = NULL, task, allow_stale = FALSE)
fm_plot(object, run_id = NULL, task, plot = NULL, allow_stale = FALSE)
}
\arguments{
\item{object}{An fmobj, Seurat v4/v5 object, or legacy fmverse_result.}
\item{run_id}{Run identifier. NULL selects the active run for extraction or generates a run ID during migration.}
\item{embedding}{The original or integrated embedding of the selected run.}
\item{task}{Analysis name: cluster, annotation, integration, markers, enrichment, trajectory, reconstruct, perturbation, grn, or drug_response.}
\item{plot}{A full plot name or its task-specific suffix, for example umap for cluster_umap. NULL selects the only plot if the task has one.}
\item{allow_stale}{Allow explicit inspection of an analysis made obsolete by a changed upstream result.}
}
\value{Constructors and fm_use return an fmobj. fm_runs returns a data.frame, fm_seurat a Seurat object, fm_run and fm_analysis a list, fm_embeddings a matrix, and fm_plot an editable ggplot.}
\details{An fmobj has five slots: seurat, runs, active_run, history, version. Embeddings are stored once in Seurat reductions and referenced by each run. Changing an active run does not rerun computation. Normal saving uses saveRDS and readRDS; checkpoint-dependent tasks reload their original models as needed. The class validity method checks run identifiers, cell order and artifact references.}
