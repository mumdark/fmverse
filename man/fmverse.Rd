% Generated documentation for the folder-based inference API.
\name{fmverse}
\alias{fmverse}
\alias{fmverse_geneformer}
\alias{fmverse_scbert}
\alias{fmverse_scfoundation}
\alias{fmverse_uce}
\title{Run a pretrained single-cell model from its model folder}
\description{
Run one of the supported R-native single-cell models using a directory that
contains an \file{fmverse.json} manifest and all files required by that model.
The manifest selects the model and checkpoint; callers do not provide either
argument.
}
\usage{
fmverse(object, model_dir, ..., run_id = NULL, overwrite = FALSE)

fmverse_geneformer(object, model_dir, ..., run_id = NULL, overwrite = FALSE)

fmverse_scbert(object, model_dir, ..., run_id = NULL, overwrite = FALSE)

fmverse_scfoundation(object, model_dir, ..., run_id = NULL,
  overwrite = FALSE)

fmverse_uce(object, model_dir, ..., run_id = NULL, overwrite = FALSE)
}
\arguments{
\item{object}{A Seurat v4/v5 object or an existing \code{fmobj}.}
\item{model_dir}{Path to a self-contained model version folder. It must contain
\file{fmverse.json}.}
\item{...}{Named model settings, including \code{batch_size}, \code{cells},
\code{assay}, \code{device}, and supported model-specific preprocessing
options. Do not pass \code{model}, \code{checkpoint}, \code{checkpoint_dir},
or resource-file overrides such as \code{vocab_file}; model identity and
resource paths are read from the folder.}
\item{run_id}{Optional run name. If omitted, an ID is generated from the model
and checkpoint name and made unique within the \code{fmobj}.}
\item{overwrite}{Whether to replace an existing run with the same explicit
\code{run_id}.}
}
\value{
An S4 \code{fmobj} with the new run as its active run.
}
\details{
The JSON manifest uses schema version 1 and declares the model, checkpoint,
and relative paths to weights and required supporting files. Paths must stay
inside the selected model folder. UCE folders also declare a default species;
the \code{species} argument may override it. Batch-size defaults remain
model-specific (16 for scGPT, 8 for Geneformer, and 1 for scBERT, scFoundation,
and UCE).

The \code{fmverse_geneformer}, \code{fmverse_scbert},
\code{fmverse_scfoundation}, and \code{fmverse_uce} functions are typed
convenience aliases. They accept the same folder-based options and check that
the manifest identifies the named model family.
}
\seealso{\code{fm_runs()}, \code{fm_embeddings()}, \code{fm_cluster()}}
