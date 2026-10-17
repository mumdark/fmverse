\name{fmverse-model-adapters}
\alias{fmverse_scbert_checkpoints}
\alias{fmverse_scfoundation_checkpoints}
\alias{fmverse_load_scbert_weights}
\alias{fmverse_load_scfoundation_weights}
\alias{fmverse_scbert_build_model}
\alias{fmverse_scfoundation_build_model}
\title{scBERT and scFoundation adapters}
\description{
List supported checkpoints, construct the R-native scFoundation cell model,
and strictly audit flat checkpoint tensors.
}
\usage{
fmverse_scbert_checkpoints()
fmverse_scfoundation_checkpoints()
fmverse_load_scbert_weights(model, path)
fmverse_load_scfoundation_weights(model, path)
fmverse_scbert_build_model()
fmverse_scfoundation_build_model(config)
}
\arguments{
  \item{model}{An R torch model module matching the selected checkpoint.}
  \item{path}{A flat named state dictionary or safetensors file. scBERT and
  scFoundation helpers in the workspace convert their official wrapped files
  into this format.}
  \item{config}{The JSON-decoded scFoundation cell-model configuration.}
}
\value{
Checkpoint functions return data frames. Model builders return torch modules.
Weight loaders return the loaded module and a complete tensor key/shape audit.
}
\details{
scBERT supports the human Panglao-pretrained checkpoint and reads the ordered
Panglao gene index from its model folder. scFoundation supports the human
\code{01B-resolution} cell checkpoint and reads its 19,264-gene index from its
model folder. The model-specific helper scripts perform a one-time
checkpoint-format conversion; inference executes in R torch. The common
embedding downstream functions work with either run.
}
\seealso{\code{fmverse()}, \code{fm_runs()}, \code{fm_cluster()}}
