#' Convert a legacy Seurat object to Seurat v5 assay representation
#'
#' @param object A Seurat v4 or v5 object.
#' @return A Seurat object whose legacy assays have been cast to `Assay5`.
#' @export
fmverse_to_seurat_v5 <- function(object) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must inherit from `Seurat`.", call. = FALSE)
  }

  if (utils::packageVersion("SeuratObject") < "5.0.0") {
    stop("SeuratObject >= 5.0.0 is required to create Assay5 assays.", call. = FALSE)
  }

  assay_names <- SeuratObject::Assays(object)
  already_v5 <- length(assay_names) > 0L && all(vapply(
    assay_names, function(assay_name) inherits(object[[assay_name]], "Assay5"),
    logical(1)
  ))
  if (already_v5) return(object)

  object <- SeuratObject::UpdateSeuratObject(object)
  assay_names <- SeuratObject::Assays(object)
  for (assay_name in assay_names) {
    assay <- object[[assay_name]]
    if (!inherits(assay, "Assay5")) {
      object[[assay_name]] <- methods::as(assay, Class = "Assay5")
    }
  }
  object
}
