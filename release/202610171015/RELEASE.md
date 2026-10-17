# fmverse 0.0.0.9005 release snapshot

This folder is the frozen release snapshot created from the development package at `../fmverse` on 2026-10-17. The package R code, documentation, installed data, tests, and package metadata were compared by path and SHA256; all 80 package files matched. The local `renv` library and activation files are development-environment files and are not part of the installable source package.

The release is independent of later changes to `../fmverse`. Update this folder only when explicitly requested.

## Artifacts

- `fmverse/`: package source tree for version `0.0.0.9005`.
- `fmverse_0.0.0.9005.tar.gz`: source package built with R 4.4.0 using `R CMD build`.
- `validation/`: R CMD check and test output from the source archive.

Install the source archive locally with:

```r
install.packages("fmverse_0.0.0.9005.tar.gz", repos = NULL, type = "source")
```

## Validation status

The package built and passed installation/loading checks. `testthat` completed with 0 failures and 272 passing expectations (35 warnings). The complete `R CMD check` finished with **3 warnings and 3 notes**, so this snapshot is not yet clean for a CRAN submission. The warnings concern duplicate or incomplete Rd help entries; notes include missing namespace/global declarations and the S4 `show` export. The current DESCRIPTION also uses the placeholder maintainer address `maintainers@example.org`; replace it with the actual maintainer contact before a public repository submission.

See `validation/R_CMD_check.log` and `validation/testthat.Rout` for the captured output.
