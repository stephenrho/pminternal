
# pminternal 0.1.1

- loaded packages are now provided to parallel clusters via `parallel::clusterEvalQ`, meaning that users no longer have to use `::` when desiring `cores > 1` for bootstrap
- optional argument `export` added so that other objects from environment can be supplied to parallel clusters via `parallel::clusterExport`

# pminternal 0.1.0

* Initial CRAN submission.
