# Zhang accepted-mask RCTD full / scIBD major comparison

Authorized scope: the 20 archived configurations, each on full CD1 and CD2,
without segmentation reruns, additional normalization, clustering or new labels.
Use original scIBD `RNAraw` integer counts and author `major_cluster` labels.
Native `createRctd` preprocessing/default filters and `runRctd(mode='full')`
are the only expression-model steps. The seed is fixed and identical for all jobs.

Accepted masks and bin assignments must match their original checksums and
each other exactly. Aggregate entire original 8um-bin count vectors using the
accepted bin-centre labels; verify exact per-gene assigned + unassigned counts.
No per-molecule localization is inferred. Zero-RNA regions remain in the input
container, with native-filter exclusion explicit. A zero-region mask is a valid
benchmark scope with no possible RCTD fit, not a fabricated zero-weight fit.

One final bundle per configuration/section contains the unchanged full input
SpatialExperiment, native full-mode result, original region IDs and metadata,
per-gene conservation, native numeric diagnostics and source provenance. All
native coefficients are preserved; never clip or renormalize the result.
Full mode has no doublet/singlet confidence classification: none is invented.
Native numeric output is read back and checksummed before acceptance. The final
summary describes coverage, native filter retention, numeric domains and cell
types; none of these is a ground-truth measure of segmentation accuracy.

The new task has its own frozen code, lock, PID, child identity, phase logs and
actual subprocess exit receipts. It never stops or alters the independently
running selected-method/minor-cluster pipeline. Run one full fit at a time with
8 cores and single-threaded BLAS/OpenMP to avoid nested oversubscription.
Results remain local under the benchmark's `rctd_major_full` subdirectory.

Recovery: confirm the failed process has ended, fix only its cause, archive the
failed phase's output/log/receipt, validate in the remote runtime and rerun the
entire failed scope via `--start-phase`. Completed independent scopes are reused
only after their clean exit receipts and final checksums have been verified.
Do not pass a failed partial result off as complete. Final native bundles retain
all input counts; only redundant task-created conversion inputs are removed
after verified readback. Scientific interpretation/visual review follows the
all-scope computational verification; completion JSON does not claim accuracy.
