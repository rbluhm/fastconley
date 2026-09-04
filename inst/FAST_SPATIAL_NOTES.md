# Fast spatial HAC branch

This branch replaces the dense-distance-matrix spatial meat with score
accumulation and output-sensitive neighbor search.

Implemented changes:

1. Score accumulation: rows are converted conceptually to `S_i = e_i X_i`, and the meat is computed as `M + t(M)` where `M = sum_i S_i' C_i` and `C_i = 0.5 S_i + sum_{j > i} w_ij S_j`.
2. Candidate pruning: the default `neighbor = "grid"` buckets 3D unit-sphere coordinates into cubic cells and searches only the current and forward-neighbor cells. The older latitude-band scan is retained as deprecated `neighbor = "band"` compatibility behavior.
3. Balanced CSR and unbalanced streaming: with `balanced_pnl = TRUE`, accepted pairs are stored as `row_ptr`, `col_idx`, and optional `weight`, then reused for every period. The general path streams pairs directly into score accumulators and retains no neighbor list.

The serial HAC component runs as one panel-wide C++ call, grouped by unit and
using deterministic chunked reduction.

Manual tests live in `tests/manual/test-fast-spatial.R`.
