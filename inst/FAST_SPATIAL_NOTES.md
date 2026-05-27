# Fast spatial HAC branch

This branch replaces the dense-distance-matrix spatial meat with a fast CSR path.

Implemented changes:

1. Score accumulation: rows are converted conceptually to `S_i = e_i X_i`, and the meat is computed as `M + t(M)` where `M = sum_i S_i' C_i` and `C_i = 0.5 S_i + sum_{j > i} w_ij S_j`.
2. Candidate pruning: observations are sorted by latitude within each time block; the pair loop breaks once latitude distance exceeds the cutoff and uses a conservative longitude screen before evaluating the exact distance.
3. CSR neighbor lists: neighbor pairs within the cutoff are stored as `row_ptr`, `col_idx`, and `weight`, avoiding dense `n x n` distance matrices. With `balanced_pnl = TRUE`, the CSR graph is computed once from the first period and reused for each time block.

The serial HAC component is intentionally kept close to the original implementation, although its internal C++ loop also avoids dense temporary matrices.

Manual tests live in `tests/manual/test-fast-spatial.R`.
