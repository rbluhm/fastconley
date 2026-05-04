# Fast spatial HAC branch

This branch replaces the dense-distance-matrix spatial meat with a fast CSR path.

Implemented changes:

1. Score accumulation: rows are converted conceptually to `S_i = e_i X_i`, and the meat is computed as `M + t(M)` where `M = sum_i S_i' C_i` and `C_i = 0.5 S_i + sum_{j > i} w_ij S_j`.
2. Candidate pruning: observations are sorted by latitude within each time block; the pair loop breaks once latitude distance exceeds the cutoff and uses a conservative longitude screen before evaluating the exact distance.
3. CSR neighbor lists: neighbor pairs within the cutoff are stored as `row_ptr`, `col_idx`, and `weight`, avoiding dense `n x n` distance matrices. With `balanced_pnl = TRUE`, the CSR graph is computed once from the first period and reused for each time block.

The serial HAC component is intentionally kept close to the original implementation, although its internal C++ loop also avoids dense temporary matrices.

Manual tests live in `tests/manual/test-fast-spatial.R`.

## Known divergence from upstream: `dist_fn = "flatearth"`

Upstream's `flatearth` distance is `sqrt((Δlat·111)² + (cos(lat₁)·Δlon·111)²)` — asymmetric in `(i,j)`. Upstream walks pairs in input-row order, so the cosine factor uses the lat of the lower-input-index point. The fast path sorts by latitude before iterating, so the cosine factor uses the lat of the lower-latitude point.

For pair-distances near the cutoff this changes which neighbors are kept and what weight they get. On the validation panels (lat ∈ [35°, 55°], `dist_cutoff = 1800`), VCOV entries differ from upstream by `O(1e-4)` to `O(3e-3)` in absolute terms. Other distance functions (`haversine`, `spherical`, `chord`) match upstream to machine or near-machine precision because they are symmetric in `(i,j)`.

This is an intentional consequence of the sort-and-screen design and was not fixed in this branch.
