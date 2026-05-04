import math
import random
import numpy as np

EARTH = 6371.0
DEG = math.pi / 180.0

def sqr(x):
    return x * x

def lon_abs_diff_deg(a, b):
    diff = abs(a - b) % 360.0
    if diff > 180.0:
        diff = 360.0 - diff
    return diff

def haversine(lat1, lon1, lat2, lon2):
    lat1 *= DEG; lon1 *= DEG; lat2 *= DEG; lon2 *= DEG
    a = sqr(math.sin((lat2 - lat1) / 2.0)) + math.cos(lat1) * math.cos(lat2) * sqr(math.sin((lon2 - lon1) / 2.0))
    return 2 * EARTH * math.asin(min(1.0, math.sqrt(a)))

def spherical(lat1, lon1, lat2, lon2):
    lat1 *= DEG; lon1 *= DEG; lat2 *= DEG; lon2 *= DEG
    inner = math.sin(lat1) * math.sin(lat2) + math.cos(lat1) * math.cos(lat2) * math.cos(lon1 - lon2)
    inner = min(1.0, max(-1.0, inner))
    return EARTH * math.acos(inner)

def flatearth(lat1, lon1, lat2, lon2):
    return math.sqrt((111 * (lat1 - lat2)) ** 2 + (math.cos(lat1 * DEG) * 111 * (lon1 - lon2)) ** 2)

def dist(fn, a, b):
    if fn == 'haversine':
        return haversine(a[0], a[1], b[0], b[1])
    if fn == 'spherical':
        return spherical(a[0], a[1], b[0], b[1])
    if fn == 'flatearth':
        return flatearth(a[0], a[1], b[0], b[1])
    raise ValueError(fn)

def fast_edges(coords, cutoff, kernel, fn):
    n = coords.shape[0]
    order = sorted(range(n), key=lambda i: coords[i, 0])
    lat_rad = coords[:, 0] * DEG
    cos_lat = np.cos(lat_rad)
    edges = []
    lat_cutoff = cutoff / 110.0 + 1e-12
    for pos_i, obs_i in enumerate(order):
        lat_i = coords[obs_i, 0]
        lon_i = coords[obs_i, 1]
        for pos_j in range(pos_i + 1, n):
            obs_j = order[pos_j]
            if coords[obs_j, 0] - lat_i > lat_cutoff:
                break
            mean_lat = 0.5 * (lat_rad[obs_i] + lat_rad[obs_j])
            min_abs_cos = max(1e-8, min(abs(cos_lat[obs_i]), abs(cos_lat[obs_j]), abs(math.cos(mean_lat))))
            lon_cutoff = cutoff / (110.0 * min_abs_cos) + 1e-12
            if lon_abs_diff_deg(lon_i, coords[obs_j, 1]) > lon_cutoff:
                continue
            a, b = obs_i, obs_j
            if fn == 'flatearth' and a > b:
                a, b = b, a
            d = dist(fn, coords[a], coords[b])
            if d <= cutoff:
                w = 1 - d / cutoff if kernel == 'bartlett' else 1.0
                if w != 0:
                    edges.append((obs_i, obs_j, w))
    return edges

def fast_meat(coords, X, e, cutoff, kernel, fn):
    S = X * e[:, None]
    C = 0.5 * S.copy()
    for i, j, w in fast_edges(coords, cutoff, kernel, fn):
        C[i, :] += w * S[j, :]
    M = S.T @ C
    return M + M.T

def dense_meat(coords, X, e, cutoff, kernel, fn):
    n = len(e)
    S = X * e[:, None]
    W = np.eye(n)
    for i in range(n):
        for j in range(i + 1, n):
            d = dist(fn, coords[i], coords[j])
            if d <= cutoff:
                w = 1 - d / cutoff if kernel == 'bartlett' else 1.0
                W[i, j] = w
                W[j, i] = w
    return S.T @ W @ S

random.seed(123)
np.random.seed(123)
for fn in ['haversine', 'spherical', 'flatearth']:
    for kernel in ['bartlett', 'uniform']:
        for n in [8, 25, 60]:
            k = 5
            coords = np.column_stack([
                np.random.uniform(-60, 60, n),
                np.random.uniform(-170, 170, n),
            ])
            X = np.random.normal(size=(n, k))
            X[:, 0] = 1.0
            e = np.random.normal(size=n)
            cutoff = 1200.0
            got = fast_meat(coords, X, e, cutoff, kernel, fn)
            exp = dense_meat(coords, X, e, cutoff, kernel, fn)
            err = np.max(np.abs(got - exp))
            if err > 1e-8:
                raise SystemExit(f'mismatch fn={fn} kernel={kernel} n={n} err={err}')
            # Also ensure screening did not drop true neighbors.
            edge_pairs = {tuple(sorted((i, j))) for i, j, _ in fast_edges(coords, cutoff, kernel, fn)}
            dense_pairs = set()
            for i in range(n):
                for j in range(i + 1, n):
                    if dist(fn, coords[i], coords[j]) <= cutoff:
                        dense_pairs.add((i, j))
            if edge_pairs != dense_pairs:
                raise SystemExit(f'edge mismatch fn={fn} kernel={kernel} n={n}')
print('python reference tests passed')
