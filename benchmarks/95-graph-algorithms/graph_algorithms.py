#!/usr/bin/env python3
"""Experiment 95 -- graph algorithms (networkx column).

Same cases, same graphs (the identical LCG construction as graph_algorithms.m),
same order.  networkx is pure Python.  Idioms: maximum_flow_value (preflow
push), stoer_wagner, edge_connectivity, max_weight_matching(maxcardinality)
for general graphs, min_edge_cover, max_weight_clique(weight=None) for the
maximum clique and -- on the complement -- the maximum independent set /
minimum vertex cover, vf2pp_is_isomorphic / vf2pp_isomorphism, and
check_planarity.  networkx has no Hamiltonian-cycle search: SKIP.

Cold and warm are the same for networkx (it caches nothing), so the cold
cases rebuild the graph object, untimed, exactly like the CAS side.
"""

import sys, os; sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import time
import networkx as nx
from harness import bench, bench_once, check, require

require(["networkx", "networkx:maximum_flow_value", "networkx:stoer_wagner",
         "networkx:max_weight_matching", "networkx:max_weight_clique",
         "networkx:vf2pp_is_isomorphic", "networkx:check_planarity"])


def lcg(length, s):
    out, x = [], s
    for _ in range(length):
        x = (1103515245 * x + 12345) % 2147483648
        out.append(x)
    return out


def rnd_pairs(n, m, s):
    r = lcg(2 * m + m + 1000, s)
    a = [(x // 65536) % n + 1 for x in r[0::2]]
    b = [(x // 65536) % n + 1 for x in r[1::2]]
    seen, p = set(), []
    for x, y in zip(a, b):
        e = (min(x, y), max(x, y))
        if e[0] < e[1] and e not in seen:
            seen.add(e); p.append(e)
    return p[:m]


def perm(n, s):
    v = lcg(n, s)
    return [i + 1 for i in sorted(range(n), key=lambda i: v[i])]


def grid_pairs(w, h):
    e = [((i - 1) * w + j, (i - 1) * w + j + 1) for i in range(1, h + 1) for j in range(1, w)]
    e += [((i - 1) * w + j, i * w + j) for i in range(1, h) for j in range(1, w + 1)]
    return e


def tri_pairs(w, h):
    return grid_pairs(w, h) + [((i - 1) * w + j, i * w + j + 1)
                               for i in range(1, h) for j in range(1, w)]


def cubic_pairs(n, s):
    p, q = perm(n, s), perm(n, s + 1)
    c = list(zip(p[:-1], p[1:])) + [(p[-1], p[0])]
    mt = [(q[2 * k], q[2 * k + 1]) for k in range(n // 2)]
    seen, out = set(), []
    for x, y in c + mt:
        e = (min(x, y), max(x, y))
        if e not in seen:
            seen.add(e); out.append(e)
    return out


def mk(n, pairs, caps=None):
    G = nx.Graph()
    G.add_nodes_from(range(1, n + 1))
    if caps is None:
        G.add_edges_from(pairs)
    else:
        G.add_edges_from((u, v, {"capacity": c}) for (u, v), c in zip(pairs, caps))
    return G


nR = 20000; rP = rnd_pairs(nR, 100000, 7)
nC = 2000;  cP = rnd_pairs(nC, 20000, 11)
dP = rnd_pairs(200, 10000, 13)
sP = rnd_pairs(100, 300, 17)
gw, gh = 250, 200; nG = gw * gh; gP = grid_pairs(gw, gh)
fw, fh = 500, 400; nF = fw * fh; fP = grid_pairs(fw, fh)
capF = [(7919 * k) % 100 + 1 for k in range(1, len(fP) + 1)]
nM = 200000; mP = rnd_pairs(nM, 1000000, 7)
capM = [(104729 * k) % 1000 + 1 for k in range(1, len(mP) + 1)]
tw = 316; nT = tw * tw; tP = tri_pairs(tw, tw)
nI = 10000; iP = cubic_pairs(nI, 28); rl = perm(nI, 29)

rnd1 = mk(nR, rP); rnd2 = mk(nC, cP); dense = mk(200, dP); sparse = mk(100, sP)
grid = mk(nG, gP); gridF = mk(nF, fP, capF); rndM = mk(nM, mP, capM); tri = mk(nT, tP)
cub = mk(nI, iP); cubP = mk(nI, [(rl[u - 1], rl[v - 1]) for u, v in iP])


def mis_size(G):
    return len(nx.max_weight_clique(nx.complement(G), weight=None)[0])


class _Timeout(Exception):
    pass


def capped(label, fn, secs, value_fn=len):
    """One timed run, abandoned after `secs` seconds (SIGALRM): some networkx
    routines (edge_connectivity, the pure-Python blossom matching) run for
    hours at these sizes.  A capped case prints SKIP with the cap, so the
    table records only that networkx took longer than the cap."""
    import signal
    def on_alarm(sig, frame):
        raise _Timeout()
    old = signal.signal(signal.SIGALRM, on_alarm)
    signal.alarm(secs)
    try:
        t0 = time.perf_counter()
        r = fn()
        dt = time.perf_counter() - t0
    except _Timeout:
        print("SKIP\t%s\tnetworkx:timeout>%ds" % (label, secs), flush=True)
        return
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)
    print("BENCH\t%s\t%.3f" % (label, 1000 * dt), flush=True)
    check(label, value_fn(r))


def skip(label):
    print("SKIP\t%s\t%s" % (label, "networkx:hamiltonian_cycle"), flush=True)


L = "FindMaximumFlow, 2*10^5-vertex grid, capacities"
capped(L, lambda: nx.maximum_flow_value(gridF, 1, nF), 300, lambda r: r)
L = "FindMaximumFlow, 10^6-edge random graph, capacities"
capped(L, lambda: nx.maximum_flow_value(rndM, 1, 2), 300, lambda r: r)
L = "FindMinimumCut, 2000 vertices 20000 edges"
bench_once(L, lambda: nx.stoer_wagner(rnd2)); check(L, nx.stoer_wagner(rnd2)[0])
L = "EdgeConnectivity, 50000-vertex grid"
capped(L, lambda: nx.edge_connectivity(grid), 120, lambda r: r)
L = "FindIndependentEdgeSet, 10^5-edge random graph"
capped(L, lambda: nx.max_weight_matching(rnd1, maxcardinality=True), 300)
L = "FindIndependentEdgeSet, 10^5-vertex triangulated grid"
capped(L, lambda: nx.max_weight_matching(tri, maxcardinality=True), 300)
L = "FindEdgeCover, 50000-vertex grid"
capped(L, lambda: nx.min_edge_cover(grid), 300)
L = "FindClique, dense 200-vertex random graph"
bench(L, lambda: nx.max_weight_clique(dense, weight=None)); check(L, len(nx.max_weight_clique(dense, weight=None)[0]))
L = "FindIndependentVertexSet, 100 vertices 300 edges"
bench_once(L, lambda: mis_size(sparse)); ms = mis_size(sparse); check(L, ms)
L = "FindVertexCover, 100 vertices 300 edges"
bench_once(L, lambda: mis_size(sparse)); check(L, 100 - ms)
L = "HamiltonianGraphQ, 400-vertex random cubic graph"; skip(L)
L = "IsomorphicGraphQ, 10^4-vertex 3-regular random graph, permuted"
capped(L, lambda: nx.vf2pp_is_isomorphic(cub, cubP), 300, lambda r: r)
L = "FindGraphIsomorphism, 10^4-vertex 3-regular random graph"
capped(L, lambda: nx.vf2pp_isomorphism(cub, cubP), 300, lambda r: 1 if r else 0)
L = "PlanarGraphQ, 10^5-vertex triangulated grid"
bench(L, lambda: nx.check_planarity(tri)[0]); check(L, nx.check_planarity(tri)[0])

# cold: identical for networkx, rebuilt untimed
def cold(label, build, f, value):
    ts = []
    for _ in range(3):
        G = build()
        t0 = time.perf_counter(); f(G); ts.append(time.perf_counter() - t0)
    print("BENCH\t%s\t%.3f" % (label, 1000 * min(ts)), flush=True)
    check(label, value)

G = mk(nF, fP, capF)
capped("FindMaximumFlow, 2*10^5-vertex grid, capacities (cold)", lambda: nx.maximum_flow_value(G, 1, nF), 300, lambda r: r)
G = mk(nM, mP, capM)
capped("FindMaximumFlow, 10^6-edge random graph, capacities (cold)", lambda: nx.maximum_flow_value(G, 1, 2), 300, lambda r: r)
L = "FindMinimumCut, 2000 vertices 20000 edges (cold)"
bench_once(L, lambda: nx.stoer_wagner(mk(nC, cP))); check(L, nx.stoer_wagner(rnd2)[0])
L = "FindIndependentEdgeSet, 10^5-edge random graph (cold)"
G = mk(nR, rP); capped(L, lambda: nx.max_weight_matching(G, maxcardinality=True), 300)
L = "FindIndependentEdgeSet, 10^5-vertex triangulated grid (cold)"
G = mk(nT, tP); capped(L, lambda: nx.max_weight_matching(G, maxcardinality=True), 300)
cold("FindClique, dense 200-vertex random graph (cold)", lambda: mk(200, dP),
     lambda G: nx.max_weight_clique(G, weight=None), len(nx.max_weight_clique(dense, weight=None)[0]))
L = "FindVertexCover, 100 vertices 300 edges (cold)"
bench_once(L, lambda: mis_size(mk(100, sP))); check(L, 100 - ms)
skip("HamiltonianGraphQ, 400-vertex random cubic graph (cold)")
G2 = (mk(nI, iP), mk(nI, [(rl[u - 1], rl[v - 1]) for u, v in iP]))
capped("IsomorphicGraphQ, 10^4-vertex 3-regular random graph, permuted (cold)",
       lambda: nx.vf2pp_is_isomorphic(G2[0], G2[1]), 300, lambda r: r)
cold("PlanarGraphQ, 10^5-vertex triangulated grid (cold)", lambda: mk(nT, tP),
     lambda G: nx.check_planarity(G)[0], nx.check_planarity(tri)[0])
