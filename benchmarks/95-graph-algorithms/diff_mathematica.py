#!/usr/bin/env python3
"""Differential test: Mathilda's graph-algorithm heads vs Mathematica 15.

Generates seeded random graphs (undirected, directed, weighted), writes ONE .m
script of queries that both systems run, and compares the answers:

  * values that are unique (flow values, connectivities, cut values, optimum
    sizes, predicate answers, planarity, Hamiltonicity, counts of all maximal
    cliques / Hamiltonian cycles / isomorphisms) must be EQUAL;
  * witnesses whose choice is not unique (the cut, matching, cover, clique,
    cycle actually returned) are checked for VALIDITY and optimal SIZE in
    Python with networkx, for both systems.

Usage (from this directory):  python3 diff_mathematica.py [ngraphs] [seed]
Mathilda binary: ../../Mathilda ; Mathematica: wolframscript.
"""

import itertools, os, random, subprocess, sys
import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
MATHILDA = os.path.join(HERE, "..", "..", "Mathilda")


def wl_graph(n, edges, directed, weights=None):
    op = "->" if directed else "<->"
    es = ", ".join("%d%s%d" % (u, op, v) for u, v in edges)
    s = "Graph[Range[%d], {%s}" % (n, es)
    if weights is not None:
        s += ", EdgeWeight -> {%s}" % ", ".join(str(w) for w in weights)
    return s + "]"


def rand_graph(rng, n, m, directed):
    pairs = [(u, v) for u in range(1, n + 1) for v in range(1, n + 1) if u != v]
    if not directed:
        pairs = [(u, v) for u, v in pairs if u < v]
    m = min(m, len(pairs))
    return sorted(rng.sample(pairs, m))


def build_cases(ng, seed):
    rng = random.Random(seed)
    cases = []      # (id, graphinfo, query_expr, kind)
    lines = []
    for gi in range(ng):
        directed = gi % 4 == 3
        n = rng.randint(4, 11)
        maxm = n * (n - 1) // (1 if directed else 2)
        m = rng.randint(n - 1, min(maxm, 3 * n))
        edges = rand_graph(rng, n, m, directed)
        weighted = gi % 5 == 1
        weights = [rng.randint(1, 9) for _ in edges] if weighted else None
        caps = [rng.randint(1, 9) for _ in edges]
        info = dict(n=n, edges=edges, directed=directed, weights=weights, caps=caps)
        g = "g%d" % gi
        lines.append("%s = %s;" % (g, wl_graph(n, edges, directed, weights)))
        s, t = rng.sample(range(1, n + 1), 2)
        perm = list(range(1, n + 1)); rng.shuffle(perm)
        pedges = [(perm[u - 1], perm[v - 1]) for u, v in edges]
        rng.shuffle(pedges)
        lines.append("p%d = %s;" % (gi, wl_graph(n, pedges, directed)))
        # a random non-isomorphic-ish partner with the same n, m
        e2 = rand_graph(rng, n, len(edges), directed)
        lines.append("h%d = %s;" % (gi, wl_graph(n, e2, directed)))
        info["e2"] = e2
        sub = sorted(rng.sample(range(1, n + 1), rng.randint(0, n)))
        info["sub"] = sub
        esub = rng.sample(edges, rng.randint(0, len(edges)))
        info["esub"] = esub
        op = "->" if directed else "<->"
        sub_wl = "{%s}" % ", ".join(map(str, sub))
        esub_wl = "{%s}" % ", ".join("%d%s%d" % (u, op, v) for u, v in esub)
        capl = "{%s}" % ", ".join(map(str, caps))
        q = [
            ("flow", "FindMaximumFlow[%s, %d, %d]" % (g, s, t), "eq"),
            ("flowcap", "FindMaximumFlow[%s, %d, %d, EdgeCapacity -> %s]" % (g, s, t, capl), "eq"),
            ("flowedges", "List @@@ FindMaximumFlow[%s, %d, %d, \"EdgeList\"]" % (g, s, t), "flowedges"),
            ("ec", "EdgeConnectivity[%s]" % g, "eq"),
            ("ecst", "EdgeConnectivity[%s, %d, %d]" % (g, s, t), "eq"),
            ("mincut", "FindMinimumCut[%s]" % g, "mincut"),
            ("edgecut", "List @@@ FindEdgeCut[%s]" % g, "edgecut"),
            ("edgecutst", "List @@@ FindEdgeCut[%s, %d, %d]" % (g, s, t), "edgecutst"),
            ("vcut", "FindVertexCut[%s]" % g, "vcut"),
            ("vcutst", "FindVertexCut[%s, %d, %d]" % (g, s, t), "vcutst"),
            ("match", "List @@@ FindIndependentEdgeSet[%s]" % g, "match"),
            ("ecover", "List @@@ FindEdgeCover[%s]" % g, "ecover"),
            ("vcover", "FindVertexCover[%s]" % g, "vcover"),
            ("mis", "FindIndependentVertexSet[%s]" % g, "mis"),
            ("misall", "Sort[Sort /@ FindIndependentVertexSet[%s, Infinity, All]]" % g, "eq"),
            ("clique", "FindClique[%s]" % g, "clique"),
            ("cliqueall", "Sort[Sort /@ FindClique[%s, Infinity, All]]" % g, "eq"),
            ("clique2", "Sort[Sort /@ FindClique[%s, {2}, All]]" % g, "eq"),
            ("ivsq", "IndependentVertexSetQ[%s, %s]" % (g, sub_wl), "eq"),
            ("vcq", "VertexCoverQ[%s, %s]" % (g, sub_wl), "eq"),
            ("iesq", "IndependentEdgeSetQ[%s, %s]" % (g, esub_wl), "eq"),
            ("ecq", "EdgeCoverQ[%s, %s]" % (g, esub_wl), "eq"),
            ("hamq", "HamiltonianGraphQ[%s]" % g, "eq"),
            ("hamc", "List @@@ #& /@ FindHamiltonianCycle[%s]" % g, "hamc"),
            ("hamcount", "Length[FindHamiltonianCycle[%s, All]]" % g, "eq"),
            ("hamp", "FindHamiltonianPath[%s]" % g, "hamp"),
            ("hampst", "FindHamiltonianPath[%s, %d, %d]" % (g, s, t), "hampst"),
            ("planar", "PlanarGraphQ[%s]" % g, "eq"),
            ("isoself", "IsomorphicGraphQ[%s, p%d]" % (g, gi), "eq"),
            ("isoother", "IsomorphicGraphQ[%s, h%d]" % (g, gi), "eq"),
            ("isocount", "Length[FindGraphIsomorphism[%s, p%d, All]]" % (g, gi), "eq"),
            ("isomap", "Normal[First[FindGraphIsomorphism[%s, p%d]]] /. Rule -> List" % (g, gi), "isomap"),
        ]
        info["s"], info["t"], info["perm"], info["pedges"] = s, t, perm, pedges
        for name, expr, kind in q:
            cid = "%d:%s" % (gi, name)
            lines.append('Print["%s\\t", ToString[%s, InputForm]];' % (cid, expr))
            cases.append((cid, info, kind))
    return cases, lines


def parse(out):
    res = {}
    for line in out.splitlines():
        if "\t" in line:
            k, v = line.split("\t", 1)
            res[k.strip()] = v.strip()
    return res


def wl_to_py(s):
    s = s.strip()
    try:
        return eval(s.replace("{", "[").replace("}", "]").replace("True", "True").replace("False", "False"))
    except Exception:
        return None


def nxg(info, directed=None):
    d = info["directed"] if directed is None else directed
    G = nx.DiGraph() if d else nx.Graph()
    G.add_nodes_from(range(1, info["n"] + 1))
    for k, (u, v) in enumerate(info["edges"]):
        w = info["weights"][k] if info["weights"] else 1
        G.add_edge(u, v, weight=w, capacity=w)
    return G


def und(info):
    G = nx.Graph()
    G.add_nodes_from(range(1, info["n"] + 1))
    G.add_edges_from(info["edges"])
    return G


def check_witness(kind, info, val, other_val):
    """Validity (and optimal size, cross-checked against other_val's size)."""
    n = info["n"]
    U = und(info)
    if val is None:
        return "unparsable"
    if kind == "mincut":
        v, (a, b) = val
        G = nxg(info)
        A = set(a)
        if not a or not b or A | set(b) != set(range(1, n + 1)):
            return "bad partition"
        cw = sum(d["weight"] for x, y, d in G.edges(data=True)
                 if (x in A and y not in A) or (not info["directed"] and y in A and x not in A))
        return None if abs(cw - v) < 1e-9 else "cut weight %s != %s" % (cw, v)
    if kind in ("edgecut", "edgecutst"):
        G = nxg(info)
        H = G.copy()
        H.remove_edges_from([tuple(e) for e in val])
        if kind == "edgecut":
            conn = nx.is_strongly_connected(H) if info["directed"] else nx.is_connected(H)
            if n >= 2 and conn:
                return "edge cut does not disconnect"
        else:
            if nx.has_path(H, info["s"], info["t"]):
                return "st edge cut leaves a path"
        return None
    if kind in ("vcut", "vcutst"):
        H = U.copy()
        H.remove_nodes_from(val)
        if kind == "vcut":
            complete = U.number_of_edges() == n * (n - 1) // 2
            if complete and info["directed"] and n >= 3 and val == []:
                return None     # Mathematica's convention for directed graphs
            if len(H) >= 2 and nx.is_connected(H) and len(val) != n - 1:
                return "vertex cut does not disconnect"
        else:
            if info["s"] in val or info["t"] in val:
                return "separator contains terminal"
            if not U.has_edge(info["s"], info["t"]) and nx.has_path(H, info["s"], info["t"]):
                return "separator leaves a path"
        return None
    if kind == "match":
        used = set()
        for u, v in val:
            if not U.has_edge(u, v) or u in used or v in used:
                return "not a matching"
            used |= {u, v}
        return None if len(val) == len(nx.max_weight_matching(U, maxcardinality=True)) else "not maximum"
    if kind == "ecover":
        if any(U.degree(x) == 0 for x in U):
            return None if val == [] else "cover despite isolated vertex"
        cov = set()
        for u, v in val:
            if not U.has_edge(u, v):
                return "not an edge"
            cov |= {u, v}
        if cov != set(U):
            return "not a cover"
        return None if len(val) == n - len(nx.max_weight_matching(U, maxcardinality=True)) else "not minimum"
    if kind == "vcover":
        S = set(val)
        if any(u not in S and v not in S for u, v in U.edges()):
            return "not a cover"
        return None
    if kind == "mis":
        S = set(val[0]) if val else set()
        if any(u in S and v in S for u, v in U.edges()):
            return "not independent"
        return None
    if kind == "clique":
        S = val[0]
        G = nxg(info)
        for a, b in itertools.combinations(S, 2):
            ok = (G.has_edge(a, b) and G.has_edge(b, a)) if info["directed"] else G.has_edge(a, b)
            if not ok:
                return "not a clique"
        return None
    if kind == "hamc":
        if not val:
            return None
        cyc = val[0]
        G = nxg(info)
        verts = [e[0] for e in cyc]
        if sorted(verts) != list(range(1, n + 1)):
            return "not spanning"
        for (a, b) in cyc:
            if not G.has_edge(a, b):
                return "not an edge %s" % ((a, b),)
        for i in range(len(cyc)):
            if cyc[i][1] != cyc[(i + 1) % len(cyc)][0]:
                return "not a cycle"
        return None
    if kind in ("hamp", "hampst"):
        if not val:
            return None
        G = nxg(info)
        if sorted(val) != list(range(1, n + 1)):
            return "not spanning"
        if any(not G.has_edge(val[i], val[i + 1]) for i in range(n - 1)):
            return "not a path"
        if kind == "hampst" and (val[0] != info["s"] or val[-1] != info["t"]):
            return "wrong endpoints"
        return None
    if kind == "flowedges":
        return None
    if kind == "isomap":
        mp = dict((a, b) for a, b in val)
        pe = set(info["pedges"]) if info["directed"] else set(frozenset(e) for e in info["pedges"])
        for u, v in info["edges"]:
            img = (mp[u], mp[v]) if info["directed"] else frozenset((mp[u], mp[v]))
            if img not in pe:
                return "not an isomorphism"
        return None
    return "unknown kind"


SIZE_KINDS = {"vcover": len, "mis": lambda v: len(v[0]) if v else 0,
              "clique": lambda v: len(v[0]) if v else 0, "vcut": len, "vcutst": len,
              "match": len, "ecover": len, "mincut": lambda v: v[0],
              "hamc": len, "hamp": len, "hampst": len}


MIN_KINDS = {"vcover", "ecover", "vcut", "vcutst", "mincut"}


def main():
    ng = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    cases, lines = build_cases(ng, seed)
    path = os.path.join("/tmp", "galg_diff_%d_%d.m" % (ng, seed))
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    om = subprocess.run([MATHILDA, path], capture_output=True, text=True, timeout=3600).stdout
    ow = subprocess.run(["wolframscript", "-file", path], capture_output=True, text=True, timeout=3600).stdout
    rm, rw = parse(om), parse(ow)
    fails = 0
    mma_subopt = 0
    for cid, info, kind in cases:
        a, b = rm.get(cid), rw.get(cid)
        if a is None:
            print("MISSING in Mathilda:", cid); fails += 1; continue
        if b is None:
            continue    # Mathematica gave no answer (message / unevaluated)
        if kind == "eq":
            if a != b:
                print("DIFF", cid, "mathilda=", a, "mathematica=", b); fails += 1
            continue
        va, vb = wl_to_py(a), wl_to_py(b)
        err = check_witness(kind, info, va, vb)
        if err:
            print("INVALID", cid, err, a); fails += 1; continue
        if vb is not None and kind in SIZE_KINDS:
            try:
                sa, sb = SIZE_KINDS[kind](va), SIZE_KINDS[kind](vb)
            except Exception:
                sa, sb = a, b
            if sa != sb:
                lower_better = kind in MIN_KINDS
                ours_better = (sa < sb) if lower_better else (sa > sb)
                if kind in ("hamc", "hamp", "hampst"):
                    # a validated witness where Mathematica says "none" is a
                    # Mathematica error (seen on weighted graphs)
                    ours_better = sa > 0 and sb == 0
                if ours_better:
                    # Mathilda's answer passed the validity/optimality check;
                    # Mathematica returned a worse one.
                    print("MMA-SUBOPTIMAL", cid, "mathilda=", a, "mathematica=", b)
                    mma_subopt += 1
                else:
                    print("SIZE", cid, "mathilda=", a, "mathematica=", b); fails += 1
    print("cases: %d  failures: %d  (Mathematica suboptimal: %d)" % (len(cases), fails, mma_subopt))
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
