# Experimental context

`` Experimental` `` holds native C kernels for workloads that otherwise run as Wolfram-language code. Each head is written to match a reference definition **exactly** (SameQ, term order included). Its contract can change between releases. `` Experimental`FeynCalc` `` covers FeynCalc-shaped computations; the sources are in `src/experimental/feyncalc/`.

The paired `.m` definitions are timed against these kernels in `benchmarks/123-native-vs-m/`.

## `` Experimental`FeynCalc`DiracTr ``

`` Experimental`FeynCalc`DiracTr[a1, ..., an] `` is the trace of n Dirac gamma matrices in 4 dimensions, with metric `g` (`` Global`g ``):

```
Tr[a1 ... an] = 4 Sum_{perfect matchings M} sign(M) Prod_{(i<j) in M} g[ai, aj]
```

The sum has (n-1)!! terms. In each `g`, the first index is the argument that comes earlier. `DiracTr[]` is `4`, and an odd n gives `0`.

The value is exactly that of this recursion (Rolf Mertig, Woxi issue #2142):

```
diracTr[] := 4;
diracTr[a_, b_] := 4 g[a, b];
diracTr[a_, rest___] := Module[{list = {rest}, n},
  Sum[(-1)^(m + 1) * Expand[g[a, list[[m]]] * diracTr @@ Delete[list, m]],
    {m, 1, Length[list]}]];
```

```
In[1]:= Experimental`FeynCalc`DiracTr[a, b, c, d]
Out[1]= 4 g[a, d] g[b, c] - 4 g[a, c] g[b, d] + 4 g[a, b] g[c, d]

In[2]:= Experimental`FeynCalc`DiracTr[a, b, a, b]
Out[2]= 4 g[a, b]^2 + 4 g[a, b] g[b, a] - 4 g[a, a] g[b, b]

In[3]:= Length[Experimental`FeynCalc`DiracTr @@ Array[x, 14]]
Out[3]= 135135
```

- **Repeated arguments.** Equal factors fold into powers, equal terms add, and zero terms drop out, as `Expand` would do. Indices are not contracted: `g[mu, mu]` stays.
- **Evaluator-free.** The answer is built directly in canonical order, never by calling the evaluator. The `n(n-1)/2` nodes `g[ai, aj]` are shared by reference across all terms, and every node is stamped as an evaluated fixed point.
  - Distinct arguments: the terms are ordered by an LSD radix sort on their factor ranks, with no comparisons.
  - Repeated arguments: factor and term order come from the canonical comparator.
- **Rules on g.** When `g`, an argument's head (through UpValues), or `Plus`/`Times`/`Power` carries rules or attributes, the raw sum is returned wrapped in `Expand`, so the evaluator applies those rules as the recursion would.
- **Size guard.** More than 2^32 terms (n ≥ 22) gives the message `DiracTr::toobig` and the call is left unevaluated.
- **Attributes:** `{Protected}`. The head is symbolic, so it has no packed-array or `Compile` surface.
