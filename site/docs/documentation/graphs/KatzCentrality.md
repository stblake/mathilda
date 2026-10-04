# KatzCentrality

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`KatzCentrality[g, a] gives the Katz centrality x = a A^T x + 1 of each vertex; KatzCentrality[g, a, b] uses b (a number or a list) in place of 1. Edge weights are ignored.`**

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

```mathematica
In[1]:= KatzCentrality[PathGraph[3], 0.1]
Out[1]= {1.12245, 1.22449, 1.12245}

In[2]:= KatzCentrality[CompleteGraph[4], 1/2]
Out[2]= {-2.0, -2.0, -2.0, -2.0}

In[3]:= KatzCentrality[PathGraph[3], 0.1, 2]
Out[3]= {2.2449, 2.44898, 2.2449}

In[4]:= KatzCentrality[PathGraph[3], 0.1, {1, 2, 3}]
Out[4]= {1.2449, 2.44898, 3.2449}

In[5]:= KatzCentrality[PathGraph[3], 0]
Out[5]= {1, 1, 1}

In[6]:= KatzCentrality[Graph[{1,2,3},{}], 0.3]
Out[6]= {1, 1, 1}

In[7]:= KatzCentrality[CompleteGraph[4], 1/3]
Out[7]= KatzCentrality[Graph[<4 vertices, 6 edges>], 1/3]
```

### Applications (4)

A vertex-transitive graph: all scores equal

```mathematica
In[8]:= KatzCentrality[CycleGraph[4], 1/4]
Out[8]= {2.0, 2.0, 2.0, 2.0}
```

With attenuation 1/2, each sink inherits half its source's score

```mathematica
In[9]:= KatzCentrality[Graph[{1 -> 2, 1 -> 3}], 1/2]
Out[9]= {1.0, 1.5, 1.5}
```

A directed triangle, small attenuation

```mathematica
In[10]:= KatzCentrality[Graph[{1 -> 2, 2 -> 3, 3 -> 1}], 1/10]
Out[10]= {1.11111, 1.11111, 1.11111}
```

A uniform baseline b = 2 in place of 1

```mathematica
In[11]:= KatzCentrality[PathGraph[{1, 2, 3}], 1/2, 2]
Out[11]= {6.0, 8.0, 6.0}
```

## Implementation notes

**Algorithm.** `builtin_katz_centrality` solves the Katz fixed point
`x = b + α Aᵀx` (Katz 1953), where `α` is the attenuation factor (second argument)
and `b` the baseline (the third argument — a scalar or a length-`n` list —
defaulting to the all-ones vector). The transpose `Aᵀ` is realised directly: the
graph is built as an in-arc CSR (`GMET_IN`), so each vertex accumulates over its
in-neighbours. The primary path is a fixed-point iteration
`xₙₑₓₜ = b + α Aᵀx`, up to 20000 sweeps, declared converged when the L1 change
falls to `≤ 1e-15` of the L1 norm; a growth guard aborts after 50 consecutive
increases (a divergent `α`). On non-convergence with `n ≤ 4000` it falls back to a
dense LAPACK solve (`mat_lapack_dgesv`) of `(I − α Aᵀ) x = b` in column-major form.

**Data structures.** A `GMET_IN` `GmetCSR`, two `double` work vectors for the
iteration, and — only on the dense fallback — an `n × n` column-major matrix `M`
holding `I − α Aᵀ` (then its LU factors) plus a pivot array. Results are memoised
through `gmet_cache`, keyed on the whole call.

**Complexity / limits.** `O(iter · E)` for the iteration, `O(n³)` for the dense
fallback (hence the `n ≤ 4000` cap). Edge weights are ignored. Two cases short-
circuit to `b` exactly (with exact entries preserved): a zero `α`, and a graph with
no edges. The fallback carries a singularity guard — if the smallest LU pivot is
below `1e-12` of the largest, the system is treated as singular (an `α` at a
reciprocal eigenvalue of `A`) and the call is left unevaluated rather than
returning a `~1e16` garbage answer.

- Solves `x = a A^T x + b`, with `b` = 1 by default, a number, or a list.
  Weights ignored.
- Answered even beyond the convergence radius, like Wolfram (`K4` with
  `a = 1/2` gives `-2`); uses a dense LAPACK solve, `n <= 4000`. A system
  singular up to rounding (LU pivot ratio below `1e-12`, e.g. `K4` with
  `a = 1/3`) is left unevaluated. Mathematica 15 returns rounding noise there
  (about `1.8*10^16` per vertex for `K4`, `a = 1/3`).
- No edges, or an exact `a = 0`, return `b` exactly.

**Attributes:** `Protected`.

## References

- L. Katz, *A new status index derived from sociometric analysis*, Psychometrika **18** (1953) 39-43.
- Source: [`src/graph/gmet_spectral.c`](https://github.com/stblake/mathilda/blob/main/src/graph/gmet_spectral.c)
- Specification: [`docs/spec/builtins/graphs.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/graphs.md)
- Tests: [`tests/test_graph_metrics.c`](https://github.com/stblake/mathilda/blob/main/tests/test_graph_metrics.c)

## Notes & additional examples

### Notes

The score vector solves `x = α Aᵀ x + b`: each vertex gets a baseline `b`
(the third argument, defaulting to `1`) plus `α` times the summed scores of the
vertices pointing at it. The attenuation factor `α` must be smaller than the
reciprocal of the largest eigenvalue of the adjacency matrix for the sum to
converge; at that boundary the system is singular and the call is left
unevaluated.

Because the propagation uses `Aᵀ`, a directed edge `u -> v` lets `u`'s status flow
into `v`. Edge weights are ignored — only the graph's adjacency enters. The
baseline may also be a length-`n` list to weight the vertices differently.
