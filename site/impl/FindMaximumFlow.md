---
references:
  - "E. A. Dinic, *Algorithm for solution of a problem of maximum flow in a network with power estimation*, Soviet Math. Dokl. **11** (1970) 1277-1280."
source: src/graph/galg_flow.c
---
**Algorithm.** `builtin_find_maximum_flow` computes the value of a maximum flow
from `s` to `t` (either may be a list of sources/sinks) with **Dinic's
algorithm**. A fourth argument selects the property: `"FlowValue"` (default),
`"FlowMatrix"` (a dense `n x n` matrix of per-edge flows), or `"EdgeList"` (the
edges carrying flow, oriented along it, in flow-matrix row-major order).
Capacities come from `EdgeCapacity -> {c1, ...}` (else the graph's own
`EdgeCapacity`, else `1`); `EdgeWeight` is ignored. A `VertexCapacity` option
caps the flow through each vertex, modelled by the split-vertex construction
(`v_in -> v_out`). Several sources/sinks are wired to a super-source/super-sink.
An undirected edge carries flow either way up to its capacity; a directed edge
only forwards.

**Data structures.** Capacities are parsed into a scaled-`int64` `GfCap` (exact
integers; reals scaled by a shared power of two; `Infinity` as a finite stand-in
above every finite sum) — including a packed-buffer fast path that reads an
`NDArray` capacity list without unpacking. The network is a `GfNet` CSR residual
graph (adjacent forward/reverse arc pairs) with level-BFS truncated at the sink's
level and an iterative current-arc blocking-flow search. The flow matrix is
offered to the packer (`pack_offer`).

**Complexity / limits.** `O(E sqrt(E))` on unit capacities; a handful of phases
in practice. Exact `Integer` for integer capacities, `Real` otherwise, `Infinity`
when the flow reaches the `Infinity` stand-in. A negative or symbolic capacity
leaves the call unevaluated.
