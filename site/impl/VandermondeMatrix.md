---
source: src/linalg/vandermondemat.c
---
**Algorithm.** `builtin_vandermondematrix` constructs the matrix of successive powers of a node list: entry `(i,j) = x_i^(j−1)`. `VandermondeMatrix[{x1,...,xn}]` is `n × n`; `VandermondeMatrix[{x}, k]` is `n × k`. The builder `vm_build` emits each entry via `vm_entry`: the exponent-0 column is the literal Integer `1` (so `0^0` reads as 1, matching interpolation semantics), and every other entry is a `Power[x_i, j]` node which the evaluator later folds (numeric powers to their value, `Power[x,1]` to `x`), leaving symbolic nodes as clean `Power` expressions. Nodes are deep-copied and need not be numeric or distinct.

**Limits.** Zero arguments emit `VandermondeMatrix::argt`. The single-matrix structured-array conversion form `VandermondeMatrix[vmat]` is unsupported (Mathilda has no structured-array representation), so a list-of-lists argument (`vm_is_matrix`) is left unevaluated.
