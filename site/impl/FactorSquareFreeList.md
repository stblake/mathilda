---
source: src/poly/facpoly_list.c
---
**Algorithm.** `builtin_factorsquarefreelist` is the square-free sibling of
`FactorList`: it forwards every argument verbatim to
`FactorSquareFree[poly, opts...]` (so all `Extension` handling lives in
`FactorSquareFree`, which runs the Yun/Musser decomposition via GCDs of the
polynomial with its derivative), then splits the resulting product into
`{factor, exponent}` pairs through the same shared `product_to_pair_list` that
`FactorList` uses. If `FactorSquareFree` declines — coming back as an
unevaluated `FactorSquareFree[...]` expression on a bad arity or an unrecognised
option — the decline is propagated (the head is detected and `NULL` returned),
so the call stays unevaluated rather than wrapping the inert form as a spurious
factor.

**Data structures.** Identical to `FactorList`: `absorb_factor` folds number
literals into a running overall numerical factor `c`, reads an integer-exponent
`Power[base, e]` as `{base, e}`, and treats everything else as `{factor, 1}`,
accumulating parallel `Expr**` bases/exponents assembled into
`{{c, 1}, {base_i, exp_i}, ...}`. The leading pair is the overall numerical
factor (`{1, 1}` when absent); the `exp_i` here are the multiplicities grouping
repeated factors.

**Complexity / limits.** Cost is `FactorSquareFree`'s (polynomial GCDs — cheaper
than full `Factor`, and sufficient when only multiplicities are needed); the
pair-split is one linear pass. A symbolic structural head — no packed/NDArray or
`Compile[]` path. `Listable` and `Protected`.
