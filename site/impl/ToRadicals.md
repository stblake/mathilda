---
references:
  - "G. Cardano, *Ars Magna*, 1545 (cubic); L. Ferrari (quartic resolvent, via Cardano)."
source: src/radicals.c
---
**Algorithm.** `builtin_to_radicals` (`builtin_to_radicals`) converts held `Root[Function[poly], k]` objects into closed-form radical expressions. The top-level walker is a structural recurrence that rebuilds every `EXPR_FUNCTION` node, so a `Root` buried inside `List`/`Equal`/`Less`/`And`/`Or`/... is handled identically.

Per `Root` node: (1) extract the polynomial body, accepting both the `Slot[1]` form `Function[expr]` and the bound-variable form `Function[t, expr]`; (2) substitute the slot/variable with a fresh symbol `x$` so the standard `get_coeff`/`get_degree_poly` univariate machinery applies; (3) dispatch on degree d — `d=1` linear (`-c0/c1`), `d=2` quadratic formula, `d=3` Cardano, `d=4` Ferrari (depressed quartic + resolvent cubic), `d≥5` only the binomial fast-path `a·x^n + b`, otherwise the `Root` is left untouched — each path producing all d radical roots as a fresh `Expr**`; (4) select the k-th root in Mathilda's canonical `Root` ordering by computing `N[Root[poly, k]]` at machine precision (`root_numericalize`) and picking the radical root closest in the complex plane, falling back to the natural per-formula order (index `k-1`) when coefficients are parametric and numeric evaluation is unavailable.

**Data structures.** `Expr*` trees; degree dispatch reuses the polynomial coefficient extractors from `src/poly/poly.c`. Intermediate radical-expression bookkeeping rides on `eval_and_free` and the Plus/Times/Power normalisation. Inputs are borrowed and deep-copied into the output.

**Complexity / limits.** Closed forms exist only up to degree 4 (Abel–Ruffini); degree ≥ 5 is supported solely for binomials. Root selection costs one numeric `Root` evaluation per node.
