---
source: src/core.c
---
**Algorithm.** `builtin_inexactnumberq` delegates to `core_is_inexact_number`,
which is the mirror of `core_is_exact_number`:

- a machine `EXPR_REAL` — `True`;
- an arbitrary-precision `EXPR_MPFR` number (when built `USE_MPFR`) — `True`;
- a `Complex[re, im]` — `True` when **either** part is inexact (the test recurses
  into the two parts and ORs them);
- anything else — `False`.

So an exact number (integer, rational, or an all-exact `Complex`) is `False`, and a
`Complex` with even one floating-point part is `True`. Among numbers it is exactly
the complement of `ExactNumberQ`.

**Data structures.** None beyond the argument; the recursion descends only the two
parts of a `Complex`.

**Complexity / limits.** `O(1)`. A non-number (symbol or symbolic expression) is
`False`. A call that is not single-argument declines and stays symbolic.
