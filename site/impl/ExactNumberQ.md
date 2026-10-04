---
source: src/core.c
---
**Algorithm.** `builtin_exactnumberq` delegates to `core_is_exact_number`, which is
exact iff the argument is built entirely from exact pieces:

- an `EXPR_INTEGER` or `EXPR_BIGINT` — `True`;
- a `Rational[p, q]` — `True`;
- a `Complex[re, im]` — `True` only when **both** `re` and `im` are themselves
  exact (the test recurses into the two parts and ANDs them);
- anything else — `False`.

The complement among numbers is `InexactNumberQ`: a machine `Real` or MPFR number
is inexact, so `ExactNumberQ` is `False` for it.

**Data structures.** None beyond the argument tree; the recursion only descends the
two parts of a `Complex`, so it is bounded.

**Complexity / limits.** `O(1)` (a `Complex` adds two leaf checks). A non-number —
a symbol or symbolic expression such as `Pi` or `x` — is `False`, not left
symbolic, since it is a definite negative answer. A call that is not
single-argument declines and stays symbolic.
