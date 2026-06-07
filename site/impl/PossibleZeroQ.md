---
references:
  - "J. T. Schwartz, \"Fast probabilistic algorithms for verification of polynomial identities\", JACM 27 (1980)."
  - "R. Zippel, \"Probabilistic algorithms for sparse polynomials\", EUROSAM 1979."
source: src/zero_test.c
---
**Algorithm.** `builtin_possible_zero_q` (`src/zero_test.c`) calls `zero_test_decide`, a staged hybrid symbolic-numeric pipeline that early-exits on the first definite verdict:

- *Stage 0 — structural:* O(1) shortcuts for literal `Integer`/`Real`/`BigInt`/`MPFR` zero, `Complex[0,0]`, lists of zeros, and unbound symbols (`decide_structural`).
- *Stage 1 — rational normalisation:* `Together ∘ Cancel` plus `Expand`, then a polynomial zero test, deciding every identity in `Q(x_1,...,x_n)` (`decide_rational`). A `True` here is trusted; a `False` is not trusted alone.
- *Stage 2 — numeric:* for symbol-free inputs, numericalize at machine precision and compare `|z|` against an IEEE catastrophic-cancellation threshold, climbing a precision ladder (53 -> 200 -> 500 -> 1000 bits) while ambiguous; a true zero shrinks geometrically across precisions (`decide_numeric`).
- *Stage 3 — Schwartz–Zippel:* for inputs with free symbols, substitute each free symbol with a random rational drawn from `Q[i]` (to probe branch cuts) and require several independent Stage-2 confirmations (`decide_schwartz_zippel`).

**Result mapping.** A definite `False` returns `False`; both `True` and `UNKNOWN` return `True`, matching Mathematica's documented "assume zero when uncertain" behaviour (the accompanying `PossibleZeroQ::ztest1` message is not emitted). The symbol is `Listable`.
