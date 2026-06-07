---
source: src/numbertheory.c
---
**Algorithm.** `builtin_binomial` dispatches on argument kind. (1) Integer/integer: it coerces both to `mpz_t` and calls GMP's `mpz_bin_ui` (requires the lower index to fit a `ulong`), returning `0` when the lower index is negative or exceeds a non-negative `n`, and handling negative `n` via the Pascal/upper-negation extension `C(n,k) = (-1)^k C(k-n-1, k)`. (2) Machine reals: the Gamma form `tgamma(n+1)/(tgamma(m+1) tgamma(n-m+1))`. (3) Symmetry/polynomial reduction: if `n - m` evaluates to a small non-negative integer, or `m` itself is a small concrete non-negative integer (`<= 32`), it expands the falling-factorial polynomial via `binomial_polynomial` so that symbolic `n` produces a degree-`m` polynomial that downstream `Expand`/`D` can act on. The `n - m` Subtract is evaluated with arithmetic warnings muted (the exploratory difference may hit a spurious `Power::infy`).

**Data structures.** GMP `mpz_t` for the exact path (results normalised by `expr_bigint_normalize`); symbolic expansions build `Times`/`Plus` trees through `eval_and_free`.
