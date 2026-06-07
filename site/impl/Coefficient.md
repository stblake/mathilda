---
source: src/poly/poly.c
---
**Algorithm.** `builtin_coefficient` (in `src/poly/poly.c`) extracts the coefficient of `form^n` (default n = 1) from `expr`. It first runs `expr_expand` to a flat sum of monomials, then decomposes both the target `form` and each summand into base–exponent pairs via `decompose_to_bp`. The helper `get_k` computes the integer power at which the target's base(s) divide each term; terms with `k == n` contribute, with the matching base factors stripped out (`get_k` handles multi-factor monomial forms like `x y`). For `n == 0` the whole term is taken. Surviving residual factors are reassembled with `internal_times`/`internal_plus`.

**Data structures.** `BPList` — a list of `{base, exp}` pairs (initialised with `bp_init`, freed with `bp_free`) — is the core representation for both the target form and each term.
