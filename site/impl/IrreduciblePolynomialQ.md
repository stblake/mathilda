---
source: src/poly/irrpolyq.c
---
**Algorithm.** `builtin_irreduciblepolynomialq` always returns `True`/`False` on a structurally valid call. It parses the `GaussianIntegers` and `Extension` options, resolves the factoring field by precedence (`Extension -> All` = absolute irreducibility; explicit `α`/`{α_i}` = `Q(α)`/compositum; `Automatic` = `extension_autodetect`; `GaussianIntegers -> True` or a complex coefficient = `Q(i)`; else `Q`), then `Factor`s the polynomial over that field and counts non-constant factors with multiplicity (`irr_dispatch`): `0` → False (constant), `1` → True, `>= 2` → False. `Extension -> All` treats degree-1 univariate as absolutely irreducible and approximates the multivariate case by factoring over `Q(i)`. For multivariate inputs that name an extension, a Hilbert-style specialisation probe (`irr_multivariate_specialize_probe`) can flip a True verdict to False (covering cases the cheap univariate-only extension path misses).

**Data structures.** Delegates to the polynomial factoring subsystem (`facpoly`/`qafactor`); options are read off `Rule` heads in the argument list. Wrong arity emits `IrreduciblePolynomialQ::argx`; malformed options emit `::nonopt` (both return `NULL`). Registered `Listable`, so list inputs thread before this handler runs.

**Complexity / limits.** Dominated by the underlying `Factor`. The `Extension -> All` multivariate approximation is incomplete — it detects `Q(i)`-conjugate splits like `x^2+y^2` but not reducibility over a general real quadratic field.
