---
source: src/simp/simp_builtins.c
---
**Algorithm.** `builtin_element` decides domain membership `Element[x, dom]`,
returning `True`/`False`, or staying unevaluated (NULL) when undetermined.
`Element[{x1,...}, dom]` and `Element[x1|x2|..., dom]` are treated as the
conjunction over components: collapse to `True`/`False` only when every component
decides, otherwise leave the original in place so the assumption framework
(`ctx_walk` in `simp_assume.c`) can still split it into per-variable facts.

For a symbol domain, it reads the current `$Assumptions` (`read_dollar_assumptions`),
builds an `AssumeCtx`, and calls `element_decide`, which returns 1/0/-1. The
decision first checks the fact set directly (`fact_in_domain`), then applies
literal rules per domain: `Integers` (Integer/BigInt yes, integer-valued Real
yes, Rational/Complex no), `Rationals`/`Algebraics`/`Reals`/`Complexes` (with
canonical `Complex` always carrying a non-zero imaginary part, so `Element[I,
Reals]` is False), `Booleans` (only literal `True`/`False`), and
`Primes`/`Composites` (delegating to `PrimeQ`). Provable facts from assumptions
are consulted via `prov_int`/`prov_re`. Unknown or symbolic cases return -1 →
the call stays unevaluated.

**Data structures.** `AssumeCtx` (flat fact array) built per call from
`$Assumptions`; `element_decide` is a `strcmp`-dispatched cascade with no
persistent state.
