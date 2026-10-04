---
source: src/poly/exponent.c
---
**Algorithm.** `builtin_exponent` first expands its argument (`expr_expand`) and
splits the result into additive terms — the arguments of a top-level `Plus`, or
the whole expression as a single term. The genuine zero polynomial expands to a
numeric `0`, which has *no* terms, so its exponent set is empty and the wrapper
fires on nothing (`Max[]` = `-Infinity`). `form` is decomposed into `(base, fe)`
pairs by `form_pairs` — a symbol or kernel gives `(form, 1)`, a `Power[b, e]`
gives `(b, e)`, and a `Times` gives one pair per non-numeric factor. For each
monomial term `exp_in_term` reads the power of each base (`base_exp_in_monomial`:
1 if the base *is* the term, the exponent of a matching `Power`, the sum over a
`Times`, else `0`); for a single-base form the term's exponent is that power
divided by `fe`, and for a product form it is the `Min` over bases of
`base-exponent / fe` — the largest `k` with `form^k` dividing the term. The
exponents are insertion-sorted by `expr_compare`, de-duplicated with `expr_eq`,
and handed to `h` (default `Max`); `Exponent[expr, form, h]` substitutes any
other `h`.

**Data structures.** Everything is `Expr` trees. A growable `Expr**` holds the
per-term exponent set (one entry per additive term before dedup), and
`form_pairs` fills parallel `Expr**` arrays of bases and form-exponents. Symbolic
or rational exponents are kept as expressions and flow through `Max`/`Min`
unevaluated (so `Exponent[x^(1/2) + x^(1+n), x]` returns a `Max[1/2, 1 + n]`).

**Complexity / limits.** Dominated by the initial `Expand`; the exponent-set
sort is an `O(n^2)` insertion sort over the small set of distinct term degrees.
`Exponent` is **purely syntactic** — no zero-coefficient recognition, so a term
with a coefficient that is zero but not in normal form still counts. It is a
symbolic structural head with no packed/NDArray or `Compile[]` path; `Listable`
makes `Exponent[expr, {f1, f2, ...}]` thread into a per-form list for free.
