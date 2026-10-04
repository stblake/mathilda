---
source: src/poly/facpoly_list.c
---
**Algorithm.** `builtin_factorlist` is a thin wrapper over `Factor`. It forwards
every argument verbatim — the polynomial and any options (`GaussianIntegers`,
`Extension -> {a1, ...}`, ...) — to `Factor[poly, opts...]`, then calls the
shared `product_to_pair_list` to split the evaluated product into
`{factor, exponent}` pairs. Positions after the first must be option rules
(`Rule`/`RuleDelayed`); a non-rule there raises `FactorList::nonopt` and leaves
the call unevaluated.

**Data structures.** `product_to_pair_list` walks the `Times` (or the single
bare factor) returned by `Factor`. `absorb_factor` classifies each multiplicative
factor: a number literal (`Integer`, `Rational`, `Real`, `Complex`, MPFR) is
multiplied into a running overall numerical factor `c`; a `Power[base, e]` with
an *integer* `e` becomes the pair `{base, e}` (a multiplicity, positive or
negative — denominator factors of a rational function carry negative exponents);
anything else, including `Power[base, 1/2]` = `Sqrt[base]`, is an irreducible
factor `{factor, 1}`. Bases and exponents accumulate in parallel growable
`Expr**` arrays, assembled into `{{c, 1}, {base_i, exp_i}, ...}` — the leading
pair is always the numerical factor (`{1, 1}` when there is none).

**Complexity / limits.** All real cost is in `Factor`; the pair-split is a single
linear pass over its output. A symbolic structural head — no packed/NDArray or
`Compile[]` path. `Listable` and `Protected`.
