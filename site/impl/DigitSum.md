---
source: src/int.c
---
**Algorithm.** `builtin_digitsum` takes an integer `n` (a machine `int64` or a
GMP bignum) and an optional base `b >= 2` (default 10). Each argument is
validated first: a concrete-but-non-integer argument (Real, Rational, Complex)
raises `DigitSum::int` / `DigitSum::base` through the `mth_message` funnel and
returns `NULL`, while a symbolic argument is left unevaluated silently so a
downstream rewrite can still rewire the call. The sign of `n` is discarded
(`mpz_abs`), then the base-`b` digits are summed by repeated Euclidean division
with no intermediate digit list: a base that fits in `unsigned long` takes the
fast path `mpz_tdiv_q_ui`, which returns the remainder directly, and a bignum
base takes `mpz_tdiv_qr`. `DigitSum[0]` is `0`.

**Data structures.** All arithmetic is done in GMP `mpz_t`, so machine integers
and arbitrary-precision bignums are handled uniformly on one path. The running
sum accumulates in a single `mpz_t` in one pass over the digits and is demoted
back to an `EXPR_INTEGER` when it fits a signed long, otherwise returned as a
bigint via `expr_new_bigint_from_mpz`.

**Complexity / limits.** `O(d)` divisions for a `d`-digit base-`b`
representation. `DigitSum` is `Listable | NumericFunction | Protected`, so it
threads element-wise over ordinary lists (`DigitSum[{...}]`), but — unlike its
sibling `IntegerDigits` — it has **no** NDArray buffer kernel and no `Compile[]`
lowering: a visible `NDArray` integer argument is left unevaluated, and
`CompileDiagnostics` reports `Compiled -> False`. Semantically it is
`Total[IntegerDigits[n, b]]`.
