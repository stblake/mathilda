---
source: src/bitwise/bitlength.c
---
**Algorithm.** `builtin_bitlength` gives the number of binary bits needed to represent an
integer `n`. All arithmetic is done in GMP, so machine integers (`EXPR_INTEGER`) and bignums
(`EXPR_BIGINT`) are handled uniformly. The magnitude whose base-2 size is wanted is `n` itself
for `n >= 0` and `BitNot[n] = -n - 1` for `n < 0` (computed as `-(n + 1)`, which is `>= 0`
there); its bit length is `mpz_sizeinbase(m, 2)` — exact, because base 2 is a power of two —
and `0` is special-cased to `0`. So `BitLength[n] = Floor[Log2[n]] + 1` for `n > 0` without
ever going through floating point, `BitLength[0] = 0`, `BitLength[-1] = 0`,
`BitLength[-2] = 1`, and `BitLength[-2^k] = k`.

A concrete non-integer numeric argument (`Real`, `Rational`, `Complex`, …) emits
`BitLength::int` and leaves the call unevaluated; a symbolic argument flows through silently
(returns `NULL`); a wrong arity emits `BitLength::argx`. Both diagnostics route through the
`mth_message` funnel, so they honour `Quiet[]`/`Check[]`.

**Data structures.** Two short-lived GMP `mpz_t` values (`n` and the magnitude `m`); the
result is a machine `Integer`, since a bit length is bounded by the operand's own size and
always lands within `EXPR_INTEGER` range.

**Complexity / limits.** `O(1)` on a machine integer, `O(size of n)` on a bignum (the cost of
`mpz_sizeinbase`). `BitLength` is `Protected` and `Listable`, so it threads element-wise over
a list. The packed/`NDArray` fast path is the `int64` kernel `ndk_BitLength_ii`
(`src/ndinteger.c`) and the `Compile[]` lowering is `OP_BLEN_I`; both cover the full `int64`
range including `INT64_MIN` (`BitLength[-2^63]` is `63`), which the complement-based magnitude
computes without overflow.
