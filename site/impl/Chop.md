---
source: src/core.c
---
**Algorithm.** `builtin_chop` reads an optional `delta` (default `1.0e-10`) via
`chop_extract_delta` (accepting Integer, BigInt, Real, MPFR, and `Rational[n,d]`,
returning the absolute value; a non-coercible argument makes the whole call stay
unevaluated). It then recurses with `chop_recursive`. For a "machine real" leaf
(detected by `chop_is_machine_real`), it replaces the value by the exact integer
`0` when `fabs(value) < delta`. `Complex` leaves are handled specially: when both
parts are machine reals it independently chops each, returning `0` if both are
small, the real part alone if only the imaginary part chops, or
`Complex[0., im]` (built with `make_complex` to avoid re-evaluation) if only the
real part chops. Non-machine `Complex` nodes recurse into each part and let
`builtin_complex` collapse `Complex[r, 0]` on the next evaluator pass. Function
nodes are rebuilt with each argument chopped; all other atoms (integers, bigints,
symbols, strings, rational components) are copied unchanged.

**Data structures.** Pure `Expr*` tree rewrite — no auxiliary structure. The
comparison is `double`-based (`chop_to_double`), so only inexact numeric leaves
are ever affected; exact integers and rationals are never chopped.
