---
source: src/core.c
---
**Definition.** `I` is the imaginary unit √(−1), the number with `I^2 == -1`. Unlike
the real constants, `I` is **not** a row of the `kConstants[]` table: `core_init`
(`src/core.c`) gives it an *OwnValue* mapping the symbol `I` to `Complex[0, 1]`
(`symtab_add_own_value("I", sym_I, make_complex(0,1))`), so the evaluator rewrites `I`
to that `Complex` atom on sight. Its only attribute is `Protected`, set in the
default-attributes table in `src/attr.c` (`{"I", ATTR_PROTECTED}`); it is *not*
`Constant`, so `Attributes[I]` is `{Protected}`. Interned name `SYM_I`
(`src/sym_names.c`); docstring in `src/info.c`. Because it resolves to a number,
`I^2 -> -1`, `Sqrt[-1] -> I`, `Conjugate[I] -> -I`, `Abs[I] -> 1` and `Arg[I] -> Pi/2`
all fall out of the generic complex arithmetic (`src/complex.c`, `src/times.c`,
`src/power.c`), and `D[I, x] -> 0` because the result is constant.

**Representation & numeric value.** The surface symbol `I` evaluates to the two-field
`Complex[0, 1]` atom — `FullForm[I]` is `Complex[0, 1]` and `Head[I]` is `Complex`.
`N[I]` numericalises the two exact integer parts to machine reals, giving
`0. + 1. I`. There is **no** arbitrary-precision MPFR filler for `I` and none is
needed: because the real and imaginary parts are the exact integers `0` and `1`,
`N[I, k]` leaves them machine-valued (`0.0 + 1.0 I`) rather than inflating them to
`k` digits — precision enters only when `I` is combined with an inexact quantity.

**Usage & limits.** `I` is how all complex literals are written (`a + b I` is
`Complex` once the parts are numeric); `NumericQ[I]` is `True` (whitelisted in
`is_numeric_quantity`, `src/core.c`) and it threads through every complex-aware head —
`ComplexExpand`, `Re`/`Im`/`Conjugate`/`Abs`/`Arg`, and the Euler identities
`Exp[I Pi] -> -1`, `Exp[I Pi/2] -> I`. The only subtlety is the one above: `I` is a
fixed exact number, so asking for extra precision on `I` alone changes nothing.
