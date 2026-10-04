---
source: src/solve/solve.c
---
**Definition.** `Modulus` is an **option symbol**, not a function. It names the integer
`p` that moves an operation from the integers (or the rationals) into the finite ring
`Z/pZ` (or the field `GF(p)` when `p` is prime). It has no builtin of its own and no
value — it is read out of an option sequence by the heads that honour it. `Modulus -> p`
is recognised by `Solve`, `Factor`, `PolynomialGCD`, `PolynomialReduce`,
`GroebnerBasis`, `Reduce` and related polynomial heads.

**Representation.** A bare `EXPR_SYMBOL` (interned `SYM_Modulus`). Each consumer scans
its trailing arguments for a `Rule[Modulus, p]` and branches on `p`:
`Solve[poly == 0, x, Modulus -> p]` solves a single-variable polynomial over `Z/pZ` by
residue enumeration, returning `{{x -> r}, ...}` with `r` ascending in `[0, p)`;
`Factor[..., Modulus -> p]` factors in `GF(p)[x]`; the Gröbner/reduction heads route
through the `gbmod.c` `gfp_divmod` engine.

**Usage & limits.** `Protected`. For `Solve`, supported for `2 <= p <= 100000`; a
system, a multivariable spec, a non-polynomial equation, or an out-of-range modulus
leaves the call unevaluated. For the polynomial heads, `Modulus -> p` generally
requires a prime `p` (a composite is reported as unsupported rather than silently
mis-factored), and the default `Modulus -> 0` means the ordinary integer ring. The
option is inert outside a head that reads it: `Modulus` on its own just evaluates to
itself.
