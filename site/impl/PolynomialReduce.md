---
references:
  - "D. Cox, J. Little and D. O'Shea, *Ideals, Varieties, and Algorithms*, 4th ed. (Springer, 2015), ch. 2 §3 — the multivariate division algorithm."
source: src/poly/polynomialreduce.c
---
**Algorithm.** `builtin_polynomialreduce` is the multivariate-division sibling of
`GroebnerBasis` and shares its options. It gives `{{a1, ..., an}, b}` with
`a1 p1 + ... + an pn + b == poly` and `b` fully reduced — no term of `b` is
divisible by the leading term of any `pi` under the chosen `MonomialOrder`
(default `Lexicographic`). It extracts options (`MonomialOrder`,
`CoefficientDomain`, `Modulus`, `ParameterVariables`), resolves the main
variables (explicit argument, else `Variables` minus any `ParameterVariables`),
treats any remaining free symbols as coefficient-field parameters, normalises the
inputs (`pr_normalise`), and dispatches to one of three engines:

- **pure `Q`** (`pr_reduce_rational`): the exact `GBPoly` divisor engine
  `gb_divmod` (`groebner.c`), with a FLINT `fmpq_mpoly` fast path for
  `Lexicographic`;
- **rational-function field `Q(params)`** (`pr_reduce_field`): division in the
  coefficient field via `Together`/`Cancel` field arithmetic;
- **`GF(p)`** (`pr_reduce_modular`): the `gbmod.c` `gfp_divmod` engine
  (`Modulus -> p` prime; supports `Lexicographic` and
  `DegreeReverseLexicographic`).

If the `pi` form a Gröbner basis, `b` is the unique normal form.

**Data structures.** Polynomials cross into the engine as `GBPoly`
(sorted-monomial representation, built by `gb_from_expr` against the main-variable
array and a resolved `GBOrder`/weight matrix) or `GFpPoly` for the modular path;
the field engine uses `RPoly` terms with `Expr`-valued field coefficients. The
quotients and remainder come back through `gb_to_expr`/`rpoly_to_expr` into the
`{List{a1..an}, b}` result.

**Complexity / limits.** The reduction loop is worst-case exponential in the
divisor count and degrees, as multivariate division is; the FLINT path dominates
for large pure-`Q` inputs. `CoefficientDomain -> Integers`/`InexactNumbers`,
nonzero `Tolerance`, `Modulus` with parameters, and unsupported order/modulus
combinations decline with a note. A symbolic head — no packed/NDArray or
`Compile[]` path. `Protected`.
