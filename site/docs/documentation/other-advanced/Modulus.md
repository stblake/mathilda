# Modulus

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

Modulus is an option for Solve.  Solve\[poly == 0, x, Modulus -\> p\] solves a single-variable polynomial equation over the finite ring Z/pZ by residue enumeration, returning {{x -\> r}, ...} with r ascending in \[0, p).  Supported for 2 \<= p \<= 100000; systems, multivariable specs, non-polynomial equations, or an out-of-range modulus leave Solve unevaluated.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

Roots in Z/5Z

```mathematica
In[1]:= Solve[x^2 + 1 == 0, x, Modulus -> 5]
Out[1]= {{x -> 2}, {x -> 3}}
```

Reducible over GF(5): (x+2)(x+3)

```mathematica
In[2]:= Factor[x^2 + 1, Modulus -> 5]
Out[2]= (2 + x) (3 + x)
```

2 is a quadratic residue mod 7

```mathematica
In[3]:= Solve[x^2 == 2, x, Modulus -> 7]
Out[3]= {{x -> 3}, {x -> 4}}
```

## Implementation notes

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

**Attributes:** none registered.

## References

- Source: [`src/solve/solve.c`](https://github.com/stblake/mathilda/blob/main/src/solve/solve.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`Modulus -> p` is an **option** that moves an operation into the finite ring `Z/pZ`
(or the field `GF(p)` for prime `p`). It is read by `Solve`, `Factor`, `PolynomialGCD`,
`PolynomialReduce`, `GroebnerBasis`, `Reduce` and related polynomial heads; `Modulus`
itself is an inert, `Protected` symbol with no value of its own.

For `Solve`, a single-variable polynomial equation is solved by residue enumeration,
returning `{{x -> r}, ...}` with `r` ascending in `[0, p)`, for `2 <= p <= 100000`;
systems, multivariable specs, non-polynomial equations, or an out-of-range modulus
leave `Solve` unevaluated. The polynomial heads generally require `p` prime (a composite
is reported unsupported rather than silently mis-factored); the default `Modulus -> 0`
is the ordinary integer ring.
