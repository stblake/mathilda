# PossibleZeroQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PossibleZeroQ[expr] gives True if symbolic and numerical methods suggest that expr has value zero, and False otherwise.
The general problem of deciding whether an expression is zero is undecidable; PossibleZeroQ is a quick but not always accurate test.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= PossibleZeroQ[E^(I Pi/4) - (-1)^(1/4)]
Out[1]= True

In[2]:= PossibleZeroQ[(x + 1)(x - 1) - x^2 + 1]
Out[2]= True

In[3]:= PossibleZeroQ[(E + Pi)^2 - E^2 - Pi^2 - 2 E Pi]
Out[3]= True

In[4]:= PossibleZeroQ[E^Pi - Pi^E]
Out[4]= False

In[5]:= PossibleZeroQ[2^(2 I) - 2^(-2 I) - 2 I Sin[Log[4]]]
Out[5]= True

In[6]:= PossibleZeroQ[Sqrt[x^2] - x]
Out[6]= False

In[7]:= PossibleZeroQ[Sin[x]^2 + Cos[x]^2 - 1]
Out[7]= True
```

## Implementation notes

**Algorithm.** `builtin_possible_zero_q` (`src/zero_test.c`) calls `zero_test_decide`, a staged hybrid symbolic-numeric pipeline that early-exits on the first definite verdict:

- *Stage 0 — structural:* O(1) shortcuts for literal `Integer`/`Real`/`BigInt`/`MPFR` zero, `Complex[0,0]`, lists of zeros, and unbound symbols (`decide_structural`).
- *Stage 1 — rational normalisation:* `Together ∘ Cancel` plus `Expand`, then a polynomial zero test, deciding every identity in `Q(x_1,...,x_n)` (`decide_rational`). A `True` here is trusted; a `False` is not trusted alone.
- *Stage 2 — numeric:* for symbol-free inputs, numericalize at machine precision and compare `|z|` against an IEEE catastrophic-cancellation threshold, climbing a precision ladder (53 -> 200 -> 500 -> 1000 bits) while ambiguous; a true zero shrinks geometrically across precisions (`decide_numeric`).
- *Stage 3 — Schwartz–Zippel:* for inputs with free symbols, substitute each free symbol with a random rational drawn from `Q[i]` (to probe branch cuts) and require several independent Stage-2 confirmations (`decide_schwartz_zippel`).

**Result mapping.** A definite `False` returns `False`; both `True` and `UNKNOWN` return `True`, matching Mathematica's documented "assume zero when uncertain" behaviour (the accompanying `PossibleZeroQ::ztest1` message is not emitted). The symbol is `Listable`.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- J. T. Schwartz, "Fast probabilistic algorithms for verification of polynomial identities", JACM 27 (1980).
- R. Zippel, "Probabilistic algorithms for sparse polynomials", EUROSAM 1979.
- Source: [`src/zero_test.c`](https://github.com/stblake/mathilda/blob/main/src/zero_test.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= PossibleZeroQ[(x - 1) (x + 1) - (x^2 - 1)]
Out[1]= True
```

```mathematica
In[1]:= PossibleZeroQ[x^2 + 1]
Out[1]= False
```

```mathematica
In[1]:= PossibleZeroQ[Sin[x]^2 + Cos[x]^2 - 1]
Out[1]= True
```

```mathematica
In[1]:= PossibleZeroQ[Sqrt[2] + Sqrt[3] - Sqrt[5 + 2 Sqrt[6]]]
Out[1]= True
```

```mathematica
In[1]:= PossibleZeroQ[Log[2] + Log[3] - Log[6]]
Out[1]= True
```

### Notes

`PossibleZeroQ[expr]` uses combined symbolic and numerical heuristics to decide
whether `expr` is identically zero. It sees through polynomial cancellation
(`(x-1)(x+1) - (x^2-1) = 0`), the Pythagorean identity
`Sin[x]^2 + Cos[x]^2 - 1`, the nested-radical denesting
`Sqrt[2] + Sqrt[3] = Sqrt[5 + 2 Sqrt[6]]`, and the logarithm law
`Log[2] + Log[3] = Log[6]`. A nonzero expression like `x^2 + 1` returns
`False`. Because deciding whether a closed-form expression is exactly zero is
undecidable in general, `PossibleZeroQ` is a fast but not infallible test:
`True` strongly suggests a zero and `False` rules one out for the cases it can
analyse.
