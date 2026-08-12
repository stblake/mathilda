# PossibleZeroQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PossibleZeroQ[expr] gives True if symbolic and numerical methods suggest that expr has value zero, and False otherwise.`**

**`PossibleZeroQ[expr, Assumptions -> assum] tests under the assumptions assum; it also respects an ambient Assuming[] / $Assumptions scope. Assumptions restrict the numeric sampler to the assumed region (integer, real or complex domain, sign, range, and Re/Im-part constraints), so identities that hold only there are recognised.`**

<details>
<summary>Notes</summary>

The general problem of deciding whether an expression is zero is undecidable; PossibleZeroQ is a quick but not always accurate test.

</details>

## Examples (13)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

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

### Options (1)

```mathematica
In[8]:= PossibleZeroQ[Sqrt[x^2] - x, Assumptions -> x >= 0]
Out[8]= True
```

### Applications (5)

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

## Algorithm

zero_test.c — PossibleZeroQ: hybrid symbolic-numeric zero recognition.

Pipeline (early exit at any stage that yields a definite verdict):

```text
  Stage 0 — O(1) structural shortcuts: literal Integer/Real/BigInt/MPFR
            zero, Complex[0, 0], List of zeros, unbound symbol, …

  Stage 1 — Rational normalisation via Together ∘ Cancel + Expand,
            then is_zero_poly. Decides every identity in Q(x_1,…,x_n).

  Stage 2 — Closed-form numeric path: numericalize at machine precision,
            compare |z| against the IEEE catastrophic-cancellation
            ambiguity threshold scale * 2^(-p/2 + 4). If ambiguous, bump
            precision (machine → 200 → 500 → 1000 bits) and retry.
            A non-zero result stays roughly constant across precisions;
            a true zero shrinks geometrically. Surviving the full ladder
            implies "True".

  Stage 3 — Schwartz–Zippel. For inputs with free symbols, substitute
            each free symbol with a random REAL rational of moderate
            magnitude, recurse into Stage 2, and require independent
            confirmations. Sampling is real-line only: an analytic identity
            holding on a real interval holds on a complex neighbourhood
            (identity theorem), so real points confirm it, while complex
            samples needlessly cross branch cuts (Log/ArcTan/Sqrt) and blow
            up special functions (Gamma), manufacturing false negatives.
            The draw stream is seeded deterministically from the input's
            structural hash, so the verdict is a pure function of the input
            (no run-to-run flakiness) and the user's RNG stream is untouched.
```

See ZERO_RECOGNISE_PLAN.md for design notes and references.

## Implementation notes

**Algorithm.** `builtin_possible_zero_q` (`src/zero_test.c`) calls `zero_test_decide`, a staged hybrid symbolic-numeric pipeline that early-exits on the first definite verdict:

- *Stage 0 — structural:* O(1) shortcuts for literal `Integer`/`Real`/`BigInt`/`MPFR` zero, `Complex[0,0]`, lists of zeros, and unbound symbols (`decide_structural`).
- *Stage 1 — rational normalisation:* `Together ∘ Cancel` plus `Expand`, then a polynomial zero test, deciding every identity in `Q(x_1,...,x_n)` (`decide_rational`). A `True` here is trusted; a `False` is not trusted alone.
- *Stage 2 — numeric:* for symbol-free inputs, numericalize at machine precision and compare `|z|` against an IEEE catastrophic-cancellation threshold, climbing a precision ladder (53 -> 200 -> 500 -> 1000 bits) while ambiguous; a true zero shrinks geometrically across precisions (`decide_numeric`).
- *Stage 3 — Schwartz–Zippel:* for inputs with free symbols, substitute each free symbol with a random rational drawn from `Q[i]` (to probe branch cuts) and require several independent Stage-2 confirmations (`decide_schwartz_zippel`).

**Result mapping.** A definite `False` returns `False`; both `True` and `UNKNOWN` return `True`, following the documented "assume zero when uncertain" behaviour (the accompanying `PossibleZeroQ::ztest1` message is not emitted). The symbol is `Listable`.

**Attributes:** `Protected`.

## See also

[Pi](../../mathematical-constants/Pi/), [Together](../../algebra/Together/), [Cancel](../../algebra/Cancel/), [Expand](../../algebra/Expand/), [Plus](../../arithmetic/Plus/), [Times](../../arithmetic/Times/), [Power](../../arithmetic/Power/), [CoefficientList](../../algebra/CoefficientList/)

## References

- J. T. Schwartz, "Fast probabilistic algorithms for verification of polynomial identities", JACM 27 (1980).
- R. Zippel, "Probabilistic algorithms for sparse polynomials", EUROSAM 1979.
- Source: [`src/zero_test.c`](https://github.com/stblake/mathilda/blob/main/src/zero_test.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_divisorsigma.c`](https://github.com/stblake/mathilda/blob/main/tests/test_divisorsigma.c)
- Tests: [`tests/test_extension_auto_builtins.c`](https://github.com/stblake/mathilda/blob/main/tests/test_extension_auto_builtins.c)
- Tests: [`tests/test_findintegernullvector.c`](https://github.com/stblake/mathilda/blob/main/tests/test_findintegernullvector.c)
- Tests: [`tests/test_flint_bridge.c`](https://github.com/stblake/mathilda/blob/main/tests/test_flint_bridge.c)

## Notes & additional examples

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
