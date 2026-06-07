# Cancel

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Cancel[expr] cancels out common factors in the numerator and denominator of expr.
Option Extension -> alpha cancels factors over Q(alpha) (e.g. simplifies
(x^2 - 2)/(x - Sqrt[2]) to x + Sqrt[2] when Extension -> Sqrt[2]).
Default Extension -> None treats algebraic numbers as opaque.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Cancel[(x^2 - 1) / (x - 1)]
Out[1]= 1 + x

In[2]:= Cancel[(x - y)/(x^2 - y^2) + (x^3 - 27)/(x^2 - 9)]
Out[2]= (9 + 3 x + x^2)/(3 + x) + 1/(x + y)

In[3]:= Cancel[(y - 1)/(Sqrt[y] - 1)]
Out[3]= 1 + Sqrt[y]

In[4]:= Cancel[(y - 1)/(y^(1/3) - 1)]
Out[4]= 1 + y^(1/3) + y^(2/3)

In[5]:= Cancel[1/(y^(2/3) - 1/y^(1/3))]
Out[5]= y^(1/3)/(-1 + y)

In[6]:= Cancel[(x^2 - 2)/(x - Sqrt[2]), Extension -> Sqrt[2]]
Out[6]= Sqrt[2] + x

In[7]:= Cancel[(x^3 - 2)/(x - 2^(1/3)), Extension -> 2^(1/3)]
Out[7]= 2^(2/3) + 2^(1/3) x + x^2

In[8]:= Cancel[Sqrt[2]/(Sqrt[2] + Sqrt[2] x^4)]
Out[8]= 1/(1 + x^4)
```

## Implementation notes

- `Protected`, `Listable`.
- Threads over equations, inequalities, logic functions, and sums dynamically.
- Evaluates greatest common divisors via polynomial GCD derivations avoiding extraneous expansions.
- Handles a single symbolic base appearing with rational fractional exponents (e.g. `Sqrt[y]`, `y^(1/3)`) by treating it as an algebraic generator: substitutes `y -> g^m` where `m` is the LCM of denominators, runs the polynomial cancellation in `g`, then substitutes back.
- The algebraic-generator pass runs `Together` on the substituted form (not just GCD-cancellation), so inputs whose `g`-substituted denominator is a Plus of terms with different `g`-denominators (e.g. `1/(g^2 - 1/g)` from `1/(y^(2/3) - 1/y^(1/3))`) are handled correctly.
- Extracts algebraic-constant atoms (`Sqrt[2]`, `Sqrt[3]`, `CubeRoot[5]`, `2^(2/3)`, ...) that appear in every summand of *both* numerator *and* denominator and divides them out before the polynomial GCD step. Closes a long-standing gap where `PolynomialGCD[Sqrt[2], Sqrt[2] + Sqrt[2] x^4]` returns 1 (the integer-content recursion treats `Sqrt[2]` as having content 1), so cancellations whose only shared factor was an algebraic constant survived as-is. The pass is intentionally narrow — only `Power[integer, rational/non-integer]` factors are eligible — to avoid disturbing the rational-function intermediates that the integration dispatcher pattern-matches against.
- **Option `Extension -> alpha`** (Phase 0 of the Integrate plan) cancels common factors over `Q(alpha)` instead of `Q`. Implementation: lifts numerator and denominator into `Q(alpha)[x]` via the QAUPoly machinery, runs `qaupoly_gcd`, divides both sides by `g`, and re-renders. Works for single-fraction inputs; `Plus` inputs (sums of fractions) currently fall back to the no-extension path because `PolynomialQuotient` does not yet accept `Extension` (Phase 0.5 follow-up).

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- von zur Gathen & Gerhard, "Modern Computer Algebra", on polynomial GCD computation.
- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), on rational function simplification.
- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Cancel[(x^2 - 1)/(x - 1)]
Out[1]= 1 + x
```

```mathematica
In[1]:= Cancel[(x^2 + 2 x + 1)/(x + 1)]
Out[1]= 1 + x
```

```mathematica
In[1]:= Cancel[(x^2 - 1)/(x^2 - 2 x + 1)]
Out[1]= (1 + x)/(-1 + x)
```

### Notes

Cancel removes the common polynomial factors between numerator and denominator by
dividing out their GCD, without performing a partial-fraction split. It reduces
`(x^2-1)/(x-1)` to `1 + x` and recognises perfect-square numerators such as
`(x^2+2x+1)/(x+1)`. When a common factor remains after cancellation, the result
stays a reduced quotient: `(x^2-1)/(x^2-2x+1)` becomes `(1+x)/(-1+x)` because both
share the factor `(x-1)` but the leftover `(x+1)/(x-1)` is already in lowest
terms. Cancel does not expand the surviving factors back out.
