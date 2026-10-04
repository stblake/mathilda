# AlgebraicNumber

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AlgebraicNumber[theta, {c0, c1, ..., cn}]`**

represents the algebraic number c0 + c1 theta + ... + cn theta^n in the field Q(theta).

<details>
<summary>Notes</summary>

The generator theta may be given as a radical, a Root object, or another AlgebraicNumber; the coefficients ci must be integers or rationals. The object is automatically reduced so that theta is an algebraic integer and the coefficient list has length equal to the degree of the minimal polynomial of theta.  AlgebraicNumber objects in the same field are combined by arithmetic; those representing a rational number reduce to explicit rational form.  They are treated as numeric quantities: N gives their value to any precision and RootReduce converts them to Root objects.  Attributes: NHoldAll, Protected.

</details>

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= AlgebraicNumber[Root[#^3 + # + 1 &, 3], {1, 2, 1}]
Out[1]= AlgebraicNumber[Root[1 + #1 + #1^3 &, 3], {1, 2, 1}]

In[2]:= 1 + %^2
Out[2]= 1 + Out[-1]^2
```

Generator -> algebraic integer

```mathematica
In[3]:= AlgebraicNumber[(1 + I)/2, {1, 3}]
Out[3]= AlgebraicNumber[1 + I, {1, 3/2}]
```

Fold over the minimal polynomial

```mathematica
In[4]:= AlgebraicNumber[3^(1/5), {1, 2, 1, 3, 3, 1}]
Out[4]= AlgebraicNumber[Root[-3 + #1^5 &, 1], {4, 2, 1, 3, 3}]
```

```mathematica
In[5]:= AlgebraicNumber[Sqrt[2], {1, 1/2}] + AlgebraicNumber[Sqrt[2], {1, 2}]
Out[5]= AlgebraicNumber[Sqrt[2], {2, 5/2}]

In[6]:= N[AlgebraicNumber[Sqrt[2] I, {1, -1}], 50] 1.4142135623730950488016887242096980785696718753769 I
Out[6]= 2.0 + 1.4142135623730950488016887242096980785696718753769*I
```

### Applications (5)

Represents 1 + 2 Sqrt[2] in Q(Sqrt[2])

```mathematica
In[7]:= AlgebraicNumber[Sqrt[2], {1, 2}]
Out[7]= AlgebraicNumber[Sqrt[2], {1, 2}]
```

Only the constant term: a rational, reduced out

```mathematica
In[8]:= AlgebraicNumber[Sqrt[2], {3, 0}]
Out[8]= 3
```

Treated as a numeric quantity

```mathematica
In[9]:= N[AlgebraicNumber[Sqrt[2], {1, 2}]]
Out[9]= 3.82843
```

Arithmetic in the same field

```mathematica
In[10]:= AlgebraicNumber[Sqrt[2], {1, 2}] + AlgebraicNumber[Sqrt[2], {3, 4}]
Out[10]= AlgebraicNumber[Sqrt[2], {4, 6}]
```

Recover the defining polynomial

```mathematica
In[11]:= AlgebraicNumberPolynomial[AlgebraicNumber[Sqrt[2], {1, 2}], x]
Out[11]= 1 + 2 x
```

## Implementation notes

**Algorithm.** `builtin_algebraicnumber` validates `AlgebraicNumber[theta,
{c0..cn}]` (arity 2) and delegates all field computation to
`flint_qqbar_algebraic_number`; the builtin itself only canonicalises and guards
against re-evaluation churn. The engine:

1. Builds the **algebraic-integer generator** `phi = lc·alpha` of `Q(theta)`,
   where `alpha = to_qqbar(theta)` and `lc` is the positive leading coefficient
   of `alpha`'s primitive integer minimal polynomial (`algint_generator`); `phi`
   is monic of degree `n`.
2. Forms `p(x) = Σ (c_i / lc^i) x^i` over `Q` and reduces it modulo `phi`'s monic
   minimal polynomial `M` (`fmpq_poly_rem`), giving a representative of degree
   `< n` in the power basis of `phi`.
3. Renders it with `poly_to_algnum`: a plain `Integer`/`Rational` when the value
   is rational (only the constant coefficient survives), otherwise the reduced
   `AlgebraicNumber[g, {d0..d_{n-1}}]` with the coefficient list padded to length
   `n` (the degree of the minimal polynomial).

A **fixpoint guard** (`expr_eq(cand, res)` → return `NULL`) leaves an input that
is already canonical untouched, so the evaluator does not loop — the
canonicalisation is idempotent. A **fast path** (via the `GENFD` generator-field
cache) short-circuits when `theta` is already the canonical integer generator
and the coefficients are already rational and reduced, returning the padded form
directly without re-running `to_qqbar`.

**Data structures.** The `AlgebraicNumber[theta, {c0..cn}]` representation
itself; FLINT `qqbar_t` for `alpha`/`phi`; `fmpq_poly` for `p` and `M`; and the
per-generator `GENFD` cache keyed by the generator expression (the minimal
polynomial is a pure function of it, so the cache never goes stale). The head
carries `NHoldAll` so `N` reaches the dedicated AlgebraicNumber branch rather
than threading into the generator and coefficient list.

**Complexity / limits.** Dominated by `to_qqbar` of the generator and the
`fmpq_poly` reduction modulo `M`; bounded by the degree cap
`QQBAR_DEGREE_CAP = 120`. Declines (stays unevaluated) for a non-algebraic
generator, malformed (non-integer/rational) coefficients, a degree-cap overflow,
or FLINT compiled out.

**Attributes:** `NHoldAll`, `Protected`.

## References

**See also:** [Root](../../solutions-of-equations/Root/), [N](../../arithmetic/N/), [RootReduce](../../algebra/RootReduce/), [Re](../../arithmetic/Re/), [Im](../../arithmetic/Im/), [Abs](../../arithmetic/Abs/), [Round](../../arithmetic/Round/), [Less](../../comparisons/Less/)

- H. Cohen, *A Course in Computational Algebraic Number Theory*, GTM 138 (Springer, 1993), §4.2–4.3 (representation of number-field elements in a power basis).
- The FLINT library (https://flintlib.org), `qqbar` module (exact algebraic numbers) and `nf`/`nf_elem` modules (number-field element arithmetic).
- Source: [`src/poly/algebraicnumber.c`](https://github.com/stblake/mathilda/blob/main/src/poly/algebraicnumber.c)
- Specification: [`docs/spec/builtins/algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/algebra.md)
- Tests: [`tests/test_algebraicnumber.c`](https://github.com/stblake/mathilda/blob/main/tests/test_algebraicnumber.c)
- Tests: [`tests/test_algebraicnumberdenominator.c`](https://github.com/stblake/mathilda/blob/main/tests/test_algebraicnumberdenominator.c)
- Tests: [`tests/test_algebraicnumbernorm.c`](https://github.com/stblake/mathilda/blob/main/tests/test_algebraicnumbernorm.c)
- Tests: [`tests/test_algebraicnumberpolynomial.c`](https://github.com/stblake/mathilda/blob/main/tests/test_algebraicnumberpolynomial.c)

## Notes & additional examples

### Notes

`AlgebraicNumber[theta, {c0, c1, …, cn}]` denotes `c0 + c1 theta + … + cn theta^n`
in the field `Q(theta)`. The object is automatically reduced so that `theta` is
an algebraic integer and the coefficient list has length equal to the degree of
`theta`'s minimal polynomial; an object representing a rational number collapses
to explicit rational form.

Objects in the same field combine under `+`, `*` and integer powers. The head
carries `NHoldAll`, so `N` reaches the dedicated numeric branch (giving the
value to any precision) rather than numericalising the stored generator and
coefficients in place. Requires FLINT for the canonicalisation.
