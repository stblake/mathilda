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

## Examples (6)

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

## Implementation notes

**Attributes:** `NHoldAll`, `Protected`.

## References

**See also:** [Root](../../solutions-of-equations/Root/), [N](../../arithmetic/N/), [RootReduce](../../algebra/RootReduce/), [Re](../../arithmetic/Re/), [Im](../../arithmetic/Im/), [Abs](../../arithmetic/Abs/), [Round](../../arithmetic/Round/), [Less](../../comparisons/Less/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/algebra.md)
- Tests: [`tests/test_algebraicnumber.c`](https://github.com/stblake/mathilda/blob/main/tests/test_algebraicnumber.c)
- Tests: [`tests/test_algebraicnumberdenominator.c`](https://github.com/stblake/mathilda/blob/main/tests/test_algebraicnumberdenominator.c)
- Tests: [`tests/test_algebraicnumbernorm.c`](https://github.com/stblake/mathilda/blob/main/tests/test_algebraicnumbernorm.c)
- Tests: [`tests/test_algebraicnumberpolynomial.c`](https://github.com/stblake/mathilda/blob/main/tests/test_algebraicnumberpolynomial.c)
