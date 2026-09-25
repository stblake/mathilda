# AlgebraicIntegerQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AlgebraicIntegerQ[x]`**

gives True if x is an algebraic integer (a root of a monic polynomial with integer coefficients) and False otherwise.  Rational integers are algebraic integers; non-integer rationals and non-algebraic quantities are not.  Requires FLINT.  Attribute: Protected.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= AlgebraicIntegerQ[Sqrt[2]]
Out[1]= True

In[2]:= AlgebraicIntegerQ[1/2]
Out[2]= False

In[3]:= AlgebraicIntegerQ[(1 + Sqrt[5])/2]
Out[3]= True

In[4]:= AlgebraicIntegerQ[2^(1/3)/2]
Out[4]= False
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Pi](../../mathematical-constants/Pi/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/algebra.md)
- Tests: [`tests/test_algebraicnumberdenominator.c`](https://github.com/stblake/mathilda/blob/main/tests/test_algebraicnumberdenominator.c)
- Tests: [`tests/test_numberfieldintegralbasis.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numberfieldintegralbasis.c)
