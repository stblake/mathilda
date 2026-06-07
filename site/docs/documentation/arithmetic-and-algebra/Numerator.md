# Numerator

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Numerator[expr]
    gives the numerator of expr regarded as a rational expression.
    Picks out factors of expr that do not carry a superficially negative
    exponent; constants and symbols pass through as-is.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Numerator[(x-1)(x-2)/(x-3)^2]
Out[1]= (-2 + x) (-1 + x)

In[2]:= Numerator[3/7 + I/11]
Out[2]= 33 + 7*I
```

## Implementation notes

- `Protected`, `Listable`.
- Picks out terms which do not have superficially negative exponents.
- Can be used on rational and complex numbers.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), on rational normal forms.
- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Numerator[6/8]
Out[1]= 3
```

```mathematica
In[1]:= Numerator[(x+1)/(x-1)]
Out[1]= 1 + x
```

```mathematica
In[1]:= Numerator[x^(-2)]
Out[1]= 1
```

```mathematica
In[1]:= Numerator[a/b + c/d]
Out[1]= a/b + c/d
```

### Notes

Numerator extracts the numerator of the structural rational form of its argument.
A rational constant is first reduced to lowest terms, so `Numerator[6/8] = 3`. For
symbolic quotients it returns the literal top of the `expr/expr` form, giving
`1 + x` for `(x+1)/(x-1)`. Factors with negative exponents are treated as
denominators, so `Numerator[x^(-2)] = 1`. Note that Mathilda does not auto-combine
a sum into a single fraction first: `Numerator[a/b + c/d]` returns the unevaluated
sum, since the expression is a `Plus`, not a single quotient. Apply `Together`
first if you want the combined numerator.
