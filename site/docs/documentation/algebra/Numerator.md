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

`builtin_numerator` calls the shared `extract_num_den` splitter and returns the numerator (freeing the denominator); `Denominator` is the mirror. `extract_num_den` handles literal rationals (`n/d`), complex numbers (clearing the common denominator of the real/imaginary parts), `Power[b, e]`/`Exp[e]` (a negative integer or rational exponent — or a `Plus` exponent with superficially-negative terms — moves the factor into the denominator), and `Times[...]` (recurse on each factor, partition the results into numerator and denominator products). Anything else is its own numerator over denominator 1. It does *not* combine a `Plus` over a common denominator — that is `Together`'s job — so `Numerator[a/b + c/d]` returns the input's surface numerator, not the combined one. `Numerator` carries `ATTR_LISTABLE | ATTR_PROTECTED`.

- `Protected`, `Listable`.
- Picks out terms which do not have superficially negative exponents.
- Can be used on rational and complex numbers.

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), on rational normal forms.
- Source: [`src/rat.c`](https://github.com/stblake/mathilda/blob/main/src/rat.c)
- Specification: [`docs/spec/builtins/algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/algebra.md)

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
