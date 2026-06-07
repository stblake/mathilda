# Simplify

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Simplify[expr]
    performs a sequence of algebraic and other transformations on expr and returns the simplest form it finds.
Simplify[expr, assum]
    does simplification using assumptions assum.

Options:
  Assumptions (default $Assumptions) -- facts assumed while simplifying.
  ComplexityFunction (default: leaf count plus integer-digit count, matching Mathematica) -- ranks candidate forms; the lowest-scoring form is returned.
  TransformationFunctions (default Automatic) -- the functions applied to try to transform parts of expr. Automatic uses the built-in collection; {f1, f2, ...} uses only the fi; {Automatic, f1, ...} uses the built-in functions together with the fi.

The built-in collection tries Together, Cancel, Expand, Factor, FactorSquareFree, Apart, TrigExpand, TrigFactor, and a TrigToExp/ExpToTrig roundtrip, keeping the smallest result.
Under positivity / reality assumptions Simplify also applies Log/Power identities -- Log[a b] -> Log[a] + Log[b], (a b)^c -> a^c b^c, (a^p)^q -> a^(p q), Log[a^p] -> p Log[a] and the like -- whenever the operand-domain conditions are provable from the assumption set.
Assumptions can be equations, inequalities, domain specifications such as Element[x, Integers], or logical combinations of these. Lists of assumptions are converted to conjunctions.
Simplify automatically threads over lists, equations, inequalities, and logic functions.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Simplify[(x - 1)(x + 1)(x^2 + 1) + 1]
Out[1]= x^4

In[2]:= Simplify[3/(x + 3) + x/(x + 3)]
Out[2]= 1

In[3]:= Simplify[a x + b x + c]
Out[3]= c + (a + b) x

In[4]:= Simplify[Sin[x]^2 + Cos[x]^2]
Out[4]= 1

In[5]:= Simplify[2 Tan[x]/(1 + Tan[x]^2)]
Out[5]= Sin[2 x]

In[6]:= Simplify[(E^x - E^(-x))/Sinh[x]]
Out[6]= 2

In[7]:= Simplify[{Sin[x]^2 + Cos[x]^2, 3/(x + 3) + x/(x + 3)}]
Out[7]= {1, 1}
```

```mathematica
In[1]:= Simplify[Sqrt[2] Sqrt[3]]
Out[1]= Sqrt[6]
```

## Implementation notes

- `Protected`. **Not** `Listable`: a `List` in the assumption position is a

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/simplification.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/simplification.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Simplify[(x^2 - 1)/(x - 1)]
Out[1]= 1 + x
```

```mathematica
In[1]:= Simplify[Sin[x]^2 + Cos[x]^2]
Out[1]= 1
```

```mathematica
In[1]:= Simplify[x + x + x]
Out[1]= 3 x
```

```mathematica
In[1]:= Simplify[Sqrt[x^2], x > 0]
Out[1]= x
```

### Notes

`Simplify` tries a collection of transformations — `Together`, `Cancel`,
`Expand`, `Factor`, `Apart`, `TrigExpand`, `TrigFactor`, and a `TrigToExp` round
trip — and keeps the smallest result, so it can cancel `(x^2-1)/(x-1)` and reduce
the Pythagorean identity to `1`. A second argument supplies assumptions:
`Simplify[Sqrt[x^2], x > 0]` uses `x > 0` to drop the absolute value and return
`x`. Assumptions may be equations, inequalities, or domain statements like
`Element[x, Integers]`. `Simplify` threads automatically over lists, equations,
and logical combinations.
