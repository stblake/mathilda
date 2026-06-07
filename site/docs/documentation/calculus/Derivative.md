# Derivative

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
f' represents the derivative of a function f of one argument.
Derivative[n1, n2, ...][f] is the general form, representing a function
obtained from f by differentiating n1 times with respect to the first
argument, n2 times with respect to the second argument, and so on.

f' is equivalent to Derivative[1][f]; f'' evaluates to Derivative[2][f].
Derivative is a functional operator acting on functions to give derivative
functions. Derivative is generated when D is applied to functions whose
derivatives the system does not know.

Mathilda attempts to convert Derivative[n1,...,nm][f] to a pure function.
When f is a symbol carrying DownValues, the evaluator rewrites the head
as Function[{t1,...,tm}, f[t1,...,tm]] with the rule expanded into the
body, then differentiates that pure function. If no DownValue matches,
the original Derivative form is returned.

Attributes: Protected, ReadProtected.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Derivative[2][f][x]
Out[1]= Derivative[2][f][x]

In[2]:= D[%, x]
Out[2]= Derivative[3][f][x]

In[3]:= D[f[x, y], y]
Out[3]= Derivative[0, 1][f][x, y]

In[4]:= f'[x]
Out[4]= 18 x^2 + 5 x^4

In[5]:= f'[5]
Out[5]= 3575

In[6]:= Derivative[1, 1][g][a, b]
Out[6]= 6 a b^2
```

## Implementation notes

- `Protected`, `ReadProtected`.
- Acts primarily as a tag carried through the differentiation

**Attributes:** `Protected`, `ReadProtected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (Kluwer, 1992), ch. 2.
- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= f'[x]
Out[1]= Derivative[1][f][x]
```

```mathematica
In[1]:= D[f[x], x]
Out[1]= Derivative[1][f][x]
```

```mathematica
In[1]:= Derivative[2][Cos]
Out[1]= Derivative[2][Cos]
```

```mathematica
In[1]:= D[f[g[x]], x]
Out[1]= Derivative[1][g][x] Derivative[1][f][g[x]]
```

### Notes

`Derivative[n][f]` is the functional operator representing `f` differentiated `n` times; the surface forms `f'` and `f''` parse to `Derivative[1][f]` and `Derivative[2][f]`. It is the object `D` generates whenever it differentiates an unknown function head, which is why `D[f[x], x]` returns `Derivative[1][f][x]` and the chain rule on `f[g[x]]` yields a product of `Derivative[1]` operators. Note that `Derivative` does not auto-resolve against the known elementary table here: `Derivative[1][Sin]` and `Derivative[2][Cos]` stay in operator form rather than collapsing to `Cos` or `-Cos`. Apply the operator to an explicit argument (e.g. `Derivative[1][f][a]`) to obtain the evaluated-at form.
