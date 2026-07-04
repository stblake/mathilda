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
Out[2]= 0

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

**Algorithm.** `Derivative[n1, ..., nm]` is primarily a *tag* head: the actual
differentiation work is done inside the `D` dispatch (src/calculus/deriv.c).
`builtin_derivative` itself returns NULL, leaving `Derivative[n]` in canonical
unevaluated form; the builtin exists chiefly so attributes can be registered on
the symbol. Two pieces of real logic apply the tag:

- `derivative_of_pure_function(deriv_head, pure_fn)` differentiates
  `Derivative[n1,...,nm][Function[{t1,...,tm}, body]]` by partial-differentiating
  the body `ni` times in each slot `ti` via the shared `compute_deriv` core.
- `derivative_of_symbol(deriv_head, fsym)` reduces `Derivative[...][f]` when the
  symbol `f` carries DownValues: it mints fresh temporary slot symbols, builds
  and evaluates `f[t1,...,tm]` (triggering the DownValue rewrite), wraps the
  substituted body in a synthetic `Function`, and delegates to
  `derivative_of_pure_function`. If the call did not rewrite (no matching
  DownValue) it aborts to NULL. All `ni` must be nonnegative integers.

For unknown functions, `compute_deriv`'s chain rule emits `Derivative[...][f]`
factors, so the tag composes naturally through the rest of differentiation.

**Data structures.** `Expr*` trees; uses a static counter to generate
collision-free temporary slot-variable names (`Derivative$<id>$<k>`) without
registering them in the symbol table.

**Complexity / limits.** Only nonnegative integer derivative orders are
reduced; symbolic or negative orders stay unevaluated.

- `Protected`, `ReadProtected`.
- Acts primarily as a tag carried through the differentiation

**Attributes:** `Protected`, `ReadProtected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (Kluwer, 1992), ch. 2.
- Source: [`src/calculus/deriv.c`](https://github.com/stblake/mathilda/blob/main/src/calculus/deriv.c)
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

```mathematica
In[1]:= D[f[g[x], h[x]], x]
Out[1]= Derivative[1][g][x] Derivative[1, 0][f][g[x], h[x]] + Derivative[1][h][x] Derivative[0, 1][f][g[x], h[x]]
```

```mathematica
In[1]:= D[(f[x])^2, x]
Out[1]= 2 f[x] Derivative[1][f][x]
```

```mathematica
In[1]:= Derivative[1][#^2 + 1 &]
Out[1]= 2 #1 &
```

### Notes

`Derivative[n][f]` is the functional operator representing `f` differentiated `n` times; the surface forms `f'` and `f''` parse to `Derivative[1][f]` and `Derivative[2][f]`. It is the object `D` generates whenever it differentiates an unknown function head, which is why `D[f[x], x]` returns `Derivative[1][f][x]` and the chain rule on `f[g[x]]` yields a product of `Derivative[1]` operators. Note that `Derivative` does not auto-resolve against the known elementary table here: `Derivative[1][Sin]` and `Derivative[2][Cos]` stay in operator form rather than collapsing to `Cos` or `-Cos`. Apply the operator to an explicit argument (e.g. `Derivative[1][f][a]`) to obtain the evaluated-at form.
