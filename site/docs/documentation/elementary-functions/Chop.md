# Chop

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Chop[expr]
    replaces approximate real numbers in expr that are close to zero
    by the exact integer 0.
Chop[expr, delta]
    replaces numbers smaller in absolute magnitude than delta by 0.

The default tolerance is 10^-10. Chop walks the entire expression
tree, so small real-valued subterms inside arbitrary heads, lists,
and held forms are all chopped. Exact numbers -- integers, bigints,
rationals, and symbolic constants -- pass through untouched.

For machine complex numbers Complex[re, im] whose real and imaginary
parts are both machine reals: if only the imaginary part is below
tolerance the whole Complex wrapper is dropped and the real part
is returned; if only the real part is below tolerance the result
is Complex[0., im], preserving the machine-complex shape with a
machine zero. If both parts are below tolerance the result is the
exact integer 0.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Chop[Exp[N[Range[4] Pi I]]]
Out[1]= {-1.0, 1.0, -1.0, 1.0}

In[2]:= Chop[N[Pi] - Rationalize[N[Pi], 10^-12]] === 0
Out[2]= True

In[3]:= Chop[N[Pi] - Rationalize[N[Pi], 10^-12], 10^-14] === 0
Out[3]= False

In[4]:= Chop[10.^-12 + 2. I]
Out[4]= 0.0 + 2.0*I

In[5]:= Chop[2. + 10.^-12 I]
Out[5]= 2.0
```

## Implementation notes

- `Protected`.
- Walks the entire expression tree, so small real-valued subterms inside

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/elementary-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/elementary-functions.md)
