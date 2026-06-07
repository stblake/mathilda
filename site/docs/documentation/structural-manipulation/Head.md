# Head

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Head[expr]
    gives the head of expr.
Head[expr, h]
    wraps the result with h, i.e. returns h[Head[expr]].

For atoms, Head returns Integer, Real, BigInt, Rational, Complex,
Symbol, or String; for a compound expression f[...], Head returns f.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Head[f[x]]
Out[1]= f

In[2]:= Head[3/4]
Out[2]= Rational

In[3]:= Head[a + b, f]
Out[3]= f[Plus]
```

## Implementation notes

- For functions, returns the symbol or expression acting as the head.
- For atoms, returns the symbolic type name: `Integer`, `Real`, `Rational`, `Complex`, `Symbol`, or `String`.

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
