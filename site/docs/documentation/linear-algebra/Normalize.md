# Normalize

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Normalize[v]
    gives the normalized form of a vector v (effectively v / Norm[v]).
Normalize[z]
    gives the normalized form of a scalar (incl. complex) z, namely z / Abs[z].
Normalize[expr, f]
    normalizes with respect to the norm function f, i.e. expr / f[expr].
Zero input is returned unchanged.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Normalize[{1, 5, 1}]
Out[1]= {1/3/Sqrt[3], 5/3/Sqrt[3], 1/3/Sqrt[3]}

In[2]:= Normalize[{3, 4}]
Out[2]= {3/5, 4/5}

In[3]:= Normalize[3 + 4 I]
Out[3]= 3/5 + 4/5*I

In[4]:= Normalize[{x, y}]
Out[4]= {x/Sqrt[Abs[x]^2 + Abs[y]^2], y/Sqrt[Abs[x]^2 + Abs[y]^2]}

In[5]:= Normalize[{x, y}, f]
Out[5]= {x/f[{x, y}], y/f[{x, y}]}

In[6]:= Normalize[{0, 0, 0}]
Out[6]= {0, 0, 0}
```

## Implementation notes

- `Protected` (not `Listable` — it acts on the whole vector, not element-wise).
- `Normalize[v]` is `v / Norm[v]` when `v` is a vector or tensor; the empty list `{}` and any all-zero input is returned unchanged.
- `Normalize[z]` for a scalar (possibly complex) `z` is `z / Abs[z]`, with `Normalize[0]` returning `0`.
- `Normalize[expr, f]` is `expr / f[expr]`, again with the zero short-circuit.
- The short-circuit is triggered only by a literal numeric zero (Integer, Real, BigInt, MPFR) in the evaluated norm — a symbolic norm that happens to be zero stays in the symbolic division so the user can see what they wrote.
- Wrong arity (`0` or `≥ 3` arguments) prints `Normalize::argt` and leaves the call unevaluated.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)
