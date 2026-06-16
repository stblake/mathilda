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

**Algorithm.** `builtin_normalize` returns `expr / f[expr]`, where `f` defaults to `Norm` (which itself reduces to `Abs` for scalars). It builds `f[expr]`, evaluates it (`eval_and_free`), then returns `Times[expr, Power[norm_val, -1]]` evaluated. Because `Times` is `Listable`, the single reciprocal threads across every leaf of a vector / matrix / higher-rank tensor, and the same path handles scalars (including complex `z / Abs[z]`).

**Limits.** The zero short-circuit uses *exact* numeric-zero detection (`norm_is_numeric_zero`: literal Integer/Real/BigInt/MPFR `0`) — a zero vector is returned unchanged, but a symbolic norm that merely happens to vanish is left as a symbolic division so the input stays visible. Arity other than 1 or 2 emits `Normalize::argt`.

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

- Source: [`src/linalg/normalize.c`](https://github.com/stblake/mathilda/blob/main/src/linalg/normalize.c)
- Specification: [`docs/spec/builtins/linear-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/linear-algebra.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Normalize[{3, 4}]
Out[1]= {3/5, 4/5}
```

Normalizing an integer vector whose norm is rational stays exact:

```mathematica
In[1]:= Normalize[{1, 2, 2}]
Out[1]= {1/3, 2/3, 2/3}
```

A symbolic vector normalizes to its general unit-vector formula:

```mathematica
In[1]:= Normalize[{a, b}]
Out[1]= {a/Sqrt[Abs[a]^2 + Abs[b]^2], b/Sqrt[Abs[a]^2 + Abs[b]^2]}
```

A complex scalar is divided by its modulus, giving a unit-modulus phase:

```mathematica
In[1]:= Normalize[3 + 4 I]
Out[1]= 3/5 + 4/5*I
```

`Normalize[expr, f]` normalizes with respect to any norm function — here `Total`
turns counts into a probability distribution:

```mathematica
In[1]:= Normalize[{1, 1}, Total]
Out[1]= {1/2, 1/2}
```

### Notes

`Normalize[v]` returns `v / Norm[v]`, the unit vector in the direction of `v`;
for a scalar (including a complex number) it returns `z / Abs[z]`. The two-argument
form `Normalize[expr, f]` divides by `f[expr]` instead, so any norm or aggregating
function may be supplied. Zero input is returned unchanged.
