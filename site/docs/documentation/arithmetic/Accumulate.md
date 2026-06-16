# Accumulate

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Accumulate[list]
    gives a list of the successive accumulated totals of elements in
    list. The result has the same length as list.

Accumulate[list] is effectively equivalent to FoldList[Plus, list].
Accumulate works with integers, arbitrary-precision bignums, machine
doubles, and symbolic expressions, and threads naturally over rows
(so for a matrix it accumulates within columns). The head of the
input is preserved:
    Accumulate[{a, b, c, d}]    ->  {a, a + b, a + b + c, a + b + c + d}
    Accumulate[f[a, b, c, d]]   ->  f[a, a + b, a + b + c, a + b + c + d]

Accumulate[list, Method -> "CompensatedSummation"] uses Kahan
compensated summation to reduce numerical error when every element
of list is a machine number. For symbolic or mixed input the option
is ignored and the standard symbolic accumulation is returned.

Accumulate has the attribute Protected.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Accumulate[{a, b, c, d}]
Out[1]= {a, a + b, a + b + c, a + b + c + d}

In[2]:= Accumulate[{{a, b}, {c, d}, {e, f}}]
Out[2]= {{a, b}, {a + c, b + d}, {a + c + e, b + d + f}}

In[3]:= Accumulate[f[a, b, c, d]]
Out[3]= f[a, a + b, a + b + c, a + b + c + d]

In[4]:= Accumulate[{1, 2, 3, 4, 5}]
Out[4]= {1, 3, 6, 10, 15}

In[5]:= Accumulate[{1.0, 2.0, 3.0}, Method -> "CompensatedSummation"]
Out[5]= {1.0, 3.0, 6.0}
```

## Implementation notes

`builtin_accumulate` returns the list of cumulative sums (prefix sums), keeping the input expression's head. The default path folds `running = Plus[running, elem]` left to right via the evaluator, so it works on any addable elements (integers, rationals, symbolics, matrix rows). When the optional `Method -> "CompensatedSummation"` is supplied and every element is a machine number, it instead runs **Kahan compensated summation** in `double` precision (tracking a running correction term `c`), emitting `Real` partial sums. An empty list returns a copy unchanged.

- `Protected`.
- `Accumulate[list]` has the same length as `list`, and is effectively equivalent to `FoldList[Plus, list]`.
- The head of the input is preserved, so `Accumulate[f[a, b, c]]` returns `f[a, a + b, a + b + c]`.
- Threads naturally over rows via `Listable` `Plus`, so `Accumulate` of a matrix accumulates within columns.
- Works on machine integers, GMP arbitrary-precision integers, machine-precision doubles, and symbolic expressions.
- With `Method -> "CompensatedSummation"`, Kahan compensated summation is used in double precision when every element is a machine number, reducing floating-point round-off error. For symbolic or mixed input the option is silently ignored and the standard symbolic accumulation is returned.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Accumulate[{1, 2, 3, 4, 5}]
Out[1]= {1, 3, 6, 10, 15}
```

```mathematica
In[1]:= Accumulate[{a, b, c, d}]
Out[1]= {a, a + b, a + b + c, a + b + c + d}
```

```mathematica
In[1]:= Accumulate[Table[1/k, {k, 1, 5}]]
Out[1]= {1, 3/2, 11/6, 25/12, 137/60}
```

```mathematica
In[1]:= Accumulate[{{1, 2}, {3, 4}, {5, 6}}]
Out[1]= {{1, 2}, {4, 6}, {9, 12}}
```

### Notes

`Accumulate[list]` gives the running (prefix) sums, equivalent to
`FoldList[Plus, list]`, and the result has the same length as the input. It
works on exact rationals — `Accumulate[Table[1/k, {k, 1, 5}]]` returns the
partial sums of the harmonic series — as well as symbolic terms. For a matrix
(a list of rows), accumulation is performed row-wise, so column totals build up
down the matrix. The head of the input is preserved.
