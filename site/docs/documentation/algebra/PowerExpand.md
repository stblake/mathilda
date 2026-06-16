# PowerExpand

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PowerExpand[expr]
    expands (a b)^c to a^c b^c and (a^b)^c to a^(b c), and expands
    Log and Arg of products and powers.
PowerExpand[expr, {x1, x2, ...}]
    expands only with respect to the listed variables.
The transformations are correct in general only when c is an integer or a and b
    are positive reals; PowerExpand otherwise disregards branch cuts.
With the Assumptions option, Assumptions->True gives a universally-correct result
    and Assumptions->assum a result valid under assum.
PowerExpand threads over lists, equations, inequalities, and logic functions.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= PowerExpand[Sqrt[x y]]
Out[1]= Sqrt[x] Sqrt[y]

In[2]:= PowerExpand[Log[(a b)^c]]
Out[2]= (Log[a] + Log[b]) c

In[3]:= PowerExpand[Sqrt[a b] + Sqrt[c d], {a, b}]
Out[3]= Sqrt[a] Sqrt[b] + Sqrt[c d]

In[4]:= PowerExpand[Sqrt[z^2], Assumptions -> z < 0]
Out[4]= -z

In[5]:= PowerExpand[Log[x y], Assumptions -> True]
Out[5]= Log[x] + Log[y] + (2*I) Pi Floor[1/2 - 1/2 (Arg[x] + Arg[y])/Pi]
```

## Implementation notes

**Algorithm.** `builtin_powerexpand` distributes powers over products and collapses nested
powers and logarithms, applying the rewrites `(a b)^c -> a^c b^c`, `(a^b)^c -> a^(b c)`,
`Log[a b] -> Log[a] + Log[b]`, `Log[a^b] -> b Log[a]`, and `Arg[a b] -> Arg[a] + Arg[b]`. Since
`Sqrt[x]` is stored as `Power[x, 1/2]` and `Log[1/z]` as `Log[Power[z, -1]]`, these cover
`Sqrt[x y]`, `Sqrt[z^2] -> z`, `Log[1/z] -> -Log[z]` with no special-casing. The transform
(`pe_rec`) is applied top-down to a fixed point, so an outermost rule fires first and the result
is reprocessed. Three modes are selected by the `Assumptions` option: **Automatic** (default,
the textbook transforms, valid for positive-real bases / integer exponents), **`-> True`** (emit
universally-correct formulas with a branch-correction term built from `Floor`/`Arg`/`Im`/`E`/
`I`/`Pi`), and **`-> assum`** (emit the True-mode formula then refine the correction terms under
the assumptions via `pe_refine`, degrading gracefully to the symbolic form where the reasoning
runs out). It threads over `List`, equations, inequalities, and logic heads, and supports the
variable-restricted form `PowerExpand[expr, {x1, …}]`.

- `Protected`.
- Applies `(a b)^c -> a^c b^c`, `(a^b)^c -> a^(b c)`, `Log[a b] -> Log[a] + Log[b]`,

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/expand_power.c`](https://github.com/stblake/mathilda/blob/main/src/expand_power.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= PowerExpand[Sqrt[a^2 b^2]]
Out[1]= a b
```

`PowerExpand` turns the logarithm of a product into a sum of logarithms:

```mathematica
In[1]:= PowerExpand[Log[a b c]]
Out[1]= Log[a] + Log[b] + Log[c]
```

It pulls exponents out of `Log` of a power:

```mathematica
In[1]:= PowerExpand[Log[x^n]]
Out[1]= n Log[x]
```

Restricting to a variable list expands only with respect to those variables (here `Sqrt[x^2] -> x`):

```mathematica
In[1]:= PowerExpand[Sqrt[x^2], {x}]
Out[1]= x
```

With `Assumptions -> True` the result is universally correct: the omitted branch-cut term reappears as an explicit `Floor[...]` correction:

```mathematica
In[1]:= PowerExpand[(a b)^n, Assumptions -> True]
Out[1]= a^n b^n E^((2*I) Pi Floor[1/2 - 1/2 (Arg[a] + Arg[b])/Pi] n)
```

### Notes

`PowerExpand[expr]` rewrites `(a b)^c` as `a^c b^c`, `(a^b)^c` as `a^(b c)`, and
expands `Log` and `Arg` of products and powers. These transformations are valid
in general only when `c` is an integer or `a, b` are positive reals; by default
`PowerExpand` disregards branch cuts. Use the variable-list form
`PowerExpand[expr, {x1, ...}]` to expand selectively, or
`Assumptions -> True` to obtain a branch-cut-correct result in which the
discarded `E^(2 Pi I Floor[...])` factor is made explicit.
