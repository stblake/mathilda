# PowerExpand

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PowerExpand[expr]`**

expands (a b)^c to a^c b^c and (a^b)^c to a^(b c), and expands Log and Arg of products and powers.

**`PowerExpand[expr, {x1, x2, ...}]`**

expands only with respect to the listed variables.

<details>
<summary>Notes</summary>

The transformations are correct in general only when c is an integer or a and b are positive reals; PowerExpand otherwise disregards branch cuts. With the Assumptions option, Assumptions-\>True gives a universally-correct result and Assumptions-\>assum a result valid under assum. PowerExpand threads over lists, equations, inequalities, and logic functions.

</details>

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= PowerExpand[Sqrt[x y]]
Out[1]= Sqrt[x] Sqrt[y]

In[2]:= PowerExpand[Log[(a b)^c]]
Out[2]= (Log[a] + Log[b]) c

In[3]:= PowerExpand[Sqrt[a b] + Sqrt[c d], {a, b}]
Out[3]= Sqrt[a] Sqrt[b] + Sqrt[c d]
```

### Options (2)

```mathematica
In[4]:= PowerExpand[Sqrt[z^2], Assumptions -> z < 0]
Out[4]= -z

In[5]:= PowerExpand[Log[x y], Assumptions -> True]
Out[5]= Log[x] + Log[y] + (2*I) Pi Floor[1/2 - 1/2 (Arg[x] + Arg[y])/Pi]
```

### Applications (5)

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

## Options & behaviour

Limitations: in `Assumptions -> True` mode the nested-power correction is left
in `Im[...]` form rather than rewritten to `Arg[...]` (mathematically equal).

## Algorithm

expand_power.c -- PowerExpand.

PowerExpand distributes powers over products and collapses nested powers and logarithms:

```text
    (a b)^c   -> a^c b^c
    (a^b)^c   -> a^(b c)
    Log[a b]  -> Log[a] + Log[b]
    Log[a^b]  -> b Log[a]
    Arg[a b]  -> Arg[a] + Arg[b]
```

Sqrt is stored as Power[x, 1/2] and Log[1/z] as Log[Power[z, -1]], so the rules above cover Sqrt[x y] -> Sqrt[x] Sqrt[y], Sqrt[z^2] -> z and Log[1/z] -> -Log[z] without any special-casing.

The rewrites are applied top-down to a fixed point (mirroring Mathematica's ReplaceRepeated semantics): a rule fires at the outermost

```text
matching node and the transformed result is reprocessed.  This is why, e.g.
```

Log[(a b)^c] becomes c (Log[a] + Log[b]) rather than expanding the inner power first.

Three modes, selected by the Assumptions option:

```text
  - Automatic (default): the textbook transforms above.  Correct when the
    bases are positive reals (or the exponents integers); branch issues are
    ignored.
  - Assumptions -> True: emit the universally-correct formulas, attaching a
    branch-correction term built from Floor / Arg / Im / E / I / Pi.
  - Assumptions -> assum: emit the True-mode formula and then refine the
    correction terms under the assumptions (Arg / Im of known-sign reals,
    Floor over assumption-bounded intervals).  Faithful on the documented
    examples; where the assumptions fall outside this reasoning it degrades
    gracefully to the symbolic True-mode form rather than a wrong value.
```

PowerExpand threads over List, equations, inequalities and logic functions, and supports the variable-restricted form PowerExpand[expr, {x1, ...}].

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
  `Log[a^b] -> b Log[a]`, and `Arg[a b] -> Arg[a] + Arg[b]`. Because `Sqrt[x]` is
  `Power[x, 1/2]` and `Log[1/z]` is `Log[z^-1]`, these also give
  `Sqrt[x y] -> Sqrt[x] Sqrt[y]`, `Sqrt[z^2] -> z`, and `Log[1/z] -> -Log[z]`.
- Rules are applied top-down to a fixed point, so `Log[(a b)^c]` becomes
  `c (Log[a] + Log[b])` rather than expanding the inner power first.
- When a product raised to a non-integer (root) exponent carries a negative
  numeric coefficient together with a `Plus` factor, the sign is folded into
  that factor so the root stays real, e.g.
  `Sqrt[-4 Dt[u]^2 (-1 + u)] -> 2 Dt[u] Sqrt[1 - u]` rather than
  `2 I Dt[u] Sqrt[-1 + u]`. This is value-preserving since `(-1)^c p^c = (-p)^c`.
- `f^-1[f[x]] -> x` for the inverse-trig / inverse-hyperbolic pairs
  (`ArcTan[Tan[x]] -> x`, `ArcSin[Sin[x]] -> x`, …).
- The default `Assumptions -> Automatic` makes the textbook transforms, which
  are correct when the bases are positive reals or the exponents are integers;
  branch cuts are otherwise ignored.
- `Assumptions -> True` emits the universally-correct result, attaching a
  branch-correction term built from `Floor`, `Arg`, `Im`, `E`, `I`, and `Pi`.
- `Assumptions -> assum` emits the `True`-mode formula and then refines the
  correction terms under `assum` (resolving `Arg`/`Im` of known-sign reals and
  evaluating `Floor` over assumption-bounded intervals). This is faithful on
  the documented examples; for assumptions outside this reasoning it degrades
  to the symbolic `True`-mode form.
- `PowerExpand[expr, {x1, ...}]` expands only subexpressions that mention one
  of the listed variables.
- Threads over `List`, equations, inequalities, and logic functions.

**Attributes:** `Protected`.

## See also

[Plus](../../arithmetic/Plus/), [Floor](../../arithmetic/Floor/), [Arg](../../arithmetic/Arg/), [Im](../../arithmetic/Im/), [E](../../mathematical-constants/E/), [I](../../mathematical-constants/I/), [Pi](../../mathematical-constants/Pi/), [List](../../other-advanced/List/)

## References

- Source: [`src/expand_power.c`](https://github.com/stblake/mathilda/blob/main/src/expand_power.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
- Tests: [`tests/test_cherry_dilog.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_dilog.c)
- Tests: [`tests/test_cherry_li.c`](https://github.com/stblake/mathilda/blob/main/tests/test_cherry_li.c)
- Tests: [`tests/test_powerexpand.c`](https://github.com/stblake/mathilda/blob/main/tests/test_powerexpand.c)

## Notes & additional examples

### Notes

`PowerExpand[expr]` rewrites `(a b)^c` as `a^c b^c`, `(a^b)^c` as `a^(b c)`, and
expands `Log` and `Arg` of products and powers. These transformations are valid
in general only when `c` is an integer or `a, b` are positive reals; by default
`PowerExpand` disregards branch cuts. Use the variable-list form
`PowerExpand[expr, {x1, ...}]` to expand selectively, or
`Assumptions -> True` to obtain a branch-cut-correct result in which the
discarded `E^(2 Pi I Floor[...])` factor is made explicit.
