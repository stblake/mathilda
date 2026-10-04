# Rescale

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Rescale[x, {min, max}]`**

gives x rescaled to run from 0 to 1 over the range min to max, equivalent to (x - min)/(max - min).

**`Rescale[x, {min, max}, {ymin, ymax}]`**

gives x rescaled to run from ymin to ymax over the range min to max.

**`Rescale[list]`**

rescales each element of list to run from 0 to 1 over the range Min\[list\] to Max\[list\].

<details>
<summary>Notes</summary>

Rescale threads over a list first argument and works with exact, real, complex, and symbolic quantities.

</details>

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= Rescale[2.5, {-10, 10}]
Out[1]= 0.625

In[2]:= Rescale[-3/2, {-2, 2}]
Out[2]= 1/8

In[3]:= Rescale[3, {-9, 7}, {11, 28}]
Out[3]= 95/4

In[4]:= Rescale[{-2, 0, 2}]
Out[4]= {0, 1/2, 1}

In[5]:= Rescale[1 + 2 I, {0, 1 + I}]
Out[5]= 3/2 + 1/2*I
```

### Applications (3)

The midpoint maps to 1/2

```mathematica
In[6]:= Rescale[5, {0, 10}]
Out[6]= 1/2
```

One argument uses the data's own Min and Max

```mathematica
In[7]:= Rescale[Range[5]]
Out[7]= {0, 1/4, 1/2, 3/4, 1}
```

Rescale into an explicit target range

```mathematica
In[8]:= Rescale[{1, 2, 3, 4}, {0, 10}, {0, 100}]
Out[8]= {10, 20, 30, 40}
```

## Implementation notes

**Algorithm.** `builtin_rescale` maps a value from one range to another:
`Rescale[x, {min, max}]` is `(x - min)/(max - min)`, the three-argument form adds
the affine target range `y0 + (y1 - y0)(x - min)/(max - min)`, and the
one-argument `Rescale[list]` reduces to `Rescale[list, {Min[list], Max[list]}]`
and re-evaluates (for a lone scalar `Min == Max`, giving `Indeterminate`, as in
Mathematica). An argument count outside 1–3 raises `Rescale::argb`.

**No per-element threading.** The affine formula is built *once* with the whole
`x` substituted in, as a `Plus`/`Times`/`Power` tree, and handed to the
evaluator. Because `Plus`, `Times` and `Power` are `Listable`, they thread at
every level — including nested lists — so no explicit per-element recursion is
needed and the arithmetic is identical to it. This is the difference between an
interpreted loop and a buffer op: the old per-element `Rescale[el, range]`
evaluation cost 2.17 s over 10⁶ packed reals; handing the whole array to the
Listable heads reaches the threaded, vectorised ND kernels instead.

**Data structures / limits.** Pure `Expr`-tree construction; nothing in this
file reads an element. `Rescale` is on `pack.c`'s `AWARE` list (it only rewrites
the call) but deliberately not `int64_ok` — `Rescale[Range[10]]` is a list of
exact rationals. Exact input stays exact and symbolic elements pass through. The
range arguments must be two-element `List`s. `ATTR_NUMERICFUNCTION | ATTR_PROTECTED`.

**Attributes:** `NumericFunction`, `Protected`.

## References

**See also:** [List](../../other-advanced/List/), [Plus](../../arithmetic/Plus/), [Times](../../arithmetic/Times/), [Power](../../arithmetic/Power/)

- Source: [`src/list/rescale.c`](https://github.com/stblake/mathilda/blob/main/src/list/rescale.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compile_coverage.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_coverage.c)
- Tests: [`tests/test_packed_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_packed_list.c)
- Tests: [`tests/test_rescale.c`](https://github.com/stblake/mathilda/blob/main/tests/test_rescale.c)

## Notes & additional examples

### Notes

`Rescale[x, {min, max}]` maps `x` linearly so that `min -> 0` and `max -> 1`:
`(x - min)/(max - min)`. The three-argument form
`Rescale[x, {min, max}, {y0, y1}]` maps into an arbitrary target interval, and
the one-argument `Rescale[list]` uses `{Min[list], Max[list]}` as the source
range.

`Rescale` threads over lists at every level and keeps input exact — `Rescale[Range[5]]`
is `{0, 1/4, 1/2, 3/4, 1}`, not floats — because it rewrites its call into a
single `Plus`/`Times` expression over the whole argument and lets the `Listable`
arithmetic heads do the work, which also routes a packed array through the fast
vectorised kernels.
