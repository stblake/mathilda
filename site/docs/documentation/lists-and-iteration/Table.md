# Table

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Table[expr, n]`**

generates a list of n copies of expr.

**`Table[expr, {i, imax}]`**

generates a list of the values of expr with i running from 1 to imax.

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= Table[i^2, {i, 4}]
Out[1]= {1, 4, 9, 16}
```

### Applications (7)

```mathematica
In[2]:= Table[i^2, {i, 1, 5}]
Out[2]= {1, 4, 9, 16, 25}

In[3]:= Table[i + j, {i, 1, 2}, {j, 1, 3}]
Out[3]= {{2, 3, 4}, {3, 4, 5}}

In[4]:= Table[x, 4]
Out[4]= {x, x, x, x}

In[5]:= Table[i, {i, 0, 1, 1/2}]
Out[5]= {0, 1/2, 1}

In[6]:= Table[Fibonacci[n], {n, 1, 12}]
Out[6]= {1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144}

In[7]:= Table[Sum[1/k, {k, 1, n}], {n, 1, 5}]
Out[7]= {1, 3/2, 11/6, 25/12, 137/60}

In[8]:= Table[If[i == j, 1, 0], {i, 1, 3}, {j, 1, 3}]
Out[8]= {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}
```

## Performance

Against other systems, from the benchmark suite (same input, results cross-checked for agreement):

| case | Mathilda | Wolfram | Python |
|---|---:|---:|---:|
| Interpolation evaluate, 20000 points | 18.9 s | 18.1 s | 0.073 s |
| Interpolation over 10^5 array | 1.26 s | 4.57 s | 4.44 s |
| Fit degree-8 polynomial, 5000 points | 0.999 s | 1.29 s | 0.327 s |
| LeastSquares 2000x20 | 0.496 s | 0.612 s | 0.293 s |
| Fit trig basis, 5000 points | 0.224 s | 0.409 s | 0.092 s |
| Fit quadratic, 5000 points | 0.223 s | 0.413 s | 0.141 s |

## Implementation notes

**Algorithm.** `builtin_table` (in `src/list.c`) implements iterator-driven list construction. `Table` carries `HoldAll | Protected` (set in `src/attr.c`), so the body and iterator specs arrive unevaluated. Multi-iterator forms `Table[expr, spec1, ..., specN]` are handled by *rewriting*: the handler peels off the outermost iterator, wraps the remaining iterators in a fresh `Table[Table[expr, ..., specN], spec1]`, and recurses through the evaluator — so nesting depth equals iterator count and the rightmost spec varies fastest.

For a single iterator the spec is parsed by the shared `iter_spec_parse` (`src/iter.c`) into an `IterSpec` discriminated by `kind`: `ITER_KIND_COUNT` (`{n}` or bare `n` — repeat the body), `ITER_KIND_LIST` (`{i, {v1, v2, ...}}` — iterate over explicit values), or `ITER_KIND_RANGE` (`{i, imax}`, `{i, imin, imax}`, `{i, imin, imax, di}`). Range bounds are resolved to doubles via `iter_spec_resolve_numeric` (with `allow_inf=false`, since `Table` never iterates to `Infinity`).

**Localization.** The iteration variable is dynamically scoped, not renamed. `iter_spec_shadow(var)` saves and clears the symbol's `own_values` chain; each step calls `symtab_add_own_value` to bind the current index, `evaluate`s the body, then the next index is computed symbolically (`evaluate(Plus[curr, di])`) so exact integers/rationals are preserved while a parallel `double` accumulator drives loop termination (with a `1e-14` tolerance and a 1,000,000-step safety cap). `iter_spec_restore` reinstalls the saved `own_values` afterward, so the variable's outer value is untouched.

**Data structures.** Results accumulate in a geometrically grown `Expr**` buffer, finally wrapped as `List[...]`.

- `HoldAll`: `expr` is evaluated once for each step.
- Supports nested iterators to create matrices.
- **Machine-real iterators take a compiled fast path.** When the iterator is
  inexact (`{x, 0., 1., 0.01}`) and `expr` is in the compilable subset, the body
  is compiled once and run per point over machine numbers instead of through the
  evaluator. Exact (Integer / BigInt / Rational) iterators are untouched and
  stay bit-for-bit identical. Each element keeps its own type, so an
  integer-valued body still yields Integers. A point where the compiled program
  cannot produce a value (non-finite, or complex where a real was promised)
  falls back to the interpreter for that element alone.
  Nested iterators compile too: the outer variables are folded in as the
  constants they currently hold. Calls to a `CompiledFunction` inside the body
  are inlined rather than dispatched per point.
- **Large machine-number results pack.** A result at or above the packing
  threshold that is a rectangular block of uniformly exact or uniformly inexact
  machine numbers is returned as a [packed list](../packed-arrays/index.md): an ordinary `List` held as a
  dense buffer, distinguishable only by `NDArrayQ`. A machine-real compiled body writes the
  buffer directly, without building the elements first; every other branch is
  offered for packing once built. `Table[i j, {i, 300}, {j, 300}]` is one rank-2
  packed array, not a list of packed rows.
- **Exact termination and a real element backstop.** An integer range terminates
  on an exact comparison of the running value against the bound, so it is correct
  anywhere in the `int64` range — `Table[i, {i, 9223372036854775805, 9223372036854775807}]`
  is the three-element list, not a runaway. A range longer than
  `100000000` elements (matching `Sum`/`Product`) returns `Table[…]` unevaluated
  rather than silently truncating; an exact-integer range is rejected up front,
  before any element is allocated.

**Attributes:** `HoldAll`, `Protected`.

## References

**See also:** [HoldAll](../../expression-information/HoldAll/), [List](../../other-advanced/List/), [NDArrayQ](../../other-advanced/NDArrayQ/), [Sum](../../calculus/Sum/), [Product](../../calculus/Product/)

- Source: [`src/list.c`](https://github.com/stblake/mathilda/blob/main/src/list.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_bernoullib.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bernoullib.c)
- Tests: [`tests/test_catch_throw.c`](https://github.com/stblake/mathilda/blob/main/tests/test_catch_throw.c)

## Notes & additional examples

### Notes

The single-argument iterator `{i, imax}` runs `i` from 1 to `imax`; `{i, imin,
imax}` and `{i, imin, imax, step}` give explicit bounds and a step. The step may
be an exact rational, so the values stay exact (`{0, 1/2, 1}`). Multiple iterator
specifications nest: the leftmost varies slowest, producing a list of lists.
`Table[expr, n]` with a plain count simply repeats `expr` `n` times. `Table` holds
its arguments, so the body is only evaluated as each iterator value is assigned —
the body may itself be a `Sum`, `D`, or any other computation (giving exact
harmonic numbers, identity matrices, and the like).
