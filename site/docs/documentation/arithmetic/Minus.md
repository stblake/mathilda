# Minus

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Minus[x] is the arithmetic negation of x, equivalent to -x (Times[-1, x]).`**

Listable; Minus\[x, y\] (any count other than one) is left unevaluated with Minus::argx.

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Minus[3]
Out[1]= -3

In[2]:= Minus[{1, 2.5, 1/2}]
Out[2]= {-1, -2.5, -1/2}

In[3]:= SortBy[{3, 1, 2}, Minus]
Out[3]= {3, 2, 1}
```

### Applications (6)

Unary negation, -x

```mathematica
In[4]:= Minus[5]
Out[4]= -5
```

-1 distributes through the sum

```mathematica
In[5]:= Minus[a + b]
Out[5]= -(a + b)
```

Negates a complex number

```mathematica
In[6]:= Minus[2 + 3 I]
Out[6]= -2 - 3*I
```

Listable: threads over the list

```mathematica
In[7]:= Minus[{1, -2, 3}]
Out[7]= {-1, 2, -3}
```

As a sort key, orders descending

```mathematica
In[8]:= SortBy[{3, -1, 2, -5}, Minus]
Out[8]= {3, 2, -1, -5}
```

The single rewrite Minus[x] -> Times[-1, x]

```mathematica
In[9]:= FullForm[Minus[x]]
Out[9]= Times[-1, x]
```

## Algorithm

minus.c -- Minus[x], the functional form of unary negation.

The parser already reads `-x` as Times[-1, x], so Minus only has to exist as a callable head: Minus[x] rewrites to Times[-1, x] and lets Times do the arithmetic (numbers, Rationals, Complex, Plus distribution, packed arrays). This is what makes SortBy[list, Minus] and KeySortBy[a, Minus] work.

```text
Attributes match Mathematica: Listable, NumericFunction, Protected.  Any
```

other argument count emits Minus::argx and stays unevaluated (WL prints Minus[x, y] as x − y but does not evaluate it).

## Implementation notes

**Algorithm.** The parser already reads `-x` as `Times[-1, x]`, so `Minus` only
has to exist as a *callable* head. `builtin_minus` returns `Times[-1,
expr_copy(arg)]` and lets `Times` do all the arithmetic — integers, Rationals,
Complex, `Plus` distribution, packed/NDArray buffers. This single-rewrite form
is precisely what makes `Minus` usable as a sort key (`SortBy[list, Minus]`,
`KeySortBy[a, Minus]`). Any argument count other than one emits `Minus::argx`
(through `builtin_arg_error`) and the call is left unevaluated, matching
Wolfram Language, which prints `Minus[x, y]` as `x - y` but does not evaluate it.

**Data structures.** Pure `Expr`: one freshly-built `Times[-1, x]` node, with
the argument deep-copied so the evaluator can still free the original call. Since
the result is simply a `Times`, every surface `Times` already supports is
inherited for free — `Minus` is on `pack.c`'s `AWARE` list for exactly this
reason (a packed or visible `NDArray` argument stays on the buffer once `Times`
sees it), and `Compile[]` lowers the negation at both scalar and rank-1 array
shapes (`compile_emit_arith.c`, with the machine- and GMP-integer domains in
`compile_mgd.c` and the type rule in `compile_infer.c`).

**Complexity / limits.** `O(1)` to build the node; the real cost is whatever
`Times[-1, x]` then does. Attributes are `Listable | NumericFunction |
Protected`, so `Minus` threads over a list and `CompileDiagnostics` reports
`Compiled -> True` at scalar and rank-1 shapes. The same module also defines
`` Internal`SyntacticNegativeQ ``, the leading-minus-sign predicate the ported
integrator's term ordering is built on.

- `Listable`, `NumericFunction`, `Protected`, as in Mathematica.
- Rewrites to `Times[-1, x]`, so numbers, rationals, complex numbers, `Plus`
  distribution and symbolic arguments all behave exactly like `-x`.
- Usable as a sort key: `SortBy[list, Minus]`, `KeySortBy[a, Minus]`.
- Packed arrays and visible `NDArray`s stay on the buffer (`Minus` is on the
  packed-aware list and hands the buffer to `Times`); `Compile[]` lowers it at
  scalar and array shapes.
- Any other argument count emits `Minus::argx` and stays unevaluated.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Plus](../../arithmetic/Plus/), [NDArray](../../linear-algebra/NDArray/), [Times](../../arithmetic/Times/)

- Source: [`src/minus.c`](https://github.com/stblake/mathilda/blob/main/src/minus.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_assoc_forms.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_forms.c)

## Notes & additional examples

### Notes

`Minus[x]` is the functional form of unary negation, equivalent to `-x`. The
parser already turns the `-` operator into `Times[-1, x]`, so `Minus` exists
mainly as a head you can pass by name: `SortBy[list, Minus]` and
`KeySortBy[assoc, Minus]` sort by the negated value (i.e. descending), and
`Map[Minus, list]` negates each element. The rewrite to `Times[-1, x]` means
`Minus` inherits everything `Times` can do — distributing over a `Plus`,
negating a `Complex`, and operating on packed or `NDArray` arguments directly on
the buffer. It also lowers inside `Compile[]`.

`Minus` is strictly unary. `Minus[x, y]` (or any count other than one) emits the
message `Minus::argx` and is left unevaluated; use `Subtract[x, y]` (or `x - y`)
for a binary difference.
