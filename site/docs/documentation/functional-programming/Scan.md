# Scan

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Scan[f, expr]`**

Applies f to each element of expr for its side effects and returns Null, discarding the results.

**`Scan[f, expr, levelspec]`**

Applies f to the parts of expr selected by levelspec (default {1}), depth-first with leaves before roots. Option Heads-\>True also visits heads. Throw exits to an enclosing Catch; Return\[ret\] makes the final value ret. Over an association, applies f to each value.

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= s = 0; Scan[(s = s + #) &, <|"a" -> 1, "b" -> 2, "c" -> 3|>]; s
Out[1]= 6

In[2]:= Scan[Print, {{{a}}}, Infinity] a {a} {{a}}
Out[2]= {{Null a^3}}

In[3]:= Catch[Scan[If[# > 5, Throw[#]] &, {2, 4, 6, 8}]]
Out[3]= 6
```

### Options (1)

```mathematica
In[4]:= Scan[Print, {a, b}, Heads -> True] List a b
Out[4]= List Null a b
```

### Applications (4)

Scan visits each element; Sow records the order, Scan itself returns Null

```mathematica
In[5]:= Reap[Scan[Sow, {1, 2, 3}]]
Out[5]= {Null, {{1, 2, 3}}}
```

A level spec reaches the inner elements

```mathematica
In[6]:= Reap[Scan[Sow, {{1, 2}, {3, 4}}, {2}]]
Out[6]= {Null, {{1, 2, 3, 4}}}
```

Depth-first, leaves before roots: the atoms a, b, c

```mathematica
In[7]:= Reap[Scan[Sow, a + b c, {-1}]]
Out[7]= {Null, {{a, b, c}}}
```

Applied purely for the side effect; the value is Null

```mathematica
In[8]:= Scan[Print, {a, b}]
```

## Implementation notes

**Algorithm.** `builtin_scan` is `Map` run for effect: it applies `f` to each
selected part of `expr`, discards every result, and returns `Null`. It parses an
optional level spec (default `{1}`) and a `Heads -> True` option, then walks the
tree in `scan_at_level` depth-first, **leaves before roots** — sub-parts (and,
with `Heads -> True`, the head) are visited before the node itself. Each visit
goes through `scan_apply`, which builds `f[part]`, evaluates it, and throws the
result away unless control flow intervenes: an in-flight `Throw` is propagated to
an enclosing `Catch`, and a `Return[ret]` is classified against the `Scan`
boundary (`eval_classify_return`) and, when consumed there, becomes the call's
return value.

An association at the default level scans its values (mirroring `Map`). A
numeric-closed body at the default level takes the `numloop_scan` fast path — it
still runs the body so a non-finite element hands control back to the interpreter.
A visible `NDArray` at the default level iterates its leading axis directly
(scalar leaf for rank 1, sub-array row otherwise); any other spec materialises to
a nested list first.

**Data structures.** `scan_apply` adopts the part copy into the `f[...]` call,
evaluates, and frees. The traversal returns `NULL` to continue and a non-`NULL`
sentinel (`Throw`/`Return`) to stop and hand back, so no result list is ever
accumulated — the point of `Scan` over `Map`.

**Complexity / limits.** `O(n)` over the parts; a non-negative upper bound stops
the descent once past the maximum level, so `Scan[f, list]` does not walk each
element's whole subtree. `get_depth` is consulted only for a negative bound.

**Attributes:** `Protected`.

## References

**See also:** [Map](../../data-structures/Map/), [Throw](../../control-flow/Throw/), [Catch](../../control-flow/Catch/), [NDArray](../../linear-algebra/NDArray/)

- Source: [`src/funcprog.c`](https://github.com/stblake/mathilda/blob/main/src/funcprog.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_assoc_atomicity.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_atomicity.c)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_catch_throw.c`](https://github.com/stblake/mathilda/blob/main/tests/test_catch_throw.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)

## Notes & additional examples

### Notes

`Scan[f, expr]` applies `f` to each element the way `Map` does but keeps nothing:
it is for side effects (printing, sowing, logging) and always returns `Null`. The
traversal is depth-first with **leaves before roots**, so a level spec such as
`{-1}` visits the atoms of an expression before the subexpressions that contain
them.

A `Throw` inside `f` exits to an enclosing `Catch`; a `Return[ret]` evaluated
directly as `f` makes `Scan` itself return `ret`. Over an association, `f` is
applied to each value.
