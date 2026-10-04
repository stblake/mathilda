# FirstPosition

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FirstPosition[expr, pattern] gives the position of the first element in expr matching pattern (in depth-first order), or Missing["NotFound"] if no such element is found.`**

**`FirstPosition[expr, pattern, default] gives default if no element matching pattern is found; default is evaluated only when it is returned.`**

**`FirstPosition[expr, pattern, default, levelspec] finds only objects on the levels specified by levelspec.`**

<details>
<summary>Notes</summary>

Over an association, the position is a key, e.g. {Key\[k\], ...}. The default level specification is {0, Infinity} with Heads -\> True; a position of {} represents the whole of expr.

</details>

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= FirstPosition[{a, b, c, b}, b]
Out[1]= {2}

In[2]:= FirstPosition[{{1, 2}, {3, 4}}, 4]
Out[2]= {2, 2}

In[3]:= FirstPosition[{1, 2, 3}, 5, None]
Out[3]= None
```

### Applications (6)

```mathematica
In[4]:= FirstPosition[{a, b, c, b}, b]
Out[4]= {2}
```

A nested position, depth first

```mathematica
In[5]:= FirstPosition[{{1, 2}, {3, 4}}, 4]
Out[5]= {2, 2}
```

```mathematica
In[6]:= FirstPosition[{1, 2, 3}, 5]
Out[6]= Missing["NotFound"]
```

The held default is returned when nothing matches

```mathematica
In[7]:= FirstPosition[{1, 2, 3}, 5, None]
Out[7]= None
```

Over an association the position is a key

```mathematica
In[8]:= FirstPosition[<|"a" -> 1, "b" -> 2|>, 2]
Out[8]= {Key["b"]}
```

Restricted to level 1

```mathematica
In[9]:= FirstPosition[{1, 2, 3}, _Integer, None, {1}]
Out[9]= {1}
```

## Implementation notes

**Algorithm.** `builtin_first_position` (`src/patterns.c`) gives the position of
the first subexpression matching a pattern, in depth-first order. Rather than
re-implement traversal it **delegates to `Position` with a first-match cap** and
takes the first result, inheriting `Position`'s levelspec parsing, `Heads`
handling, association value → `Key[...]` remapping, and depth-first ordering
unchanged. A visible `NDArray` first argument is an atom to the matcher, so
`patterns_delist_visible` materialises it before matching.

The handler splits the trailing `Heads -> True|False` option from the positional
arguments `expr`, `pattern`, `default`, `levelspec`. It then builds the delegated
call: for a 2-argument association it uses the 2-argument `Position` form (the
only one that does the value→`Key` remap); otherwise it passes the supplied
levelspec — or the `{0, Infinity}` default, which includes level 0 (the whole
expression, position `{}`) — plus the match count `1` so `Position` stops at the
first hit, forwarding any `Heads` option. If `Position` returns a non-empty list
the first element is the answer; otherwise the held `default` is returned
(evaluated by the fixed-point loop only when returned), falling back to
`Missing["NotFound"]`.

**Attributes & limits.** `FirstPosition` carries `HoldRest | Protected`, which is
how `default` stays unevaluated until it is actually produced. It takes two to
four positional arguments.

- Attributes `{HoldRest, Protected}`; `default` is evaluated only when it is returned.
- Same defaults as `Position`: levels `{0, Infinity}` with `Heads -> True`; a position of `{}` is the whole of `expr`.
- Over an association the position is a key, e.g. `{Key[k], ...}`.
- Delegates to `Position` (with the first-match cap), so traversal, level specs, the `Heads` option, and association handling match `Position` exactly.

**Attributes:** `HoldRest`, `Protected`.

## References

**See also:** [Position](../../data-structures/Position/)

- Source: [`src/patterns.c`](https://github.com/stblake/mathilda/blob/main/src/patterns.c)
- Specification: [`docs/spec/builtins/pattern-matching.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/pattern-matching.md)
- Tests: [`tests/test_assoc_atomicity.c`](https://github.com/stblake/mathilda/blob/main/tests/test_assoc_atomicity.c)
- Tests: [`tests/test_firstposition.c`](https://github.com/stblake/mathilda/blob/main/tests/test_firstposition.c)

## Notes & additional examples

### Notes

`FirstPosition[expr, pattern]` gives the position — a list of indices — of the
first subexpression matching `pattern` in depth-first order, or
`Missing["NotFound"]` if there is none. `FirstPosition[expr, pattern, default]`
returns `default` instead; `default` is held and evaluated only when it is
actually returned. A fourth argument is a level specification.

It delegates to `Position` (with the match count fixed at 1), so it inherits
`Position`'s traversal, default level spec `{0, Infinity}`, `Heads -> True`
handling, and association value → `Key[...]` remapping exactly. Because the
default spec includes level 0, a pattern that matches the whole of `expr` gives
the position `{}`.
