# FirstCase

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FirstCase[expr, patt]`**

Gives the first element of expr matching patt, or Missing\["NotFound"\]. FirstCase\[expr, patt, default\] uses default. Over an association, matches values and returns the first match.

## Examples (8)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= SelectFirst[{1, 3, 4, 5, 6}, EvenQ]
Out[1]= 4

In[2]:= FirstCase[{1, 2, 3, 4}, _?EvenQ]
Out[2]= 2

In[3]:= SelectFirst[{1, 3, 5}, EvenQ, None]
Out[3]= None
```

### Applications (5)

The first element matching the pattern

```mathematica
In[4]:= FirstCase[{1, "a", 2, "b"}, _String]
Out[4]= "a"
```

A conditional pattern

```mathematica
In[5]:= FirstCase[{1, 2, 3, 4}, x_ /; x > 2]
Out[5]= 3
```

A transformation rule returns the rewritten match

```mathematica
In[6]:= FirstCase[{1, 2, 3, 4}, x_?EvenQ -> x^2]
Out[6]= 4
```

No match, so the supplied default

```mathematica
In[7]:= FirstCase[{1, 2, 3}, _String, None]
Out[7]= None
```

No match and no default: Missing["NotFound"]

```mathematica
In[8]:= FirstCase[{1, 2, 3}, _String]
Out[8]= Missing["NotFound"]
```

## Implementation notes

**Algorithm.** `builtin_first_case` is a thin front-end over `Cases`: it builds
`Cases[expr, pattern]`, evaluates it, and returns a copy of the first element of
the resulting match list. Delegating this way means `FirstCase` inherits `Cases`'s
whole surface for free — the pattern may be a plain pattern or a transformation
rule `patt -> rhs` / `patt :> rhs`, in which case the returned value is the
transformed first match, not the matching element. A visible `NDArray` is an atom
to the matcher, so it is materialised to a list first (`patterns_delist_visible`).

If `Cases` yields no match, the three-argument form `FirstCase[expr, pattern,
default]` returns `default` and the two-argument form returns
`Missing["NotFound"]`.

**Data structures.** One synthesised `Cases[...]` call `Expr`, evaluated and
freed; the result is a copy of the match list's first argument (or of the
supplied default).

**Complexity / limits.** Because it delegates to the two-argument `Cases`, the
full match list is materialised before the first element is taken — there is no
first-match early stop (contrast `FirstPosition`, which passes `n = 1` to
`Position`). A level spec is not accepted: the third argument is the default, not
a level.

**Attributes:** `Protected`.

## References

**See also:** [SelectFirst](../../functional-programming/SelectFirst/)

- Source: [`src/patterns.c`](https://github.com/stblake/mathilda/blob/main/src/patterns.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_ndarray_functions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndarray_functions.c)

## Notes & additional examples

### Notes

`FirstCase[expr, pattern]` gives the first element of `expr` matching `pattern`,
the single-element companion to `Cases`. When `pattern` is a transformation rule
`patt -> rhs`, the result is the rewritten first match rather than the raw element
(here `4` is `2^2`). With no match it returns `Missing["NotFound"]`, or the
`default` passed as a third argument.
