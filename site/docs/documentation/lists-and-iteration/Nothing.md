# Nothing

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Nothing`**

is a symbol that is automatically removed from any list in which it appears as an element: {a, Nothing, b} evaluates to {a, b}. It is the identity element for list construction, so Table\[If\[test, val, Nothing\], ...\] builds a list of just the values for which test held. Any Nothing\[...\] form is removed likewise. Non-list heads treat Nothing as an ordinary symbol.

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= {1, Nothing, 2}
Out[1]= {1, 2}

In[2]:= Table[If[EvenQ[i], i, Nothing], {i, 6}]
Out[2]= {2, 4, 6}
```

### Applications (3)

Every Nothing simply disappears

```mathematica
In[3]:= {a, Nothing, b, Nothing, c}
Out[3]= {a, b, c}
```

A Nothing[...] form is removed too

```mathematica
In[4]:= {1, 2, Nothing[x, y], 3}
Out[4]= {1, 2, 3}
```

Keep only the elements that pass a test

```mathematica
In[5]:= Table[If[PrimeQ[i], i, Nothing], {i, 10}]
Out[5]= {2, 3, 5, 7}
```

## Implementation notes

**Algorithm.** `Nothing` is not a builtin function but a special symbol handled
by the evaluator: it is the identity element for list construction and vanishes
from any `List` it appears in. `strip_nothing` (`src/eval.c`), called during
`evaluate_step`, scans a `List`'s arguments and removes every element that is
either the symbol `Nothing` or any `Nothing[...]` form (an `EXPR_FUNCTION` whose
head is `Nothing`). It makes one counting pass to find how many survive, and only
if some element is a `Nothing` does it rebuild the argument array — freeing the
removed nodes, keeping the rest, and invalidating the expression's cached hash.

**List-specific.** The rewrite fires only for the head `List`; for any other
head `Nothing` is an ordinary symbol. This is what makes
`Table[If[test, val, Nothing], ...]` the idiom for conditionally building a list:
the iterator produces a plain `List` whose unwanted slots hold `Nothing`, and the
evaluator strips them on the next step.

**Data structures / limits.** O(n) over the list's arguments, one reallocation
when a `Nothing` is present (none otherwise). `Nothing` carries `ATTR_PROTECTED`
(`src/attr.c`). No numeric surface — it is a structural identity element.

**Attributes:** `Protected`.

## References

- Source: [`src/eval.c`](https://github.com/stblake/mathilda/blob/main/src/eval.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
- Tests: [`tests/test_sequence.c`](https://github.com/stblake/mathilda/blob/main/tests/test_sequence.c)

## Notes & additional examples

### Notes

`Nothing` is a symbol that is automatically removed from any list in which it
appears. It is the identity element for list construction, which makes
`Table[If[test, val, Nothing], ...]` the standard idiom for building a list of
just the values for which a test holds — no `Select` or `DeleteCases` pass is
needed. Any `Nothing[...]` form is stripped likewise.

The removal is list-specific: for a non-`List` head, `Nothing` is an ordinary
symbol and is left in place.
