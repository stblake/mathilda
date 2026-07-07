# Count

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Count[list, pattern] gives the number of elements in list that match pattern.
Count[expr, pattern, levelspec] gives the total number of subexpressions matching pattern that appear at the levels in expr specified by levelspec.
Count[pattern] represents an operator form of Count that can be applied to an expression.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Cases[<|"a" -> 1, "b" -> 2, "c" -> 3|>, x_ /; x > 1]
Out[1]= {2, 3}

In[2]:= Count[<|"a" -> 1, "b" -> 2, "c" -> 3|>, x_ /; x > 1]
Out[2]= 2

In[3]:= DeleteCases[<|"a" -> 1, "b" -> 2, "c" -> 3|>, x_ /; x > 1]
Out[3]= <|"a" -> 1|>
```

## Implementation notes

**Algorithm.** `builtin_count` (`src/patterns.c`) tallies the subexpressions matching a pattern. Level-spec (default `{1,1}` — immediate elements only), the `Heads -> True|False` option, and the argument shapes mirror `Cases`. The worker `do_count_at_level` recurses depth-first into the head (when `heads`) and every argument, and at each in-range node calls `match(e, pattern, env)` from `src/match.c`, incrementing a `size_t` counter on success; level membership for negative specs is resolved via `get_expr_depth_patterns`. Unlike `Cases`/`Position` it stores nothing and has no result limit — it just returns the final count as an `EXPR_INTEGER`. `Count[pat]` with one argument returns the operator form `Function[Count[#1, pat]]`.

**Data structures.** A single `size_t` accumulator; one `MatchEnv` per node tested. No results buffer.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/patterns.c`](https://github.com/stblake/mathilda/blob/main/src/patterns.c)
- Specification: [`docs/spec/builtins/data-structures.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/data-structures.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= Count[{1,2,1,3,1},1]
Out[1]= 3
```

```mathematica
In[1]:= Count[{1,2,3,4,5,6},_?EvenQ]
Out[1]= 3
```

```mathematica
In[1]:= Count[Range[100], _?PrimeQ]
Out[1]= 25
```

```mathematica
In[1]:= Count[{a, b, {c, a}, a, {a, {a}}}, a, Infinity]
Out[1]= 5
```

```mathematica
In[1]:= Count[IntegerDigits[2^100], _?(# > 5 &)]
Out[1]= 11
```

### Notes

The second argument is a pattern, so `Count[list, _?EvenQ]` counts elements satisfying a predicate, not just literal matches. A level specification such as `Infinity` makes `Count` recurse into subexpressions, so it counts every matching leaf at any depth — here `Count[Range[100], _?PrimeQ]` recovers the prime-counting value 25, and the last example tallies how many decimal digits of `2^100` exceed 5.
