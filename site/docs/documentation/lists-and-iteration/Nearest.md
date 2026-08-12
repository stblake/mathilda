# Nearest

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Nearest[list, x]`**

Gives the element of list closest to x, as a list. All elements tied at the minimum distance Abs\[element - x\] are returned, in their original order; an empty list gives {}. Returns unevaluated unless every distance is a real number, so a symbolic element or target leaves the expression unchanged rather than dropping it from the result.

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= Nearest[{1, 5, 10}, 3]
Out[1]= {1, 5}

In[2]:= Nearest[{10, 20, 30}, 100]
Out[2]= {30}

In[3]:= Nearest[{1, a, 3}, 2]
Out[3]= Nearest[{1, a, 3}, 2]
```

## Algorithm

Nearest[list, x] -- the element(s) of `list` closest to the target `x`.

```text
  Nearest[{1, 5, 10}, 3]    -> {1, 5}
  Nearest[{10, 20, 30}, 100] -> {30}
  Nearest[{}, 3]             -> {}
```

Distance is Abs[element - x], composed from the existing internal_subtract and internal_abs rather than a bespoke helper, so a complex element uses its modulus for free: Nearest[{3 + 4 I, 1}, 0] is {1} because the distances are 5 and 1.

ALL elements tied at the minimum distance are returned, in their original order. That requirement is what fixes the algorithm. The obvious neighbour, RankedMin (sort.c), selects with a quickselect whose comparator carries an original-index tiebreak (sort.c:932) existing precisely to make ties IMPOSSIBLE, so that exactly one element can win -- the opposite of what is needed here. The shape used instead is MinimalBy's (sort.c:663-716): one pass to find the minimum, a second to collect every distance equal to it. Input order among ties then falls out of the ascending collect pass, so there is no tie logic in this file at all.

Nearest diverges from MinimalBy in exactly one respect, and deliberately. Every distance must be a real number, or the whole call stays unevaluated. MinimalBy[{1, a, 3}, Abs[# - 2] &] answers {1, 3}: expr_compare orders symbols after all numbers (sort.c:379-380), so the symbolic element is never minimal and silently vanishes from a result that still looks plausible. Gating on the DISTANCE rather than on the element covers a symbolic element, a symbolic target, and a non-real complex in a single check. A symbolic real such as Pi is rejected too -- Abs[Pi - 3] stays as Abs[-3 + Pi] -- rather than being numericalized the way RankedMin's ranked_numeric_key would.

Cost: O(n) distance evaluations, O(n) comparisons, O(n) peak extra memory. The two evaluate passes per element dominate, so this is an interpreter-speed path and not a buffer path. There is no exact-hit short circuit: a later element can tie at distance 0, and dropping it would break the tie contract.

Only the two-argument form lives here. The n-nearest, radius, rule, all-pairs, and NearestTo operator forms, and the DistanceFunction option, are separate follow-ups.

## Performance

Measured on arm64 Darwin at commit `2dea9cc05`.

| case | n | time |
|---|---:|---:|
| nearest to a target | 1,000 | 367 us |
| nearest to a target | 10,000 | 3.5 ms |
| nearest to a target | 100,000 | 34.3 ms |

## Implementation notes

- `Protected`.
- Distance is `Abs[element - x]`, so a complex element uses its modulus.
- **All** elements tied at the minimum distance are returned, in their original
  order: `Nearest[{1, 5, 10}, 3]` gives `{1, 5}`, not `{1}`.
- Distances are compared by numeric value, so exact and inexact distances of
  equal value tie: `Nearest[{0, 2.0}, 1]` gives `{0, 2.0}`.
- An empty list gives `{}`.
- Returns unevaluated unless every distance is a real number. A symbolic element
  or target leaves the expression unchanged rather than being dropped from the
  result — unlike
  [`MinimalBy`](../functional-programming/index.md), which orders
  symbolic keys after all numbers and so silently omits them.
- A rational with a bigint component declines, because `Abs` does not evaluate
  one: `Nearest[{1/10^25, 1}, 0]` stays unevaluated.
- Only the two-argument form is implemented. The `n`-nearest, radius, rule, and
  all-pairs forms, the `NearestTo` operator form, and the `DistanceFunction`
  option are not yet available.

**Attributes:** `Protected`.

## See also

[MinimalBy](../../functional-programming/MinimalBy/), [Abs](../../arithmetic/Abs/)

## References

- Source: [`src/list/list_init.c`](https://github.com/stblake/mathilda/blob/main/src/list/list_init.c)
- Specification: [`docs/spec/builtins/lists-and-iteration.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/lists-and-iteration.md)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
