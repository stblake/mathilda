# ComplexQ

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`ComplexQ[expr]`**

gives True if expr is a Complex number, False otherwise.

<details>
<summary>Notes</summary>

Tests for the Complex\[re, im\] head; a purely real number is not Complex.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

Stored as Complex[2, 3], so the head test passes

```mathematica
In[1]:= ComplexQ[2 + 3 I]
Out[1]= True
```

Inexact parts are still a Complex

```mathematica
In[2]:= ComplexQ[2.0 + 3.0 I]
Out[2]= True
```

A bare real is not Complex

```mathematica
In[3]:= ComplexQ[2.0]
Out[3]= False
```

Nor is an integer

```mathematica
In[4]:= ComplexQ[5]
Out[4]= False
```

## Implementation notes

**Algorithm.** `builtin_complexq` is a one-argument syntactic predicate. It tests
whether its (already-evaluated) argument has head `Complex` — the single check
`core_head_is(arg, SYM_Complex)` on an interned-pointer comparison — and returns
the symbol `True` or `False` accordingly. A purely real number is *not* `Complex`:
Mathilda stores `2 + 3 I` as `Complex[2, 3]` but keeps a real as a bare `Real`, so
`ComplexQ[2.0]` is `False`. The test is on the head alone, so it does not care
whether the parts are exact or inexact.

**Data structures.** None beyond the argument node. The comparison is a pointer
equality against the interned `SYM_Complex` name; the result is a fresh
`EXPR_SYMBOL` (`True`/`False`).

**Complexity / limits.** `O(1)`. With fewer or more than one argument it declines
(returns `NULL`) and stays symbolic. Being a syntactic head test, it reports on
representation, not mathematical identity — a value that happens to be real-valued
but is written with a `Complex` head is a matter for the arithmetic that built it,
not for `ComplexQ`.

**Attributes:** `Protected`.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`ComplexQ[expr]` tests one thing: whether `expr` has head `Complex`. It is a
*syntactic* test on the representation, not a check of mathematical type — a real
number is stored with a real head, so `ComplexQ` is `False` for it regardless of
value. The parts may be exact or inexact; only the head matters. Anything that is
not a single-argument call is left unevaluated.
