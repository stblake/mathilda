# PrimitiveRootList

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PrimitiveRootList[n]
    gives the sorted list of all primitive roots of n in the canonical residues {1, ..., n-1}.

Returns an empty list unless n is 2, 4, an odd prime power p^k, or twice an odd prime power 2 p^k.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= PrimitiveRootList[9]
Out[1]= {2, 5}

In[2]:= PrimitiveRootList[19]
Out[2]= {2, 3, 10, 13, 14, 15}

In[3]:= PrimitiveRootList[12]
Out[3]= {}

In[4]:= Union[Table[PowerMod[2, i, 9], {i, 6}]]
Out[4]= {1, 2, 4, 5, 7, 8}
```

## Implementation notes

`builtin_primitiverootlist` returns the sorted list of all primitive roots of `n` in `[1, n-1]`. It classifies `n` for cyclicity (`pr_classify`; non-cyclic or `n ≤ 1` gives `{}`), computes `φ(n)` and its distinct prime divisors, and finds the smallest primitive root `g` of `n` (`pr_smallest_primitive_root`, same test as `PrimitiveRoot`). The full set is then enumerated as `{g^i mod n : 1 ≤ i ≤ φ(n), gcd(i, φ(n)) = 1}` — there are exactly `φ(φ(n))` of them — and the residues are sorted. Wrong arg count emits `PrimitiveRootList::argx`; non-integer input returns unevaluated with no diagnostic. The enumeration bails (NULL) if `φ(n)` does not fit in `unsigned long`. GMP `mpz_t` throughout.

- `Protected`, `Listable`.
- Returns `{}` if `n` is not 2, 4, an odd prime power, or twice an odd

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/numbertheory.c`](https://github.com/stblake/mathilda/blob/main/src/numbertheory.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)
