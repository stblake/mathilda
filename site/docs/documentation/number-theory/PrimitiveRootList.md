# PrimitiveRootList

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PrimitiveRootList[n] gives the sorted list of all primitive roots of n in the canonical residues {1, ..., n-1}; when n admits any, there are exactly EulerPhi[EulerPhi[n]] of them.`**

<details>
<summary>Notes</summary>

Returns an empty list unless n is 2, 4, an odd prime power p^k, or twice an odd prime power 2 p^k.

</details>

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

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

### Applications (5)

```mathematica
In[5]:= PrimitiveRootList[7]
Out[5]= {3, 5}

In[6]:= PrimitiveRootList[18]
Out[6]= {5, 11}

In[7]:= PrimitiveRootList[15]
Out[7]= {}

In[8]:= Length[PrimitiveRootList[101]]
Out[8]= 40

In[9]:= EulerPhi[EulerPhi[101]]
Out[9]= 40
```

## Implementation notes

`builtin_primitiverootlist` returns the sorted list of all primitive roots of `n` in `[1, n-1]`. It classifies `n` for cyclicity (`pr_classify`; non-cyclic or `n ≤ 1` gives `{}`), computes `φ(n)` and its distinct prime divisors, and finds the smallest primitive root `g` of `n` (`pr_smallest_primitive_root`, same test as `PrimitiveRoot`). The full set is then enumerated as `{g^i mod n : 1 ≤ i ≤ φ(n), gcd(i, φ(n)) = 1}` — there are exactly `φ(φ(n))` of them — and the residues are sorted. Wrong arg count emits `PrimitiveRootList::argx`; non-integer input returns unevaluated with no diagnostic. The enumeration bails (NULL) if `φ(n)` does not fit in `unsigned long`. GMP `mpz_t` throughout.

- `Protected`, `Listable`.
- Returns `{}` if `n` is not 2, 4, an odd prime power, or twice an odd
  prime power.
- Enumerates the $\varphi(\varphi(n))$ primitive roots as $g^i \bmod n$
  for $i \in \{1, \ldots, \varphi(n)\}$ with $\gcd(i, \varphi(n)) = 1$,
  where $g$ is the smallest primitive root of `n`. The list is sorted
  ascending.
- Falls back to unevaluated when $\varphi(n)$ exceeds `unsigned long`,
  since the enumeration cannot be represented.
- Non-integer numeric inputs (e.g. `11.0`, `11 + I`) flow through
  unevaluated with no diagnostic, matching Mathematica.
- Diagnostic: `PrimitiveRootList::argx` if not called with exactly 1
  argument.

**Attributes:** `Listable`, `Protected`.

## References

- K. Ireland and M. Rosen, *A Classical Introduction to Modern Number Theory*, 2nd ed., Springer, 1990 — primitive roots and the structure of `(Z/nZ)*` (Chapter 4).
- Source: [`src/numbertheory.c`](https://github.com/stblake/mathilda/blob/main/src/numbertheory.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_primitive_root.c`](https://github.com/stblake/mathilda/blob/main/tests/test_primitive_root.c)

## Notes & additional examples

### How many primitive roots?

When `n` admits a primitive root at all (`n = 2, 4, p^k, 2p^k`), the group `(Z/nZ)*` is
cyclic of order `EulerPhi[n]`, and a cyclic group of order `m` has exactly `EulerPhi[m]`
generators. So the number of primitive roots is `EulerPhi[EulerPhi[n]]` — e.g. two for
`n = 7`, since `EulerPhi[EulerPhi[7]] = EulerPhi[6] = 2`.

### Notes

`PrimitiveRootList[n]` returns every primitive root of `n` in canonical
residues `{1, ..., n-1}`, sorted. The list is non-empty only when `n` is `2`,
`4`, an odd prime power `p^k`, or `2 p^k`; otherwise it is `{}`. When primitive
roots exist there are exactly `EulerPhi[EulerPhi[n]]` of them, since the cyclic
group of order `EulerPhi[n]` has that many generators.
