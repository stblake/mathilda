# PrimitiveRoot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PrimitiveRoot[n]
    gives a primitive root of n.
PrimitiveRoot[n, k]
    gives the smallest primitive root of n greater than or equal to k.

A primitive root of n is a generator of the multiplicative group of integers modulo n relatively prime to n.  PrimitiveRoot returns unevaluated unless n is 2, 4, an odd prime power p^k, or twice an odd prime power 2 p^k.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= PrimitiveRoot[9]
Out[1]= 2

In[2]:= PrimitiveRoot[10]
Out[2]= 7

In[3]:= PrimitiveRoot[10, 1]
Out[3]= 3

In[4]:= PrimitiveRoot[10, 4]
Out[4]= 7

In[5]:= PrimitiveRoot[{9, 7, 19}]
Out[5]= {2, 3, 2}

In[6]:= PrimitiveRoot[12]
Out[6]= PrimitiveRoot[12]
```

## Implementation notes

**Algorithm.** `builtin_primitiveroot` returns the smallest primitive root of `n` ≥ an optional second-argument start (`PrimitiveRoot[n]` / `PrimitiveRoot[n, k]`). It first classifies `n` with `pr_classify` to confirm the unit group `(Z/nZ)*` is cyclic (i.e. `n ∈ {1, 2, 4, p^e, 2p^e}` for odd prime p), then computes `φ(n)` and its distinct prime divisors. `pr_smallest_primitive_root` scans candidates `g`, testing each with `pr_is_primitive_root`: `g` is a primitive root iff `gcd(g, n) = 1` and `g^(φ(n)/q) ≢ 1 (mod n)` for every prime `q | φ(n)` (via `mpz_powm`). Non-integer numeric input emits `PrimitiveRoot::intg`; `n < 2` likewise; wrong arg count emits `PrimitiveRoot::argt`; symbolic input returns unevaluated.

**Data structures.** GMP `mpz_t`; distinct primes of `φ(n)` in a fixed `mpz_t[]` array.

**Complexity / limits.** Primitive-root density is `φ(φ(n))/φ(n)`, so the scan finds one in roughly `O(log log p)` candidates on average; each test is `ω(φ(n))` modular exponentiations.

- `Protected`, `Listable`.
- Returns unevaluated unless `n` is 2, 4, an odd prime power $p^k$, or

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/numbertheory.c`](https://github.com/stblake/mathilda/blob/main/src/numbertheory.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)
