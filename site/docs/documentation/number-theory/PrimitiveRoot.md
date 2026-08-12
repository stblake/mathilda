# PrimitiveRoot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`PrimitiveRoot[n]`**

gives a primitive root of n.

**`PrimitiveRoot[n, k]`**

gives the smallest primitive root of n greater than or equal to k.

<details>
<summary>Notes</summary>

A primitive root of n is a generator of the multiplicative group of integers modulo n relatively prime to n.  PrimitiveRoot returns unevaluated unless n is 2, 4, an odd prime power p^k, or twice an odd prime power 2 p^k.

</details>

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

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

### Applications (5)

```mathematica
In[1]:= PrimitiveRoot[7]
Out[1]= 3
```

It finds generators even for large primes; 5 generates the multiplicative group modulo the prime `10^9 + 7`:

```mathematica
In[1]:= PrimitiveRoot[10^9 + 7]
Out[1]= 5
```

Prime powers are supported directly:

```mathematica
In[1]:= PrimitiveRoot[3^5]
Out[1]= 2
```

The two-argument form returns the smallest primitive root not below `k`:

```mathematica
In[1]:= PrimitiveRoot[7, 5]
Out[1]= 5
```

When the multiplicative group is non-cyclic (e.g. `n = 8`), no primitive root exists and the call stays unevaluated:

```mathematica
In[1]:= PrimitiveRoot[8]
Out[1]= PrimitiveRoot[8]
```

## Implementation notes

**Algorithm.** `builtin_primitiveroot` returns the smallest primitive root of `n` ≥ an optional second-argument start (`PrimitiveRoot[n]` / `PrimitiveRoot[n, k]`). It first classifies `n` with `pr_classify` to confirm the unit group `(Z/nZ)*` is cyclic (i.e. `n ∈ {1, 2, 4, p^e, 2p^e}` for odd prime p), then computes `φ(n)` and its distinct prime divisors. `pr_smallest_primitive_root` scans candidates `g`, testing each with `pr_is_primitive_root`: `g` is a primitive root iff `gcd(g, n) = 1` and `g^(φ(n)/q) ≢ 1 (mod n)` for every prime `q | φ(n)` (via `mpz_powm`). Non-integer numeric input emits `PrimitiveRoot::intg`; `n < 2` likewise; wrong arg count emits `PrimitiveRoot::argt`; symbolic input returns unevaluated.

**Data structures.** GMP `mpz_t`; distinct primes of `φ(n)` in a fixed `mpz_t[]` array.

**Complexity / limits.** Primitive-root density is `φ(φ(n))/φ(n)`, so the scan finds one in roughly `O(log log p)` candidates on average; each test is `ω(φ(n))` modular exponentiations.

- `Protected`, `Listable`.
- Returns unevaluated unless `n` is 2, 4, an odd prime power $p^k$, or
  twice an odd prime power $2 p^k$ (the moduli for which $(\mathbb{Z}/n\mathbb{Z})^*$
  is cyclic). For all other `n`, the call is left unevaluated.
- The 1-argument form returns a canonical primitive root: smallest for
  $n \in \{2, 4\}$ and odd prime powers; for $n = 2 p^k$ the formula
  $g$ if $g$ is odd else $g + p^k$ is applied, where $g$ is the smallest
  primitive root of $p^k$. This matches Mathematica's convention so that,
  e.g. `PrimitiveRoot[10] == 7` while `PrimitiveRoot[10, 1] == 3`.
- The 2-argument form walks forward from `k`; if `k > n - 1` the call is
  left unevaluated.
- All arithmetic uses GMP `mpz_t`, so machine integers, bignums, and
  symbolic bignum products like `Prime[1000000]^1000000` are handled
  uniformly. The prime-power detection iteratively strips prime exponents
  via `mpz_root`, which runs in $O(\omega(k))$ root extractions.
- Diagnostics:
  - `PrimitiveRoot::argt` if not called with 1 or 2 arguments.
  - `PrimitiveRoot::intg` if `n` (or the 2nd-arg `k` when numeric) is not
    an integer greater than 1.

**Attributes:** `Listable`, `Protected`.

## References

- Source: [`src/numbertheory.c`](https://github.com/stblake/mathilda/blob/main/src/numbertheory.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_multiplicative_order.c`](https://github.com/stblake/mathilda/blob/main/tests/test_multiplicative_order.c)
- Tests: [`tests/test_primitive_root.c`](https://github.com/stblake/mathilda/blob/main/tests/test_primitive_root.c)

## Notes & additional examples

### Notes

`PrimitiveRoot[n]` returns a generator of the multiplicative group of integers
coprime to `n`. Such a generator exists only when `n` is `2`, `4`, an odd prime
power `p^k`, or twice one (`2 p^k`); for all other `n` the group is non-cyclic
and the call is left unevaluated. `PrimitiveRoot[n, k]` returns the smallest
primitive root that is `>= k`.
