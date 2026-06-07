---
source: src/numbertheory.c
---
`builtin_primitiverootlist` returns the sorted list of all primitive roots of `n` in `[1, n-1]`. It classifies `n` for cyclicity (`pr_classify`; non-cyclic or `n ≤ 1` gives `{}`), computes `φ(n)` and its distinct prime divisors, and finds the smallest primitive root `g` of `n` (`pr_smallest_primitive_root`, same test as `PrimitiveRoot`). The full set is then enumerated as `{g^i mod n : 1 ≤ i ≤ φ(n), gcd(i, φ(n)) = 1}` — there are exactly `φ(φ(n))` of them — and the residues are sorted. Wrong arg count emits `PrimitiveRootList::argx`; non-integer input returns unevaluated with no diagnostic. The enumeration bails (NULL) if `φ(n)` does not fit in `unsigned long`. GMP `mpz_t` throughout.
