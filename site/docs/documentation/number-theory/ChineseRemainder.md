# ChineseRemainder

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ChineseRemainder[{r1, r2, ...}, {m1, m2, ...}] gives the smallest x >= 0 satisfying the integer congruences x mod mi == ri mod mi for all i, lying in 0 <= x < LCM[m1, m2, ...].`**

**`ChineseRemainder[{r1, r2, ...}, {m1, m2, ...}, d]`**

gives the smallest such x \>= d, lying in d \<= x \< d + LCM\[m1, m2, ...\].

<details>
<summary>Notes</summary>

The moduli need not be pairwise coprime; a solution exists iff every pair of congruences agrees modulo the gcd of their moduli, and ChineseRemainder returns unevaluated when the system is inconsistent (e.g. ChineseRemainder\[{1, 2}, {6, 10}\]).  Solved by a streaming extended-Euclidean CRT fold, exact via GMP, so residues and moduli may be arbitrary-precision integers -- the recovery step of a residue number system.  Protected.

</details>

## Examples (12)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= ChineseRemainder[{3, 4}, {4, 5}]
Out[1]= 19

In[2]:= ChineseRemainder[{2, 3, 5}, {3, 5, 7}]
Out[2]= 68

In[3]:= ChineseRemainder[{2, 3, 5}, {3, 5, 7}, 100]
Out[3]= 173

In[4]:= ChineseRemainder[{1, 2}, {6, 10}]
Out[4]= ChineseRemainder[{1, 2}, {6, 10}]

In[5]:= Mod[ChineseRemainder[{123, 456, 789}, {1000000007, 1000000009, 1000000021}], {1000000007, 1000000009, 1000000021}]
Out[5]= {123, 456, 789}
```

### Worked examples (1)

```mathematica
In[6]:= ChineseRemainder[{}, {}]
Out[6]= 0
```

### Applications (6)

```mathematica
In[7]:= ChineseRemainder[{2, 3}, {3, 5}]
Out[7]= 8

In[8]:= ChineseRemainder[{1, 2, 3}, {2, 3, 5}]
Out[8]= 23

In[9]:= ChineseRemainder[{3, 4}, {5, 7}, 2]
Out[9]= 18

In[10]:= ChineseRemainder[{2, 8}, {6, 10}]
Out[10]= 8

In[11]:= ChineseRemainder[{1, 2}, {6, 10}]
Out[11]= ChineseRemainder[{1, 2}, {6, 10}]

In[12]:= ChineseRemainder[{123, 456}, {10^9 + 7, 10^9 + 9}]
Out[12]= 499999841499998989
```

## Implementation notes

**Algorithm.** `builtin_chineseremainder` reconstructs the integer `x` from a list of
residues `{r1, ..., rk}` and moduli `{m1, ..., mk}` by a **streaming
extended-Euclidean CRT fold**. It carries a running pair `(x, M)` -- the current
solution and the LCM of the moduli seen so far -- initialised from the first
congruence, and folds in each subsequent `(ri, mi)` in turn.

*The fold step.* To merge `x ≡ a (mod M)` with `x ≡ ri (mod mi)`, it runs the extended
Euclidean algorithm on `M` and `mi` to obtain `g = GCD[M, mi]` and Bezout cofactors.
The merged congruence is solvable iff `g` divides `ri - a`; if it does not, the system
is inconsistent and the builtin returns `NULL`, leaving the `ChineseRemainder[...]`
surface unevaluated. When it is solvable, the update sets the new modulus to
`LCM[M, mi] = M*mi/g` and lifts `x` to the unique residue modulo that LCM. Handling the
general `g > 1` case here is exactly what lets the moduli be non-coprime.

*Range normalisation.* The final `x` is reduced into `0 <= x < M` (or, given the
optional offset `d`, into `d <= x < d + M` by adding the appropriate multiple of `M`).

**Data structures.** All arithmetic is on `mpz_t` GMP integers, so residues and moduli
of arbitrary size are exact; results are normalised back to `EXPR_INTEGER` /
`EXPR_BIGINT` via `expr_bigint_normalize`. Inputs are validated to be two equal-length
lists of integer-like values (plus an optional integer offset) before the fold begins.

**Complexity.** Linear in the number of congruences, each fold step costing one extended
GCD -- `O(k * M(n) log n)` for `k` congruences of `n`-bit moduli, dominated by GMP's
sub-quadratic GCD. No factorisation of the moduli is required.

- Integer-only; machine integers and GMP bigints share one path, so a residue-number-system recovery over large coprime moduli returns the full bignum (`ChineseRemainder[data, keys]` with prime keys is the CRT-encryption idiom).
- The moduli need **not** be pairwise coprime: a solution exists iff every pair of congruences agrees modulo the gcd of their moduli. When the system is inconsistent, `ChineseRemainder` is left unevaluated, e.g. `ChineseRemainder[{1, 2}, {6, 10}]` (since `1 ≢ 2 (mod gcd(6, 10) = 2)`).
- Computed by a streaming pairwise CRT fold: an accumulator congruence `x ≡ x_acc (mod m_acc)` is merged with each `(ri, |mi|)` via `mpz_gcdext` (`x_acc += ((ri − x_acc)/g)·u·m_acc`, `m_acc = lcm`), reducing modulo the running lcm each step. The final `x = d + ((x_acc − d) mod L)` gives the smallest `x ≥ d` and handles a negative `d` and empty lists (`L = 1`).
- If all `0 ≤ ri < mi` the result satisfies `x mod mi == ri`; unreduced or negative residues are handled (the result gives `x mod mi == ri mod mi`).
- Edge cases: `ChineseRemainder[{}, {}]` → `0`; a zero modulus, a length mismatch between the two lists, and non-integer or symbolic arguments leave the call unevaluated.
- Diagnostic: `ChineseRemainder::argt` when called with other than 2 or 3 arguments.

**Attributes:** `Protected`.

## References

- Knuth, "The Art of Computer Programming, Vol. 2: Seminumerical Algorithms", §4.3.2 (modular arithmetic and the CRT).
- Ding, Pei & Salomaa, "Chinese Remainder Theorem: Applications in Computing, Coding, Cryptography" (1996).
- von zur Gathen & Gerhard, "Modern Computer Algebra", §5.4 (the extended Euclidean algorithm and CRT).
- von zur Gathen & Gerhard, "Modern Computer Algebra", §5.4.
- Source: [`src/numbertheory/chineseremainder.c`](https://github.com/stblake/mathilda/blob/main/src/numbertheory/chineseremainder.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_chinese_remainder.c`](https://github.com/stblake/mathilda/blob/main/tests/test_chinese_remainder.c)

## Notes & additional examples

### Notes

`ChineseRemainder[{r1, ..., rk}, {m1, ..., mk}]` inverts the residue map: it recovers
the unique `x` in `0 <= x < LCM[m1, ..., mk]` whose remainder modulo each `mi` is `ri`.
This is the constructive Chinese Remainder Theorem, and it is what makes a *residue
number system* usable -- a large integer can be carried as its residues against several
moduli, added and multiplied component-wise in parallel, and reassembled by
`ChineseRemainder` at the end.

Unlike the textbook statement, the implementation does not require the moduli to be
pairwise coprime. A solution to a pair of congruences `x ≡ ri (mod mi)` and
`x ≡ rj (mod mj)` exists if and only if `ri ≡ rj (mod GCD[mi, mj])`; when that holds
for every pair the congruences fuse into one modulo `LCM`, and when it fails the system
has no solution and the call is returned unevaluated. The optional offset `d` reports
the least solution at or above `d`, in `d <= x < d + LCM`, which is convenient when the
natural representative is wanted in a particular range.
