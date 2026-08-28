---
status: Stable
references:
  - "Knuth, \"The Art of Computer Programming, Vol. 2: Seminumerical Algorithms\", §4.3.2 (modular arithmetic and the CRT)."
  - "Ding, Pei & Salomaa, \"Chinese Remainder Theorem: Applications in Computing, Coding, Cryptography\" (1996)."
  - "von zur Gathen & Gerhard, \"Modern Computer Algebra\", §5.4 (the extended Euclidean algorithm and CRT)."
---
### Worked examples

The two-congruence classic: the smallest non-negative `x` with `x ≡ 2 (mod 3)` and `x ≡ 3 (mod 5)`.

```mathematica
In[1]:= ChineseRemainder[{2, 3}, {3, 5}]
Out[1]= 8
```

Any number of congruences at once; the answer lies in `0 <= x < LCM` of the moduli.

```mathematica
In[1]:= ChineseRemainder[{1, 2, 3}, {2, 3, 5}]
Out[1]= 23
```

A third argument `d` shifts the search window, returning the smallest solution `>= d`.

```mathematica
In[1]:= ChineseRemainder[{3, 4}, {5, 7}, 2]
Out[1]= 18
```

The moduli need **not** be pairwise coprime. When they overlap, a solution exists exactly when the congruences agree on the shared part:

```mathematica
In[1]:= ChineseRemainder[{2, 8}, {6, 10}]
Out[1]= 8
```

and the call stays unevaluated when the system is inconsistent, rather than inventing a false answer:

```mathematica
In[1]:= ChineseRemainder[{1, 2}, {6, 10}]
Out[1]= ChineseRemainder[{1, 2}, {6, 10}]
```

Everything is exact through GMP, so the residues and moduli may be arbitrary-precision integers -- the reconstruction step of a residue number system:

```mathematica
In[1]:= ChineseRemainder[{123, 456}, {10^9 + 7, 10^9 + 9}]
Out[1]= 499999841499998989
```

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
