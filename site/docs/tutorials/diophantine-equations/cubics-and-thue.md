# Cubics: Mordell, three cubes, and Thue

At degree three the elementary tricks run out and Diophantine solving becomes
algebraic number theory in earnest. This page covers three landmark families:
integral points on **Mordell curves** \(y^2 = x^3 + k\), the **sum of three
cubes** \(x^3 + y^3 + z^3 = k\), and **Thue equations** \(F(x, y) = m\) for an
irreducible binary form \(F\) of degree \(\ge 3\). Each has finitely many
solutions, and each is decided by a genuine theorem rather than a scan.

```mathematica
In[1]:= Solve[y^2 == x^3 - 2, {x, y}, Integers]
Out[1]= {{x -> 3, y -> -5}, {x -> 3, y -> 5}}
```

Fermat claimed, and Euler proved, that \(y^2 = x^3 - 2\) has only the solutions
\((3, \pm 5)\) — \(27 - 2 = 25\). Mathilda reproduces the theorem.

---

## 1. Mordell curves \(y^2 = x^3 + k\)

For \(k < 0\) with \(|k|\) squarefree and a class-number condition met, the
equation factors in the quadratic ring \(\mathbb{Z}[\sqrt{k}]\): writing \(y^2 - k
= x^3\) as \((y - \sqrt{k})(y + \sqrt{k}) = x^3\) and using unique factorisation
of ideals forces \(y \pm \sqrt{k}\) to be a cube, which pins down \(x\) and \(y\)
completely. The result is the *complete* integral-point set, and a failed descent
is a proof of emptiness:

```mathematica
In[1]:= Solve[y^2 == x^3 - 5, {x, y}, Integers]
Out[1]= {}

In[2]:= Solve[y^2 == x^3 - 13, {x, y}, Integers]
Out[2]= {{x -> 17, y -> -70}, {x -> 17, y -> 70}}
```

\(y^2 = x^3 - 5\) has **no** integer points — the `{}` is the descent argument in
\(\mathbb{Z}[\sqrt{-5}]\), not an exhausted search — whereas \(y^2 = x^3 - 13\)
has exactly \((17, \pm 70)\), since \(17^3 - 13 = 4913 - 13 = 4900 = 70^2\).

### References

1. L. J. Mordell, *Diophantine Equations*, Academic Press, 1969 — the Mordell
   curve and the finiteness of integral points.
2. J. Gebel, A. Pethő, H. G. Zimmer, *Computing integral points on elliptic
   curves*, Acta Arith. **68** (1994), 171–192.

---

## 2. Sum of three cubes \(x^3 + y^3 + z^3 = k\)

The three-cubes problem is famously delicate: some \(k\) have no representation
at all, while others have one whose entries are astronomically large. Mathilda
attacks it from both ends — a global impossibility proof, and an efficient
bounded search.

### 2.1 A global obstruction: the mod-9 impossibility

Every cube is congruent to \(-1\), \(0\), or \(1\) modulo \(9\) (check
\(0^3,\dots,8^3\)). A sum of three cubes is therefore a sum of three values from
\(\{-1, 0, 1\}\), which lands in \(\{-3,\dots,3\} \equiv \{0,1,2,3,6,7,8\}
\pmod 9\) and can **never** be \(\equiv 4\) or \(5 \pmod 9\). So
\(x^3 + y^3 + z^3 = k\) with \(k \equiv \pm 4 \pmod 9\) has no integer solution
whatsoever — and Mathilda returns that as a proof, with **no bound required**:

```mathematica
In[1]:= Solve[x^3 + y^3 + z^3 == 4, {x, y, z}, Integers]
Out[1]= {}
```

The empty set here is the congruence argument, not an exhausted search: \(4
\equiv 4 \pmod 9\), so no triple of cubes can reach it. The same holds for every
\(k \equiv 4, 5 \pmod 9\) (and, since \(-k^3\) is also a cube residue, for the
signed forms \(x^3 + y^3 - z^3\) too).

### 2.2 Booker's bounded search

For a reachable residue, there is no shortcut to *finding* the solutions, and
Mathilda uses **Booker's method**: for each small value of the divisor \(d = |x
+ y|\), solve the cube-root congruence \(z^3 \equiv k \pmod d\) and recover the
candidates — far more efficient over a box than a triple loop.

```mathematica
In[1]:= Solve[x^3 + y^3 + z^3 == 3 && Abs[x] < 8000 && Abs[y] < 8000 && Abs[z] < 8000, {x, y, z}, Integers]
Out[1]= {{x -> -5, y -> 4, z -> 4}, {x -> 1, y -> 1, z -> 1}, {x -> 4, y -> -5, z -> 4}, {x -> 4, y -> 4, z -> -5}}
```

Within the box, \(k = 3\) has the well-known \((1,1,1)\) and the \((4,4,-5)\)
family. This is a *bounded* search, so here an empty answer is a proof about the
box, not the whole of \(\mathbb{Z}\):

```mathematica
In[2]:= Solve[x^3 + y^3 + z^3 == 42 && Abs[x] < 8000 && Abs[y] < 8000 && Abs[z] < 8000, {x, y, z}, Integers]
Out[2]= {}
```

There is no representation of \(42\) with all coordinates below \(8000\) — the
smallest, found by Booker and Sutherland in 2019, has 17-digit entries. Because
this search is intrinsically bounded, an *unbounded* query for a **reachable**
residue (one the mod-9 proof cannot dismiss) is declined rather than answered:

```mathematica
In[3]:= Solve[x^3 + y^3 + z^3 == 3, {x, y, z}, Integers]
Out[3]= Solve[x^3 + y^3 + z^3 == 3, {x, y, z}, Integers]
```

\(k = 3\) is \(\equiv 3 \pmod 9\), so no congruence rules it out, and without a
bound Mathilda has no complete procedure — it returns the input unchanged, an
honest gap exactly as the [second guarantee](index.md) promises, rather than a
possibly-incomplete list.

### References

1. A. R. Booker, *Cracking the problem with 33*, Research in Number Theory **5**
   (2019), 26.
2. A. R. Booker & A. V. Sutherland, *On a question of Mordell*, PNAS **118**
   (2021), e2022377118.

---

## 3. Thue equations \(F(x, y) = m\)

### 3.1 The Tzanakis–de Weger method

A **Thue equation** sets an irreducible homogeneous binary form of degree
\(\ge 3\) equal to a constant, e.g. \(x^3 - 2y^3 = 1\). Axel Thue proved in 1909
that such an equation has only finitely many solutions, but his proof was
ineffective. The **Tzanakis–de Weger** algorithm makes it effective:

1. work in the number field \(K = \mathbb{Q}(\theta)\), where \(\theta\) is a root
   of \(F(x, 1)\), and factor \(F(x, y) = a_0 \prod_i (x - \theta_i y)\);
2. reduce the equation to **unit equations** in the ring of integers \(O_K\),
   parameterising \(x - \theta y\) by the fundamental units;
3. bound \(\max(|x|, |y|)\) with **Baker's theory of linear forms in
   logarithms** — an astronomically large but *explicit* bound;
4. shrink that bound to a searchable range by **LLL lattice reduction**;
5. enumerate the reduced range and verify.

Mathilda implements the whole chain (`src/solve/solvethue.c`, with a real
number-field layer in `src/numbertheory/`), so an unbounded Thue equation returns
its complete, finite solution set:

```mathematica
In[1]:= Solve[x^3 - 2 y^3 == 1, {x, y}, Integers]
Out[1]= {{x -> -1, y -> -1}, {x -> 1, y -> 0}}

In[2]:= Solve[x^3 - 7 y^3 == 1, {x, y}, Integers]
Out[2]= {{x -> 1, y -> 0}, {x -> 2, y -> 1}}
```

Both are complete — \(x^3 - 2y^3 = 1\) has only \((1,0)\) and \((-1,-1)\), and
\(x^3 - 7y^3 = 1\) adds \((2, 1)\) since \(8 - 7 = 1\). When the form represents
\(m\) in no way at all, the empty set is again a proof:

```mathematica
In[3]:= Solve[x^3 - 2 y^3 == 5, {x, y}, Integers]
Out[3]= {}
```

### 3.2 Quartic and totally-complex forms

The engine is not limited to cubics. A quartic Thue form works the same way, and
when the field is **totally complex** (no real embeddings) Mathilda uses an
elementary rigorous bound \(|y| \le (|m| / \prod_i |\mathrm{Im}\,\theta_i|)^{1/n}\)
in place of Baker's machinery — faster, and enough to close the solution set:

```mathematica
In[1]:= Length[Solve[x^4 + y^4 == 17, {x, y}, Integers]]
Out[1]= 8

In[2]:= Solve[x^4 - 2 y^4 == -1, {x, y}, Integers]
Out[2]= {{x -> -1, y -> -1}, {x -> -1, y -> 1}, {x -> 1, y -> -1}, {x -> 1, y -> 1}}
```

\(x^4 + y^4 = 17\) has the eight sign-and-swap variants of \((1, 2)\) and
\((2, 1)\) — \(1 + 16 = 17\) — over the field \(\mathbb{Q}(\zeta_8)\); the
Ljunggren-type \(x^4 - 2y^4 = -1\) has the four \((\pm 1, \pm 1)\). Thue equations
Mathilda cannot yet reach (very large regulators, some non-monogenic fields, some
\(|m| \ne 1\)) are **declined**, never guessed — a comparison against PARI/GP's
`thue()` is on the [performance page](performance.md).

### References

1. A. Thue, *Über Annäherungswerte algebraischer Zahlen*, J. Reine Angew. Math.
   **135** (1909), 284–305 — the finiteness theorem.
2. N. Tzanakis & B. M. M. de Weger, *On the practical solution of the Thue
   equation*, J. Number Theory **31** (1989), 99–132.
3. A. Baker, *Linear forms in the logarithms of algebraic numbers*,
   Mathematika **13** (1966), 204–216.

---

Next: **[Exponential Diophantine equations](exponential.md)** — equations where
the unknown sits in the exponent, decided by Mihăilescu's theorem and the
Ramanujan–Nagell analysis.
