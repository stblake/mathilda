# Quadratic forms and conics

A degree-two Diophantine equation in two or three variables is a **quadratic
form**, and its integer solutions live on a conic. Pell (the previous page) is
the special hyperbola \(x^2 - Dy^2 = 1\); this page covers the rest — general
binary forms with cross and linear terms, ellipses, and the three-variable
(ternary) case where the decisive tool is a theorem of Legendre rather than a
search.

```mathematica
In[1]:= Solve[n^2 + n + 41 == y^2 && n > 0 && y > 0, {n, y}, Integers]
Out[1]= {{n -> 40, y -> 41}}
```

That is Euler's famous prime-generating polynomial \(n^2 + n + 41\), asked in
reverse: for which \(n\) is it a perfect square? Only \(n = 40\), where it equals
\(41^2\).

---

## 1. Binary forms: complete the square, then factor

### 1.1 Conics with a square leading coefficient

When \(y^2 = A x^2 + B x + C\) has a perfect-square \(A\), completing the square
turns it into a **difference of two squares** \((2pY)^2 - (2Ax + B)^2 = \Delta\),
and every solution corresponds to a factorisation of the constant \(\Delta\).
Mathilda enumerates the divisors of \(|\Delta|\) — a finite, exhaustive process,
so an empty result is a proof. Euler's conic above is the archetype: \(y^2 - (n^2
+ n + 41)\) rearranges to \((2y)^2 - (2n+1)^2 = 163\), and the only positive
factorisation of \(163\) (a prime) gives \(n = 40\).

### 1.2 Factorable forms — Runge's method

If the quadratic part \(A x^2 + B x y + C y^2\) has a discriminant \(B^2 - 4AC\)
that is a positive perfect square, the form **factors into two rational lines**,
and the equation \(= m\) becomes a product of two integer linear factors equal to
\(m\). Enumerating the divisor pairs of \(m\) solves it completely:

```mathematica
In[1]:= Solve[2 x^2 + 3 x y - 2 y^2 == 7, {x, y}, Integers]
Out[1]= {{x -> -3, y -> 1}, {x -> 3, y -> -1}}
```

Here \(2x^2 + 3xy - 2y^2 = (2x - y)(x + 2y)\), so the equation is \((2x - y)(x +
2y) = 7\); the four divisor pairs of \(7\) yield these two integer points.

### 1.3 Ellipses — a finite interval to scan

When the discriminant is *negative* the conic is an **ellipse**, which encloses
only finitely many lattice points. Mathilda treats the equation as a quadratic in
\(x\) for each fixed \(y\), over exactly the finite interval of \(y\) where the
\(x\)-discriminant stays non-negative, and solves that integer quadratic
exactly — rotations and linear terms included:

```mathematica
In[1]:= Solve[x^2 + x y + y^2 == 3, {x, y}, Integers]
Out[1]= {{x -> -2, y -> 1}, {x -> -1, y -> -1}, {x -> -1, y -> 2}, {x -> 1, y -> -2}, {x -> 1, y -> 1}, {x -> 2, y -> -1}}
```

All six lattice points on the tilted ellipse \(x^2 + xy + y^2 = 3\). Because the
scan of the \(y\)-interval is exhaustive, an ellipse with no lattice points
returns `{}` as a proof.

### References

1. T. Andreescu & D. Andrica, *Quadratic Diophantine Equations*, Springer, 2015.
2. D. Alpern, *Methods to solve \(ax^2 + bxy + cy^2 + dx + ey + f = 0\)* — the
   discriminant classification Mathilda's binary-form router follows.

---

## 2. Ternary quadratics and Legendre's theorem

### 2.1 Deciding solvability without searching

A homogeneous ternary quadratic such as \(x^2 + y^2 = c\,z^2\) either has only the
trivial solution \((0,0,0)\) or infinitely many — and **Legendre's theorem**
decides which by congruence conditions (quadratic residues / Hilbert symbols),
with no search at all. When the conditions fail, the *only* integer solution is
the trivial one, and Mathilda reports exactly that:

```mathematica
In[1]:= Solve[x^2 + y^2 == 3 z^2, {x, y, z}, Integers]
Out[1]= {{x -> 0, y -> 0, z -> 0}}
```

There is no nonzero way to write \(3z^2\) as a sum of two squares, because a
number \(\equiv 3 \pmod 4\) is never a sum of two squares; the single trivial
tuple is Legendre's obstruction made explicit — a proof, not an empty guess.

### 2.2 When it *is* solvable — the complete parametric family

Contrast the Pythagorean form \(x^2 + y^2 = z^2\), which is solvable. Here
Mathilda returns the **entire** solution set as a two-parameter family — the
classical \((m^2 - n^2,\, 2mn,\, m^2 + n^2)\) parametrisation, together with its
sign and coordinate-swap variants and the degenerate axis solutions:

```mathematica
In[1]:= Solve[x^2 + y^2 == z^2, {x, y, z}, Integers]
Out[1]= {{x -> -2 C[1] C[2] C[3], y -> C[3] (C[1]^2 - C[2]^2), z -> C[3] (C[1]^2 + C[2]^2)}, {x -> C[3] (C[1]^2 - C[2]^2), y -> -2 C[1] C[2] C[3], z -> C[3] (C[1]^2 + C[2]^2)}, {x -> -2 C[1] C[2] C[3], y -> C[3] (C[2]^2 - C[1]^2), z -> C[3] (C[1]^2 + C[2]^2)}, {x -> C[3] (C[1]^2 - C[2]^2), y -> 2 C[1] C[2] C[3], z -> C[3] (C[1]^2 + C[2]^2)}, {x -> 2 C[1] C[2] C[3], y -> C[3] (C[1]^2 - C[2]^2), z -> C[3] (C[1]^2 + C[2]^2)}, {x -> C[3] (C[2]^2 - C[1]^2), y -> -2 C[1] C[2] C[3], z -> C[3] (C[1]^2 + C[2]^2)}, {x -> 2 C[1] C[2] C[3], y -> C[3] (C[2]^2 - C[1]^2), z -> C[3] (C[1]^2 + C[2]^2)}, {x -> C[3] (C[2]^2 - C[1]^2), y -> 2 C[1] C[2] C[3], z -> C[3] (C[1]^2 + C[2]^2)}, {x -> 0, y -> C[1], z -> C[1]}, {x -> C[1], y -> 0, z -> C[1]}, {x -> 0, y -> -C[1], z -> C[1]}, {x -> -C[1], y -> 0, z -> C[1]}}
```

The three parameters are the two generating integers \(m = C[1]\), \(n = C[2]\)
and an overall scale \(C[3]\). Picking values from the first non-degenerate
branch produces an ordinary triple:

```mathematica
In[2]:= Simplify[Solve[x^2 + y^2 == z^2, {x, y, z}, Integers][[5]] /. {C[1] -> 2, C[2] -> 1, C[3] -> 1}]
Out[2]= {x -> 4, y -> 3, z -> 5}
```

\(m = 2\), \(n = 1\) gives the \((3, 4, 5)\) triangle; scaling \(C[3]\) and
varying \(m, n\) reaches every Pythagorean triple exactly once per branch.

### References

1. A.-M. Legendre, *Théorie des nombres*, 1798 — the solvability theorem for
   ternary quadratic forms.
2. J. E. Cremona & D. Rusin, *Efficient solution of rational conics*, Math. Comp.
   **72** (2003), 1417–1441.
3. H. Cohen, *A Course in Computational Algebraic Number Theory*, Springer,
   1993 — §5.3 (conics and ternary forms).

---

## 3. Systems that reduce to a conic

A quadratic combined with a *linear* constraint can often be collapsed to a
one-variable problem. The **Pythagorean triangle of fixed perimeter** is the
model case: \(x^2 + y^2 = z^2\) with \(x + y + z = p\). Eliminating \(z\) leaves a
bilinear equation whose solutions come from the divisors of a single number:

```mathematica
In[1]:= Solve[x^2 + y^2 == z^2 && x + y + z == 12 && 0 < x < y && z > 0, {x, y, z}, Integers]
Out[1]= {{x -> 3, y -> 4, z -> 5}}

In[2]:= Solve[x^2 + y^2 == z^2 && x + y + z == 1000 && 0 < x < y && z > 0, {x, y, z}, Integers]
Out[2]= {{x -> 200, y -> 375, z -> 425}}
```

The only right triangle with integer sides and perimeter \(12\) is \((3,4,5)\);
the only one with perimeter \(1000\) is \((200, 375, 425)\). No enumeration of
triangles is needed — the linear equation turns the search into a divisor
listing.

---

Next: **[Cubics: Mordell, three cubes, and Thue](cubics-and-thue.md)** — degree
three and above, where the methods become genuine algebraic number theory:
factorisation in quadratic rings, cube roots modulo a divisor, and the
Tzanakis–de Weger algorithm for Thue equations.
