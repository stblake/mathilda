# Famous results and impossibility proofs

This page is the showcase. The additive equations below — sums of like powers,
taxicab numbers, the Euler brick — are the problems whose *individual* solutions
made history, and several were first cracked by famous computer searches.
Mathilda solves each with a purpose-built method (meet-in-the-middle, modular
descent, divisor enumeration), and closes the page with equations whose answer is
the empty set: a proof that no solution exists at all.

```mathematica
In[1]:= Solve[x^5 + y^5 + z^5 + w^5 == r^5 && 0 < x < y < z < w < r < 700, {x, y, z, w, r}, Integers]
Out[1]= {{x -> 27, y -> 84, z -> 110, w -> 133, r -> 144}, {x -> 54, y -> 168, z -> 220, w -> 266, r -> 288}, {x -> 81, y -> 252, z -> 330, w -> 399, r -> 432}, {x -> 108, y -> 336, z -> 440, w -> 532, r -> 576}}
```

---

## 1. Lander and Parkin — Euler's conjecture falls

In 1769 Euler conjectured that a \(k\)-th power needs at least \(k\) other
\(k\)-th powers to sum to it — no three fourth powers add to a fourth power, no
four fifth powers to a fifth power. In 1966, using a CDC 6600, **Lander and
Parkin** found the counterexample that ended the fifth-power case:
\[
27^5 + 84^5 + 110^5 + 133^5 = 144^5.
\]
That is the first tuple in `In[1]` above; the other three are its \(2\times\),
\(3\times\), and \(4\times\) multiples inside the search box. Mathilda finds them
with an **ordering-aware meet-in-the-middle** search: it tabulates the partial
sums of one half of the variables and binary-searches for the complement,
carrying 128-bit integers so the fifth powers do not overflow. Tightening the box
to isolate the primitive solution returns it alone:

```mathematica
In[1]:= Solve[x^5 + y^5 + z^5 + w^5 == r^5 && 0 < x <= y <= z <= w < 150 && r < 200, {x, y, z, w, r}, Integers]
Out[1]= {{x -> 27, y -> 84, z -> 110, w -> 133, r -> 144}}
```

## 2. Frye's fourth-power counterexample

The fourth-power case of Euler's conjecture held out until 1988, when **Roger
Frye**, using Elkies' elliptic-curve analysis to localise the search, found the
smallest solution on a Connection Machine:
\[
95800^4 + 217519^4 + 414560^4 = 422481^4.
\]
Mathilda reaches it directly via **Frye's own mod-5 descent** — the divisibility
\(625 \mid (w^4 - C^4)\) plus secondary moduli prune the \((w, C)\) pairs — over a
narrow window around the known \(w\):

```mathematica
In[1]:= Solve[x^4 + y^4 + z^4 == w^4 && 0 < x < y < z < w && 422400 < w < 422500, {x, y, z, w}, Integers]
Out[1]= {{x -> 95800, y -> 217519, z -> 414560, w -> 422481}}
```

!!! note "Search windows"
    The narrow window keeps this example under a second. The full cold search over
    \(0 < w < 10^6\) is a single-core computation of roughly ten minutes — Frye's
    original run used a 65 536-processor machine. The [performance
    page](performance.md) puts the timing in context.

## 3. Taxicab numbers — two cubes, two ways

Ramanujan's taxicab number \(1729\) is the smallest integer expressible as a sum
of two cubes in two different ways: \(1729 = 1^3 + 12^3 = 9^3 + 10^3\). For a
*fixed* target, Mathilda uses the **divisor method** — for a sum of two cubes,
\(x + y\) divides the target, so the answer comes from the divisors of \(1729\),
not a two-dimensional scan:

```mathematica
In[1]:= Solve[x^3 + y^3 == 1729 && 0 < x <= y && x < 10^5, {x, y}, Integers]
Out[1]= {{x -> 1, y -> 12}, {x -> 9, y -> 10}}
```

Asked directly for two *distinct* representations, the same number falls out of a
meet-in-the-middle search:

```mathematica
In[2]:= Solve[x^3 + y^3 == z^3 + w^3 && 0 < x < y && 0 < z < w && x != z && x^3 + y^3 < 1730, {x, y, z, w}, Integers]
Out[2]= {{x -> 1, y -> 12, z -> 9, w -> 10}, {x -> 9, y -> 10, z -> 1, w -> 12}}
```

## 4. The Euler brick

An **Euler brick** is a rectangular box whose three edges *and* all three face
diagonals are integers: \(x^2 + y^2 = a^2\), \(x^2 + z^2 = b^2\), \(y^2 + z^2 =
c^2\). Mathilda solves the coupled system by **staged elimination** — only the
free edge variables are enumerated, and each diagonal is closed as an exact square
root per candidate:

```mathematica
In[1]:= Solve[x^2 + y^2 == a^2 && x^2 + z^2 == b^2 && y^2 + z^2 == c^2 && 0 < x < y < z < 250 && a > 0 && b > 0 && c > 0, {x, y, z, a, b, c}, Integers]
Out[1]= {{x -> 44, y -> 117, z -> 240, a -> 125, b -> 244, c -> 267}}
```

The smallest Euler brick, discovered by Paul Halcke in 1719, has edges
\((44, 117, 240)\).

## 5. Two more classics

The **Markov equation** \(x^2 + y^2 + z^2 = 3xyz\) generates the Markov triples,
and Mathilda counts thirteen with all entries \(\le 1000\); **Brocard's problem**
\(n! + 1 = m^2\) has only the three known *Brown numbers*:

```mathematica
In[1]:= Length[Solve[x^2 + y^2 + z^2 == 3 x y z && 0 < x <= y <= z <= 1000, {x, y, z}, Integers]]
Out[1]= 13

In[2]:= Solve[Factorial[n] + 1 == m^2 && n > 0 && m > 0 && n < 100, {n, m}, Integers]
Out[2]= {{n -> 4, m -> 5}, {n -> 5, m -> 11}, {n -> 7, m -> 71}}
```

\(4! + 1 = 25\), \(5! + 1 = 121\), \(7! + 1 = 5041 = 71^2\) — whether any fourth
Brown number exists is still open, so the bounded query is the honest statement:
these are all with \(n < 100\).

---

## 6. When the answer is a proof of impossibility

The deepest results are the ones where the solution set is *empty*, and Mathilda
returns `{}` only when it can prove it. The supreme example is **Fermat's Last
Theorem** — that \(x^n + y^n = z^n\) has no positive-integer solution for \(n \ge
3\), proved by Andrew Wiles in 1995. Mathilda recognises the form and returns the
proof instantly:

```mathematica
In[1]:= Solve[x^3 + y^3 == z^3 && x > 0 && y > 0 && z > 0, {x, y, z}, Integers]
Out[1]= {}
```

That single `{}` stands for one of the most celebrated theorems in mathematics.
The **Erdős–Straus / Egyptian-fraction** equation, by contrast, is solvable, and
Mathilda's reciprocal recursion returns every ordered unit-fraction decomposition
of \(1\) into three parts:

```mathematica
In[2]:= Solve[1 == 1/x + 1/y + 1/z && 0 < x <= y <= z, {x, y, z}, Integers]
Out[2]= {{x -> 2, y -> 3, z -> 6}, {x -> 2, y -> 4, z -> 4}, {x -> 3, y -> 3, z -> 3}}
```

The three ways to write \(1\) as a sum of three unit fractions:
\(\tfrac12 + \tfrac13 + \tfrac16\), \(\tfrac12 + \tfrac14 + \tfrac14\), and
\(\tfrac13 + \tfrac13 + \tfrac13\).

### References

1. L. J. Lander & T. R. Parkin, *Counterexample to Euler's conjecture on sums of
   like powers*, Bull. Amer. Math. Soc. **72** (1966), 1079.
2. R. Frye, *Finding \(95800^4 + 217519^4 + 414560^4 = 422481^4\) on the
   Connection Machine*, Proc. Supercomputing '88, IEEE, 1988, 106–116.
3. N. D. Elkies, *On \(A^4 + B^4 + C^4 = D^4\)*, Math. Comp. **51** (1988),
   825–835.
4. A. Wiles, *Modular elliptic curves and Fermat's Last Theorem*, Ann. of Math.
   **141** (1995), 443–551.

---

Next: **[How Mathilda compares](performance.md)** — the head-to-head against
sympy and PARI/GP, with the coverage table and the honest gaps.
