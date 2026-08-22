# Exponential Diophantine equations

In an **exponential Diophantine equation** the unknown appears in the exponent,
so the equation is not a polynomial at all — \(3^m - 2^n = 1\), \(2^n - 7 = x^2\).
These are exactly the equations that sympy's `diophantine` cannot express (see
the [performance page](performance.md)), and they are decided by some of the
deepest theorems in the subject.

```mathematica
In[1]:= Solve[3^m - 2^n == 1 && m > 0 && n > 0, {m, n}, Integers]
Out[1]= {{m -> 1, n -> 1}, {m -> 2, n -> 3}}
```

The only powers of \(3\) and \(2\) that differ by \(1\) are \(3^1 - 2^1 = 1\) and
\(3^2 - 2^3 = 9 - 8 = 1\).

---

## 1. Catalan's equation and Mihăilescu's theorem

Catalan conjectured in 1844 that \(8\) and \(9\) are the *only* consecutive
perfect powers — that \(x^a - y^b = 1\) with \(x, y, a, b > 1\) has the unique
solution \(3^2 - 2^3 = 1\). This resisted proof for 158 years until Preda
Mihăilescu settled it in 2002 using the theory of cyclotomic fields. Mathilda
knows the theorem:

```mathematica
In[1]:= Solve[x^a - y^b == 1 && x > 1 && y > 1 && a > 1 && b > 1, {x, y, a, b}, Integers]
Out[1]= {{x -> 3, y -> 2, a -> 2, b -> 3}}
```

The single tuple \((x, y, a, b) = (3, 2, 2, 3)\) is the entire solution set of the
Catalan equation — one of the rare Diophantine problems whose answer is a
*theorem with a name*. The same machinery handles fixed-base variants like the
\(3^m - 2^n = 1\) shown above, where Mihăilescu's result plus the small linear
cases give the complete list.

### References

1. E. Catalan, *Note extraite d'une lettre adressée à l'éditeur*, J. Reine Angew.
   Math. **27** (1844), 192.
2. P. Mihăilescu, *Primary cyclotomic units and a proof of Catalan's
   conjecture*, J. Reine Angew. Math. **572** (2004), 167–195.

---

## 2. The Ramanujan–Nagell equation \(2^n - 7 = x^2\)

In 1913 Ramanujan asked which \(2^n - 7\) are perfect squares; Nagell proved in
1948 that there are **exactly five**. It is the prototypical equation of the shape
\(x^2 + D = 2^n\). Mathilda solves it by factoring in the ring
\(\mathbb{Z}[(1 + \sqrt{-7})/2]\) (which has class number one), reducing the
problem to a Lucas-sequence term equal to \(\pm 1\); the **Bilu–Hanrot–Voutier**
primitive-divisor theorem then bounds \(n \le 32\), so a short finite scan returns
the whole set:

```mathematica
In[1]:= Solve[2^n - 7 == x^2 && n > 0 && x > 0, {n, x}, Integers]
Out[1]= {{n -> 3, x -> 1}, {n -> 4, x -> 3}, {n -> 5, x -> 5}, {n -> 7, x -> 11}, {n -> 15, x -> 181}}
```

The five solutions — the largest being \(2^{15} - 7 = 32761 = 181^2\) — are the
famous **Ramanujan–Nagell numbers**. The \(n = 15\) case is the reason the problem
is celebrated: it is far enough out that no short search would think to look
there, yet the primitive-divisor bound proves nothing lies beyond it.

Outside the class-number-one gate the general analysis (linear forms in
logarithms for an arbitrary base or \(D\)) is not built, so such inputs are
declined rather than half-answered:

```mathematica
In[2]:= Solve[3^n - 7 == x^2 && n > 0 && x > 0, {n, x}, Integers]
Out[2]= Solve[-7 + 3^n == x^2 && n > 0 && x > 0, {n, x}, Integers]
```

The equation returns unchanged — a reported gap, not a wrong or empty answer.

### References

1. S. Ramanujan, *Question 464*, J. Indian Math. Soc. **5** (1913), 130.
2. T. Nagell, *The Diophantine equation \(x^2 + 7 = 2^n\)*, Ark. Mat. **4**
   (1961), 185–187.
3. Y. Bilu, G. Hanrot, P. M. Voutier, *Existence of primitive divisors of Lucas
   and Lehmer numbers*, J. Reine Angew. Math. **539** (2001), 75–122.

---

Next: **[Famous results and impossibility proofs](famous-power-sums.md)** — the
showpieces the whole engine was built to reach: the Lander–Parkin and Frye
counterexamples, taxicab numbers, the Euler brick, and Fermat's Last Theorem.
