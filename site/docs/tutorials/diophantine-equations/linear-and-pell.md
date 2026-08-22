# Linear equations and Pell

The two oldest strands of Diophantine analysis are the **linear** equation,
solved by the extended Euclidean algorithm, and the **Pell** equation
\(x^2 - Dy^2 = 1\), solved by continued fractions. Both have infinitely many
solutions that fall on a regular pattern, and Mathilda returns that whole pattern
as a **parametric family** in an integer parameter `C[1]`, rather than a finite
list it could never finish printing.

```mathematica
In[1]:= Solve[x + 2 y == 5, {x, y}, Integers]
Out[1]= {{x -> 5 + 2 C[1], y -> -C[1]}}
```

Every integer solution of \(x + 2y = 5\) is obtained by choosing an integer for
`C[1]`. That single rule *is* the complete answer.

---

## 1. Linear equations

### 1.1 A single equation — the Euclidean staircase

A linear Diophantine equation \(a x + b y = c\) is solvable if and only if
\(\gcd(a, b)\) divides \(c\); when it does, the extended Euclidean algorithm
produces one particular solution and the general solution is that particular
solution plus integer multiples of \((b, -a)/\gcd(a,b)\). Substituting integers
for `C[1]` walks along the family:

```mathematica
In[1]:= Solve[x + 2 y == 5, {x, y}, Integers] /. C[1] -> 0
Out[1]= {{x -> 5, y -> 0}}

In[2]:= Solve[x + 2 y == 5, {x, y}, Integers] /. C[1] -> 3
Out[2]= {{x -> 11, y -> -3}}
```

When \(\gcd(a,b) \nmid c\), there is **no** solution — and Mathilda returns the
empty set as a *proof*, from the divisibility test alone:

```mathematica
In[3]:= Solve[6 x + 9 y == 5, {x, y}, Integers]
Out[3]= {}
```

Here \(\gcd(6, 9) = 3\) does not divide \(5\), so \(6x + 9y\) is always a
multiple of \(3\) and can never equal \(5\). The `{}` is that congruence
argument, not a failed search.

### 1.2 Systems — Hermite Normal Form

Several linear equations at once become a matrix equation \(A\mathbf{x} =
\mathbf{b}\) over \(\mathbb{Z}\). Mathilda reduces \(A\) to **Hermite Normal
Form** (the integer analogue of row echelon form, `src/linalg/hnf.c`): a
particular solution comes from back-substitution with an exact divisibility test
at each step, and the kernel lattice supplies the free parameters. An
underdetermined \(2\times 3\) system leaves a one-parameter family:

```mathematica
In[1]:= Solve[{x + 2 y + 3 z == 10, x - y + z == 2}, {x, y, z}, Integers]
Out[1]= {{x -> 18 + 5 C[1], y -> 8 + 2 C[1], z -> -8 - 3 C[1]}}

In[2]:= Solve[{x + 2 y + 3 z == 10, x - y + z == 2}, {x, y, z}, Integers] /. C[1] -> -3
Out[2]= {{x -> 3, y -> 2, z -> 1}}
```

`In[2]` picks the small representative \((3, 2, 1)\) out of the family — you can
check \(3 + 4 + 3 = 10\) and \(3 - 2 + 1 = 2\). If the divisibility test fails at
any HNF pivot the system has no integer solution and `{}` is returned as a proof;
this is why Mathilda declines to *integer-filter* a rational family, which could
manufacture a spurious empty set.

### 1.3 Positivity rays and huge coefficients

Two specialised linear paths are worth knowing. A homogeneous system with a
positivity constraint has solutions forming a **ray** — a single primitive
direction scaled by `C[1] >= 1`:

```mathematica
In[1]:= Solve[2 x == 3 y && 4 z == 3 y && x > 0 && y > 0 && z > 0, {x, y, z}, Integers]
Out[1]= {{x -> ConditionalExpression[6 C[1], C[1] >= 1], y -> ConditionalExpression[4 C[1], C[1] >= 1], z -> ConditionalExpression[3 C[1], C[1] >= 1]}}
```

The smallest positive solution is \((6, 4, 3)\), and every positive solution is a
positive multiple of it. Second, when a single linear equation is confined to a
large box but has *huge* coefficients, a naive scan is hopeless. Mathilda instead
enumerates the intersection of the solution line with the box using an
**LLL-reduced lattice**, so the work is proportional to the number of answers,
not the box volume:

```mathematica
In[2]:= Length[Solve[1000003 x + 999983 y == 7 && Abs[x] < 10^9 && Abs[y] < 10^9, {x, y}, Integers]]
Out[2]= 2000
```

All \(2000\) lattice points in a box of \(4\times10^{18}\) cells, found without
touching more than a handful of candidates.

### References

1. G. H. Hardy & E. M. Wright, *An Introduction to the Theory of Numbers*, 6th
   ed., Oxford, 2008 — Ch. 2 (continued fractions) and the linear congruence.
2. H. Cohen, *A Course in Computational Algebraic Number Theory*, Springer,
   1993 — §2.4 (Hermite Normal Form), §1.3 (extended gcd).
3. A. K. Lenstra, H. W. Lenstra, L. Lovász, *Factoring polynomials with rational
   coefficients*, Math. Ann. **261** (1982), 515–534 — the LLL algorithm.

---

## 2. The Pell equation

### 2.1 Bounded and unbounded solutions

Pell's equation \(x^2 - D y^2 = 1\) (for \(D > 0\) non-square) has infinitely many
solutions generated from a single **fundamental solution** \((x_1, y_1)\), which
Mathilda finds from the continued-fraction expansion of \(\sqrt D\) (the PQa
algorithm). With an upper bound present, it lists the orbit inside the box:

```mathematica
In[1]:= Solve[x^2 - 2 y^2 == 1 && x > 0 && y > 0 && x < 100, {x, y}, Integers]
Out[1]= {{x -> 3, y -> 2}, {x -> 17, y -> 12}, {x -> 99, y -> 70}}
```

The fundamental solution is \((3, 2)\); the next two, \((17, 12)\) and
\((99, 70)\), are its powers under \((3 + 2\sqrt2)^k\). Drop the bound and
Mathilda returns the *entire* orbit as one closed-form family:

```mathematica
In[2]:= Solve[x^2 - 2 y^2 == 1 && x > 0 && y > 0, {x, y}, Integers]
Out[2]= {{x -> ConditionalExpression[1/2 ((3 - 2 Sqrt[2])^C[1] + (3 + 2 Sqrt[2])^C[1]), C[1] >= 1], y -> ConditionalExpression[(1/2 ((3 + 2 Sqrt[2])^C[1] - (3 - 2 Sqrt[2])^C[1]))/Sqrt[2], C[1] >= 1]}}
```

That is the exact expression for the \(k\)-th solution, with `C[1]` playing the
role of \(k\). Substituting a value and simplifying recovers a concrete pair:

```mathematica
In[3]:= Simplify[Solve[x^2 - 2 y^2 == 1 && x > 0 && y > 0, {x, y}, Integers] /. C[1] -> 2]
Out[3]= {{x -> 17, y -> 12}}
```

### 2.2 Why continued fractions matter — a large fundamental solution

The reason continued fractions are indispensable, rather than a scan, is that the
fundamental solution can be *enormous* even for a small \(D\). The classic
example is \(D = 61\):

```mathematica
In[1]:= Solve[x^2 - 61 y^2 == 1 && x > 0 && y > 0 && x < 10^10, {x, y}, Integers]
Out[1]= {{x -> 1766319049, y -> 226153980}}
```

No search bounded by anything reasonable would reach \(x = 1\,766\,319\,049\); the
continued fraction produces it directly. (This is precisely the case where a
naive integer-search fallback times out — see the [performance
page](performance.md).)

### 2.3 Negative Pell and general \(N\)

The **negative Pell** equation \(x^2 - D y^2 = -1\) is solvable only for special
\(D\) — exactly when the continued-fraction period of \(\sqrt D\) is odd. Mathilda
decides this and, when solvable, returns the family:

```mathematica
In[1]:= Length[Solve[x^2 - 13 y^2 == -1 && x > 0 && y > 0, {x, y}, Integers]]
Out[1]= 1
```

For \(D = 13\) there is one solution class, built on the fundamental
\((18, 5)\): indeed \(18^2 - 13\cdot 5^2 = 324 - 325 = -1\). The **generalised
Pell** equation \(x^2 - D y^2 = N\) for arbitrary \(N\) is handled the same way,
by the Nagell–Lagrange–Matthews–Mollin method — one fundamental solution per
solution class, each grown by the \(D\)-unit:

```mathematica
In[2]:= Solve[x^2 - 2 y^2 == 7 && x > 0 && y > 0 && x < 100, {x, y}, Integers]
Out[2]= {{x -> 3, y -> 1}, {x -> 5, y -> 3}, {x -> 13, y -> 9}, {x -> 27, y -> 19}, {x -> 75, y -> 53}}
```

When the congruence conditions rule out every class, the equation has no solution
and `{}` is returned as a proof.

### References

1. H. W. Lenstra Jr., *Solving the Pell equation*, Notices of the AMS **49**
   (2002), 182–192.
2. T. Nagell, *Introduction to Number Theory*, Wiley, 1951 — Ch. VI, the
   generalised Pell equation.
3. K. Matthews, *The Diophantine equation \(x^2 - Dy^2 = N\), \(D > 0\)*, Expo.
   Math. **18** (2000), 323–331.

---

Next: **[Quadratic forms and conics](quadratic-forms.md)** — what happens when
the equation is a genuine binary or ternary quadratic, with cross terms, an
ellipse instead of a hyperbola, or three variables and Legendre's theorem to
decide solvability.
