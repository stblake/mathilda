# Integration along a contour (line integrals)

`Method -> "LineIntegral"` evaluates a **complex contour integral** — the
integral of an analytic function along a path in the complex plane made of one or
more straight segments:

\[
\int_\gamma f(z)\,dz \;=\; \sum_k \int_{z_k}^{z_{k+1}} f(z)\,dz.
\]

You give the path as a list of complex points. Two points is a single segment;
more points make a **piecewise-linear contour**, and repeating the first point at
the end closes it into a loop.

```mathematica
In[1]:= Integrate[z, {z, 0, 1 + I}]
Out[1]= I
```

That integrates \(f(z) = z\) along the straight segment from \(0\) to \(1 + i\);
the antiderivative \(z^2/2\) gives \((1+i)^2/2 = i\). A `{z, z0, z1, …}` spec with
a complex endpoint auto-dispatches to this method, so you rarely need to name it —
but `Method -> "LineIntegral"` forces it, and `Integrate`LineIntegral[…]` is the
explicit entry point.

---

## 1. The mathematical idea

### Contour integrals

A contour integral sums the integrand along a path \(\gamma\). Mathilda
parametrises each straight segment \(z_k \to z_{k+1}\) by a real parameter
\(t \in [0, 1]\),

\[
\gamma(t) = z_k + t\,(z_{k+1} - z_k), \qquad \gamma'(t) = z_{k+1} - z_k,
\]

so the segment integral is \(\int_0^1 f(\gamma(t))\,\gamma'(t)\,dt\). This real
parametrisation is the key idea in the implementation: it turns questions about
the complex path (where are the singularities? what are the endpoint values?)
into questions about a real variable \(t\).

### Cauchy's theorem — path independence

If \(f\) is **analytic** (holomorphic) in a region, its contour integral depends
only on the endpoints, not the path between them. Two consequences drive most
examples:

- the integral over any **closed** loop of an analytic \(f\) is \(0\);
- along an open path, you may deform the path freely without changing the value.

### The Cauchy integral formula and residues

Analyticity fails at poles, and there the loop integral is *not* zero. The
prototype is

\[
\oint \frac{dz}{z} = 2\pi i
\]

around any loop that encircles the origin once counter-clockwise. More generally,
the **residue theorem** makes a closed contour integral equal to \(2\pi i\) times
the sum of the residues of \(f\) at the poles the contour encloses. Contour
integration is thus the constructive, path-based counterpart to the
[residue method](residue.md).

### Branch cuts

The antiderivative of a rational function carries logarithms and inverse
tangents, whose branch cuts the naive difference \(F(b) - F(a)\) can cross.
Mathilda computes the *continuous* change of \(F\) along a straight segment
branch-correctly: because a straight segment subtends an angle less than \(\pi\)
at any point off it, each logarithm of an affine argument can be combined into a
single principal \(\operatorname{Log}\) of a ratio
\(\operatorname{Log}[u(b)/u(a)]\). This is exactly what makes \(\oint dz/z\) come
out to \(2\pi i\) rather than \(0\).

### References

1. E. T. Whittaker, G. N. Watson, *A Course of Modern Analysis*, 4th ed.,
   Cambridge, 1927 — Ch. V–VI (Cauchy's theorem, the residue calculus).
2. L. V. Ahlfors, *Complex Analysis*, 3rd ed., McGraw-Hill, 1979 — Ch. 4.
3. J. E. Marsden, M. J. Hoffman, *Basic Complex Analysis*, 3rd ed., Freeman,
   1999 — Ch. 2 (contour integration and path independence).
4. NIST *Digital Library of Mathematical Functions*, §1.10 (functions of a
   complex variable), [dlmf.nist.gov/1.10](https://dlmf.nist.gov/1.10).

---

## 2. How Mathilda realises it

Each segment \(z_k \to z_{k+1}\) is handled in four steps:

1. **Antiderivative.** Compute `F = Integrate[f, z]` in the friendly real
   variable \(z\) (rational/elementary coefficients). If no antiderivative is
   found, the method declines.

2. **On-path singularities.** Substitute the parametrisation and locate the real
   roots \(t^* \in (0, 1)\) of the denominator of \(f(\gamma(t))\). A root is a
   singularity sitting *on* the segment, so the contour integral diverges —
   Mathilda reports `Integrate::idiv` and leaves it unevaluated.

3. **Continuous change of \(F\).** The segment value is \(F(b) - F(a)\), with
   endpoint values taken by the `Limit` engine as real one-sided limits in \(t\)
   when direct substitution is singular. When the segment crosses a branch cut,
   the value is recovered branch-correctly from the affine-\(\operatorname{Log}\)-
   ratio rule above.

4. **Numerical crosscheck.** Every symbolic segment value is cross-checked against
   a complex quadrature of \(f(\gamma(t))\,\gamma'(t)\). A symbolic candidate is
   accepted only when it agrees; an uncorrectable branch crossing leaves the
   integral unevaluated rather than return a wrong branch.

Step 4 makes this method **verified**, not merely correct-by-construction: unlike
the [Mellin](mellin.md) and [residue](residue.md) methods (which never call
`NIntegrate`), the line integrator confirms each segment numerically before
trusting it.

---

## 3. Invoking the method

There are three surface forms, all equivalent:

```mathematica
Integrate[f, {z, a, b}]                        (* auto: a or b non-real *)
Integrate[f, {z, z0, z1, ..., zn}]             (* piecewise-linear contour *)
Integrate`LineIntegral[f, {z, z0, ..., zn}]    (* explicit entry point *)
```

You can inspect where the path meets a singularity without evaluating the
integral, with `Integrate`PathSingularPoints`:

```mathematica
In[1]:= Integrate`PathSingularPoints[1/(z^2 + 1), {z, -2 + I, 2 + I}]
Out[1]= {I}
```

The horizontal segment from \(-2 + i\) to \(2 + i\) passes straight through the
pole at \(z = i\).

---

## 4. Introductory examples

### 4.1 Straight segments

Integrate \(z^2\) from the origin to \(1 + i\):

```mathematica
In[1]:= Integrate[z^2, {z, 0, 1 + I}]
Out[1]= -2/3 + 2/3*I
```

### 4.2 Path independence and Cauchy's theorem

Because \(z^2\) is analytic everywhere, bending the path into the two-segment
polyline \(0 \to 1 \to 1 + i\) gives the **same** value:

```mathematica
In[1]:= Integrate[z^2, {z, 0, 1, 1 + I}]
Out[1]= -2/3 + 2/3*I
```

And any *closed* loop of an analytic integrand integrates to zero — here a square,
and a triangle:

```mathematica
In[1]:= Integrate[z^2, {z, 1, I, -1, -I, 1}]
Out[1]= 0

In[2]:= Integrate[Exp[z], {z, 0, I, 1 + I, 1, 0}]
Out[2]= 0

In[3]:= Integrate[z^3 - 2 z, {z, 0, 2, 1 + I, 0}]
Out[3]= 0
```

---

## 5. Worked examples

The set below moves from Cauchy's integral formula through winding numbers,
branch-dependence, and the residue theorem, to a divergent contour.

**The Cauchy integral formula.** A counter-clockwise square around the origin
picks up the pole of \(1/z\), giving the archetypal \(2\pi i\):

```mathematica
In[1]:= Integrate[1/z, {z, 1, I, -1, -I, 1}]
Out[1]= (2*I) Pi
```

**Winding number.** Traverse the same square **twice** and the value doubles —
the contour now winds around the pole twice:

```mathematica
In[1]:= Integrate[1/z, {z, 1, I, -1, -I, 1, I, -1, -I, 1}]
Out[1]= (4*I) Pi
```

**Branch dependence.** Two *open* paths with the same endpoints \(1 \to -1\) can
give different values when they pass on opposite sides of the singularity at
\(0\). Going through the upper half-plane (via \(i\)) versus the lower (via
\(-i\)) flips the sign — exactly the branch jump of \(\operatorname{Log}\),
resolved correctly:

```mathematica
In[1]:= Integrate[1/z, {z, 1, I, -1}]
Out[1]= I Pi

In[2]:= Integrate[1/z, {z, 1, -I, -1}]
Out[2]= -I Pi
```

**The residue theorem.** A contour enclosing a single pole of \(1/(z^2+1)\)
returns \(2\pi i\) times that residue. The upper pole \(z = i\) has residue
\(1/(2i)\), so a square around it yields \(\pi\) (shown via `N` — the symbolic
form is a sum of logarithms equal to \(\pi\)):

```mathematica
In[1]:= N[Integrate[1/(z^2 + 1), {z, 1 + I/2, 1 + 3 I/2, -1 + 3 I/2, -1 + I/2, 1 + I/2}]]
Out[1]= 3.14159
```

Enlarge the square so it encloses **both** poles \(\pm i\) and their residues
cancel — the loop integral is zero:

```mathematica
In[1]:= N[Integrate[1/(z^2 + 1), {z, 2 - 2 I, 2 + 2 I, -2 + 2 I, -2 - 2 I, 2 - 2 I}]]
Out[1]= 0.0
```

A pole away from the origin behaves the same way — a rectangle around \(z = i\)
returns \(2\pi i\):

```mathematica
In[1]:= N[Integrate[1/(z - I), {z, 1, 1 + 2 I, -1 + 2 I, -1, 1}]]
Out[1]= 0.0 + 6.28319*I
```

**A divergent contour.** When a pole lies *on* the path — here \(1/z\) integrated
straight across the origin from \(-1\) to \(1\) — the integral genuinely diverges.
Mathilda detects the on-path singularity (via `Integrate`PathSingularPoints`,
which reports `{0}`) and leaves the integral unevaluated rather than return a
finite number:

```mathematica
In[1]:= Integrate`PathSingularPoints[1/z, {z, -1, 1}]
Out[1]= {0}

In[2]:= Integrate[1/z, {z, -1, 1}]
Out[2]= Integrate[1/z, {z, -1, 1}]
```

---

## 6. Divergence and verification

- **On-path poles diverge.** A singularity strictly inside a segment makes the
  contour integral divergent; the method emits `Integrate::idiv` and returns the
  integral unevaluated. Endpoint singularities, by contrast, are taken as real
  one-sided limits in the parameter \(t\).

- **Unverifiable branches are declined.** Every segment value is confirmed by a
  numerical quadrature. If a branch crossing cannot be corrected and the symbolic
  candidate disagrees with the quadrature, the integral is left unevaluated — the
  method never returns a value on the wrong branch.

---

## 7. Surface forms

| Form | Meaning |
|---|---|
| `Integrate[f, {z, a, b}]` | single segment \(a \to b\); auto-dispatches here when \(a\) or \(b\) is non-real |
| `Integrate[f, {z, z0, z1, …, zn}]` | piecewise-linear contour; repeat `z0` at the end to close it |
| `Integrate`LineIntegral[f, {z, …}]` | explicit entry point (same evaluation) |
| `Integrate`PathSingularPoints[f, {z, …}]` | list the singularities of \(f\) lying on the contour |

---

## 8. Limitations

- **An antiderivative must exist.** The method integrates \(f\) symbolically in
  \(z\) first, so it handles rational and elementary integrands. A non-analytic
  integrand such as `Conjugate[z]` has no complex antiderivative, and the integral
  is returned unevaluated:

  ```mathematica
  In[1]:= Integrate[Conjugate[z], {z, 0, 1 + I}]
  Out[1]= Integrate[Conjugate[z], {z, 0, 1 + I}]
  ```

- **Straight segments only.** Contours are piecewise-linear. A curved arc must be
  approximated by a polyline, or handled analytically through the
  [residue method](residue.md), which closes standard arcs symbolically.

- **Symbolic simplification of log sums.** A residue value at a pole off the
  origin is returned as a correct sum of logarithms; it may print un-reduced.
  Apply `N` for the numerical value — but note that `Simplify` can mis-combine
  these principal-branch logarithms, so trust `N` (which matches the internal
  quadrature) over a symbolic re-simplification here.

---

## See also

- [Residue theorem](residue.md) — closes the standard *arcs* (semicircle, wedge,
  keyhole) symbolically; the constructive counterpart to this path-based method.
- [Mellin transforms](mellin.md) — the third contour method, for half-line
  integrals.
- [Calculus tutorial](../08-calculus.md) — `Integrate`, `D`, `Series`, `Limit`.
- `?Integrate`LineIntegral` in the REPL for the built-in help string.
