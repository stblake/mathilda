# Zeta

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Zeta[s]`**

is the Riemann zeta function zeta(s) = Sum\_{k\>=1} k^-s.

**`Zeta[s, a]`**

is the Hurwitz zeta function zeta(s, a) = Sum\_{k\>=0} (k+a)^-s.

<details>
<summary>Notes</summary>

Even positive integers give rational multiples of Pi^(2n), negative integers give rationals, Zeta\[0\] is -1/2, and Zeta\[1\] is ComplexInfinity; odd positive integers stay symbolic. Hurwitz zeta at a positive integer a reduces to Zeta\[s\] minus a finite power sum, and Zeta\[s, 1/2\] is (2^s - 1) Zeta\[s\]. Real, complex, machine and arbitrary-precision (MPFR) numeric arguments evaluate numerically via mpfr\_zeta (real Riemann) or an Euler-Maclaurin kernel. Listable.

</details>

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= Zeta[2]
Out[1]= 1/6 Pi^2

In[2]:= Series[Zeta[x], {x, 1, 2}] // Normal
Out[2]= EulerGamma + 1/(-1 + x) - StieltjesGamma[1] (-1 + x) + 1/2 StieltjesGamma[2] (-1 + x)^2
```

### Applications (9)

```mathematica
In[3]:= Zeta[2]
Out[3]= 1/6 Pi^2

In[4]:= Zeta[6]
Out[4]= 1/945 Pi^6

In[5]:= Zeta[-1]
Out[5]= -1/12

In[6]:= Zeta[-3]
Out[6]= 1/120

In[7]:= Table[Zeta[-2 n], {n, 1, 4}]
Out[7]= {0, 0, 0, 0}

In[8]:= N[Zeta[3], 40]
Out[8]= 1.2020569031595942853997381615114499907651

In[9]:= N[Zeta[1/2 + 14.134725 I], 10]
Out[9]= 1.76743e-08 - 1.1102e-07*I

In[10]:= Series[Zeta[s], {s, 1, 2}]
Out[10]= 1/(s - 1) + EulerGamma + -StieltjesGamma[1] (s - 1) + 1/2 StieltjesGamma[2] (s - 1)^2 + O[s - 1]^3

In[11]:= Zeta[4, 5]
Out[11]= -22369/20736 + 1/90 Pi^4
```

## Algorithm

Mathilda -- the Riemann and Hurwitz zeta functions.

```text
  Zeta[s]      Riemann zeta      zeta(s)   = Sum_{k>=1} k^-s          (Re s > 1)
  Zeta[s, a]   Hurwitz zeta      zeta(s,a) = Sum_{k>=0} (k+a)^-s      (Re s > 1)
```

Both are defined elsewhere by analytic continuation; the evaluator routes each kind of argument to the cheapest exact or fastest numeric path:

```text
  exact integer s         ->  closed form:
                                s = 1        : ComplexInfinity (pole)
                                s = 0        : -1/2
                                s = 2n > 0   : rational * Pi^(2n)   (Bernoulli)
                                s = -m < 0   : rational            (Bernoulli)
                                s = 2n+1 > 0 : stays symbolic (no closed form)
  exact Hurwitz, integer a -> Zeta[s] - Sum_{k=1}^{a-1} k^-s
  machine / MPFR real s    -> MPFR mpfr_zeta (Riemann only)
  complex s, or any a != 1 -> Euler-Maclaurin complex-MPFR kernel
  everything else          -> stays symbolic (return NULL)
```

MPFR provides mpfr_zeta for real Riemann zeta only -- it has no Hurwitz and no complex zeta -- so the Hurwitz / complex kernel is implemented here from the Euler-Maclaurin summation formula (DLMF 25.11.5):

```text
  zeta(s,a) = Sum_{k=0}^{N-1} (a+k)^-s
            + (a+N)^(1-s)/(s-1)
            + 1/2 (a+N)^-s
            + Sum_{j>=1} B_{2j}/(2j)! (s)_{2j-1} (a+N)^(-s-2j+1)
```

with (s)_{2j-1} the rising factorial. N is chosen from the working precision and |s|; the correction series is truncated at its optimal (smallest) term.

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Algorithm.** `builtin_zeta` evaluates `Zeta[s]` (Riemann) and `Zeta[s, a]`
(Hurwitz). Exact integer `s`: `1 -> ComplexInfinity` (pole), `0 -> -1/2`, even
`2n > 0 -> rational · Pi^(2n)` (via exact Bernoulli), negative `-m -> rational`
(Bernoulli; `0` at even `m`), odd `2n+1 > 0` left symbolic; `s = Infinity -> 1`.
Hurwitz exact: `a = 1 -> Zeta[s]`; `a = 1/2 -> (2^s - 1) Zeta[s]`; a positive
integer `a -> Zeta[s] - Sum_{k=1}^{a-1} k^-s`. Numeric: real Riemann zeta via
`mpfr_zeta`; complex `s`, or any `a != 1`, via a Euler-Maclaurin complex-MPFR
kernel (DLMF 25.11.5) — head sum plus the `(a+N)^(1-s)/(s-1)`, `(1/2)(a+N)^-s`
and Bernoulli correction terms, with `N` chosen from the working precision and
`|s|` and the correction series truncated at its optimal (smallest) term; the
two-argument form uses the *symmetric* power `((a+k)^2)^(-s/2)` for `Re a < 0`.
Intervals route through `interval_apply_function`.

**Data structures.** `Expr`; an exact `mpq` Bernoulli cache; a local `zcx`
(`mpfr_t` re/im) toolkit for the Euler-Maclaurin kernel; GMP; `mpfr_zeta` for
real Riemann. ND: unary kernel `NDKU_Zeta = { NULL, ndk_Zeta_r, ... }` (real
Riemann, via `sf_machine_zeta`), registered `REG_U`, so `packed_aware`. The
two-argument Hurwitz form has **no** ND buffer kernel on `Zeta` — the separate
`HurwitzZeta` head carries the binary kernel `NDKB_HurwitzZeta`. Attributes:
`Listable`, `NumericFunction`, `Protected`.

**Complexity / limits.** Exact-integer cap `ZETA_EXACT_INT_CAP = 10000`;
Hurwitz-integer-`a` cap `100000`. The Euler-Maclaurin head term count is
`N ~ digits + |s|`, the correction series optimally truncated. Real zeta goes
only through `mpfr_zeta`. `Compile[]` lowers `Zeta[s]` at both scalar and rank-1
array shapes (`Compiled -> True`).

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [HurwitzZeta](../../special-functions/HurwitzZeta/)

- DLMF §25.2 — the Riemann zeta function; §25.11 — the Hurwitz zeta function.
- DLMF §25.11.5 — the Euler-Maclaurin summation formula used for the Hurwitz / complex kernel.
- Source: [`src/special_functions/zeta.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/zeta.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_elliptic.c`](https://github.com/stblake/mathilda/blob/main/tests/test_elliptic.c)
- Tests: [`tests/test_findroot.c`](https://github.com/stblake/mathilda/blob/main/tests/test_findroot.c)
