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

A trivial value — Euler's solution of the Basel problem, ζ(2) = π²/6:

```mathematica
In[1]:= Zeta[2]
Out[1]= 1/6 Pi^2
```

Even positive integers are rational multiples of powers of π (via Bernoulli numbers):

```mathematica
In[1]:= Zeta[6]
Out[1]= 1/945 Pi^6
```

Negative-integer values are rational, exposing the Bernoulli-number connection ζ(−n) = −Bₙ₊₁/(n+1); the negative even integers are the trivial zeros:

```mathematica
In[1]:= Zeta[-1]
Out[1]= -1/12

In[1]:= Zeta[-3]
Out[1]= 1/120

In[1]:= Table[Zeta[-2 n], {n, 1, 4}]
Out[1]= {0, 0, 0, 0}
```

Apéry's constant ζ(3) to 40 digits, evaluated through the MPFR kernel:

```mathematica
In[1]:= N[Zeta[3], 40]
Out[1]= 1.2020569031595942853997381615114499907651
```

A value on the critical line at the first nontrivial zero ρ ≈ 1/2 + 14.134725 i — the result is numerically zero to within the input precision, confirming the zero of the Riemann hypothesis:

```mathematica
In[1]:= N[Zeta[1/2 + 14.134725 I], 10]
Out[1]= 1.7674298414e-08 - 1.1102028931e-07*I
```

The Laurent expansion about the pole s = 1 defines the Stieltjes constants, with residue 1 and constant term the Euler–Mascheroni constant γ:

```mathematica
In[1]:= Series[Zeta[s], {s, 1, 2}]
Out[1]= 1/(s - 1) + EulerGamma + -StieltjesGamma[1] (s - 1) + 1/2 StieltjesGamma[2] (s - 1)^2 + O[s - 1]^3
```

The Hurwitz zeta at an integer second argument reduces to ζ(s) minus a finite power sum, mixing exact π-power and rational parts:

```mathematica
In[1]:= Zeta[4, 5]
Out[1]= -22369/20736 + 1/90 Pi^4
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

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## See also

[HurwitzZeta](../../special-functions/HurwitzZeta/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_autocompile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_autocompile.c)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_findroot.c`](https://github.com/stblake/mathilda/blob/main/tests/test_findroot.c)
- Tests: [`tests/test_flint_bridge.c`](https://github.com/stblake/mathilda/blob/main/tests/test_flint_bridge.c)
