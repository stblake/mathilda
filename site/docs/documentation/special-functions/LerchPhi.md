# LerchPhi

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`LerchPhi[z, s, a]`**

is the Lerch transcendent Phi(z, s, a) = Sum\_{k\>=0} z^k/(k + a)^s.

**`Zeta[s, a] and z LerchPhi[z, s, 1] is PolyLog[s, z]. Exact reductions`**

<details>
<summary>Notes</summary>

It generalizes Zeta, HurwitzZeta and PolyLog: LerchPhi\[1, s, a\] is cover z = 0 (a^-s), s = 0 (1/(1-z)), z = +-1, positive integer a (a PolyLog form) and negative integer s (a rational function of z). The options DoublyInfinite -\> True (sum k from -Infinity to Infinity) and IncludeSingularTerm -\> True (keep the k + a = 0 term) are supported. Inexact arguments with |z| \< 1 evaluate numerically at machine or arbitrary (MPFR) precision; |z| \> 1 stays symbolic. Listable.

</details>

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= LerchPhi[z, s, 1]
Out[1]= PolyLog[s, z]/z

In[2]:= LerchPhi[0.5, 3, 2.5]
Out[2]= 0.0794983
```

### Applications (7)

At z = 1 the Lerch transcendent is a Hurwitz zeta

```mathematica
In[3]:= LerchPhi[1, 2, 1]
Out[3]= 1/6 Pi^2
```

Only the k = 0 term survives

```mathematica
In[4]:= LerchPhi[0, s, a]
Out[4]= a^(-s)
```

A geometric sum, independent of a

```mathematica
In[5]:= LerchPhi[z, 0, a]
Out[5]= 1/(1 - z)
```

The alternating harmonic series

```mathematica
In[6]:= LerchPhi[-1, 1, 1]
Out[6]= Log[2]
```

A negative integer s gives a rational function

```mathematica
In[7]:= LerchPhi[z, -1, 1]
Out[7]= 1/(1 - z)^2
```

Reduces to a polylogarithm

```mathematica
In[8]:= LerchPhi[2, 3, 1]
Out[8]= 1/2 PolyLog[3, 2]
```

```mathematica
In[9]:= N[LerchPhi[1/2, 2, 1], 20]
Out[9]= 1.16448105293002501181
```

## Algorithm

Mathilda -- the Lerch transcendent LerchPhi.

```text
  LerchPhi[z, s, a]   Phi(z, s, a) = Sum_{k>=0} z^k / (k + a)^s
                      (|z| < 1; analytic continuation elsewhere, branch cut
                       z in [1, Infinity))
```

For Re(a) < 0 the principal value uses the symmetric power ((k+a)^2)^(-s/2),

```text
and any term with k + a = 0 is excluded.  LerchPhi is the common
```

generalization of Zeta, HurwitzZeta and PolyLog:

```text
  Phi(1, s, a) = Zeta[s, a],      z Phi(z, s, 1) = PolyLog[s, z].
```

The evaluator routes each kind of argument to the cheapest exact or fastest numeric path (mirroring src/special_functions/{zeta,hurwitzzeta,polylog}.c):

```text
  exact reductions (any z, s, a):
      z = 0                  ->  a^-s                        (the k = 0 term)
      s = 0                  ->  1/(1 - z)                   (geometric sum)
      z = 1                  ->  Zeta[s, a]
      z = -1                 ->  2^-s (Zeta[s,a/2] - Zeta[s,(a+1)/2])
      a positive integer m   ->  z^-m (PolyLog[s,z] - Sum_{j<m} z^j j^-s)
      s negative integer -n  ->  (z d/dz + a)^n [1/(1-z)]    (rational in z)
  options:
      IncludeSingularTerm->True at a non-positive integer a  ->  ComplexInfinity
      DoublyInfinite->True   ->  Phi(z,s,a) + z^-1 Phi(1/z, s, 1-a)
  numeric (>= 1 inexact operand):
      |z| < 1                ->  complex-MPFR power series
      |z| > 1                ->  diverges, stays symbolic (no continuation)
  everything else            ->  stays symbolic (return NULL)
```

Attributes: Listable, NumericFunction, Protected.

## Implementation notes

**Algorithm.** `builtin_lerchphi` evaluates `LerchPhi[z, s, a] = Sum_{k>=0}
z^k/(k+a)^s`, the common generalization of `Zeta`, `HurwitzZeta` and `PolyLog`.
Exact reductions: `z = 0 -> a^-s`; `s = 0 -> 1/(1-z)`; `z = 1 -> Zeta[s, a]`;
`z = -1 -> 2^-s (Zeta[s, a/2] - Zeta[s, (a+1)/2])` (with special handling of
`s = 1` via a digamma difference, and of half-integer `a` via a Dirichlet-beta
closed form); `a` a positive integer `m -> z^-m (PolyLog[s, z] - Sum_{j<m}
z^j j^-s)`; `a` a non-positive integer reduces onto `PolyLog`; `s` a negative
integer `-n -> (z d/dz + a)^n [1/(1-z)]`, a rational function built with the
Euler operator and `Together`. Options `IncludeSingularTerm -> True` (at a
non-positive-integer `a` gives `ComplexInfinity`) and `DoublyInfinite -> True`
(`Phi(z,s,a) + z^-1 Phi(1/z,s,1-a)`). Numeric (at least one inexact operand):
`|z| < 1` (or `|z| = 1` with `Re s > 1`) sums the complex-MPFR power series
(symmetric power `((k+a)^2)^(-s/2)` for `Re a < 0`); `|z| > 1` off the cut with
`|Log z| < 2 pi` uses the Erdelyi large-`z` continuation through `HurwitzZeta`
and `Gamma`, else stays symbolic.

**Data structures.** `Expr`; `lcx` (`mpfr_t` re/im) MPFR-complex toolkit; the
exact and continuation paths reuse the `PolyLog`/`Zeta`/`HurwitzZeta`/`Gamma`
builtins. ND: N-ary kernel `NDKN_LerchPhi` (arity 3), element-wise over the `z`
buffer; `packed_aware` from kernel registration. Attributes: `Listable`,
`NumericFunction`, `Protected`.

**Complexity / limits.** Series term cap `wp*64 + 100000`; the large-`z`
continuation truncates at `NMAX = 2000` with an optimal-truncation safeguard.
`|z| > 1` is generally left symbolic outside the Erdelyi domain (integer `s` on
the cut is not implemented). `Compile[]` lowers at scalar shape
(`Compiled -> True`) but **not** at rank-1 array shape (`Compiled -> False`); the
ND kernel still serves the packed / visible-`NDArray` fast path at the REPL.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Zeta](../../special-functions/Zeta/), [HurwitzZeta](../../special-functions/HurwitzZeta/), [PolyLog](../../special-functions/PolyLog/)

- DLMF §25.14 — the Lerch transcendent.
- A. Erdelyi et al., Higher Transcendental Functions, Vol. I, §1.11(8) — the large-z continuation.
- Source: [`src/special_functions/lerchphi.c`](https://github.com/stblake/mathilda/blob/main/src/special_functions/lerchphi.c)
- Specification: [`docs/spec/builtins/special-functions.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/special-functions.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_lerchphi.c`](https://github.com/stblake/mathilda/blob/main/tests/test_lerchphi.c)
- Tests: [`tests/test_numeric_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_numeric_stress.c)

## Notes & additional examples

### Notes

`LerchPhi[z, s, a]` is the Lerch transcendent `Sum_{k>=0} z^k/(k+a)^s`, the
common generalization of `Zeta`, `HurwitzZeta` and `PolyLog`:
`LerchPhi[1, s, a] = Zeta[s, a]` and `z LerchPhi[z, s, 1] = PolyLog[s, z]`.

Many argument shapes reduce in closed form — `z = 0` gives `a^-s`, `s = 0` gives
`1/(1-z)`, a negative integer `s` gives a rational function of `z`, and a
positive integer `a` shifts onto `PolyLog`. Otherwise a numeric value is summed
from the defining series for `|z| < 1`; `|z| > 1` uses an analytic continuation
within its domain and is otherwise left symbolic. The options
`DoublyInfinite -> True` and `IncludeSingularTerm -> True` select the two-sided
sum and the singular-term convention.
