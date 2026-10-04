# NLimit

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NLimit[expr, z -> z0]`**

numerically finds the limiting value of expr as z approaches z0.

<details>
<summary>Notes</summary>

A geometric sequence of sample points approaching z0 is constructed (z0 may be finite, complex, or an infinite point such as Infinity or I Infinity) and the limit is recovered by sequence acceleration. Method -\> Automatic (default) runs both Richardson/Romberg and Wynn's epsilon and keeps the most self-consistent estimate, so branch-point / fractional-power approaches (which Richardson alone mishandles) are resolved accurately (Levin's u-transform also participates when the samples are settling). Method -\> EulerSum forces Richardson/Romberg; Method -\> SequenceLimit forces Wynn's epsilon; Method -\> "Levin" forces Levin's transformation ("LevinU" | "LevinT" | "LevinV" select the u/t/v variant). expr must be numerical when z is numerical. Small spurious residuals are not recognised as zero -- Chop if needed. An expression whose sampled values oscillate with a non-decaying envelope has no limit: NLimit::osc is issued and the form is returned unevaluated, rather than reporting the meaningless extrapolant. Options: Method (Automatic | EulerSum | SequenceLimit | "Levin"), WorkingPrecision (default MachinePrecision), AccuracyGoal (default MachinePrecision), PrecisionGoal, Direction (Automatic == -1, or a complex approach vector), Scale (initial step / distance, default 1), Terms (starting extrapolation depth, default 13, grown adaptively up to meet AccuracyGoal -- the depth, not WorkingPrecision, sets the accuracy on branch-point / fractional-power approaches), WynnDegree (SequenceLimit iterations, default 1). Each method is also callable directly as NLimit\`m\[expr, z -\> z0\]: NLimit\`Automatic, NLimit\`EulerSum, NLimit\`SequenceLimit, NLimit\`Levin, NLimit\`LevinU, NLimit\`LevinT, NLimit\`LevinV.

</details>

## Examples (16)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= NLimit[Sin[x]/x, x -> 0]
Out[1]= 1.0

In[2]:= NLimit[(1 + 1/n)^n, n -> Infinity]
Out[2]= 2.71828

In[3]:= NLimit[(1 + I/x)^x, x -> Infinity]
Out[3]= 0.540302 + 0.841471*I

In[4]:= NLimit[Tanh[Pi x]/(1 + x^2), x -> I] // Chop
Out[4]= 0.0 - 1.5708*I

In[5]:= NLimit[x Sin[x], x -> Infinity] NLimit::osc: The sampled values oscillate with a non-decaying envelope; ...
Out[5]= RepeatedNull[Null]
```

Same oscillation, decaying envelope

```mathematica
In[6]:= NLimit[Sin[x]/x, x -> Infinity]
Out[6]= -0.000155258
```

### Options (5)

```mathematica
In[7]:= NLimit[(10^x - 1)/x, x -> 0, Terms -> 10, Method -> SequenceLimit]
Out[7]= 2.30259

In[8]:= NLimit[z + Conjugate[z]/z, z -> 0, Direction -> -I] // Chop
Out[8]= -1.0

In[9]:= NLimit[Tan[z], z -> Infinity I, Method -> SequenceLimit] // Chop
Out[9]= 0.0 + 1.0*I

In[10]:= NLimit[(2^x - 1)/x, x -> 0, WorkingPrecision -> 30, Terms -> 14]
Out[10]= 0.6931471805599453094172321284473

In[11]:= NLimit[Sin[x]/x, x -> 0, Method -> "Levin"]
Out[11]= 1.0
```

### Applications (5)

```mathematica
In[12]:= NLimit[Sin[x]/x, x -> 0]
Out[12]= 1.0

In[13]:= NLimit[(1 + 1/n)^n, n -> Infinity]
Out[13]= 2.71828

In[14]:= NLimit[n (2^(1/n) - 1), n -> Infinity]
Out[14]= 0.693147

In[15]:= NLimit[Zeta[x] - 1/(x - 1), x -> 1]
Out[15]= 0.577216

In[16]:= NLimit[Sin[x]/x, x -> 0, Method -> "Levin"]
Out[16]= 1.0
```

## Algorithm

nlimit.c — NLimit[expr, z -> z0, opts]

```text
Numerically estimates  lim_{z -> z0} expr  by sampling expr on a geometric
```

sequence of points approaching z0 and accelerating the resulting sequence.

```text
  Sample points
  -------------
  Finite z0:    z_k = z0 - d * Scale * 2^-k          (k = 0 .. Terms-1)
                d is the (unit) Direction vector; the points sit on the -d
                side of z0, i.e. one moves *along* d to reach z0.  The
                default Direction -> Automatic == -1 places the points at
                z0 + Scale*2^-k, approaching "from larger values".
  Infinite z0:  z_k = u * Scale * 2^k                (k = 0 .. Terms-1)
                u is the unit ray of the infinite limit point (the direction
                of Infinity / I Infinity / DirectedInfinity[..]); the points
                march outward to infinity along that ray.

  Acceleration
  ------------
  Method -> Automatic (default): run BOTH Richardson and Wynn's epsilon (at
    every admissible degree) and return the estimate whose internal
    convergence residual is smallest.  Richardson models an integer-power
    (analytic) error tail; its fixed 2^j-1 denominators cannot annihilate a
    geometric or fractional-power (branch-point) tail — e.g. the sqrt(step)
    imaginary part of 2 ArcTan[Sqrt[(1+x)/(1-x)]] approaching x -> 1 from the
    larger-values side.  Wynn's epsilon captures exactly those tails.
    Selecting by best self-consistency (smallest step) picks the right tool
    per problem, so smooth limits keep Richardson's accuracy while
    algebraic/branch approaches gain Wynn's.

  Method -> EulerSum: Richardson / Romberg extrapolation of the sequence S_k,
    treated as a function of the geometric step.  Following the same
    convention as ND's "EulerSum" (see nderiv.c) the tableau uses the
    all-powers denominator 2^j - 1:
        T(i,0) = S_i,
        T(i,j) = T(i,j-1) + (T(i,j-1) - T(i-1,j-1)) / (2^j - 1),
        result = T(Terms-1, Terms-1).

  Method -> SequenceLimit: Wynn's epsilon algorithm (iterated Shanks /
    Aitken).  With ε_{-1}=0, ε_0 = S_n,
        ε_{k+1}^{(n)} = ε_{k-1}^{(n+1)} + 1/(ε_k^{(n+1)} - ε_k^{(n)});
    the even columns ε_{2d} are the limit estimates.  WynnDegree -> d selects
    column 2d and requires at least 2(d+1) terms for a convergence check.

  Robustness
  ----------
  Two independent gates, both returning the form unevaluated when they fire.

  Oscillatory divergence (NLimit::osc).  An extrapolator returns a number for
  any input, including a sequence with no limit, so before trusting one we
  check that the oscillation envelope actually decays: samples whose
  increments keep reversing direction get a wide auxiliary ladder, and a
  non-decaying envelope is refused.  See the NL_OSC_* block below.

  Noise (NLimit::noise).  The last two extrapolates are compared; if they
  fail to agree to a loose tolerance relative to the result magnitude (or the
  result is non-finite), the form is returned unevaluated.  As with
  Mathematica, spurious tiny residuals are *not* recognised as zero — apply
  Chop when needed.
```

Options: Method (Automatic | EulerSum | SequenceLimit), WorkingPrecision (MachinePrecision

```text
| digits -> MPFR), Direction (Automatic == -1, or a complex approach vector),
```

Scale (initial step / distance, default 1), Terms (default 13), WynnDegree

```text
(default 1).  The default of 13 (not the historical 7) is deliberate: the
```

accuracy of a branch-point / fractional-power approach is limited by the *depth of the extrapolation tableau*, not by the arithmetic precision, so a short sequence starves the extrapolators regardless of WorkingPrecision. Thirteen samples resolve such tails to ~12 machine digits while leaving smooth limits at their existing roundoff floor; the extra six evaluations are sub-millisecond, and the Automatic best-of selector discards any deep tableau that only adds roundoff, so the larger default never degrades an easy case.

Memory: receives `res` owned by the evaluator; returns a fresh Expr* on

```text
success or NULL (unevaluated).  Never frees `res`.  Every temporary OwnValue
```

created for the sampler is removed on all return paths.

## Implementation notes

**Algorithm.** `builtin_nlimit` parses `NLimit[expr, z -> z0]` and constructs a
**geometric sequence of sample points** approaching `z0`: for a finite target,
`z_k = z0 − d·Scale·2^{−k}` (the default `Direction -> Automatic == −1`
approaches from larger values); for an infinite target (`Infinity`,
`ComplexInfinity`, `DirectedInfinity[d]`, `I Infinity`, …), `z_k = u·Scale·2^k`
marches outward along the ray. Each sample binds `z` Block-style and bumps the
evaluation clock. The limit is recovered by **sequence acceleration**; `Method
-> Automatic` runs three accelerators and keeps the one with the smallest
internal convergence residual:

- **Richardson/Romberg** (`EulerSum`, the `2^j − 1` tableau shared with `ND`) —
  trusted even at an exact-zero residual (integer-power tails converge exactly);
- **Wynn's epsilon** (`SequenceLimit`, the iterated Shanks transform) at degrees
  `1 … (terms−1)/2` — admitted only with a strictly positive residual;
- **Levin's u-transform** — admitted only when the sample increments are
  contracting (otherwise it can collapse to a spurious finite value with a
  deceptively small residual).

`Method -> "Levin"` (`"LevinU"`/`"LevinT"`/`"LevinV"`) forces the chosen Levin
variant; the kernels live in `seqaccel.c` (shared with `NSum`). Two robustness
gates both return the form unevaluated rather than a meaningless extrapolant: a
**non-decaying oscillation** screen (`NLimit::osc`) that counts increment
reversals, then — only if the screen fires — resamples on a 20-octave ladder at
two offsets (one scaled by the golden ratio to defeat power-of-two aliasing) and
refuses when the far/near envelope ratio stays above 0.6; and a **noise/
divergence** gate (`NLimit::noise`) that rejects a non-finite result, a result
far larger than the sample scale, or successive extrapolates that do not settle.
Spurious tiny residuals are *not* chopped. Each method is also exposed directly
as `NLimit`Automatic`/`EulerSum`/`SequenceLimit`/`Levin`/`LevinU`/`LevinT`/`LevinV`.

**Data structures.** The sample sequence is a `double _Complex*` of the sample
count on the machine path, or paired `mpfr_t*` real/imag buffers under MPFR; the
oscillation ladder is a small `double` magnitude array. `WorkingPrecision`
(`MachinePrecision` or a digit count) selects the path. The limit point is
classified and its direction normalised to a unit vector before sampling.

**Complexity / limits.** `[start, 50]` symbolic samples with an `O(terms²)`
Richardson tableau, plus Wynn over `~terms/2` degrees and one Levin pass; the
oscillation stage 2 (only when triggered) adds 40 samples. Options: `Method`
(`Automatic` | `EulerSum` | `SequenceLimit` | `"Levin"` variants), `Direction`
(`Automatic == −1` or a complex vector), `Scale` (default 1), `Terms` (default
**13**, grown adaptively — the tableau depth, not the arithmetic precision, sets
the accuracy of a branch-point/fractional-power approach), `WynnDegree`,
`WorkingPrecision`, `AccuracyGoal` (default `MachinePrecision`), `PrecisionGoal`.
Diagnostics (`badmeth`, `notnum`, `osc`, `noise`, `nc_warn_goal`) route through
`mth_message`.

**Attributes:** `Protected`.

## References

**See also:** [Chop](../../elementary-functions/Chop/), [Pi](../../mathematical-constants/Pi/), [ND](../../numerical-calculus/ND/), [Information](../../expression-information/Information/), [AccuracyGoal](../../other-advanced/AccuracyGoal/), [PrecisionGoal](../../other-advanced/PrecisionGoal/)

- P. Wynn, *On a device for computing the e_m(S_n) transformation*, MTAC **10** (1956) 91–96 — the epsilon algorithm (SequenceLimit).
- D. Levin, *Development of non-linear transformations for improving convergence of sequences*, Internat. J. Comput. Math. **3** (1973) 371–388.
- Source: [`src/numerical_calculus/nlimit.c`](https://github.com/stblake/mathilda/blob/main/src/numerical_calculus/nlimit.c)
- Specification: [`docs/spec/builtins/numerical-calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/numerical-calculus.md)
- Tests: [`tests/test_accuracygoal.c`](https://github.com/stblake/mathilda/blob/main/tests/test_accuracygoal.c)
- Tests: [`tests/test_nlimit.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nlimit.c)

## Notes & additional examples

### Notes

`NLimit[expr, z -> z0]` builds a geometric sequence of sample points approaching
`z0` and recovers the limit by sequence acceleration. The first three cases give
`1`, the constant `E = 2.71828...`, and `Log[2] = 0.693147...`. The fourth is the
classic Laurent-expansion limit of the Riemann zeta function at its pole: the
constant term is the Euler–Mascheroni constant `EulerGamma = 0.577216...`. `z0`
may be finite, complex, or an infinite point such as `Infinity` or `I Infinity`.
`Method -> Automatic` (default) keeps the most self-consistent of Richardson
extrapolation (`EulerSum`), Wynn's epsilon (`SequenceLimit`) and Levin's
u-transform; `Method -> "Levin"` (`"LevinU"`/`"LevinT"`/`"LevinV"`) forces
Levin's transformation. `Chop` small spurious residuals.
