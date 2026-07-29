# Goursat Integrator — Graded Exercises

A graded battery of pseudo-elliptic integrals, each **solved by the Goursat
integrator** (`Integrate`GoursatAlgebraic`), ordered from easy to hard.  The
set deliberately spans:

* **every radical exponent** `p ∈ {1/2, 1/3, 2/3, 1/4, 3/4}`;
* **both signs of the exponent** — radical in the *denominator* (`R^(-p)`,
  negative exponent) and in the *numerator* (`R^(+p)`, positive exponent); and
* **every involution / eigenspace equation** the algorithm tests.

Each exercise is reproduced as a unit test in
[`tests/test_integrate_goursat.c`](tests/test_integrate_goursat.c)
(`test_graded`), forced through `Method -> "GoursatAlgebraic"` so the result is
guaranteed to come from this integrator and not the `Automatic` cascade.

## Form recognised

The integrand must reduce to

```
F(t) · R(t)^(n/d) ,   F rational,   R a polynomial,   n/d any rational of
                       reduced denominator d ∈ {2, 3, 4}.
```

Recognition rewrites `R^(n/d) = R^k · R^(-p)` with `p = pnum/d ∈
{1/2,1/3,2/3,1/4,3/4}` and folds the integer power `R^k` into `F`, so the
radical may **start in the numerator or the denominator** — both are handled
uniformly.

## Involution equations exercised

| Branch | Radical | Elementarity criterion |
|--------|---------|------------------------|
| Square-root V4 | `p = 1/2`, quartic `R` | trivial-character projection `P0 == 0` |
| Period-3 | `p = 1/2`, cubic `R` | order-3 character `F + F(S) + F(S²) == 0` |
| Cube-root | `p = 1/3` | `H1 == 0` |
| Cube-root | `p = 2/3` | `H0 == 0` |
| Fourth-root | `p = 1/4` | `H1 == H2 == 0` (harmonic quartic) |
| Fourth-root | `p = 3/4` | `H0 == H1 == 0` (harmonic quartic) |

Set `Integrate`GoursatDebug = True` to watch form recognition, the criterion
tests, and the recursive reductions on stderr.

---

## Tier 1 — recognise + fold (antiderivative is a single radical power)

Negative exponents; the cofactor folds so the result is just a power of `R`.

| # | ∫ f dt | Result |
|---|--------|--------|
| 1 | `t^2/(t^3-1)^(1/3)` | `(1/2) (t^3-1)^(2/3)` |
| 2 | `t^2/(t^3-1)^(2/3)` | `(t^3-1)^(1/3)` |
| 3 | `t^3/(t^4-1)^(1/4)` | `(1/3) (t^4-1)^(3/4)` |
| 4 | `t^3/(t^4-1)^(3/4)` | `(t^4-1)^(1/4)` |

## Tier 2 — positive exponent (radical in the numerator)

`R^(+p)` — the same machinery, radical promoted to the numerator.

| # | ∫ f dt | Result |
|---|--------|--------|
| 5 | `t^2 (t^3-1)^(1/3)` | `(1/4) (t^3-1)^(4/3)` |
| 6 | `t^2 (t^3-1)^(2/3)` | `(1/5) (t^3-1)^(5/3)` |
| 7 | `t^3 (t^4-1)^(1/4)` | `(1/5) (t^4-1)^(5/4)` |
| 8 | `t^3 (t^4-1)^(3/4)` | `(1/7) (t^4-1)^(7/4)` |

## Tier 3 — genuine cube-root involution (ArcTan + Log)

These trigger the `H1 == 0` (`p = 1/3`) and `H0 == 0` (`p = 2/3`) reductions
and return logarithmic / inverse-tangent closed forms.

* **9.** `∫ dt/(t^3-1)^(1/3)` &nbsp;(`p=1/3`, `H1==0`)

  ```
  -ArcTan[1/Sqrt[3] + 2 (t^3-1)^(1/3)/(Sqrt[3] t)]/Sqrt[3]
    - (1/3) Log[-1 + (t^3-1)^(1/3)/t]
    + (1/6) Log[1 + (t^3-1)^(1/3)/t + (t^3-1)^(2/3)/t^2]
  ```

* **10.** `∫ t/(t^3-1)^(2/3) dt` &nbsp;(`p=2/3`, `H0==0`, dual of #9)

  ```
  -ArcTan[1/Sqrt[3] + 2 t/(Sqrt[3] (t^3-1)^(1/3))]/Sqrt[3]
    - (1/3) Log[-1 + t/(t^3-1)^(1/3)]
    + (1/6) Log[1 + t/(t^3-1)^(1/3) + t^2/(t^3-1)^(2/3)]
  ```

* **11.** `∫ dt/(t (t^3-1)^(2/3))` &nbsp;(`p=2/3`, `H0==0`, pole in the cofactor)

  ```
  (1/3) Log[1 + (t^3-1)^(1/3)]
    + ArcTan[-1/Sqrt[3] + 2 (t^3-1)^(1/3)/Sqrt[3]]/Sqrt[3]
    - (1/6) Log[1 - (t^3-1)^(1/3) + (t^3-1)^(2/3)]
  ```

* **12.** `∫ dt/(t^3-8)^(1/3)` &nbsp;(`p=1/3`, shifted radicand) — same shape
  as #9 over `t^3-8`.

## Tier 4 — genuine fourth-root involution (ArcTan + ArcTanh)

Harmonic quartics; `H1==H2==0` (`p=1/4`) and `H0==H1==0` (`p=3/4`).

* **13.** `∫ dt/(t^4-1)^(1/4)` &nbsp;(`p=1/4`)

  ```
  -(1/2) ArcTan[(t^4-1)^(1/4)/t] + (1/2) ArcTanh[(t^4-1)^(1/4)/t]
  ```

* **14.** `∫ t^2/(t^4-1)^(3/4) dt` &nbsp;(`p=3/4`, dual of #13)

  ```
  -(1/2) ArcTan[t/(t^4-1)^(1/4)] + (1/2) ArcTanh[t/(t^4-1)^(1/4)]
  ```

* **15.** `∫ dt/(t^4-16)^(1/4)` &nbsp;— a harmonic quartic other than `t^4-1`:

  ```
  -(1/2) ArcTan[(t^4-16)^(1/4)/t] + (1/2) ArcTanh[(t^4-16)^(1/4)/t]
  ```

* **16.** `∫ t^2/(t^4-16)^(3/4) dt`

  ```
  -(1/2) ArcTan[t/(t^4-16)^(1/4)] + (1/2) ArcTanh[t/(t^4-16)^(1/4)]
  ```

## Tier 5 — Goursat's square-root V4 (ArcTanh + Log)

The original Goursat algorithm: `p = 1/2`, quartic radicand, `V4` trivial
projection `P0 == 0`.  Closed forms live over `Sqrt[(t^2-1)(t^2-4)]`.

* **17.** `∫ t/Sqrt[(t^2-1)(t^2-4)] dt` &nbsp;(Example 5.1).
* **18.** `∫ (t^2-2)/(t Sqrt[(t^2-1)(t^2-4)]) dt` &nbsp;(finite fixed point `S = 2/t`)

  ```
  ArcTanh[(Sqrt[t^4-5t^2+4] - t^2)/2] - (1/2) Log[5 + 2 (Sqrt[t^4-5t^2+4] - t^2)]
  ```

* **19.** `∫ (t^2-2)/(t Sqrt[(t^2-4)(t^2-9)]) dt` &nbsp;(different radicand).
* **20.** `∫ (t^4 + 2t^3 - 4)/(t^2 Sqrt[(t^2-1)(t^2-4)]) dt` &nbsp;— regression
  guarding the **negative leading-coefficient** descent (`Sqrt[lc·Q]` must not
  split `Sqrt[lc]` into the prefactor across a branch cut).

## Tier 6 — positive square-root exponent (numerator)

* **21.** `∫ t Sqrt[(t^2-1)(t^2-4)] dt` &nbsp;(`p=1/2`, radical in the
  numerator).  Closes to a large elementary form (a rational function of
  `Sqrt[t^4-5t^2+4]` plus an `ArcTanh` and a `Log`); the unit test verifies it
  by a numeric differentiate-back.

## Tier 7 — period-3 higher symmetry (hardest)

* **22.** `∫ (t-1)/((t+2) Sqrt[t^3-1]) dt` &nbsp;— Goursat (1887), §4.  Here the
  cubic radicand admits an **order-3** Möbius map `S` and `F` is a non-trivial
  period-3 character (`F(S) = α F`): the `V4` trivial projection is non-zero,
  but `F + F(S) + F(S²) == 0`.  The antiderivative is a complex `Log`; the unit
  test verifies it numerically on the real axis (`t > 1`, where the integrand is
  real).

---

## Negative controls — the involution gate must DECLINE

These are *not* elementary by Goursat (or not elementary at all); the forced
method must leave them unevaluated.

| ∫ f dt | Why it declines |
|--------|-----------------|
| `t/(t^3-1)^(1/3)` | cube-root `H1 ≠ 0` |
| `1/(t^3-1)^(2/3)` | cube-root `H0 ≠ 0` |
| `t/(t^4-1)^(1/4)` | fourth-root `V1` obstructive |
| `t^2 (t^4-1)^(3/4)` | non-elementary (Chebyshev binomial-differential) |
| `t^2/Sqrt[(t^2-1)(t^2-4)]` | `V4`-invariant → genuinely elliptic |
| `1/Sqrt[t^3-1]` | period-3 trivial projection ≠ 0 → elliptic |

---

## Reproducing

```wolfram
(* Solve any exercise, forcing the Goursat method: *)
Integrate[1/(t^3-1)^(1/3), t, Method -> "GoursatAlgebraic"]

(* Trace the recognition + criterion tests: *)
Integrate`GoursatDebug = True;
Integrate[t/Sqrt[(t^2-1)(t^2-4)], t, Method -> "GoursatAlgebraic"]
```

Build and run the corresponding unit tests:

```bash
cd tests/build && cmake --build . --target integrate_goursat_tests -j8
./integrate_goursat_tests
```
