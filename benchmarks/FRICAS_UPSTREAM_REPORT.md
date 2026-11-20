# Draft upstream report for FriCAS: an elementary integral left unevaluated (false "not elementary")

*Draft, not submitted. Found by the Mathilda–FriCAS comparison (`benchmarks/BRONSTEIN_VS_FRICAS.md`,
benchmark 259, the `gauss` family). Reproduced with FriCAS 1.3.6 (2022-11-13 build).*

## Minimal example

```
integrate(exp(1+2*%i*x)/(1+exp(%i*x)), x)
```

FriCAS returns the integral unevaluated:

```
   x   2 %i %A + 1
 ++  %e
 |   ------------- d%A
++      %i %A
      %e      + 1
```

An unevaluated `integrate` from FriCAS means the Risch algorithm has proved that there is no
elementary antiderivative. This integral is elementary. With `t = exp(%i*x)` and
`dt = %i*t*dx`:

```
exp(1+2*%i*x)/(1+exp(%i*x)) dx = %e*t^2/(1+t) * dt/(%i*t) = -%i*%e * t/(1+t) dt
integral = -%i*%e*(exp(%i*x) - log(1 + exp(%i*x)))
```

which `D(-%i*%e*(exp(%i*x) - log(1+exp(%i*x))), x) - exp(1+2*%i*x)/(1+exp(%i*x))`
normalizes to 0.

## Related cases

- The real-exponent analogue is integrated correctly:
  `integrate(exp(1+2*x)/(1+exp(x)), x)` gives `-%e log(%e^x + 1) + %e %e^x`.
- The difference is the offset: `exp(1+2*%i*x)` and `exp(%i*x)` are not recognised as
  algebraically dependent when the exponents differ by a non-zero constant over `Q(i)`. The
  same holds for `exp(%i*(2+x))` beside `exp(2*%i*x)`, and for `exp(1+3*%i*x)`.
- The comparison found 14 integrands of this kind. 11 are elementary by construction
  (`f = D[F]`), and 3 more have an elementary antiderivative that verifies. All have an
  unevaluated `integral` from FriCAS. Examples:
  - `integrate(D(x*exp(1+2*%i*x)/(1+exp(%i*x)), x), x)` (FriCAS returns the derivative
    unevaluated as an integral);
  - `exp(1 + 2*%i*x)*tan(x)`, `exp(2*%i*(2 + x))*tan(x)`, `exp(%i - %i*(1 + x))*tan(2 + x)`.

## Suggested cause

The structure-theorem check that decides whether a new exponential `exp(w)` is algebraic
over the field built so far looks for `w = sum k_j w_j + c` with integer `k_j`. The constant
`c` (here `1`, or a multiple of `%i`) then has to be handled as an element of the constant
field. Over `Q(i)`, with both `%e` and `exp(%i)` in play, this case appears to make the
dependency go undetected, so `exp(1+2*%i*x)` becomes an independent monomial. The integrand
is then non-elementary over that (wrong) tower.
