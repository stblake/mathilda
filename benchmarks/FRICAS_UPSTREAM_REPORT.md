# FriCAS upstream report (submitted)

Submitted 2026-11-21 as https://github.com/fricas/fricas/issues/249. Every antiderivative in it
is Mathilda's, verified in FriCAS 1.3.6 by `normalize` / `complexNormalize` of `D(r, x) - f`.
The earlier draft listed `integrate(D(x*exp(1+2*%i*x)/(1+exp(%i*x)), x), x)` as left
unevaluated; FriCAS in fact integrates that one, so it was dropped.

The text as submitted:

---

FriCAS 1.3.6 (2022-11-13 build) returns these integrals unevaluated, which normally means the Risch algorithm proved that no elementary antiderivative exists. Each one does have an elementary antiderivative. The answers below were computed by the Mathilda CAS, and each was verified in FriCAS itself: `normalize` or `complexNormalize` of `D(r, x) - f` gives `0`.

## Minimal example

```
integrate(exp(1+2*%i*x)/(1+exp(%i*x)), x)
```

FriCAS returns:

```
   x   2 %i %A + 1
 ++  %e
 |   ------------- d%A
++      %i %A
      %e      + 1
```

An elementary antiderivative, from Mathilda:

```
r := -%i*exp(1+%i*x) + %i*%e*log(1+exp(%i*x))
normalize(D(r, x) - exp(1+2*%i*x)/(1+exp(%i*x)))      -- 0
```

By hand: with `t = exp(%i*x)`, the integrand times `dx` is `-%i*%e*t/(1+t) dt`.

The real-exponent analogue is integrated correctly: `integrate(exp(1+2*x)/(1+exp(x)), x)` gives `-%e log(%e^x + 1) + %e %e^x`.

## More examples, same pattern

Each is left unevaluated by FriCAS 1.3.6. The antiderivatives are Mathilda's, verified with `complexNormalize(D(r, x) - f) = 0` in FriCAS:

```
f := exp(1+2*%i*x)*tan(x)
r := -exp(1+2*%i*x)/2 + %e*log(%e^2 + exp(2+2*%i*x))

f := exp(2*%i*(2+x))*tan(x)
r := -exp(2*%i*(2+x))/2 + exp(4*%i)*log(exp(8*%i) + exp(8*%i+2*%i*x))

f := exp(%i-%i*(1+x))*tan(2+x)
r := -exp(2*%i-%i*(2+x)) + 2*exp(2*%i)*atan(exp(-%i*(2+x)))
```

In a comparison over 640 generated transcendental integrands, 14 had this shape. 11 were elementary by construction (`f = D(F)`), and 3 more had an elementary antiderivative that verifies. FriCAS left all 14 unevaluated.

## Possible cause

In every case two exponentials have exponents that differ by a nonzero constant over Q(i): `exp(1+2*%i*x)` beside `exp(%i*x)`, and `exp(%i*(2+x))` beside `exp(2*%i*x)`. The structure-theorem test decides whether a new `exp(w)` is algebraic over the tower built so far, by looking for `w = sum k_j w_j + c` with integer `k_j` and a constant `c`. Over Q(i) it seems to miss this dependency, so `exp(1+2*%i*x)` is taken as a new independent monomial. Over that wrong tower the integrand is then (correctly) non-elementary. We have only tested 1.3.6 and do not know whether a later release fixes this.
