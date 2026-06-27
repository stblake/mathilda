# Calculus

13 built-in function(s) in this category.

- [`D`](D.md) — D[f, x] gives the partial derivative of f with respect to x.  _(Stable)_
- [`Derivative`](Derivative.md) — f' represents the derivative of a function f of one argument.  _(Stable)_
- [`DifferenceDelta`](DifferenceDelta.md) — DifferenceDelta[f, i] gives the forward difference (f /. i -> i+1) - f, the discrete analogue of D. It is the left inverse of indefinite Sum.  _(Stable)_
- [`Dt`](Dt.md) — Dt[f] gives the total derivative of f.  _(Stable)_
- [`FindMaximum`](FindMaximum.md) — FindMaximum[f, {x, x0}]  _(Stable)_
- [`FindMinimum`](FindMinimum.md) — FindMinimum[f, {x, x0}]  _(Stable)_
- [`FindRoot`](FindRoot.md) — FindRoot[f, {x, x0}]  _(Stable)_
- [`Integrate`](Integrate.md) — Integrate[f, x] gives the indefinite integral of f with respect to x.  _(Partial)_
- [`Limit`](Limit.md) — Limit[f, x -> a]  _(Stable)_
- [`PolynomialQuotientRemainder`](PolynomialQuotientRemainder.md) — PolynomialQuotientRemainder[p, q, x] returns {Quotient, Remainder}  _(Stable)_
- [`Product`](Product.md) — Product[f, {i, imax}]  _(Stable)_
- [`SubresultantPolynomialRemainders`](SubresultantPolynomialRemainders.md) — SubresultantPolynomialRemainders[a, b, x] gives the polynomial-remainder  _(Stable)_
- [`Sum`](Sum.md) — Sum[f, {i, imax}] gives the sum of f for i from 1 to imax. Sum[f, {i, imin, imax}], Sum[f, {i, imin, imax, di}] and Sum[f, {i, {i1, i2, ...}}] use the standard iterator forms; multiple iterators give nested sums. Sum[f, i] gives the indefinite sum (antidifference). Symbolic and infinite sums are evaluated in closed form via Method -> "Polynomial" | "Geometric" | "Gosper".  _(Partial)_
