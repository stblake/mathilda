# Integration methods

The [Calculus tutorial](../08-calculus.md) introduces `Integrate` and lets it
choose an algorithm for you (`Method -> Automatic`). Underneath, Mathilda's
`Integrate` is a **cascade** of specialised, correct-by-construction methods:
each one recognises a particular class of integrand, and the cascade tries them
in a deterministic order, stopping at the first that succeeds.

These advanced tutorials open up that cascade. Each page covers **one method**:
the mathematics it implements, how Mathilda realises it, the exact convergence
conditions it enforces, and a graded set of worked examples — every transcript
reproduced from the current build. You can always pin a single method with the
`Method` option:

```mathematica
Integrate[expr, {x, a, b}, Method -> "Mellin"]
```

Pinning a method is useful for two reasons: to *understand* which engine solved
an integral, and to *force* a specific engine when you know the automatic
dispatch would reach for a slower or less appropriate one.

<div class="grid cards" markdown>

-   :material-transfer: __[Mellin transforms](mellin.md)__

    Half-line integrals \(\int_0^\infty x^{s-1} f(x)\,dx\) in closed form via a
    table of base Mellin transforms, the Ramanujan Master Theorem for
    hypergeometric integrands, and rigorously enforced convergence strips.

-   :material-vector-circle: __[Residue theorem](residue.md)__

    Improper and periodic integrals via Cauchy's residue theorem — rational,
    Fourier/Jordan, and trigonometric families, plus sector, keyhole, and
    rectangular contours, principal values, and symbolic parameters.

-   :material-vector-polyline: __[Contour / line integrals](lineintegral.md)__

    Complex integrals along piecewise-linear paths — Cauchy's theorem, winding
    numbers, branch-correct logarithms, and the residue theorem worked
    constructively, with a numerical crosscheck on every segment.

</div>

!!! note "More methods coming"
    This section is being filled in one method at a time. The complete list of
    method names accepted by `Integrate[..., Method -> "..."]` is
    `"DerivativeDivides"`, `"RischNorman"`, `"BronsteinRational"`,
    `"Weierstrass"`, `"NewtonLeibniz"`, `"DiffUnderInt"`, `"LineIntegral"`,
    `"Residue"`, and `"Mellin"` (a synonym for `"RamanujanMasterTheorem"`).

## Following along

Start the REPL with `./Mathilda` and type each `In[...]` line yourself (without
the prompt). Press Return to evaluate. Every `Out[...]` shown in these pages was
produced by the actual binary; Mathilda's output form is sometimes arranged
differently from a textbook, but it is always mathematically correct.
