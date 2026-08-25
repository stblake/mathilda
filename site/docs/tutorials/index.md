# Tutorials

A guided path through Mathilda, from your first REPL session to writing your own
pattern-based rules, doing symbolic calculus, and closing infinite series and
products into constants. Every example is worked end to end and **verified
against the current Mathilda build**.

Work through them in order if you're new — each one builds on the last.

<div class="grid cards" markdown>

-   :material-rocket-launch: __[1. Getting started](01-getting-started.md)__

    Build Mathilda, launch the REPL, understand `In[]`/`Out[]`, learn the
    surface syntax, and get help on any function with `?Name`.

-   :material-file-tree: __[2. Expressions & evaluation](02-expressions-and-evaluation.md)__

    Everything is an expression. Meet `FullForm`, `Head`, the attribute system,
    the fixed-point evaluator, and how `Hold` suspends evaluation.

-   :material-vector-polyline: __[3. Pattern matching & rules](03-pattern-matching-and-rules.md)__

    Blanks and named patterns, conditions and tests, transformation rules
    (`->`, `:>`), replacement (`/.`, `//.`), and defining your own functions.

-   :material-decimal: __[4. Arithmetic](04-arithmetic.md)__

    Exact integers and rationals, fast machine-precision reals, and
    arbitrary-precision arithmetic (`N`, `Precision`); the basic operators,
    digit and radix manipulation, and combinatorial functions.

-   :material-key-variant: __[5. Number theory](05-number-theory.md)__

    `GCD`, `ExtendedGCD`, modular arithmetic and `PowerMod`, primes
    (`PrimeQ`, `FactorInteger`, `NextPrime`), `EulerPhi`, and continued
    fractions — up to RSA-style worked examples.

-   :material-function-variant: __[6. Algebra](06-algebra.md)__

    Expand and factor polynomials, dissect and divide them, reshape rational
    expressions with `Together`/`Apart`, simplify, and put the polynomial
    toolkit (`Resultant`, `GroebnerBasis`) to work on real problems.

-   :material-equal: __[7. Solutions of equations](07-solutions-of-equations.md)__

    Solve polynomial, transcendental, and simultaneous equations with `Solve`;
    `Root` objects and `ToRadicals`; eliminate variables with `Eliminate`;
    complete parametric solution sets with `Reduce`; and tackle geometry and
    optimisation problems.

-   :material-code-greater-than-or-equal: __[Solutions of inequalities](solutions-of-inequalities.md)__

    Carve out regions rather than points: `Reduce` on inequalities over the
    reals — sign diagrams, rational functions with poles, `Abs` and piecewise
    functions, integer ranges, and two- and three-variable regions solved by
    cylindrical algebraic decomposition.

-   :material-key-chain: __[Diophantine equations](diophantine-equations/index.md)__

    Solve polynomial equations over the **integers**: linear systems and Pell,
    quadratic and ternary forms, Mordell and **Thue** equations, exponential
    Diophantine equations, and the famous power-sum searches (**Lander–Parkin**,
    **Frye**, taxicabs) — with a head-to-head against sympy and PARI/GP.

-   :material-math-integral: __[8. Calculus](08-calculus.md)__

    Differentiate and integrate, expand power series, take limits, evaluate
    symbolic sums, and find roots and extrema numerically.

-   :material-transfer: __[Integration methods](integration-methods/index.md)__

    Advanced, per-method deep dives into `Integrate`'s cascade — the
    **transcendental Risch** decision procedure and **Cherry's special-function**
    extensions (`erf`, `Ei`, `li`, dilogarithm), the **Mellin transform** engine
    for half-line integrals \(\int_0^\infty x^{s-1} f\,dx\), and the **residue
    theorem** engine for improper, periodic, and contour integrals. Algorithm,
    references, and worked examples.

-   :material-calculator-variant: __[9. Numerical calculus](09-numerical-calculus.md)__

    When there is no closed form: numerical integration, differentiation,
    summation, products, limits, series, and residues — `NIntegrate`, `ND`,
    `NSum`, `NProduct`, `NLimit`, `NSeries`, `NResidue`.

-   :material-sigma-lower: __[10. Special functions](10-special-functions.md)__

    The higher transcendental functions: `Gamma`, `Zeta`, `PolyGamma`, `Erf`,
    `PolyLog`, the Bernoulli and Euler numbers, and the hypergeometric family —
    with their exact reductions and numerical values.

-   :material-arrow-expand-horizontal: __[Interval arithmetic](15-interval-arithmetic.md)__

    Rigorous enclosures that never lie: exact and symbolic endpoints, the
    dependency problem, functions with poles and discontinuities, and validated
    numerics — proving a root exists, guaranteed bounds, and rounding error made
    visible.

-   :material-sigma: __[11. Symbolic summation](11-symbolic-summation.md)__

    Close infinite series into constants with `Sum`: telescoping and Gosper's
    algorithm, the Basel problem and the zeta family, Euler sums and multiple
    zeta values, binomial sums, and the hypergeometric machines for `π`.

-   :material-pi: __[12. Infinite products](12-infinite-products.md)__

    Evaluate infinite products with `Product`: rational telescoping, the
    Wallis/Viète trigonometric factorizations, Euler prime products for `ζ`,
    and the exponential products for `e`, `γ`, and Glaisher's constant.

-   :material-cog-play: __[Compilation & auto-compilation](16-compilation.md)__

    Run numeric work at machine speed: `Compile[]` for ints, reals, arrays and
    associations; automatic compilation (`$AutoCompilation`) and packed arrays
    (`$AutoArrayPacking`) that you get for free; the compilable-subset cliff and
    how to see it — measured against Python 3.11 + NumPy at every step.

-   :material-grid: __[13. BLAS kernels](13-blas.md)__

    Call the machine-precision BLAS kernels directly through the `` BLAS` ``
    context: dot products and norms, `dgemv`, `dgemm`, the symmetric and
    triangular Level-3 routines, and their complex `z*` counterparts.

-   :material-table-large: __[14. LAPACK drivers](14-lapack.md)__

    Solve systems and least squares, factor matrices (LU, QR, Cholesky),
    compute the SVD, and solve symmetric and general eigenproblems with the
    `` LAPACK` `` context.

</div>

!!! tip "Following along"
    Start the REPL with `./Mathilda` and type each `In[...]` line yourself
    (without the prompt). Press Return to evaluate. End a line with `\` to
    continue a long expression onto the next line.
