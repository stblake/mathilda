---
hide:
  - navigation
  - toc
---

<div class="mathilda-hero" markdown>

# Mathilda

**A small Mathematica-like computer algebra system — written from scratch in C99.**

[Explore the documentation :material-book-open-variant:](documentation/index.md){ .mathilda-cta .md-button .md-button--primary }
[Start the tutorials :material-school:](tutorials/index.md){ .mathilda-cta .md-button }
[View on GitHub :fontawesome-brands-github:](https://github.com/stblake/mathilda){ .mathilda-cta .md-button }

</div>

Mathilda is a small computer algebra system that recreates the **core
architecture and evaluation semantics of the Mathematica programming language** — a recursive
expression model, attribute-driven evaluation, structural pattern matching with
backtracking, and a rewrite-rule engine — together with an extensive library of
**~695 built-in functions**.

It spans roughly **294,000 lines of C99**, uses **GMP** for arbitrary-precision
integers and **MPFR** for arbitrary-precision reals, and is licensed under
**GPLv3**.

---

## See it in action

Every example on this site — including the ones below — is run through the
actual Mathilda build and its real output captured. Nothing is transcribed by
hand. Here are some famous results, computed from scratch:

```mathematica
In[1]:= Sum[1/k^2, {k, 1, Infinity}]                    (* the Basel problem *)
Out[1]= 1/6 Pi^2

In[2]:= Integrate[Exp[-x^2], {x, -Infinity, Infinity}]  (* the Gaussian integral *)
Out[2]= Sqrt[Pi]

In[3]:= Zeta[-1]                                         (* 1 + 2 + 3 + ... "=" -1/12 *)
Out[3]= -1/12

In[4]:= Sum[HarmonicNumber[k]/k^2, {k, 1, Infinity}]    (* an Euler sum *)
Out[4]= 2 Zeta[3]

In[5]:= Integrate[Sin[x]/x, {x, 0, Infinity}]           (* the Dirichlet integral *)
Out[5]= 1/2 Pi

In[6]:= Factor[x^4 + 4]                                  (* the Sophie Germain identity *)
Out[6]= (2 - 2 x + x^2) (2 + 2 x + x^2)

In[7]:= Limit[(1 + 1/n)^n, n -> Infinity]               (* a definition of e *)
Out[7]= E

In[8]:= ContinuedFraction[GoldenRatio, 8]               (* the "most irrational" number *)
Out[8]= {1, 1, 1, 1, 1, 1, 1, 1}

In[9]:= Series[Exp[x], {x, 0, 5}]
Out[9]= 1 + x + 1/2 x^2 + 1/6 x^3 + 1/24 x^4 + 1/120 x^5 + O[x]^6

In[10]:= FactorInteger[2^67 - 1]                         (* Cole's 1903 factorization *)
Out[10]= {{193707721, 1}, {761838257287, 1}}

In[11]:= PrimePi[10^9]                                   (* primes below a billion *)
Out[11]= 50847534
```

---

## What's inside

<div class="grid cards" markdown>

-   :material-function-variant: __Evaluation engine__

    Infinite-evaluation semantics: expressions are rewritten top-down to a fixed
    point. A small generic core consults per-symbol attribute bits
    (`HoldAll`, `Flat`, `Orderless`, `Listable`, `OneIdentity`, …) to decide how
    to process each call.

-   :material-vector-polyline: __Pattern matching & rules__

    First-class `Blank` (`_`), `BlankSequence` (`__`), named bindings (`x_`),
    `Condition` (`/;`), `PatternTest` (`?`), `Optional`, and `Repeated` — with
    full sequence backtracking. Transformation rules (`->`, `:>`) and
    replacement (`/.`, `//.`).

-   :material-numeric: __Numbers__

    Arbitrary-precision integers via GMP, exact rationals and complex numbers,
    and MPFR-backed arbitrary-precision reals with precision/accuracy tracking
    (`N[expr, prec]`).

-   :material-sigma: __Symbolic mathematics__

    Differentiation, multi-method integration, `Series`, `Limit`, symbolic
    summation, polynomial factorization (over ℤ and ℚ(α)), Gröbner bases,
    dense linear algebra, and a complexity-driven `Simplify`.

-   :material-calculator-variant: __Numerical calculus__

    Machine- and arbitrary-precision numerics for the cases with no closed
    form: `NIntegrate` (adaptive, oscillatory, multidimensional), `NSum`,
    `NProduct`, `ND`, `NLimit`, `NSeries`, and `NResidue`.

-   :material-function: __Special functions__

    Gamma, log-gamma, beta and the digamma/polygamma family; the Riemann/Hurwitz
    `Zeta` and Stieltjes constants; `Erf`/`Erfc`/`Erfi`; `ExpIntegralEi`,
    `LogIntegral`, `PolyLog`, Bernoulli/Euler numbers, and the hypergeometric
    family.

-   :material-shape: __Number theory & factorization__

    `GCD`, `PowerMod`, `EulerPhi`, `PrimitiveRoot`, continued fractions, and an
    automatic integer-factorization pipeline (Pollard Rho/P−1, Williams P+1,
    Fermat, CFRAC, Dixon, ECM).

-   :material-code-braces: __Programming__

    Functional programming (`Map`, `Apply`, `Fold`, `Nest`, pure functions),
    scoping (`Module`, `Block`, `With`), control flow, and a standard library of
    lists, strings, statistics, and dates.

</div>

---

## Build & run

Mathilda builds with a C99 toolchain and links GMP, MPFR and GNU Readline.
FLINT (≥ 3.0), GMP-ECM and LAPACK/BLAS are optional (all auto-detected): FLINT
provides fast, rigorous algebraic-extension arithmetic and `acb` numerics,
GMP-ECM powers advanced integer factorization, and LAPACK/BLAS accelerates
machine-precision linear algebra.

### Install dependencies

=== "Linux (Debian / Ubuntu)"

    ```bash
    # Required libraries
    sudo apt install libgmp-dev libmpfr-dev libreadline-dev

    # Optional: FLINT (>= 3.0) for fast, rigorous algebraic-extension arithmetic
    sudo apt install libflint-dev

    # Optional: GMP-ECM for advanced integer factorization
    sudo apt install libecm-dev

    # Optional: LAPACK/BLAS (fast linear algebra) and CMake (test suite)
    sudo apt install liblapacke-dev libopenblas-dev cmake
    ```

=== "macOS (Homebrew)"

    ```bash
    brew install gmp mpfr readline cmake
    brew install flint                       # optional: FLINT-backed kernels
    brew install gmp-ecm                      # optional: advanced integer factorization
    # LAPACK/BLAS is provided by Apple's Accelerate framework — no install needed.
    ```

### Clone, build, run

```bash
git clone https://github.com/stblake/mathilda.git
cd mathilda
make -j            # builds ./Mathilda
./Mathilda         # start the interactive REPL
```

!!! tip "Building with FLINT"
    FLINT (≥ 3.0) is enabled automatically when `pkg-config` finds it — no flag
    needed. If it is missing or older than 3.0, the build prints a warning and
    falls back to `USE_FLINT=0` (the classical, still-rigorous paths). Force it
    off with `make USE_FLINT=0`, and confirm the installed version with
    `pkg-config --modversion flint`. See the
    [FLINT context](documentation/flint/index.md) for the routines it powers.

!!! tip "Building with GMP-ECM"
    GMP-ECM (the Elliptic Curve Method for integer factorization) is a plain
    system library — install `gmp-ecm` (Homebrew) or `libecm-dev` (Debian/Ubuntu)
    and the build autodetects it via a compile-link probe and links `-lecm`. When
    it is absent the build still succeeds with advanced factorization disabled;
    force that with `make USE_ECM=0`.

Then type an expression and press Return. Ask for help on any function with
`?Name`:

```mathematica
In[1]:= ?Integrate
```

---

## Where to next

- **[Documentation Center](documentation/index.md)** — every built-in function,
  grouped by category, each with a description, verified examples,
  implementation notes, status, and references.
- **[Tutorials](tutorials/index.md)** — guided, worked walkthroughs from first
  launch through pattern matching and symbolic calculus.

!!! info "About these docs"
    This site is generated from Mathilda's own docstrings and specification, and
    every code example is verified against the current build. See the
    [project README](https://github.com/stblake/mathilda) and
    [`SPEC.md`](https://github.com/stblake/mathilda/blob/main/SPEC.md) for the
    architecture in depth.

<!-- Published via GitHub Actions to https://stblake.github.io/mathilda/ -->
