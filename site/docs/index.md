---
hide:
  - navigation
  - toc
---

<div class="mathilda-hero" markdown>

# Mathilda

**A tiny, faithful symbolic computer algebra system — written from scratch in C99.**

[Explore the documentation :material-book-open-variant:](documentation/index.md){ .mathilda-cta .md-button .md-button--primary }
[Start the tutorials :material-school:](tutorials/index.md){ .mathilda-cta .md-button }
[View on GitHub :fontawesome-brands-github:](https://github.com/stblake/mathilda){ .mathilda-cta .md-button }

</div>

Mathilda is a small ("pico") computer algebra system that recreates the **core
architecture and evaluation semantics of the Wolfram Language** — a recursive
expression model, attribute-driven evaluation, structural pattern matching with
backtracking, and a rewrite-rule engine — together with an extensive library of
**~388 built-in functions**.

It spans roughly **159,000 lines of C99**, uses **GMP** for arbitrary-precision
integers and **MPFR** for arbitrary-precision reals, and is licensed under
**GPLv3**.

---

## See it in action

Every example on this site — including the ones below — is run through the
actual Mathilda build and its real output captured. Nothing is transcribed by
hand.

```mathematica
In[1]:= Factor[x^4 - 1]
Out[1]= (-1 + x) (1 + x) (1 + x^2)

In[2]:= Solve[x^2 - 5x + 6 == 0, x]
Out[2]= {{x -> 2}, {x -> 3}}

In[3]:= D[Sin[x^2], x]
Out[3]= 2 x Cos[x^2]

In[4]:= Integrate[1/(x^2 + 1), x]
Out[4]= ArcTan[x]

In[5]:= Sum[k^2, {k, 1, n}]
Out[5]= 1/6 n (1 + n) (1 + 2 n)

In[6]:= Series[Exp[x], {x, 0, 4}]
Out[6]= 1 + x + 1/2 x^2 + 1/6 x^3 + 1/24 x^4 + O[x]^5

In[7]:= Limit[Sin[x]/x, x -> 0]
Out[7]= 1

In[8]:= Eigenvalues[{{2, 0}, {0, 3}}]
Out[8]= {3, 2}

In[9]:= FactorInteger[2^67 - 1]
Out[9]= {{193707721, 1}, {761838257287, 1}}
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

```bash
git clone https://github.com/stblake/mathilda.git
cd mathilda
make -j            # builds ./Mathilda
./Mathilda         # start the interactive REPL
```

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
