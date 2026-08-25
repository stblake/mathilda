# The Mathilda Book — Roadmap

Each section below is scoped as **its own campaign**: plan it, write it, verify it,
then move on. Chapters are independent enough that they can be written in almost any
order once the foundation (Campaign 0) is in place. When you start a section, open a
plan for it, write to its `chapters/NN-slug.tex` + `examples/NN-slug/`, and update the
status here when done.

Read `CONTEXT.md` before writing any chapter.

## Status legend

- **Planned** — stub chapter exists; not yet written.
- **Drafting** — prose in progress.
- **Verified** — written, all transcripts generated, `make check-links` + `make pdf` clean.
- **Done** — Verified and reviewed.

## Foundation

| # | Section | File | Status |
|---|---------|------|--------|
| 0 | Toolchain, `CONTEXT.md`, `ROADMAP.md`, front matter | `mathilda.sty`, `tools/`, `Makefile`, `frontmatter/` | **Done** |

## Getting Started (Chapters 1–3)

| # | Section | File | Status |
|---|---------|------|--------|
| 1 | About Mathilda | `chapters/01-about.tex` | **Verified** (pilot) |
| 2 | Compiling and Running Mathilda | `chapters/02-building.tex` | **Verified** (first pass) |
| 3 | A Brief Introduction to Mathilda | `chapters/03-introduction.tex` | **Verified** (first pass) |

**Ch. 2 scope.** The practical on-ramp: obtaining the source, the toolchain (the
GCC-not-clang rule), required vs. optional libraries and the `USE_*` graceful-
degradation flags, per-OS dependency installation, building, the first interactive
session, how the REPL works (In/Out history, `%`, `;` suppression, Readline editing,
`?name`/`Names`/`Information` help, quitting), the three run modes (`-file` scripts,
the interactive REPL, the NDJSON pipe), and building/running the test suite. Every
Mathilda transcript verified; shell commands shown as hand-authored `shell` listings.

**Ch. 3 scope.** A whirlwind tour across the whole system, breadth over depth, meant
to grow over time. First pass (written) covers: numbers exact/approximate, number
theory, algebra (expand, collect, factor over Q and over algebraic extensions),
calculus (derivatives, partials, limits, integrals), numerical calculus (ND,
NIntegrate, NSum, FindRoot), and programming (pattern matching, procedural,
functional). Every example verified. Natural future additions: linear algebra, special
functions, and a first plot. **Deferred:** plotting the §3.5 NDSolve solution as a
figure — waiting on native Mathilda graphics export (coming soon), rather than
rendering sampled data through an external toolchain.

## Chapter 4 — Mathematics in Mathilda (sections 4.1–4.7)

This is one chapter (`chapters/04-mathematics.tex`) whose sections are the individual
mathematical domains, each its own file under `chapters/math/` and its own campaign.

| # | Section | File | Status |
|---|---------|------|--------|
| 4.1 | Arithmetic | `chapters/math/arithmetic.tex` | **Verified** |
| 4.2 | Algebra | `chapters/math/algebra.tex` | **Verified** |
| 4.3 | Calculus | `chapters/math/calculus.tex` | Planned |
| 4.4 | Linear Algebra | `chapters/math/linear-algebra.tex` | Planned |
| 4.5 | Numerical Calculus | `chapters/math/numerical-calculus.tex` | Planned |
| 4.6 | Number Theory | `chapters/math/number-theory.tex` | Planned |
| 4.7 | Special Functions | `chapters/math/special-functions.tex` | Planned |

**§4.1 Arithmetic scope** (broadened beyond the original three-subsection outline):
- **Types & operators** — the numeric heads (`Head`), and how `+ - * / ^` reduce to
  the `Plus`/`Times`/`Power` primitives (`FullForm`).
- **Machine-precision** — internal sizes/limits of machine ints and doubles; that
  Mathilda does **not** use significance arithmetic; the wide-exponent machine real.
- **Arbitrary precision** — how machine, GMP bignum, and MPFR computations mix and
  hand off; automatic promotion/demotion; `N[expr, prec]`; precision literals; contagion.
- **Complex numbers** — `Complex[re,im]`, exact Gaussian integers vs inexact complex,
  `Re`/`Im`/`ReIm`/`Conjugate`/`Abs`/`Arg`/`ComplexExpand`.
- **Rounding & integer division** — `Floor`/`Ceiling`/`Round` (banker's), `IntegerPart`/
  `FractionalPart`, floored `Mod`/`Quotient`/`QuotientRemainder`, `GCD`/`LCM`.
- **Rational parts & cleanup** — `Numerator`/`Denominator`/`Rationalize`/`Chop`.
- **Interval arithmetic** — overview, classical examples, advanced uses (dependency
  problem), and a comparison with the ISO standard (IEEE 1788-2015: decorations, flavors).

**§4.2 Algebra scope.** Expand/Factor/Together/Cancel/Apart; polynomials over ℚ and
over algebraic number fields; GCD/resultants; Gröbner bases; Solve for polynomial
systems. Show the algorithms (square-free decomposition, Hensel lifting) in
`underhood` boxes.

**§4.3 Calculus scope** (per outline): 4.3.1 Derivatives, 4.3.2 Limits, 4.3.3 Series,
4.3.4 Residues, 4.3.5 Indefinite integration (the Risch story), 4.3.6 Definite
integration (contour methods, Mellin, symmetry). A large section — likely split into
several campaigns.

**§4.4 Linear Algebra scope.** Vectors/matrices as expressions; Dot/Det/Inverse; the
decompositions (LU/QR/Cholesky/SVD/Schur); eigen-problems; the machine-precision
(LAPACK) vs exact/MPFR split; packed arrays for speed.

**§4.5 Numerical Calculus scope.** `ND`, `NIntegrate`, `NSum`/`NProduct`, `NDSolve`,
`NLimit`/`NSeries`, `FindRoot`; accuracy/precision goals; when and why to go numeric.

**§4.6 Number Theory scope.** Primes, factorization pipeline (Rho, P±1, ECM),
modular arithmetic, Diophantine equations (`Solve[..., Integers]`), continued
fractions.

**§4.7 Special Functions scope.** Gamma/Zeta/PolyGamma families, Bessel/Airy,
error functions, hypergeometric functions; exact values, symbolic identities, and
rigorous `acb` numerics via FLINT.

## The System (Chapters 5–9)

| # | Chapter | File | Status |
|---|---------|------|--------|
| 5 | Graphics | `chapters/05-graphics.tex` | Planned |
| 6 | Data Structures | `chapters/06-data-structures.tex` | Planned |
| 7 | Programming in Mathilda | `chapters/07-programming.tex` | Planned |
| 8 | Compilation and the Compiler | `chapters/08-compilation.tex` | Planned |
| 9 | Data I/O | `chapters/09-data-io.tex` | Planned |

**Ch. 5 Graphics scope.** `Plot`/`Plot3D`/`ListPlot`/`Graphics`/`Show` primitives,
options, the adaptive sampler. **Toolchain note:** graphics return `plot`/`image`
messages over the pipe, not text — this campaign must extend `build_examples.py` to
capture and embed rendered figures (a new deliverable, not just prose).

**Ch. 6 Data Structures scope.** Expression trees (everything is an expression);
Lists; Associations; NDArrays and packed arrays — the substrate that makes numeric
work fast. Heavy `underhood`/`performance` content.

**Ch. 7 Programming scope** (per outline): pattern matching, procedural programming,
functional programming — the three paradigms Mathilda supports, and how they compose.

**Ch. 8 Compilation scope.** `Compile[]`, auto-compilation, the bytecode VM, the
compilable subset as a cliff, and how packed arrays and compilation combine. (This is
the `Compile[]` *builtin*, distinct from Chapter 2's build-the-software material.)

**Ch. 9 Data I/O scope.** *Mostly not implemented yet in Mathilda* — write this last,
and only cover what actually exists; be explicit about what is not yet supported.

## The Project (Chapters 10–13)

| # | Chapter | File | Status |
|---|---------|------|--------|
| 10 | The Internals of Mathilda | `chapters/10-internals.tex` | Planned |
| 11 | The Development of Mathilda | `chapters/11-development.tex` | Planned |
| 12 | Contributing to Mathilda | `chapters/12-contributing.tex` | Planned |
| 13 | About the Author | `chapters/13-about-the-author.tex` | Planned |

**Ch. 10 Internals scope.** The pipeline (parser → evaluator → printer), the symbol
table, attributes, the pattern matcher, the rule engine — a reader's tour of
`SPEC.md`.

**Ch. 11 Development scope.** How Mathilda is built, including the AI-in-the-loop
workflow; the testing, benchmarking, and audit machinery (the packed-array/compile
gates).

**Ch. 12 Contributing scope.** Adapted from `docs/extending.md`: adding a builtin, a
module, a pattern, an operator; the coding standards and the C99 portability gate.

**Ch. 13 About the Author.** Short.

## Appendices

| # | Section | File | Status |
|---|---------|------|--------|
| A | System Architecture at a Glance | `chapters/appendix-architecture.tex` | **Verified** |

**Appendix A** is a high-level map of the system: the evaluation pipeline
(parser → evaluator → printer), the "everything is an expression" data model
(with verified `FullForm`/`Head` transcripts), how a call is evaluated
(attributes, the step sequence, fixed-point rule composition), the pattern
matcher and rule engine, a condensed source-tree diagram, and the startup /
module-initialization flow. It distills `SPEC.md` for a reader and forward-refers
to Ch. 10 (Internals) for depth. Further appendices can follow the same pattern.

## Back matter

Colophon (**Done** — stub in `frontmatter/colophon.tex`), Bibliography (grows with
`references.bib`), Index (populate `\index{}` entries as chapters are written).
