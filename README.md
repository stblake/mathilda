# Mathilda

Mathilda is a small, open source computer algebra system (CAS) heavily inspired by the core architecture and evaluation semantics of Mathematica. Written entirely in C99 and its own language, it implements a recursive expression model, structural pattern matching with backtracking, rewriting rules, and an extensive library of built-in mathematical functions. 

Today Mathilda spans over **370,000 lines of C99** across **577 source modules** (947 files including headers), exposing **more than 780 built-in functions** organized into **36 functional categories** — from arbitrary-precision arithmetic and symbolic calculus through polynomial factorization, dense linear algebra, a `Compile[]` bytecode compiler with transparent packed-array acceleration, numerical integration and optimization, integer and Diophantine equation solving, graph theory, machine learning, and interactive 2D/3D graphics.

## 🌟 Key Features

**Evaluation engine**
* **Infinite evaluation semantics:** expressions are repeatedly evaluated top-down until a fixed point is reached.
* **Attribute-driven evaluator:** a small generic core consults per-symbol bitflags (`HoldFirst`/`HoldAll`, `Flat`, `Orderless`, `Listable`, `OneIdentity`, `Protected`, …) to decide how to process each call.

**Pattern matching & rules**
* **First-class pattern matching:** `Blank` (`_`), `BlankSequence` (`__`), `BlankNullSequence` (`___`), named bindings (`x_`, `x_h`), `Condition` (`/;`), `PatternTest` (`?`), `Optional`, and `Repeated` — with full sequence backtracking.
* **Rule engine:** transformation rules (`->`, `:>`) and replacement operators (`/.`, `//.`, `Replace`).

**Numbers**
* **Arbitrary-precision integers** via the GNU Multiple Precision Arithmetic Library (GMP), with automatic promotion/demotion from machine integers.
* **Exact rationals and complex numbers**, plus **MPFR-backed arbitrary-precision reals** with precision/accuracy tracking (`N[expr, prec]`, precision literals such as `` 3.98`50 ``).
* **Interval arithmetic** (`Interval`) for rigorous enclosures.

**Symbolic mathematics**
* **Calculus:** symbolic differentiation (`D`, `Dt`, `Derivative`); multi-method integration (`Integrate`) cascading rational-function, Risch–Norman, Cherry/Liouvillian, radical-substitution, derivative-divides, and CRC integral-table methods, plus definite/contour integration; `Series`, `Limit` (with a Gruntz most-rapidly-varying engine), and symbolic summation and products (`Sum`, `Product` — finite and infinite families via Gosper, hypergeometric, telescoping, and zeta algorithms).
* **Polynomials:** univariate and multivariate arithmetic, factorization (square-free decomposition, Hensel lifting, irreducibility testing), algebraic-number factoring over ℚ(α), Gröbner bases (incl. over finite fields), GCD/LCM, resultants, and partial fractions — with optional FLINT acceleration.
* **Linear algebra:** `Det`, `Inverse`, `Dot`, `Cross`; LU / QR / Cholesky / SVD / Schur decompositions; eigenvalues and eigenvectors via multiple algorithm kernels; norms, rank, and condition numbers — with optional LAPACK/BLAS acceleration for machine-precision work.
* **Special functions:** `Gamma`, `LogGamma`, `PolyGamma`, `Zeta`, `HurwitzZeta`, `LerchPhi`, `PolyLog`, the Bessel and Airy families, `Erf`/`Erfc`/`Erfi`, exponential/log/sine/cosine integrals, `HypergeometricPFQ` (and `0F1`/`1F1`/`2F1`), Legendre, `BernoulliB`, `EulerE`, and more — many with rigorous `acb` numerics when built with FLINT.
* **Simplification:** `Simplify`/`FullSimplify` with a complexity-driven search, trigonometric identities and rationalization, radical denesting, and an assumptions framework (`$Assumptions`, `Assuming`, `Element`).

**Equation solving**
* **Over the reals and complexes:** `Solve` for polynomial (via radicals and `Root` objects), radical, transcendental, and nonlinear systems; `SolveAlways`, `Eliminate`, `ToRadicals`, and `GroebnerBasis`.
* **Diophantine (`Solve[…, Integers]`):** a dedicated integer-equation engine — linear systems via Hermite normal form, Pell and generalized Pell equations, Thue equations through a number-field/unit layer (Voronoi and Minkowski–LLL unit search), the sum-of-three-cubes search (Booker/Heath-Brown), Ramanujan–Nagell, ternary quadratic forms, Egyptian fractions, and power-sum searches (Lander–Parkin, Frye).

**Numerical analysis**
* **Quadrature & sums:** `NIntegrate` (adaptive, oscillatory, Levin-collocation, MPFR-precision), `NSum`/`NProduct` (Euler–Maclaurin and convergence acceleration), `NResidue`.
* **Differential equations:** `NDSolve` for ODE initial-value problems and PDE systems (with upwind schemes and optional BLAS/LAPACK acceleration).
* **Optimization:** global `NMinimize`/`NMaximize` (differential evolution, simulated annealing, basin hopping, SHGO, dual annealing) and local `FindMinimum`/`FindMaximum` (L-BFGS-B, Nelder–Mead, Powell, SLSQP, COBYLA/COBYQA, trust-region), plus mixed-integer and constrained problems.
* **Roots & local analysis:** `FindRoot`, `NSolve`, `NRoots`, `NLimit`, `NSeries`, `ND`.

**High-performance numerics**
* **Packed arrays:** ordinary `List`s of machine numbers are transparently backed by dense `NDArray` storage, and element-wise, reduction, and structural kernels operate directly on the buffer — optionally multithreaded across cores — for order-of-magnitude speedups over boxed evaluation.
* **`Compile[]`:** a bytecode compiler that lowers numeric expressions to a register VM at both scalar and rank-1 array shapes, with auto-compilation of `Table`/`Sum`/`NIntegrate` bodies, and `CompileDiagnostics` for inspecting what lowered.
* **Optional acceleration backends:** LAPACK/BLAS for dense linear algebra and FFTW for `Fourier`.

**Number theory & factorization**
* **Number theory:** `GCD`, `LCM`, `ExtendedGCD`, `PowerMod`, `Divisors`, `DivisorSigma`, `EulerPhi`, `MoebiusMu`, `PrimitiveRoot`, `MultiplicativeOrder`, continued fractions, and more.
* **Integer factorization:** a unified, automatic pipeline alongside explicit algorithms — Pollard's Rho, Pollard's $P-1$, Williams' $P+1$, Fermat, CFRAC, Dixon's Method, and the Elliptic Curve Method (ECM).

**Data structures**
* **Associations** (`<|…|>`) with O(1) hashed lookup: `Keys`, `Values`, `Lookup`, `KeySort`, `AssociationThread`, and the aggregation vocabulary `GroupBy`, `Merge`, `Counts`, `GatherBy`, `Tally`.

**Graph theory**
* `Graph` with directed/undirected/weighted edges, generators (`CompleteGraph`, `CycleGraph`, `RandomGraph`), queries (`AdjacencyMatrix`, `VertexDegree`, `ConnectedComponents`), algorithms (`FindShortestPath`, `GraphDistance`, `FindSpanningTree`), and `GraphPlot`.

**Machine learning**
* High-level `Classify`/`Predict` with `ClassifierFunction`/`PredictorFunction`, clustering (`FindClusters`, `ClusteringComponents`), dimensionality reduction (`PrincipalComponents`, `DimensionReduce`), `LinearModelFit`, `Nearest`, and `LearnDistribution`.

**Signal & image processing**
* **Fourier:** `Fourier`/`InverseFourier`, `FourierDCT`/`FourierDST` (FFTW-backed when available).
* **Images:** `Image`/`Image3D`, filtering (`GaussianFilter`, `ImageConvolve`, `EdgeDetect`), morphology (`Dilation`, `Erosion`), `Binarize`, `ColorConvert`, `ImageResize`, and histograms.

**Graphics & visualization** *(requires [Raylib](https://www.raylib.com/); gracefully omitted otherwise)*
* **2D plots:** `Plot` with adaptive sampling and re-sampling on zoom/pan; `ParametricPlot` for curves and filled regions; `StreamPlot` for vector fields with speed-gradient colouring; `ListPlot`/`ListLinePlot` for discrete data; `ContourPlot` for iso-contours with marching-squares (equation form `f == c`, list of equations `{eq1, eq2, …}`, and numeric-body shading with `ContourShading`/`ColorFunction`).
* **3D plots:** `Plot3D` for surface meshes; `ParametricPlot3D` for parametric space curves and surface patches; both rendered with per-face Lambertian shading in an interactive orbit camera (drag to rotate, scroll to zoom, right-drag to pan).
* **Graphics primitives:** hand-built `Graphics[…]` and `Graphics3D[…]` objects using `Line`, `Point`, `Arrow`, `Disk`, `Rectangle`, `Polygon`, `Text`, `RGBColor`, `Opacity`, `Thickness`, and more; combined with `Show`.
* **Legends and labels:** `PlotLegends -> Automatic` adds a colour-scale bar (contour/stream plots) or per-curve swatch box; `AxesLabel`, `PlotLabel`, `GridLines`, `Frame`, `Ticks` all pass through to the renderer.
* **Interactive window:** toolbar with close, screenshot save, and view-reset buttons; Escape or the OS close button exits cleanly.

**Programming & utilities**
* **Functional programming:** `Map`, `Apply`, `Fold`, `Nest`, `Through`, `Composition`, and pure functions (`&` / `#`).
* **Scoping & control flow:** `Module`, `Block`, `With`; `If`, `Which`, `Switch`, `Do`, `For`, `While`, `Piecewise`.
* **Standard library:** lists and iteration, strings, statistics, random numbers, date/time, and file I/O.

---

## 📚 Function Categories

The complete reference (780+ functions across 36 categories) lives in [`Mathilda_spec.md`](Mathilda_spec.md), which indexes the per-category pages under [`docs/spec/builtins/`](docs/spec/builtins/):

* Arithmetic · Algebra · Number Theory
* Calculus · Simplification · Power Series
* Solutions of Equations
* Numerical Calculus
* Elementary Functions · Special Functions · Mathematical Constants
* Linear Algebra · LAPACK · BLAS
* Packed Arrays
* Fourier Transforms
* Statistics · Machine Learning · Random Number Generation
* Graphs
* Image Processing
* Data Structures (Associations)
* Structural Manipulation · Expression Information
* Lists and Iteration · Functional Programming
* Pattern Matching · Comparisons
* Control Flow · Assignment and Rules · Scoping Constructs
* String Operations · File I/O · Time and Date
* Graphics & Visualization
* FLINT context (direct FLINT-kernel access)

Weekly change summaries are recorded under [`docs/spec/changelog/`](docs/spec/changelog/).

---

## 🚀 Getting Started

### Prerequisites

To build and run Mathilda you need:

* **GCC** (a real GCC — the build deliberately rejects Apple/LLVM clang, which hides the glibc-portability warnings the Linux CI gates on; on macOS install `gcc` via Homebrew)
* **GMP** (`libgmp` / `gmp-dev`) — arbitrary-precision integers *(required)*
* **GNU Readline** (`libreadline` / `readline-dev`) — interactive line editing *(required)*
* **MPFR** (`libmpfr` / `mpfr-dev`) — arbitrary-precision reals *(enabled by default)*
* **FLINT** ≥ 3.0 (`libflint` / `flint-dev`) — fast, rigorous polynomial arithmetic over algebraic extensions and rigorous `acb` numerics *(optional, auto-detected)*
* **GMP-ECM** (`gmp-ecm` / `libecm-dev`) — Elliptic Curve Method integer factorization *(optional, auto-detected)*
* **LAPACK / BLAS** — fast machine-precision linear algebra *(optional, auto-detected)*
* **FFTW** ≥ 3 (`libfftw3` / `fftw3-dev`) — fast `Fourier`/`FourierDCT`/`FourierDST`; falls back to a naive $O(n^2)$ transform when absent *(optional, auto-detected via `pkg-config`)*
* **Raylib** ≥ 4.0 — interactive graphics window for `Plot`, `Plot3D`, `ContourPlot`, etc. *(optional, auto-detected via `pkg-config`; falls back to a text placeholder when absent)*
* **CMake** — only required to build the test suite

POSIX threads (used to parallelize the packed-array kernels across cores) are enabled by default on macOS and Linux and need no extra package.

The optional backends are controlled by build-time flags and **degrade gracefully** when disabled or absent:

| Flag | Default | Effect when on |
|------|---------|----------------|
| `USE_MPFR`     | `1` | Arbitrary-precision reals: `N[expr, prec]`, `Precision`/`Accuracy`, precision literals. Build without it via `make USE_MPFR=0`. |
| `USE_FLINT`    | `1` | Fast, rigorous FLINT (≥ 3.0) kernels: multivariate polynomial GCD/factoring over ℚ, univariate GCD/factoring over number fields ℚ(α) (via the `gr` layer + ANTIC), the finite-field workhorse behind parametric ℚ(t)(α) work, and rigorous `acb` numerics (`Zeta`, `HurwitzZeta`, `PolyGamma`, `StieltjesGamma`). Auto-detected via `pkg-config` with a ≥ 3.0 version floor. Falls back to the classical (slower but still rigorous) path (`USE_FLINT=0`) when absent. |
| `USE_LAPACK`   | `1` | Fast machine-precision linear algebra. Auto-detected: Apple **Accelerate** on macOS, `lapacke`/`lapack`/`blas` on Linux. Falls back to the pure-C path (`USE_LAPACK=0`) if none is found. |
| `USE_ECM`      | `1` | Elliptic Curve Method factorization via the system GMP-ECM library. Auto-detected via a compile-link probe; install `gmp-ecm` / `libecm-dev`. Falls back to disabled (`USE_ECM=0`) when absent. |
| `USE_FFTW`     | `1` | FFTW-backed `Fourier`/`FourierDCT`/`FourierDST`. Auto-detected via `pkg-config fftw3`. Falls back to a naive $O(n^2)$ transform (`USE_FFTW=0`) when absent. |
| `USE_THREADS`  | `1` | POSIX-thread parallelism for the element-wise packed-array kernels (`Sin`, `Exp`, `Erf`, …) on large arrays. Enabled on macOS/Linux; build the serial path via `make USE_THREADS=0`. |
| `USE_GRAPHICS` | `1` | Interactive 2D/3D plot windows via Raylib. Auto-detected via `pkg-config raylib`. When absent, `Show`/`Plot`/`Plot3D`/`ContourPlot`/etc. print a text placeholder and return normally. Build without it via `make USE_GRAPHICS=0`. |

#### Installing dependencies

**Linux (Debian / Ubuntu):**

```bash
# Required libraries
sudo apt install libgmp-dev        # GMP — arbitrary-precision integers
sudo apt install libmpfr-dev       # MPFR — arbitrary-precision reals
sudo apt install libreadline-dev   # GNU Readline — interactive REPL

# Optional: FLINT (>= 3.0) for fast, rigorous algebraic-extension arithmetic
sudo apt install libflint-dev      # Debian Bookworm+/Ubuntu 24.04+ ship >= 3.0

# Optional: GMP-ECM for advanced integer factorization
sudo apt install libecm-dev

# Optional: LAPACK / BLAS for fast machine-precision linear algebra
sudo apt install liblapacke-dev libopenblas-dev

# Optional: FFTW for fast Fourier / FourierDCT / FourierDST
sudo apt install libfftw3-dev

# Optional: Raylib for interactive plot windows (Plot, Plot3D, ContourPlot, ...)
sudo apt install libraylib-dev      # Ubuntu 24.04+ / Debian Bookworm+
# or build from source: https://github.com/raysan5/raylib

# Optional: CMake, only needed to build the test suite
sudo apt install cmake
```

On Fedora/RHEL the equivalents are `gmp-devel`, `mpfr-devel`, `readline-devel`,
`flint-devel` (≥ 3.0), `gmp-ecm-devel`, `lapack-devel`/`openblas-devel`,
`fftw-devel`, plus `cmake`.

> **Note on FLINT versions.** Mathilda requires **FLINT ≥ 3.0** (the release that
> merged ANTIC for number-field arithmetic). Distributions that only package
> FLINT 2.x — e.g. Ubuntu 22.04 or Debian Bullseye — are detected as too old and
> the build automatically falls back to `USE_FLINT=0`. Install a newer FLINT from
> source or a backport if you want the accelerated paths on those systems.

**macOS (Homebrew):**

```bash
brew install gmp mpfr readline cmake
# Optional: FLINT (>= 3.0) for fast, rigorous algebraic-extension arithmetic:
brew install flint
# Optional: FFTW for fast Fourier / FourierDCT / FourierDST:
brew install fftw
# Optional: Raylib for interactive plot windows:
brew install raylib
# Optional: GMP-ECM for advanced integer factorization:
brew install gmp-ecm
```

LAPACK/BLAS need not be installed on macOS — the build auto-detects Apple's
**Accelerate** framework.

### Building Mathilda

The `makefile` auto-discovers `src/*.c`, configures and compiles internal dependencies, then links the main executable (`-std=c99 -O3`).

1. Clone the repository:
   ```bash
   git clone https://github.com/stblake/Mathilda.git
   cd Mathilda
   ```
   Install GMP-ECM (used for advanced integer factorization) from your package
   manager — `brew install gmp-ecm` on macOS or `sudo apt install libecm-dev`
   on Debian/Ubuntu. The build autodetects it and links `-lecm`; if it is
   absent, the build still succeeds with advanced factorization disabled
   (equivalent to `make USE_ECM=0`).
2. Build the project:
   ```bash
   make -j$(nproc)
   ```
   To build a leaner binary, disable optional backends, e.g.
   `make -j$(nproc) USE_LAPACK=0 USE_MPFR=0`.
3. Start the interactive REPL:
   ```bash
   ./Mathilda
   ```
   Or run a `.m` file as a script and exit — nothing is echoed, so the output
   is exactly what the file `Print`s:
   ```bash
   ./Mathilda -file script.m       # a bare ./Mathilda script.m works too
   ./Mathilda --help               # all options
   ```

### Running the Test Suite

Mathilda ships a comprehensive C-based unit suite — **461 `test_*.c` files**
covering evaluation, parsing, pattern matching, arithmetic, polynomials,
calculus, equation solving, numerical analysis, linear algebra, packed arrays,
the `Compile[]` compiler, and more. CMake auto-detects the same optional
backends (MPFR, FLINT, LAPACK, FFTW, ECM) during configuration.

```bash
cd tests
mkdir -p build && cd build
cmake ..
make -j$(nproc)

# Run all test binaries
for t in *_tests; do ./$t; done
```

---

## 🛠️ Developer Guide & Architecture

Everything in Mathilda is represented by an immutable-by-convention `Expr` AST
node — a tagged union over `Integer`, `Real`, `BigInt`, `Symbol`, `String`,
`Function`, and (when built with MPFR) arbitrary-precision `MPFR` reals. Compound
values are `Function` nodes: `Rational[n, d]`, `Complex[a, b]`, and lists
(`List[…]`) are all expressions, so the same generic tools (`Part`, `Map`,
`ReplaceAll`, …) operate uniformly on everything.

The system is modularized into several independent subsystems:

1. **Parser (`parse.c`)** — a Pratt parser mirroring Mathematica's operator precedences (inline lexing, no separate tokenizer).
2. **Evaluator (`eval.c`)** — the fixed-point evaluation loop; applies `Hold*`, `Flat`, `Orderless`, and `Listable` before recursively evaluating arguments.
3. **Symbol table (`symtab.c`)** — `OwnValues`, `DownValues` (user rules), attributes, docstrings, and native C built-in function pointers.
4. **Pattern matcher & rule engine (`match.c`, `replace.c`)** — structural tree unification with sequence segmenting and backtracking.

Larger mathematical domains live in dedicated subdirectories of `src/`:

| Subsystem | Responsibility |
|-----------|----------------|
| `poly/`     | Univariate/multivariate polynomial arithmetic, factorization, algebraic-number fields, Gröbner bases |
| `linalg/`   | Dense linear algebra; eigen kernels split by algorithm |
| `calculus/` | `D`/`Dt`/`Derivative`, `Series`, `Limit`, `Integrate` (incl. Risch–Norman, Cherry/Liouvillian) |
| `simp/`     | `Simplify`, trigonometric simplification, radical denesting, assumptions |
| `sum/`, `product/` | Symbolic summation and products — finite and infinite families |
| `solve/`    | `Solve` over Complexes/Reals plus the `Solve[…, Integers]` Diophantine engine (one file per method) |
| `special_functions/` | Gamma, Zeta, Bessel, Airy, error/exponential integrals, hypergeometric families |
| `numerical_calculus/`, `numerical_roots/` | `NIntegrate`, `NDSolve`, `NSum`/`NProduct`, `NMinimize`/`FindMinimum` family; `NSolve`, `NRoots`, `FindRoot` |
| `ndarray.c`, `pack.c`, `compile/` | Packed machine-precision arrays, the transparency gate, and the `Compile[]` bytecode compiler + auto-compilation |
| `graph/`    | Graph data structure, generators, and algorithms |
| `ml/`, `stats/` | Machine-learning primitives and descriptive statistics |
| `graphics/` | 2D/3D plot engine: adaptive sampler, marching-squares contours, Raylib renderer, vector font; `Plot`, `Plot3D`, `ParametricPlot`, `ParametricPlot3D`, `StreamPlot`, `ContourPlot`, `ListPlot`, `Show` |
| `internal/` | Mathematica-syntax bootstrap `.m` files (init, integral tables) loaded at startup |

A recurring design pattern is **C for performance, rules for mathematics**: hot
paths (parser, evaluator, matcher, arithmetic) are C, while higher-level
identities (integral tables, etc.) are expressed as `DownValues` in Mathilda's
own language. The full architecture is documented in [`SPEC.md`](SPEC.md), and
extension recipes in [`docs/extending.md`](docs/extending.md).

### Extending Mathilda: Adding a New Built-in Function

Adding new functionality to Mathilda is straightforward:

1. **Write the C implementation.**
   Create your evaluation logic in the appropriate `.c` module (e.g., `core.c`).
   Your function signature must be `Expr* builtin_myfunc(Expr* res)`.

   * **Memory rule (ownership contract):** the builtin **takes ownership** of
     `res`. On success, return a **new** `Expr*` — the **evaluator** frees `res`
     for you, so you must **not** call `expr_free(res)` yourself (doing so causes
     a double-free). If you cannot evaluate the input (e.g. symbolic arguments to
     a purely numeric function), return `NULL` **without freeing `res`**, and the
     evaluator retains ownership, leaving the expression unevaluated. When you
     reuse parts of `res` in your result, NULL them out first so the evaluator's
     cleanup doesn't free them twice.

   ```c
   Expr* builtin_myfunc(Expr* res) {
       if (res->data.function.arg_count != 1) return NULL;  /* leave unevaluated */
       // ... mathematical logic ...
       return expr_new_integer(42);  /* evaluator frees res — do NOT free it here */
   }
   ```

2. **Register the function.**
   In the module's initialization routine (e.g., `core_init()`), register the
   function and assign a documentation string so it is available via `?MyFunc`:
   ```c
   symtab_add_builtin("MyFunc", builtin_myfunc);
   symtab_set_docstring("MyFunc", "MyFunc[x]\n\tComputes the ultimate answer.");
   ```

3. **Assign attributes.**
   If your function threads over lists, operates symmetrically, or holds its
   arguments unevaluated, set the corresponding attributes during initialization:
   ```c
   symtab_get_def("MyFunc")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
   ```

4. **Test and document.**
   * Add test cases to the appropriate suite in `tests/` using the `TEST(...)` macro.
   * Document the function in the matching file under [`docs/spec/builtins/`](docs/spec/builtins/) and add a note to the current week's [`docs/spec/changelog/`](docs/spec/changelog/) entry.

See [`docs/extending.md`](docs/extending.md) for the full recipes (modules,
patterns, internal `.m` rules, operators) and [`CLAUDE.md`](CLAUDE.md) for the
contributor workflow.

---

## 📜 Open Source & License

Mathilda is open-source software licensed under the **GNU General Public License v3.0 (GPLv3)**.

You are heavily encouraged to explore the codebase, submit pull requests, report issues, and expand the CAS with new mathematical algorithms! Please see the `LICENSE` file for more details.
