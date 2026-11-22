# The Mathilda Book — Update Review (git history since the last book pass)

> **STATUS: IMPLEMENTED 2026-11-22** (book now current through v0.6664; 702 pages, was 564; `make pdf` clean, no undefined references, `check-links` OK). Window `a4bb2d198..541769142`, Mathilda
> **v0.507 → v0.6664**. The previous review (v0.290–v0.507, implemented 2026-10-20) is
> in git history (`git log -- book/BOOK_UPDATE_REVIEW.md`).

This is a **triage backlog**, not prose. It maps every *book-relevant* capability
that landed after the book's last substantive pass onto the chapter that must cover
it. Read `CONTEXT.md` first; this file is scoped to *the delta*. Tick `[x]` as rows land.

## Anchor and window

- **Book anchor:** `a4bb2d198` (2026-10-20, v0.507), "show the natural inputs now that
  v0.482–v0.507 fixed the kernel gaps". Later book commits were narrow: `a2b899b9c`
  (REPL features, 2026-11-11) and `6d48ea365` (relicense GPLv3 → Apache 2.0, 2026-11-20).
- **Window:** ~1,040 commits, v0.507 → v0.6664, 2026-10-20 → 2026-11-22. Sources: the
  git log and `docs/spec/changelog/2026-10-19.md` (after v0.507), `2026-10-26.md`,
  `2026-11-02.md`, `2026-11-09.md`, `2026-11-16.md` (~17.5k lines).
- **Method:** three triage passes:
  1. v0.508–v0.6119;
  2. v0.6120–v0.6664;
  3. deep dives on integration (Risch/Trager/Bronstein) and LinkedList.

  All 350 example sessions were re-run on v0.6664 and diffed against `book/generated/`.

| Tag | Meaning |
|-----|---------|
| **REVIEW** | A behaviour change that falsifies existing prose or an example. Fix first. |
| **INCLUDE** | New capability with no book coverage. |
| **UPDATE** | Extends a section that already covers the area. |

## Structural decisions

- **New Part II chapter "Symbolic integration"** (`chapters/integration.tex`). It is built from
  calculus.tex §"Indefinite integration: the Risch story", and the `calculus/integrate-*` and
  `elementary-q` examples move to `examples/integration/`.
- **New Part II chapter "Certified answers"** (`chapters/certified-answers.tex`), for Lean 4.
- **Special functions: full coverage.** Every head added since v0.507 gets a pair.

## REVIEW — things the book says that are no longer true

| # | Where | What changed (version) | Fix |
|---|-------|------------------------|-----|
| [x] R1 | all transcripts | Canonical order (v0.597); factored `Together` (v0.6073); Factor sign convention (v0.6063/6132); decimal precision marks and `*^` printing (v0.6292/6604); `Import` of .txt gives a String (v0.551) | `make examples`; spot-read every changed transcript against its caption |
| [x] R2 | `02-building.tex` :248-256, :340, :368, :400; `01-about` :103; `10-internals` :791 | The in-tree line editor replaces readline on a tty; new banner; Ctrl-C gives `$Aborted`, and a second Ctrl-C within 1 s exits 130 (v0.6538) | Rewrite the REPL section |
| [x] R3 | `algebra.tex` :142, `algebra/factor/3-4` | `FactorSquareFree` gives `x^2 (1 + x)^2` (v0.6066) | New pair |
| [x] R4 | `algebra.tex` ~:370 callout, `algebra/rational/1` | ratfun engine; `Together`/`Cancel`/`Apart` keep factored denominators (v0.6072-75) | Rewrite the callout |
| [x] R5 | `algebra.tex` :620, :675 (`reduce/6`, `reduce-transcendental-roots/5`) | `Sqrt[1-x^2]`; `Reduce[E^x > 3x]` gives `ProductLog` breakpoints (v0.6491) | Reframe as a feature, or choose a truly transcendental example |
| [x] R6 | `arithmetic.tex` :202 (`precision/2`), `10-internals` `numbers-tower/10` | `Precision[N[Pi,50]]` is `50.` (16 guard bits hidden) (v0.6292, v0.6596) | Rewrite the note |
| [x] R7 | `10-internals.tex` :723-738 (`printer-roundtrip/2`, `printer-forms/4`) | 17-digit literals are machine reals; machine InputForm has no backtick (v0.6136, v0.6604); InputForm writes `*` (v0.6133) | Rewrite |
| [x] R8 | `10-internals.tex` :275, :308-325, :558, :589 | Canonical order (v0.597); dependency-recorded fixed points (v0.6053/6119); generality ordering of DownValues (v0.6022); System`/Mathilda` contexts (v0.533/565) | Rewrite, with Under-the-hood callouts |
| [x] R9 | `07-programming.tex` :280-288, :1074, :1097-1121 | Generality ordering, not specificity (v0.6010/6022); `TimeUsed` (`timing/6`); TimeConstrained timeout is an ordinary Abort, and the Windows paragraph is out of date (v0.6623) | Rewrite |
| [x] R10 | `number-theory.tex` ~:645 | Discrete log uses Pohlig–Hellman + BSGS/ρ (v0.524), not an orbit walk | Rewrite; show a large log |
| [x] R11 | `linear-algebra.tex` :134 (`det/2`) | `Det` returns the expanded form (v0.6139) | Reword "Laplace along the first row" |
| [x] R12 | `special-functions.tex` :136 (`polygamma/5`) | `PolyGamma[0, 1/2]` stays unevaluated (v0.6131) | Wrap in `FunctionExpand` |
| [x] R13 | `calculus.tex` :1122 (`dsolve-pde/3`) | Wave IVP gives `Cos[c t] Sin[x]` (v0.6496-98) | Reword the d'Alembert note |
| [x] R14 | `optimization.tex` :443 | FindMinimum full-precision argmin; Brent's evaluation count changed (v0.6419-60) | Use the count from the transcript |
| [x] R15 | `06-data-structures.tex` :88, `expressions/3` | Rational/Complex are atoms to every structural builtin (v0.6535) | Reword |
| [x] R16 | `09-data-io.tex` :18-22, :332, `importexport/7` | `Import` reads Text/String/Lines; `DeleteFile` etc. exist (v0.529-571) | Rewrite the intro and the pair |
| [x] R17 | `03-introduction` `algebra-factor/4` (:177) | Factor sign form `-(-1 + Sqrt[2] x - x^2)(…)` | Choose a cleaner example |
| [x] R18 | calculus.tex :526-540 Theory, :556-558, :576-593 PMT | Trager and Bronstein are decision procedures that run before PMT; structure theorem, parametric RDE, hypertangent and Liouvillian levels; `Integrate::incomplete` | Superseded by the new integration chapter |
| [x] R19 | `03-introduction` `solve-numerical/4`, `inequalities/10` | FindRoot vector backtracking finds a different root (v0.6387) | Regenerate; check the prose |

## INCLUDE / UPDATE — by chapter

**Symbolic integration (new chapter)**
- [x] I1 **Rational functions.** Keep Hermite / Lazard–Rioboo–Trager.
- [x] I2 **Transcendental tower.**
  - Log, exp, hypertangent (tan/tanh by Moebius) and Liouvillian (Erf/Ei/li) levels, over Q, Q(i), number fields and parameters.
  - Structure theorem; parametric RDE; ``Risch`RischDECertified``, ``Risch`ParametricRischDE`` (v0.6517-6652).
- [x] I3 **Algebraic: Trager** (`src/internal/trager/`, v0.609x-6587).
  - Integral basis, Hermite, residues by Galois orbit, divisor torsion mod p, LogToArcTan.
  - Parameters over Q(a); `Method -> "Trager"`/`"TragerKauers"`.
- [x] I4 **Mixed: Bronstein** (`src/internal/bronstein/`, v0.6520-6652).
  - One radical over a Log/Exp tower; Puiseux at special places; Shape B.
  - Phase 5: Risch over the curve Q(x)(y).
  - Galois components and support groups.
- [x] I5 **Cascade and guarantees.**
  - Stage order, budgets, `Integrate::nonelem` only with a proof, `Integrate::incomplete`.
  - FriCAS comparison: DDI 664/803; benchmarks 256/257/259/260.
  - Open-gaps box.
- [x] I6 **ParallelMixedTower:** Kummer split, conics with parameters, kernel substitutions (v0.6099-6112).

**Ch. 6 Data structures**
- [x] I7 **SparseArray (new §)**: `EXPR_SPARSE`, CSR at any rank, the transparency gate, native kernels, `IdentityMatrix[n, SparseArray]`, scipy parity (v0.572-586).
- [x] I8 **LinkedList (new §)**:
  - mutable deque, reference semantics, chunked storage with a clz index, typed storage;
  - O(1) `Push`/`Pop`/`PopFirst`/`AppendTo`/`PrependTo`;
  - works with every head, with closure; the no-element-arrays invariant; `check-linkedlist`;
  - compiled loops (v0.6625-6664).
- [x] U1 **Strings:** named characters and raw Unicode in source, `\!\(...\)` (v0.587, v0.6001); canonical collation (v0.530).

**Ch. 7 Programming**
- [x] I9 **Messages, aborts and non-local exits (new §):**
  - `Return`, `Abort`/`CheckAbort`/`AbortProtect`/`PropagateAborts`;
  - nested `Check`/`Quiet`; `Message`/`Off`/`On`/`$MessageList`;
  - `MemoryConstrained` (v0.558-600, v0.6623).
- [x] I10 **Packages and contexts (new §):** `BeginPackage`, `Needs`, `$Packages`, `$ContextPath`, ``Mathilda` ``, `<<file`; FeynCalc 13214/13412 (v0.532-6119).
- [x] U2 **Rules and definitions:**
  - `Dispatch`, `MaxIterations`;
  - SubValues, NValues, FormatValues/`Format`, `OptionsPattern`, `OptionQ`/`FilterRules`;
  - `Replace` whole-expression (v0.6174); `Block` localises attributes (v0.6170).
- [x] U3 **Utilities:** `Echo`/`EchoTiming`, `Dialog`, `Unique` gives `x$nnn`, parallel heads with serial semantics.

**Ch. 10 Internals (Under-the-hood callouts)**
- [x] U4 Canonical order (v0.597); dependency-recorded fixed points and INERT stamps (v0.6053/6119); scoping binder summaries (v0.6653-57).
- [x] U5 Module-local reclamation and shared nodes (v0.6056-59); constraint frames (`src/memacct.h`, v0.6623); hash index for pattern-free rules (v0.6023).

**Ch. 2 Building / Ch. 1 About**
- [x] U6 **REPL:**
  - history, completion, colour, the `took…` note;
  - status line + `Monitor`;
  - slash commands, themes, mods, `SlashCommand`;
  - `MATHILDA_EDITOR=readline` (v0.6538).
- [x] U7 `Quit[n]` in scripts (v0.527); FLINT 3.0–3.2 (v0.6116); the licence argument is consistent after Apache.

**Arithmetic**
- [x] U8 **Numeric-tower callout:**
  - guard bits; printing the exact digits; `*^`; exact zeros (`0.*^-20`);
  - adaptive `N` (`N::meprec`, v0.6031); arb correct rounding;
  - `Overflow[]`/`Underflow[]`; `$MinPrecision`/`$MaxPrecision` (v0.6264-6618).
- [x] U9 `22/7 == Pi` soundness (v0.6584); comparisons with Infinity; `Unitize`; Precision/Accuracy of lists.

**Algebra**
- [x] U10 Factor `Modulus`, multivariate `Extension`/`GaussianIntegers`, `E^x` kernels (v0.6006-6071).
- [x] U11 Gröbner: F4/GF(p), certified multi-modular, FGLM, inexact, parametric, `Method -> "F4"` (v0.6077-6463).
- [x] U12 Simplify special-function engine (v0.536-540); `SeparablePolynomialDecomposition` (v0.6177); denesting/Gauss periods (v0.6494); `CountRoots` (v0.6343); `2^a 3^a -> 6^a` (v0.594-6008).

**Number theory**
- [x] I11 Bitwise family `BitXor`/`BitAnd`/`BitOr`/`BitNot`/`BitSet`/`BitGet`/shifts (v0.508-517).
- [x] I12 `Multinomial`, `Binomial` extensions (v0.598); permutation groups `Cycles`/`PermutationGroup`/`GroupOrder`/`Permute`/`Signature` (v0.549-590); `MangoldtLambda`, `CubeRoot`.
- [x] U13 BRBD timing 62.9 s → 2.8 s (v0.525-526).

**Linear algebra**
- [x] I13 **Sparse matrices (new §):**
  - SpMV/SpMM/SpGEMM;
  - direct and Krylov solvers, `LinearSolveFunction`, band solver;
  - Krylov eigensolvers, sparse SVD (v0.572-579, v0.6180-6385).
- [x] I14 `MatrixExp`/`MatrixLog`/`MatrixFunction`, non-integer `MatrixPower`, `CholeskyDecomposition`, `HessenbergDecomposition`, `SingularValueList`, `Minors`, `Adjugate`, `Orthogonalize`, `Projection`, `TensorProduct`, `KroneckerProduct`, `LUBackSubstitution`.
- [x] U14 MPFR dense linear algebra callout (v0.6242-6308); radical and parametric `RowReduce`/`NullSpace` (v0.6107-08).

**Special functions (full coverage)**
- [x] I15 Mathieu (`MathieuC/S` and primes, `MathieuCharacteristicA/B`) (v0.6165-66).
- [x] I16 Jacobi elliptic functions, `EllipticTheta`, the Weierstrass family.
- [x] I17 Appell `F1`–`F4`; `MeijerGReduce`; `BetaRegularized`.
- [x] I18 Struve L (half-integer now elementary, v0.6204), `AngerJ`, `WeberE`, Kelvin, Hankel/spherical Bessel.
- [x] I19 Zeros: `ZetaZero`, Airy zeros, `BesselYZero`; `RiemannSiegelZ`/`Theta`, `RiemannR`.
- [x] I20 `JacobiP`, `GegenbauerC`, `ParabolicCylinderD`; plus every other head in the changelog enumeration.
- [x] U15 Every special function compiles at arbitrary precision (v0.6482).

**Calculus (non-integration)**
- [x] I21 `LaplaceTransform` (v0.6491).
- [x] U16 Roach/Meijer-G definite integration, the G&R stage, elliptic definite integrals (v0.6203-6493); Bessel-product integrals (v0.544).
- [x] U17 DSolve: PDE IVP layer, RiemannP, ScalingReduction, equations polynomial in y′; corpus 843/1000 (v0.6148-6496).
- [x] U18 `D[f, {x}]`; D through Dot/Sum/Product and Bessel order; essential-singularity `Residue`; `Limit` assumptions reach Series; Sum/Product soundness, `HarmonicNumber` (v0.6044-6553).

**Numerical calculus / optimization**
- [x] I22 NDSolve: `WhenEvent`, `NDSolveValue`, BVP shooting, arbitrary-precision Taylor, FEM over regions, `DirichletCondition`, `NDEigensystem` (v0.6233-6417).
- [x] I23 `FunctionInterpolation`, `GaussianQuadratureWeights`/`Error` (v0.6243-44).
- [x] U19 NIntegrate oscillatory/Ooura–Mori/vertex singularities; NSum; NRoots/NSolve completeness; FindRoot at WorkingPrecision; FindMinimum SQP and minimax (v0.6228-6460).

**Geometry**
- [x] I24 `RotationMatrix`, `ArcLength`, `RegionMeasure`, `Annulus` (v0.6401-06).

**Ch. 9 Data I/O**
- [x] I25 `FileNames`, `Directory`/`SetDirectory`, `ParentDirectory`, `FileNameTake`, `AbsoluteFileName`, `FileType`, `DeleteFile`, `DirectoryQ`, `ReadString`, `StringToStream`/`Skip`, `Environment`.
- [x] I26 `FortranForm`, `CForm`; `$System`, `$SystemID`, `$ProcessorCount`, `$InstallationDirectory`, … (v0.6167, v0.6374-77).

**Certified answers (new chapter)**
- [x] I27 Lean 4: `Certify`, `ProofObject`, `LeanVerify`, `LeanForm`, `Refute`, `Conjecture`, prover consensus, `lean/` tactics (v0.6536-37); its relation to PossibleZeroQ and Risch certificates.

## Skipped (no book impact)

Portability (Linux/Windows/glibc), CI and test-suite infrastructure, benchmark-only work,
docs-site work, valgrind and leak fixes, and micro-speedups. Narrow bug fixes whose only
book effect is a regenerated transcript are covered by R1.

## Kernel findings from the writing pass (follow-ups)

**From the previous pass, re-checked on v0.6664**
- FIXED v0.6204: `StruveL[1/2, x]` reduces to elementary form.
- FIXED v0.6154: `Integrate[-(8 t E^(t/8) Sin[3/8 Sqrt[7] t])/(3 Sqrt[7] k), t]` integrates (Integrate splits off x-free constants).
- Possibly fixed: the rare `compute_deriv(NULL)` segfault under load (v0.6664 fixed a TimeConstrained backstop segfault).
- Still open:
  - FIXED v0.6683: `Integrate` principal value with two symbolic interior poles;
  - NOT A BUG: `Series[Gamma[s], {s,0,2}]` keeps `PolyGamma[2, 1]` by design (`FunctionExpand` gives `-2 Zeta[3]`);
  - FIXED v0.6667: `Refine[Sqrt[a^2], a > b > 0]` is `a`;
  - FIXED v0.6668: `Simplify[C[1] E^x + C[2] x^2 E^x]` gives `E^x (C[1] + x^2 C[2])`;
  - FIXED v0.6679: `D[StruveH[0, x], x]` is `2/Pi - StruveH[1, x]`;
  - FIXED v0.6683: the PV `1/(x - b)` result is `Log[1 - b] - Log[1 + b]`;
  - FIXED v0.6684: the dsolve-forcing/1 answer is `(1 - E^(1 - t)) UnitStep[-1 + t]`;
  - `g + 1/2 g` in Minimize printing (not re-checked).

**Integration, known open (written around, in the chapter's open-gaps box)**
- A radical in the middle of a tower (`Log[1 + Sqrt[Log[x]]]/x`) gives `Integrate::incomplete`; nested radicals likewise.
- Answers over Q(i) carry Root objects; transcendental constants (Log[2]) demote a correct "no".
- ``Risch`ElementaryIntegralQ`` is transcendental-only (`1/Sqrt[1+x^3]` stays unevaluated).
- `TragerKauers` floods `Thread::tdlen`.
- `Cos[Log[x]]` real form (the Pitfall stands).

**New, found in this pass**
- **FIXED v0.6669.** (ch. 7) `$MessageList` is not reset per input in pipe (NDJSON) mode: after `1/0` on an earlier input, a later `Check[1/0, c]` line and then `safeSqrt[-9]; $MessageList` lists every message of the session (`{Power::infy, Power::infy, ..., nonneg::neg}`); expected only the current input's messages (the v0.561 note says both front ends reset it per line). The chapter reads it inside `Block[{$MessageList = {}}, ...]`.
- **FIXED v0.6672.** (ch. 7) `Dispatch` gives no speedup when the literal left-hand sides share a head: `Dispatch[Table[f[i] -> i^2, {i, 4000}]]` applied with `/.` to 2000 `f[k]` takes 0.99 s against 0.94 s for the plain list; with atomic keys (`i -> i^2`) the table is ~100x faster. Expected the index to key on the whole left-hand side. The chapter shows integer keys and states the limitation in a Performance callout.
- **FIXED v0.6670 (PatternTest infix) and v0.6677 (RHS guard listed on the right).** (ch. 7) OutputForm prints a pattern test in full form: `x_?EvenQ` and `DownValues[k]` show `PatternTest[x_, EvenQ]`; expected `x_?EvenQ`. Also `q[x_, y_] := "g" /; x > y` lists as `HoldPattern[q[x_, y_] /; x > y] :> "g"`, where the condition is expected on the right-hand side. Not shown in the book (the DownValues pairs avoid tests).
- **FIXED v0.6671.** (ch. 7) A fresh session's `$ContextPath` is `{"Lean`", "SeparablePolynomialDecomposition`", "FullSimplify`", "Global`", "Mathilda`", "System`"}`: three kernel packages loaded at start-up stay on the user's path (and in `$Packages`). Shown as is in `packages/2`, and the caption says what they are; possibly intended, flagged for review.
- **FIXED v0.6666.** (algebra) `GroebnerBasis[{x^2 + y^2 - 1.5, x - 2.1 y}, {x, y}]` gives the exact `{-150 + 541 y^2, 10 x - 21 y}`; expected inexact coefficients (inexact in, inexact out, as `Factor`/`Together` do). The book shows only `CoefficientDomain -> InexactNumbers[20]`.
- **NOT A BUG: moved to FullSimplify deliberately in v0.6178 (it cost 26% of Simplify time).** (algebra) The v0.6177 Simplify stage for `SeparablePolynomialDecomposition` no longer fires: `Simplify[Expand[(3 a b - b^2)^3 (x^2 - y + 1)^2 - (a - b)^2 (x - y)]]` gives a 197-leaf collected form; the v0.6177 changelog gives `(3 a - b)^3 b^3 (1 + x^2 - y)^2 + (a - b)^2 (-x + y)` (37 leaves). Likewise `Simplify[Expand[(a + b) (x^2 + y) + (a - b)^2 (x - y)]]` is 44 leaves against the decomposition's 23. Written around as a Pitfall in §Simplification beyond polynomials.
- **FIXED v0.6688 (new builtin, every numeric surface).** (number theory) `BitClear[255, 0]` stays unevaluated: there is no `BitClear` (expected `254`). The chapter clears bits with `BitAnd[n, BitNot[mask]]`.
- **FIXED v0.6673.** (ch. 6) `\!\(x\^2\)` (linear-syntax superscript operator `\^`; also `\_` for a subscript) is a syntax error, "Invalid linear box syntax", in source and in `ToExpression`; expected `x^2` as in the reference. Only the explicit `\*SuperscriptBox[...]` forms are read (v0.6001). Written around as a Pitfall in §Characters, not bytes.
- **NOT A BUG: the expected behaviour also gives "40".** (ch. 9) `st = StringToStream["10 20 30 40"]; Skip[st, Number, 2]; Read[st, Number]; ReadString[st]` gives `"40"`; expected `" 40"` (the read point should sit just after the number read, so the separator is still unread). Shown as is in `filenames/33`; the prose makes no claim about the space.
- **FIXED v0.6687.** (ch. 9) `Export["/tmp/f.txt", "hello"]` stays unevaluated although `Import` reads `.txt` as `"Text"` (v0.551); expected the string written and the file name returned, so text round-trips like an image. The chapter states that `Export` writes no text or data format.
- **FIXED v0.6674.** (special functions) `EllipticTheta[1, 0.64`30, 0.022`30]` gives `0.459643677868141922653044722437 - 6.75632632694932196235668615652*^-49 I`, and `N[EllipticTheta[1, 64/100, 22/1000], 30]` a `*^-68 I` residue; expected a Real (real `u` and real nome, as the machine value `0.459644` is). Only theta_1 shows it (the principal `q^(1/4)` in complex arithmetic?). The chapter checks the theta identities on theta_2..4 and at u = 0 only.
- **FIXED v0.6675.** (special functions) `FunctionExpand[SphericalBesselJ[1, z]]` (likewise `SphericalBesselY`, the spherical Hankel functions, and `HankelH1[1/2, z]`-type half-integer orders) comes back unchanged; expected the elementary form `-Cos[z]/z + Sin[z]/z^2`. Written around as a Pitfall in §Bessel and Airy functions.
- **FIXED v0.6681.** (special functions) `NIntegrate[WaveletPhi[DaubechiesWavelet[2], x], {x, 0, 3}]` gives `0.999874` with no message; the piecewise-linear interpolant of the dyadic cascade integrates to exactly `1` (and the grid sums are exactly 1), so the default rule loses four digits on the kinked integrand silently. Not shown in the book.
- **FIXED v0.6685 (root cause: a lazy table re-read after another module's interrupted read).** (integration) `Integrate[Sqrt[Log[x] + Sqrt[Log[x]]], x]` in a fresh session prints about 30 `SetDelayed::write: Tag IntegrateTable in IntegrateTable[...] is Protected` messages before its `Integrate::nonelem`; so does `Integrate[Sqrt[x + Log[x]], x]` after `Integrate[(1 + 1/x)/(2 Sqrt[x + Log[x]]), x]` in the same session. Expected only the `nonelem`. A table file appears to be re-read onto the protected `IntegrateTable`, so `Check` returns the failure branch. The chapter uses `Method -> "Bronstein"` for those examples.
- **FIXED v0.6686.** (integration) `Integrate[1/(1 + E^x), x]` gives `-2 ArcTanh[1 + 2 E^x]`. Its derivative is right, but the value is complex for every real x (the argument is > 1). Expected `x - Log[1 + E^x]`. Not shown in the book.
- **NOT CHANGED: the time is real special-function work (third-kind elliptic antiderivative in ParallelMixedSpecial; DerivativeDivides/Roach), and later stages still find special-function answers after a proof; speeding up the elliptic construction is a separate task.** (integration) Proofs of non-elementarity are slow in the full cascade, because the held certificate still pays for the special-function stages: `Integrate[1/(x Sqrt[x^3 + a x + 1]), x]` takes 11.3 s against 0.9 s with `Method -> "Trager"`, and `Integrate[Sqrt[x] (x + 1)^(1/3), x]` takes 7.9 s against 0.8 s. The verdicts are the same. Written up as a Performance callout in §Budgets and messages.
- **FIXED v0.6686 (docs).** (integration, docs) `docs/spec/builtins/calculus.md` ~:4374 still places Trager/Bronstein at "15b"/"15c" after LogByParts. They have run at 08b/08c since v0.6521, and "15b" is now SpecialByParts. The open-gaps row G-mid-1 (`Log[1 + Sqrt[Log[x]]]/x`) is still incomplete under `Method -> "Bronstein"`, but Automatic now answers it (verified by differentiation).
- **FIXED v0.6682 (forward transform; InverseLaplaceTransform still absent).** (calculus) `LaplaceTransform[Cosh[2 t], t, s]` and `LaplaceTransform[Sin[b t], t, s]` (no assumptions) stay unevaluated; expected `s/(s^2 - 4)` (valid for s > 2) and `b/(b^2 + s^2)`. The defining-integral fallback assumes only s > 0, so any exponentially growing kernel is declined (`TrigToExp` first works). Also `LaplaceTransform[f'[t], t, s]` is not rewritten as `s LaplaceTransform[f[t], t, s] - f[0]`, and there is no `InverseLaplaceTransform`. Shown as `laplace/12` with a Pitfall.
- **FIXED v0.6676.** (calculus) `Sum[(k - 3)^2/k!, {k, 0, Infinity}]` gives `27/2 + 1/24 HypergeometricPFQ[{2, 2}, {5, 1}, 1]`. The value is right (13.5914...), but the v0.6530 changelog gives the closed form `5 E`; expected `5 E`. The book cites the old wrong answer `27/2` in prose only, not this output.
- **FIXED v0.6678 (Derivative[n] of a built-in symbol).** (numerical calculus) `GaussianQuadratureError[5, Exp, 0, 1]` gives `3.94497*^-13 Derivative[10][Exp]` (after `N`); expected `Derivative[10][Exp]` to evaluate to `Exp`, as applied derivatives of built-ins do since v0.6471. The book uses a symbolic `f`.
- **FIXED v0.6680.** (numerical calculus) `FindRoot[Exp[1000 x] - 2, {x, 1}]` gives `FindRoot::cvmit` and `{x -> 0.9}` (plain Newton moves 0.001 per step); the v0.6228 changelog lists it as solved by the overflow rerun. Expected `x -> 0.000693147`. The book shows the `LaguerreL[1000, x]` case instead.

**Found while removing the workarounds (v0.6688, open)**
- `Integrate[x/((x - a) (x - b)), {x, -1, 1}, PrincipalValue -> True, Assumptions -> -1 < a < b < 1]` gives `Integrate::idiv` and stays unevaluated; expected the finite principal value (numerator `x + 2`, and numeric poles, both work).
- `LaplaceTransform[y'[t] + 2 y[t] == 1, t, s]` does not thread over `Equal`; the book transforms each side.
- `DSolve[{y'[t] + 2 y[t] == 1, y[0] == 0}, y[t], t]` returns the unsimplified `(-1/2 + E^(2 t)/2)/E^(2 t)`; expected `(1 - E^(-2 t))/2`.
