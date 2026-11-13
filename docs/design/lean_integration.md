# Lean 4 integration: design note

User reference: [`docs/spec/builtins/formal-verification.md`](../spec/builtins/formal-verification.md).
Demo: [`examples/lean/`](../../examples/lean/README.md). Lean-side package:
[`lean/`](../../lean/README.md).

## Trust model: untrusted oracle, trusted checker

Mathilda is an **untrusted oracle**. The Lean 4 kernel is the **trusted checker**.

- Mathilda computes an answer (an antiderivative, a factorization, an optimum) and a
  *certificate*: cofactors, a Gram matrix, a Pratt tree, KKT multipliers, a closed form for
  an induction. Mathilda's own algorithms (Risch, Gröbner, SDP, factoring) are never trusted.
  They only have to be right often enough to be useful.
- `Certify` writes a self-contained Lean 4 file that states the result as a theorem and
  proves it from the certificate. The Lean kernel decides whether it is accepted.
- A bug in Mathilda can make `Certify` fail to find a certificate, or produce a certificate
  that Lean rejects (`"Verified" -> False`). It cannot produce a wrong `True`, except
  through the steps listed in `"TrustedSteps"` below, or through a bug in Lean itself.
- The statement is the input **as written**: `Certify` never evaluates a claim before
  translating it (arithmetic aside), so a Mathilda computation inside a claim cannot turn it
  into a tautology (third review; see below).
- In the other direction (`lean/`), the same holds: the tactics ask Mathilda for data and
  re-check every step in Lean. A wrong oracle answer makes the tactic fail
  (`lean/test/Untrusted.lean` drives the tactics with a fake oracle that lies).

### What `"Verified" -> True` means

The file is header (`set_option autoImplicit false`), prelude, theorem chunks, then a
**trusted footer** with, per theorem, `example : <statement> := mathilda_cert_k` and
`#print axioms mathilda_cert_k`. All of the following hold, as checked by
`` Lean`Internal`VerifySource `` and by the per-theorem attribution in `Lean.m`:

1. No message with severity `error` or kind `hasSorry` falls in the theorem's chunk or its
   footer lines.
2. No such message falls outside every chunk and footer (header or prelude).
3. The `#print axioms` report printed on the theorem's **own footer line** names only
   `propext`, `Classical.choice` and `Quot.sound` (a report on any other line is ignored).
   `sorryAx` and the auxiliary axiom that `native_decide` adds both fail this test.
4. Lean's exit code is 0, or every error lies in some theorem's lines (then the sibling
   theorems that own those errors are `False`, and the clean ones are re-run in a file of
   their own, because a file with errors is not compiled and so cannot be replayed).
5. **Kernel replay.** The file was compiled to an `.olean`, and a fresh `lean` process
   running `src/internal/lean/Replay.lean` (none of the checked file's code) replayed every
   declaration of the module into a new kernel environment and found only standard axioms
   under the footer theorem `mathilda_cert_k_footer : S := mathilda_cert_k`. Conditions 1-4
   read the checked run's own output; condition 5 reads only the `.olean`, which that run
   wrote. The replay is therefore defence in depth, not the foundation: the foundation is
   that the file contains **no user-defined code** (condition 0).
0. **Source check.** The generated file passes the check of `src/lean_check.c` (the same
   minimal lexical surface, word rule and grammar as `LeanVerify`, plus the generator
   auxiliaries `def mathilda_cert_<k>_<x>`, the footer lines
   `#print axioms mathilda_cert_<k>_footer` and the byte-identical, comment-free prelude). The only user-influenced
   text in it is the strictly translated statements and deny-listed `Method` scripts.

The heuristic engine never reports `True`: its answer is `Missing["Heuristic", b]`.

`Missing["NotChecked"]` means `Verify -> False` (the source was generated but not run).
`Missing["NotAvailable"]` means no Lean executable was found. A missing checker is not a
refutation, so it is never `False`.

### TrustedSteps: what Lean does not check, per kind

| Kind | Tier | `"TrustedSteps"` |
|---|---|---|
| Identity (polynomial) | 0 | none |
| Identity (with `Assumptions`, or `Domain -> Complexes` / `Reals`) | 0 | none: stated over every char-0 field (or every ordered field); `Domain -> Rationals` / `Integers` over `Rat` / `Int`, as labelled |
| Identity (rational function) | 0 | field algebra from the cross-multiplied identity (unreduced numerators and denominators) to the rational identity where the denominators are nonzero |
| Factorization | 0 | irreducibility of the factors |
| Antiderivative | 0 | the elementary derivative rules and quotient rule (applied symbolically); polynomial derivatives by `D`; the `Together` forms of the integrand, of each `u_i` and of the rational part `R = r/w`; that the cleared denominators are not identically zero (checked by Mathilda) |
| Antiderivative claim `Integrate[f, x] == F` | 0 / 1 | as the row above (Tier 0) or none (Tier 1), plus: the claim is certified as `D[F, x] == f` (`"Theorem"`); `Integrate[f, x] == F` itself holds only up to an additive locally constant term |
| Antiderivative | 1 (Mathlib) | none: `HasDerivAt (fun x => F) f x` on the stated domain; when a hypothesis is not a polynomial in `x` (a `Log`), that the hypotheses are jointly satisfiable, checked numerically by Mathilda at a witness point (otherwise `Certify::vacuous`) |
| Bezout | 0 | none |
| IdealMembership (`Implies[...]`) | 0 | none (every char-0 field) |
| IdealMembership (`PolynomialReduce`) | 0 | that the remainder is fully reduced (the division identity is proved) |
| Roots (univariate, splits over ℚ) | 0 | none: completeness `p = 0 ↔ x = r1 ∨ …` is proved |
| Roots (systems, constraints, other polynomials) | 0 | completeness (only soundness is proved: every constraint, and every denominator nonzero, at every solution) |
| MatrixInverse, Determinant, LinearSolve | 0 | none |
| NullSpace | 0 | completeness (independence and span) |
| Nonnegativity, Optimality | 0 | none: stated over every ordered field, ℝ included |
| Primality (Pratt) | 0 | Lucas's theorem (Lean evaluates the certificate checker, it does not prove that a valid certificate implies primality) |
| Primality (`PrimeQ[2]`) | 0 | none: `Mathilda.IsPrime 2` |
| Primality (composite, `n < 2`) | 0 | none: `¬ Mathilda.IsPrime n` |
| ClosedFormSum | 0 | the closed form is proved for integer `n >= a` only (Lean's statement is over `n = a + m`, `m : Nat`) |

### Tiers

- **Tier 0** uses core Lean 4 only (v4.35 or newer, for `grind`). No Mathlib, no Lake. Every
  generated file inlines `src/internal/lean/Prelude.lean` (about 90 lines: `powMod`,
  `prattCheck`, `PrattCertifies`, `IsPrime`, `matMul`, `matVec`, `HasShape`, `HasLength`, `idMat`,
  `det`, `sqnn`, `mulnn`).
  Proofs are `grind` (polynomial identities, with SOS hints), `decide +kernel` (matrices,
  Pratt) or induction (sums). A warm run takes about 0.5 s, plus about 1.1 s for the kernel
  replay.
- **Tier 1** is opt-in and needs a Lake project with Mathlib built (`$LeanProject`). It is
  used for analytic statements that core Lean cannot express: `HasDerivAt` for
  antiderivatives with `Log`/`ArcTan`. `import Mathlib` costs 10–40 s per run.

## Architecture

```
Certify[comp] ──► Lean.m generator ──► Lean source ──► lean --json ──► verdict ──► ProofObject
   (HoldFirst)      (pure .m:            (LeanForm      (system(),       (NDJSON +
                     D, Together,         strict          temp files)     #print axioms)
                     FactorInteger,       translator)
                     SOS, LinearSolve)
```

| Piece | File | Role |
|---|---|---|
| `LeanForm` | `src/print_lean.c` | display form, and the strict translator `` Lean`Internal`ToLean `` that refuses anything whose Lean text would mean something else (inexact numbers, Int division, Mathlib-only heads in Tier 0, unknown heads) |
| runner | `src/lean_run.c` | `` Lean`Internal`RunProcess ``, `` Lean`Internal`VerifySource ``: temp files, `system()`, NDJSON parsing, axiom extraction, the replay run and its report |
| replay checker | `src/internal/lean/Replay.lean` | `lean --run` script: replays the compiled module in a fresh kernel environment, prints each declaration's axioms |
| package | `src/internal/lean/Lean.m` | settings and discovery, one generator per certificate kind, batching, the session cache, `ProofObject`, `` Lean`Internal`Oracle `` |
| SOS engine | `src/solve/reduce_sos.c` | `` Reduce`SOSCertificate ``: the exact Gram blocks of `Reduce`'s Putinar stage (an additive wrapper; `Reduce` is unchanged) |
| Lean package | `lean/` | `mathilda_linear_combination`, `mathilda_sos`, `mathilda_prime`, `#mathilda`; `Mathilda/Answer.lean` parses the oracle's answer as data |

The generators are `.m` code because they only orchestrate existing builtins and assemble
strings, which follows the repo's "C for performance, rules for mathematics" split. The
package is loaded eagerly and defines rules only: it adds about 2.5 ms to startup.

**Context hygiene.** `` Lean` `` is on `$ContextPath` after load, so any symbol it exports
captures that short name in packages loaded later. The C builtins and private helpers
therefore live in `` Lean`Internal` `` and `` Lean`Private` ``, which are not on the path.
Only the public API is in `` Lean` ``: `Certify`, `ProofObject`, `LeanVerify`,
`$LeanExecutable`, `$LeanProject`, `$LeanVersion` (v0.6523). Option names are strings
(`Options[Certify]` lists `"Engine"`, `"Tier"`, `"Domain"`, `"Verify"`, `"LeanExecutable"`,
`"LeanProject"`, plus the System symbols `Method`, `TimeConstraint`, `Assumptions`), and
`OptionValue` matches a name by its text, so the user's `Engine -> v` (a `` Global` ``
symbol) and `"Engine" -> v` are the same option. Before v0.6523 the option names were
exported from `` Lean` ``: they captured the user's `Domain`, `Verify`, … (which then
survived `ClearAll["Global`*"]` and were missing from `Names["Global`*"]`). The tactic
protocol entry is `` Lean`Internal`Oracle ``, called by its full name.

## Why `system()`, and why it is safe

Strict C99 has no process API other than `system()`, and the repo forbids `popen` and
`fork` in new code. `system()` is enough because the runner only needs an exit code and two
output files. The only POSIX calls (`mkdtemp`, `open`, `fdopen`, `rmdir`, under
`_POSIX_C_SOURCE 200809L` plus `_DARWIN_C_SOURCE` on macOS, guarded by `#ifndef _WIN32`)
create the private scratch directory.

- **No user text reaches the shell.** The command line contains only the executable path,
  flags written by Mathilda, and temp-file paths. The Lean source, the claim and any
  `Method` tactic string go into the `.lean` file.
- **Quoting is complete on POSIX.** Every argument is single-quoted, with each embedded `'`
  written as `'\''`. No character is special inside POSIX single quotes.
- **No planted files.** Every run gets a fresh directory from `mkdtemp` (mode `0700`), and
  the source is created in it with `O_CREAT | O_EXCL | O_NOFOLLOW` (mode `0600`). No other
  user can create, replace or symlink anything there, so the `.out` that is parsed is the
  one the shell's redirection created for this run (if it is missing, the run is a
  failure); `set -C` (noclobber) stays as a second guard. The files and the directory are
  removed afterwards. (The earlier scheme -- predictable names in the shared temp
  directory, probed with `fopen` and then opened -- could be raced: a dangling symlink at
  the source name was followed, and a pre-created `.out` was parsed as Lean's output.)
- **Windows.** `cmd.exe` has no complete quoting, so an argument containing any of
  `" % ! ^ & | < >` is refused (`LeanVerify::badpath`). That branch compiles under MinGW but
  has not been exercised at run time.
- **Time bounds.** Each theorem carries `set_option maxHeartbeats (40000·TimeConstraint) in`,
  which is deterministic. `system()` blocks, and `TimeConstrained` measures this process's
  CPU time, so the timer does not fire during a Lean run and no child process is orphaned.
  A wall-clock kill needs a real spawn API (F7).

## Soundness hardening (v0.6518)

An adversarial review produced four false `"Verified" -> True` results and several weaker
holes. Each is closed at its root, and each reproducer is a regression test
(`test_review_*` in `tests/test_lean_verify.c`, `test_hardening` in
`tests/test_lean_certify.c`).

| Hole | Fix |
|---|---|
| **Method injection.** `Method` was pasted as `"by\n  " <> m`; a second line at column 0 became a top-level command, and a `macro_rules` turned a later `theorem mathilda_cert_2` of a batch into a `def`, whose empty `#print axioms` report passed `x == x + 1`. | Three layers. (1) A Method script is refused (`Certify::method`) if it contains a control character, tab, `«`, `»`, a backquote, `$`, `#`, `@[`, `--`, `/-` or `"`, or uses any command keyword or meta escape as a word (`theorem`, `def`, `example`, `macro`, `syntax`, `elab`, `set_option`, `attribute`, `open`, `namespace`, `end`, `axiom`, `run_tac`, `native_decide`, `sorry`, `Lean`, `IO`, …); every line is indented under `by`. (2) Each item with a user script runs in a Lean file of its own. (3) For every theorem the trusted footer, after all chunks, re-elaborates the statement (`example : S := mathilda_cert_k`) and prints its axioms, and only the report on that footer line counts -- so a chunk cannot pass by declaring something else under the theorem's name, and a report elsewhere (forged or not) is ignored. A deny-list was chosen over a tactic allow-list because lemma names and term arguments are open-ended; the deny-list blocks every way to leave the tactic block or extend the parser, and layers 2-3 bound the damage of anything it misses. |
| **Vacuous Int hypotheses.** `Certify[x == 0, Assumptions -> 2 x == 1]` proved `∀ x : Int, 2 * x = 1 → x = 0`. | The automatic domain is `Int` only for hypothesis-free integer polynomial identities (which then hold in every commutative ring). With `Assumptions` it is `Rat`; ideal membership (`Implies[...]`) is always `Rat` (its cofactor identity holds in every commutative ring anyway). `Domain -> Integers` is still honoured and labelled. Variables of the assumptions are bound. |
| **Solve and `x / 0 = 0`.** Solve's spurious root `x = 1` of `(x-1)/(x^2-1) == (x-1)/(2x-2)` was "verified", because both sides are `0` in Lean. Only `Equal` constraints were checked (`x^2 == 1 && x > 0` checked `1^2 = 1` only). | The soundness statement is every constraint (`==`, `!=`, `<`, `<=`, `>`, `>=`) at every solution, conjoined with `d ≠ 0` for every denominator `d` of the system. Completeness needs one polynomial equation and nothing else. Other heads (`||`, ...) get no certificate. Audit of the other kinds: antiderivative and identity statements are cleared polynomial identities; matrix, optimality and Pratt statements have no division; sums with a denominator in the variable now get no certificate. |
| **Pratt fuel.** `powMod` had 256 bits of fuel, so `PrattCertifies (2^256 + 1) [(2^256 + 1, 3, [(2, 256)])]` held for the composite Fermat number F8. | The fuel is `e + 1`, which bounds the `log2 e + 1` halving steps for every `e`; the kernel unfolds only as many steps as there are bits. F8 is rejected; `2^64 + 1`, `1000003`, `51·2^257 + 1` still verify (well under 1 s), and `NextPrime[2^256]` verifies (Lean 0.8 s; Mathilda's certificate search takes about 117 s, so it is not in the default suite). The checker already required the factor list to multiply out to `p - 1` and every factor to be an earlier step. |
| **LeanVerify.** `axiom bad : False` used from an `example` or a `def`, and `example … := by native_decide`, gave `"Success" -> True`. | Sources that mention `axiom`, `native_decide`, `implemented_by`, `extern`, `unsafe`, `debug`, `#eval`, `elab`, `macro`, `syntax`, `run_cmd`, ... are refused unrun (`LeanVerify::unsafe`). `example`s are renamed to named defs, and every theorem / lemma / def / abbrev / named instance (namespace-qualified) is axiom-checked through an appended, line-anchored `#print axioms` footer. |
| **Temp-file race.** Predictable names, probe-then-open, dangling symlinks followed, a pre-created `.out` parsed as Lean's output. | A private `mkdtemp` directory per run, `O_EXCL | O_NOFOLLOW` source creation, and a check that the shell created the output file. |
| **Shape truncation.** `matVec [[1, 1]] [1] = [1]` holds (`zipWith` truncates). | Matrix statements state `Mathilda.HasShape` / `Mathilda.HasLength`. |
| **Free identifiers.** A symbol not bound by the statement (e.g. `default` in `Solve[x^2 == 1 && default x == 0, x]`) was emitted as a free Lean identifier and resolved to `Inhabited.default`. | Generated files set `autoImplicit false`; theorem statements are translated in closed mode, which refuses unbound symbols, names that contain `«`/`»`/control characters, and two symbols with the same Lean text. |
| **Heuristic `True`.** `Engine -> "Mathilda"` used the same `"Verified" -> True`. | It now gives `Missing["Heuristic", True \| False]`; the summary box already reads `Mathilda Heuristic`. |
| **Cache key.** A 64-bit hash of the theorem text. | The text itself is the key. |

## Soundness hardening, second round (v0.6519)

A second adversarial review found two more false positives at the root of the trust model
and two weaker holes. Regression tests: `test_review2_*` in `tests/test_lean_verify.c`,
the hostile-answer cases in `lean/test/Untrusted.lean`.

| Hole | Fix |
|---|---|
| **LeanVerify deny-list bypass.** `import Lean` plus a user tactic elaborator (`@[tactic Lean.Parser.Tactic.decide]`, or `attribute [tactic …]`, which no word list catches) added `theorem one_eq_zero : (1:Nat) = 0` under `debug.skipKernelTC` with a name built from strings: `"Success" -> True`, "does not depend on any axioms". | **Kernel replay.** Every Lean run of `Certify` and `LeanVerify` compiles the source (`-R <dir> -o <dir>/MathildaCert.olean`, in the private temp directory) and then runs `lean --run src/internal/lean/Replay.lean <olean>` -- through `lake env` for a project. The script reads the module, imports its imports, replays every declaration with `Lean.Kernel.Environment.replay` (what the toolchain's `leanchecker` does per module) and computes each declaration's axioms from the replayed environment (imported constants through Lean's precomputed axiom table). It runs none of the checked file's code, so neither its verdict nor its axiom report can be forged by that file; a `skipKernelTC` declaration fails with `(kernel) declaration type mismatch`. `"Success"` needs the replay, standard axioms for every declaration of the module and for every declared name. The deny-list stays as an early refusal (extended with `addDecl`, `skipKernelTC`, `command_elab`, `@[tactic`, …). `leanchecker` itself was measured (0.6-0.7 s) but reports no axioms, which LeanVerify needs from a trusted process; the script costs about 1.1 s on core Lean (0.5 s of it elaborating the script, 0.45 s importing `Init`), below the 3 s threshold, so `Certify` replays every run too (defence in depth: its sources contain no user code except deny-listed `Method` scripts). |
| **Oracle answer elaborated as Lean.** The tactics parsed Mathilda's answer as an arbitrary `term` and elaborated it; a malicious oracle answered `by_elab do …`, which added a kernel-unchecked theorem, registered it as a `grind` lemma (so `mathilda_sos` closed `y ≤ x` and `(1:Rat) = 0` followed) and ran IO. | The answer is **data**. `lean/Mathilda/Answer.lean` is a hand-written tokenizer and Pratt parser (a pure function on `String`) for numerals, `v<i>`, `+ - * / ^`, unary `-`, lists and tuples, with size bounds; the tactics build the certificate's terms from the AST with fixed quotations, numerals and the goal's own fvars. `Lean`Internal`Oracle` now answers `[1, -v1 - v2]` instead of `fun (v1 v2 : _) => (… : List _)`. `LeanVerify` also refuses `import Mathilda` and `set_option mathilda.*`. |
| **Vacuous hypotheses over ℚ.** `Certify[x == 7, Assumptions -> x^2 + 4 == 0, Method -> …]` proved `∀ x : Rat, 4 + x^2 = 0 → x = 7`, false over Mathematica's complex semantics (`x = 2 I`); likewise `x^2 == 2`. | **Generic domains.** Core Lean's `grind` works over its own algebraic classes, so statements with hypotheses quantify over every field of characteristic 0 (`[Lean.Grind.Field K] [Lean.Grind.IsCharP K 0]`: ℚ, ℝ, ℂ) and, with an order relation, over every ordered field (`[LE K] [LT K] [Std.LawfulOrderLT K] [Std.IsLinearOrder K] [Lean.Grind.OrderedRing K]`: ℝ). The statement is translated over `Rat` and generalised (`Rat` → `K`); the theorem takes the binders as parameters, the footer states the full `∀ {K : Type} …`. Ideal membership, identities with assumptions, SOS and optimality use it; every showcase row still verifies, and the "Rat → Real transfer not checked" trusted step is gone. `Solve` soundness statements are closed rational facts and stay over `Rat`. |
| **Together cancels; weak primality statements.** `Certify[(x^2 - 1)/(x - 1) == x + 1]` proved `∀ x : Int, 1 + x = 1 + x`; `Certify[Not[PrimeQ[0]]]` proved `0 < 2`. | Rational identities are cross-multiplied with the unreduced numerators and denominators (`rawFraction`), so Lean checks `x ^ 2 - 1 = (x - 1) * (1 + x)`. The prelude defines `Mathilda.IsPrime n := 2 ≤ n ∧ ∀ d, d ∣ n → d = 1 ∨ d = n`; non-primes are `¬ Mathilda.IsPrime n` (from a factor, or from `2 ≤ n`), and `PrimeQ[2]` is `Mathilda.IsPrime 2`, all with no trusted step. |

(Correction, v0.6520: the claim above that the replay's verdict "cannot be forged by that
file" was wrong. The replay runs none of the file's code, but its input -- the `.olean` --
is written by the process that elaborated the file; a file that runs code at elaboration
time can write a forged module, or race to overwrite it. The third round closes this by
never running such a file.)

## Soundness hardening, third round (v0.6520)

A third adversarial review found two more false `True` verdicts and three weaker holes.
Regression tests: `test_review3_*` in `tests/test_lean_verify.c` and
`tests/test_lean_certify.c`.

| Hole | Fix |
|---|---|
| **Forged replay input.** `import Lean` plus `attribute [tactic Lean.Parser.Tactic.decide] pwn`, where `pwn` added `theorem bogus : (1 : Nat) = 0` with `addDeclCore (doCheck := false)` and wrote forged `MathildaCert.olean.server` / `.olean.private` (which `Replay.lean` read in place of the module), or spawned a shell that overwrote the `.olean` between Lean's exit and the replay: `"Success" -> True`. The root cause: the elaborating process is attacker-controlled and produces the replay's input. | **No user code runs.** LeanVerify's guarantee holds only when the checked source cannot run user-controlled code during elaboration, and that is now enforced structurally, before Lean starts, by an allowlist grammar (`` Lean`Internal`CheckSource ``, `src/lean_check.c`): (1) no `import` (core Lean only, "Tier 0"), or with a project only `import Mathlib` / `import Mathlib.<Path>` (`LeanVerify::import`; `Lean`, `Lake`, `Mathilda`, `Qq`, `Aesop`, `Batteries` are refused); (2) every command starts at column 0 (continuation lines indented; tabs and carriage returns refused) with `theorem`, `lemma` (Mathlib only), `example`, `open` (namespace names, optional `scoped` / trailing `in`), `namespace`, `section`, `noncomputable section`, `end`, `variable`, `universe`, `set_option maxHeartbeats / maxRecDepth / synthInstance.maxHeartbeats <n>` or `autoImplicit / relaxedAutoImplicit / linter.* true / false` (optional `in`), or `#print axioms <name>`; everything else (`def`, `abbrev`, `instance`, `structure`, `inductive`, `class`, `axiom`, `opaque`, `macro`, `syntax`, `notation`, `elab`, `attribute`, `@[…]`, `initialize`, `#eval` / `#check` / …, `run_cmd`, `simproc`, `mutual`, `module`, …) is `LeanVerify::cmd`; (3) inside a command, a tiny Lean lexer (comments `--` and nested `/- -/`, string, raw-string and character literals, `«»` identifiers) refuses the escapes `by_elab`, `run_tac`, `run_elab`, `run_meta`, `run_cmd`, `native_decide`, `native` (`decide +native`), `implemented_by`, `extern`, `unsafe` and `unsafe…IO`, `IO` / `EIO` / `BaseIO`, `Elab`, `MetaM` …, `sorry`, `admit`, `ofReduceBool`, `trustCompiler`, `debug`, `plausible`, `slim_check`, `polyrith`, every command keyword (Lean starts a command wherever its keyword appears, column or not), `#word`, `word%` (`eval%`), `@[`, backquotes (name literals, quotations) and interpolated strings (`LeanVerify::unsafe`); an unterminated comment / string / `«»` identifier is `LeanVerify::lex`. Without user definitions and without meta-level imports there is no user-defined code to run: built-in elaborators run only trusted code. The trusted prelude may occur once, byte-identical to `Prelude.lean`. (4) Defence in depth: `Replay.lean` reads only the main `.olean` and refuses a module-system file or any `.olean.server` / `.olean.private`; `lean_run.c` lists the private directory after elaboration and refuses to replay if it holds anything besides `MathildaCert.lean`, `.out`, `.err` and a regular-file `.olean`; and every theorem / example the source text declares (names from the checker; examples renamed as before) must be among the replayed declarations with standard axioms, so a module that drops or renames one fails. (5) `Certify`'s generated files pass the same check as an internal assertion (`Certify::unsafe`, `"Verified" -> False`): the header, the byte-identical prelude, the Tier 1 `import Mathlib`, and generator auxiliaries `def mathilda_cert_<k>_<x>` (sum preambles) are accepted; the user-influenced text (Method scripts, already deny-listed, and the strictly translated statements) must pass the grammar. |
| **Ordered field for complex variables.** `Certify[x == 7, Assumptions -> x^2 + 4 == 0 && y > 0, Method -> …]` was `True`: one order relation made *every* variable range over an ordered field, where `x^2 + 4 = 0` has no solution. | The ordered-field statement is used only when every bound variable is real: it occurs in an order relation of the assumptions, or is declared with `Element[x, Reals]` / `Element[{x, y}, Reals]` (new). Order relations with a variable that is not real give `Certify::realvars` and no certificate; without order relations the statement is over every char-0 field, as before; with every variable declared real (and no order relation) it is over every ordered field. SOS (`ForAll[vars, p >= q]`) and `Minimize` / `Maximize` keep the ordered-field form: inequalities and optimisation are over the reals by Mathematica's own semantics, so every variable there is real. |
| **Evaluated claims.** `generate` evaluated both sides of an `Equal` before translating, so `Simplify[(x^2-1)/(x-1)] == x + 1`, `Refine[Sqrt[x^2], x < 0] == -x`, `Integrate[x^2, x] == x^3/3`, `Integrate[Sin[x], {x, 0, Pi}] == 2`, `D[x^3, x] == 3 x^2`, `Factor[x^2 - 1] == (x - 1)(x + 1)`, `Exp[Log[x]] == x`, `Mod[7, 3] == 1` were `True` with tautological statements and no trusted step. | Every input is walked **held** (`rewriteInput`) before a generator sees it: numbers, symbols, `Plus` / `Times` / `Subtract` / `Divide` / `Minus`, integer powers, relations, connectives and `List` are kept; `Expand`, `Together`, `Cancel`, `Apart`, `Factor`, `FactorTerms`, `Collect` are read through (value-preserving on rational functions, so Lean checks the user's own expression; `Factor` adds the irreducibility trusted step and `"Kind" -> "Factorization"`); `D[F, x]` with `F` a polynomial is replaced by Mathilda's derivative with the trusted step "Mathilda's differentiation: D[F, x] = …"; an indefinite `Integrate[f, x]` as a whole side goes to the antiderivative generator with the **claimed** `F` (an antiderivative is certified up to a constant); a `Sum[t, {k, a, n}]` side goes to the closed-form generator with the claimed closed form; anything else inside an input gives `Certify::eval`. Identity numerators and denominators are read off the held sides (`rawH`), so `x/x == 1` is `x = x` with `x ≠ 0`, not `1 = 1`. `"Input"` is the input as written; `"Result"` and the heuristic engine use Mathilda's own evaluation. |
| **Sum bound; silent options.** `Sum[1, {k, 1, Infinity}] == Infinity` matched `n_Symbol`; `Assumptions` were silently ignored by every kind but identities. | The upper bound must be a symbol outside ``System` `` that is neither `Protected` nor `Constant`. A kind that ignores `Domain`, `Assumptions` or `Method` issues `Certify::optignored` (ignoring a hypothesis proves a stronger theorem, so it was sound, but not silent any more). |
| **Scratch leftovers.** `scratch_close` removed a fixed suffix list and missed `.ir.sig` (and anything else). | It removes the known names (now including `.ir.sig`, `.c`, `.hash`, `.trace`, …), then every remaining entry of the private directory via `readdir` (POSIX): entries are `unlink`ed (a symlink is removed, never followed), subdirectories only `rmdir`ed, never recursed into; nothing outside the directory is touched. |

## Soundness hardening, fourth round (v0.6521)

A fourth adversarial review found one more false `"Success" -> True` and two weaker holes.
Regression tests: `test_review4_*` in `tests/test_lean_certify.c`, `tests/test_lean_verify.c`
and `tests/test_lean_mathlib.c` (the reproducers' texts, with their payload paths
neutralised, are in `tests/lean_check_corpus.h`), plus the ctests `lean_check_tokens_sync`
and `lean_check_tokens_mathlib`.

| Hole | Fix |
|---|---|
| **The checker's lexer disagreed with Lean's** (three ways). (a) A `'` after a digit or a non-ASCII byte was taken for an identifier prime, where Lean starts a character literal (`↦'"'`, `1'"'`); (b) identifiers with `!` / `?` (`x!'"' "`) were split differently; (c) `/--/` closed a comment for the checker but opens a doc comment in Lean. Each let a string or comment hide, for the checker, a column-0 `#eval` that Lean ran: reproducer x2 (an `#eval` that ran a shell script forging the `.olean` and the output, then exited 0) gave `"Success" -> True` for `theorem bad : False`. | **A minimal lexical surface instead of a second lexer.** The root cause was the design: "the checker's lexer agrees with Lean's byte for byte" is a property no test can establish. It is replaced by a surface on which there is nothing to agree about. **(A)** Refused anywhere, `LeanVerify::lex` with the construct and its line: `"` (every kind of string), `/-`, `-/`, `--` (every comment), `«`, `»`, the backquote, `#` (no hash commands; `#print axioms` is dropped, the replay computes the axioms), `@[`, `$`, `%`, `\`, control characters, a `'` that is not an identifier prime (it must follow `[A-Za-z0-9_']` in a run starting with a letter or `_`: no character literal can occur), invalid or overlong UTF-8, and any non-ASCII character outside a fixed list of mathematical symbols and the Greek letters. With no string, comment or character literal there is no span whose content Lean skips, so every byte belongs to a token the checker sees. **(B)** Words are maximal runs of `[A-Za-z0-9_]`, every other byte separates them; Lean's identifiers are coarser, so every keyword Lean sees starts a word here, except after a numeric literal (`2theorem`; every suffix of a digit-led run is checked) and after a Mathlib notation token ending in a letter (`→otheorem` is `→o` then `theorem`; the list `LEAN_GLUED_TOKENS_MATHLIB` is checked) -- found by the token-table test below. Every command-leading word and every escape word is refused anywhere except an allowed starter in its starting position (at worst over-refusal, never a miss). **(C)** The grammar of v0.6520, minus `#print axioms`, minus `open ... in` and `open scoped` (namespace names only; this also avoids activating scoped notation and macros for a single command), and with each starter followed by a space or a newline (`variable?` is another command). **Derived, tested keyword list.** The command-leading words (`src/lean_check_tokens.h`) are not written by hand: `tools/lean_command_tokens.lean` reads the `command` category from `Lean.Parser.parserExtension` and prints the first tokens of every leading parser (and the token table); `tests/scripts/lean_check_tokens.sh` runs it and checks, against `` Lean`Internal`CheckSource `` itself, that every leading token is refused inside a command and (unless it is an allowed starter) at column 0, that no leading command parser is unindexed, and that a keyword glued to any token that ends in a letter is refused. With `import Lean` (core): 159 leading tokens, 120 words in the core list (the others carry `#`, `%`, `$`, `@[`, `/-`, which the surface refuses), 504 tokens, all covered; with `import Mathlib`: 248 leading tokens, 34 more words in the Mathlib supplement, 1148 tokens, 25 glued tokens inside the surface, all covered. The one exempt token is `?` (its parsers need `binder_predicate` or an antiquotation `$`, both refused). "Command list completeness" is therefore a tested property, not a trust item. **Prelude.** `Prelude.lean` (and its mirror `lean/Mathilda/Prelude.lean`) no longer has comments (its documentation moved to `Lean.m`); the checker refuses a prelude that contains a string, comment or character literal, so Lean's lexer is at the top level where the skipped prelude ends. Generated files lost their header comment. **Method** scripts satisfy (A) and (B) too: `methodProblem` runs the checker on the script as a theorem body, on top of its own list. **Greek variables** are printed as plain identifiers (`α`, not `«α»`), which the surface accepts. |
| **`Evaluate` through the held walker.** `cw`, `cwl`, `cwClaim`, `hk`, `rawH`, `rawHl` were `HoldAll`, which still evaluates an `Evaluate[...]` argument: `ProofObject[Evaluate[Simplify[(x^2 - 1)/(x - 1)]] == x + 1]` was `True` with the theorem `1 + x == 1 + x` and no trusted step. | The claim-walking helpers are `HoldAllComplete`, and `Evaluate` itself is refused (`Certify::eval`). `Certify` and `ProofObject` are `HoldAll` (were `HoldFirst`), so a top-level `Certify[Evaluate[c]]` is still the explicit way to certify a computed claim. |
| **Values substituted silently.** A symbol with an `OwnValue` was replaced by its value (`x = 5; ProofObject[x^2 == 25]` proved `(25 : Int) = 25`; `e = Simplify[...]; ProofObject[e == x + 1]` a tautology), while a `DownValue` gave `Certify::eval`. | A symbol with a value (an `OwnValue`, not `Protected`) anywhere in the claim, including a symbol holding the whole claim, is refused with `Certify::value` ("x has a value (5); Clear[x] or certify the value-substituted claim explicitly"). The `Assumptions` option is checked held for values and `Evaluate` before the options are evaluated (`a = 0; Certify[y == 7, Assumptions -> y a == 1]` would have made the hypothesis `False`). |

Kernel replay and the private-directory checks of v0.6519 / v0.6520 stay, as defence in
depth.

**Fifth review (v0.6521, follow-up fix).** In Lean 4 a plain `open X` -- not only
`open scoped X` -- and likewise `namespace X` activate X's scoped extensions, scoped
*command* syntax included (for example `Lean.PostprocessTraces` declares
`scoped syntax "postprocess_traces " term " in" command`, whose elaborator evaluates the
term). The keyword list is derived with nothing open, so such keywords were outside the
tested property. LeanVerify sources may therefore no longer use `open` or `namespace` at
all; proofs use fully qualified names. The generators never emitted either (the trusted
prelude's `namespace Mathilda` is part of the byte-identical prefix), so `Certify` is
unaffected.

## Soundness hardening, sixth review (v0.6522)

A sixth adversarial review of the `Certify` path found two false `"Verified" -> True`
results against the contract (one through `Domain`, one through `Integrate[f, x] == F`) and
five weaker holes. Regression tests: `test_review6_*` in `tests/test_lean_certify.c` (no Lean)
and `test_review6_verify` in `tests/test_lean_verify.c`.

| Hole | Fix |
|---|---|
| **`Domain` values fell back silently.** `domainType[_, ...]` mapped every unlisted value to the automatic rules, but the generic char-0 / ordered-field form was chosen only for `Domain -> Automatic`, and `Element[x, Reals]` was dropped on that path: `Certify[x == 7, Assumptions -> x^2 == -4, Domain -> Complexes, Method -> …]` proved `∀ x : Rat, x ^ 2 = -4 → x = 7` (false over ℂ). | **One meaning per value** (`$domains` in `Lean.m`): `Automatic` -- Mathematica's default (hypothesis-free: `Int` / `Rat`, equivalent for polynomial identities; with hypotheses: as `Complexes`); `Complexes` -- every char-0 field, hypotheses or not (every ordered field when every variable is declared real; `Certify::realvars` otherwise); `Reals` -- every ordered field (Tier 0; it used to be Mathlib's `Real`); `Rationals` -- `Rat`; `Integers` -- `Int`. Anything else (`Algebraics`, the string `"Complexes"`, …) is `Certify::domain` before any generator runs. `Element` assumptions follow the same rules on every path (`Element[x, Reals]` declares realness, implied memberships are dropped, any other `Element` is `Certify::element`). Variables printed `K` or `Rat` are refused in the generic form (the type variable is `K`). The identity generator is the only one that reads `Domain`; the others issue `Certify::optignored`. |
| **`Integrate[f, x] == F` certified up to a constant.** `Certify[Integrate[2 x, x] == x^2 + 5]` and `Certify[Integrate[-1/(1 + x^2), x] == ArcTan[1/x]]` were `True`; in Mathematica both claims are `False` (`Integrate` returns one antiderivative; `ArcTan[1/x]` differs from `-ArcTan[x]` by `±Pi/2`, a different constant on each side of 0). | The claim is not certifiable as written: Lean cannot know which antiderivative `Integrate` returns without evaluating it. It is now certified **as the honest claim** `D[F, x] == f`: `"Theorem"` is `ForAll[x, D[F, x] == f]` (Tier 1: the `HasDerivAt` statement), `"Input"` stays as written, `Certify::antideriv` is issued, and `"TrustedSteps"` gains "the claim Integrate[f, x] == F is certified as D[F, x] == f …; Integrate[f, x] == F itself holds only up to an additive term with derivative 0, which is locally constant (constant on each interval where F and f are defined …)". |
| **Values beyond `OwnValues`; `Protected` exempt.** `valuedQ` exempted every `Protected` symbol and looked at `OwnValues` only: `e = Simplify[...]; Protect[e]; Certify[e == x + 1]` was a tautology, `t /: t^2 = -1; Certify[t^2 == -1]` proved `(-1 : Int) = -1`, `v = 5; Protect[v]; Certify[Integrate[2 w, w] == w^2 + v]` substituted `v`; the variable slots of `Integrate` / `Sum` / `D` were not checked at all. | `symbolProblem`: a non-``System` `` symbol in a claim or its `Assumptions` must have no `OwnValues`, `DownValues`, `UpValues`, `SubValues`, `NValues` or `DefaultValues` and no attribute beyond `Protected` / `ReadProtected` / `Locked` / `Temporary`; ``System` `` symbols (the heads and constants the walker and translator understand) are accepted unless a user gave one an `OwnValue` or `UpValue` (`I :> I` excepted). The exemption is by context, not by `Protected`. The variable slots are checked like every other symbol. `Certify::value` names the kind ("t has UpValues (…)"). |
| **Private surd capture.** `Certify[Integrate[2, x] == Sqrt[2] Lean`Private`s1 x]` was `True`: the claim named the generator's surd variable, so `s1^2 - 2 = (s1^2 - 2) * 1`. | Symbols of `` Lean`Private` `` / `` Lean`Internal` `` are refused in claims and `Assumptions` (`Certify::context`), and the surd variables are fresh per call (`` Lean`Private`s<j>$<n> ``, printed `s<j>`), so capture is impossible by name or by coincidence. |
| **Minted ProofObjects.** `ProofObject[<|"Verified" -> True, …|>]["Verified"]` was `True`. | A session registry: `mkProofObject` stores the `Hash` of every association it produces; `obj["Verified"]` and the summary box report the verdict only for a registered object, `Missing["Unauthenticated"]` otherwise. Documented limits: code that deliberately writes `` Lean`Private` `` internals can forge the registry, and objects from another session are unauthenticated (re-run `Certify`, or `LeanVerify[obj["LeanSource"]]`). |
| **Unnamed trusted steps.** `Together[rpart]` for the rational part `R = r/w` of an antiderivative; a closed-form `Sum` is proved for integers `n >= a` only. | Both are now in `"TrustedSteps"` (table above). |
| **Computed `Assumptions`.** `Certify[x == 1, Assumptions -> Reduce[x^2 == 1 && x > 0, x]]` used the hypothesis `x == 1` computed by `Reduce`, with `"TrustedSteps" -> {}` (likewise `Simplify[...]`, a user function `x == f[2]`). | `heldOptionProblem` walks the `Assumptions` held with the claim walker (arithmetic, relations, `And` / `Or` / `Not` / `Implies`, `Element`, lists, the value-preserving heads); anything else, `D` included, is `Certify::eval`. A `SetOptions[Certify, Assumptions -> …]` default was evaluated by `SetOptions` before `Certify` could see it, so a default other than `{}` is refused (`Certify::optdefault`) unless the call gives its own `Assumptions`. `Domain` and `Method` defaults need no such rule: they are checked by value (`$domains`, `methodProblem`). |

## Final polish (v0.6523)

| Hole | Fix |
|---|---|
| **Vacuous Tier 1 hypotheses.** `Certify[Integrate[x, x] == x^2/2 + 1/(Log[x^2] - 2 Log[x]), Tier -> 1]` was `"Verified" -> True`: the statement `∀ x : ℝ, … → Real.log (x ^ 2) - 2 * Real.log x ≠ 0 → HasDerivAt …` is vacuous, since `Real.log x = log |x|` makes that hypothesis false at every real `x` (over ℂ it is nonzero for `x < 0`, so Mathilda's own algebra did not notice). The `x^3 + …` variant was rejected only because the proof script did not find the contradiction. | The Tier 1 builder needs a witness: a rational point (from a fixed list) at which every hypothesis is nonzero under Lean's semantics (`Log[u]` evaluated as `Log[Abs[u]]`, `Real.log 0 = 0`), decided at 40 digits with `|d| > 10^-20`. Without one the input is `Certify::vacuous` and `$Failed`; with one the point is `"Certificate"["HypothesisWitness"]`, and the numeric check is a trusted step when some hypothesis is not a polynomial in `x`. Regression test `test_tier1_vacuous_hypotheses` (`tests/test_lean_mathlib.c`). |
| **Option names captured user symbols.** `Engine`, `Tier`, `Domain`, `Verify`, `KeepFiles`, `LeanExecutable`, `LeanProject` and `Oracle` were created in `` Lean` `` on `$ContextPath`, so `Domain = 7; Verify[x_] := x + 1; ClearAll["Global`*"]` cleared nothing, and `Engine = "zzz"` broke `Certify` with "unknown Engine". | String option names matched by text (see Context hygiene); `Lean`Internal`Oracle`. A held option whose name symbol has an `OwnValue` is `Certify::optvalue` (the value would replace the name); an unknown name is `Certify::optx` / `LeanVerify::optx`. Test `test_option_names`. |
| **Inaccurate refusals.** `Certify[Sum[k^2, {k, n}]]` said "applies Sum"; `Integrate[1/Sqrt[x], x]` said "applies Power"; `Solve[x^3 + x == 0, x, Reals]` named its value `{{x -> 0}}`. | `Sum[t, {k, n}]` is read as `Sum[t, {k, 1, n}]` (`Sum`'s definition); `Certify::fracpow`; a supported computation head in an unsupported form is `Certify::nocert` for the input as written. Test `test_refusals`. |

### The guarantee, precisely

`LeanVerify[src]["Success"] -> True` means: `src` passed the source check (minimal lexical
surface, word rule, allowlist grammar); Lean elaborated it (with the trusted prelude, if
present) with no error and no `sorry`; the compiled module, read by a fresh `lean` process,
replayed declaration by declaration into a new kernel environment; every declaration of the
module and every theorem / example the source declares depends only on `propext`,
`Classical.choice` and `Quot.sound`; and the elaboration's own `#print axioms` reports
(appended by LeanVerify after the check) agree. `Certify[...]["Verified"] -> True` means the
same for the generated file and its footer theorem `mathilda_cert_k_footer`, whose statement
is the faithful translation of the input as written -- no computation in it evaluated, no
value substituted, over the `Domain` it shows -- after Mathilda's arithmetic normalisation of
the claim: the held `Plus`, `Times`, `Power` and `Rational` nodes (and those of a claimed
antiderivative `F`) are put in Mathilda's canonical form (flattened, ordered, numeric
coefficients collected, rationals reduced) before translation, and Lean checks that
canonical form. Every other unchecked step is in `"TrustedSteps"`.
The one deliberate exception is a claim `Integrate[f, x] == F`, which is certified as
`D[F, x] == f`: `"Theorem"` shows that statement, `Certify::antideriv` says so, and
`"TrustedSteps"` states that the equation itself holds only up to a locally constant term.
`"Verified"` is reported only for a `ProofObject` that `Certify` produced in this session.

**Trust base.** The Lean kernel; the Lean toolchain (elaborator, compiler, `.olean` writer
and reader) as installed; the imported compiled modules (core `Init`, and with a project
Mathlib and its dependencies) as built; and the source check, which now rests on two things
that are either trivially inspectable or tested: (1) the **minimal-surface rules** -- with
no string, comment, character literal, `#`, `$`, `%`, `\`, backquote, `@[` or French quote,
and only an allowlist of Unicode, no part of the source is skipped by Lean's lexer, so the
checker does not need to reproduce that lexer; (2) the **command-keyword list**, derived
from Lean's parser tables and tested against them for core Lean and, when a project is
configured, for Mathlib (ctest `lean_check_tokens_sync` / `lean_check_tokens_mathlib`;
re-run after a toolchain or Mathlib update). The remaining assumption, stated plainly: the
**built-in (non-user) elaborators and tactics of core Lean and, with a project, of Mathlib
do not execute user-supplied code through ordinary term and tactic syntax** (beyond the
refused escape words such as `native_decide`, `by_elab`, `run_tac`, `plausible`, `eval`).
That list of escape words is the one hand-maintained deny-list left; the replay additionally
guards against an elaboration that skips the kernel. Mathilda's own algorithms are not
trusted, except for two things: (3) the **arithmetic normalisation of the claim** described
above (Plus / Times / Power canonical form), and (4) the per-kind `"TrustedSteps"` (for
example Lucas's theorem behind a Pratt certificate). Mathilda cannot make Lean accept a false
theorem; it can make the theorem differ from the claim only through the listed
`"TrustedSteps"` and the arithmetic normalisation of the claim, and `"Theorem"` /
`"LeanTheorem"` show exactly what was proved.

## Configuration

| Setting | Meaning |
|---|---|
| `$LeanExecutable` / option `"LeanExecutable"` | `Automatic` or a path. `Automatic` tries the environment variable `MATHILDA_LEAN`, then the newest `~/.elan/toolchains/*/bin/lean` (called directly, so no `lean-toolchain` file is needed), then `lean` on `PATH`. |
| `$LeanProject` / option `"LeanProject"` | `None` or a Lake project with Mathlib, for Tier 1. Defaults to the environment variable `MATHILDA_LEAN_PROJECT`. Tier 1 runs `cd <project> && <lake> env <lean> --json FILE`, with the `lake` beside the Lean binary. |
| `MATHILDA_BIN` | Lean side: the Mathilda binary that the `lean/` tactics call (else `set_option mathilda.path`, else `./Mathilda` or `../Mathilda`, else `PATH`). |
| `MATHILDA_LAKE` | test scripts: the `lake` to use for `lean/`. |

## Limitations

- Tier 0 statements are over `Int`/`Rat`, or over a generic char-0 / ordered field (statements
  with hypotheses, `Domain -> Complexes` / `Reals`, SOS, optimality). Real-analytic content (derivatives of `Log`,
  `ArcTan`) is either reduced to a cleared polynomial identity (Tier 0, with the trusted
  steps above) or proved in Mathlib (Tier 1).
- Antiderivatives must be built from rational functions, quadratic surds, and one level of
  `Log`/`ArcTan`/`ArcTanh`. `RootSum`, `Root` and complex results have no certificate (F9).
- `Det` is limited to n ≤ 6 (Laplace expansion in the kernel; 6×6 takes about 1.5 s).
- `Domain` (`Automatic`, `Complexes`, `Reals`, `Rationals`, `Integers`; anything else is
  `Certify::domain`), `Method` and `Assumptions` apply to identity claims. Other kinds take their
  hypotheses from the claim itself (`Implies[...]`, `ForAll[vars, cond, ...]`) and pick their
  own proof script; they issue `Certify::optignored` when these options are given.
- A claim may contain arithmetic, `Expand` / `Together` / `Cancel` / `Apart` / `Factor` /
  `FactorTerms` / `Collect`, `D` of a polynomial and a whole-side `Integrate[f, x]` or `Sum`;
  any other computation inside it is refused (`Certify::eval`): certify the computation
  itself instead.
- Factor irreducibility, Risch non-elementarity, `NullSpace` completeness and Lucas's
  theorem are not proved (listed as trusted steps).
- The Lean-side oracle starts one Mathilda process per tactic call (about 0.05 s).
- `LeanVerify` renames `example`s to named defs so their axioms are checked. It runs only
  sources without user-defined code inside a minimal surface: no `def`, `instance`,
  `structure`, `macro`, `attribute`, …, no imports beyond Mathlib, and no strings,
  comments, character literals, `#` commands, `%` (write `n % 2` as `Nat.mod n 2` or
  `HMod.hMod n 2`), `$` or Unicode outside the allowlist. Command and escape words are
  refused as identifiers too (a variable named `end`, `local` or `alias`; `Polynomial.eval`
  with Mathlib). A development with its own definitions has to be checked by Lean / Lake
  directly.
- `Certify` refuses a claim that mentions a symbol with a value of any kind or an
  evaluation-changing attribute (`Certify::value`), including `c = a == b; Certify[c]`: use
  `Certify[Evaluate[c]]` to certify the value. `Assumptions` set with `SetOptions` are
  refused (`Certify::optdefault`): pass them in the call.

## Follow-ups

| ID | Item | Why deferred |
|---|---|---|
| F1 | `LinearProgramming` with an exact rational simplex | No simplex exists today. LP optimality is already certified for every LP that `Minimize`/`Maximize` solves. |
| F2 | Gosper certificates for hypergeometric sums | The induction step needs factorial/binomial facts from Mathlib. Polynomial·cᵏ terms do not. |
| F3 | Gröbner cofactor tracking (`<G> ⊆ <F>`) | A change to the hot Gröbner kernels. |
| F4 | Irreducibility certificates; Risch non-elementarity | Research items. |
| F5 | Tier 1 primality through Mathlib's `lucas_primality` | Removes the trusted Lucas step; needs `ZMod` proof assembly. |
| F6 | `DSolve` substitution certificates (Tier 1) | Needs a second-order `HasDerivAt` builder. |
| F7 | A spawn API with a wall-clock kill, and a public `RunProcess` | Portability work on three platforms. Heartbeats already bound each run deterministically. |
| F8 | A persistent NDJSON-pipe oracle for the Lean tactics | Only pays off for proofs with many tactic calls. |
| F9 | `Root`/`RootSum`/complex antiderivatives and radical `Solve` roots in Tier 0 | Needs an algebraic-number encoding beyond quadratic surds. |
| F10 | `NullSpace` completeness (rank witness) and eigenpairs | Needs a minor-based rank certificate. |
| F11 | A core-Lean field lemma from the raw rational statement to the cleared identity | Removes the remaining "field algebra" trusted step of rational identities (the cross-multiplied identity now uses the unreduced denominators); core Lean lacks `div_add_div`. |
