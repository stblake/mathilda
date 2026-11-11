# Formal verification (Lean 4)

Mathilda can hand its answers to the [Lean 4](https://lean-lang.org) theorem prover for checking. `Certify` computes a result, writes a self-contained Lean 4 file that states the result as a theorem with a proof script, and runs the Lean kernel on it. The answer is a `ProofObject`.

**Trust model.** Mathilda is an untrusted oracle and the Lean kernel is the trusted checker. `"Verified" -> True` means that Lean accepted the theorem with no `sorry` and with no axiom outside `{propext, Classical.choice, Quot.sound}`. Each `ProofObject` also lists its `"TrustedSteps"`: the facts that Lean did **not** check for that kind of certificate, such as the `Together` normal form used to clear denominators, or the irreducibility of the factors of a `Factor` result. When the theorem is the whole claim, the list is empty.

Lean is optional. With no Lean installed, `Certify[..., Verify -> False]` still generates the Lean source and the statement, and `"Verified"` is `Missing["NotChecked"]`. Tier 0 certificates need only core Lean 4 (v4.35 or newer, for `grind`) and no Mathlib. The package lives in `src/internal/lean/` (`Lean.m`, and `Prelude.lean`, which is inlined into every generated file). The C side is `src/print_lean.c` (LeanForm) and `src/lean_run.c` (the runner).

## LeanForm
- `LeanForm[expr]`: prints as Lean 4 syntax for `expr`.
- `ToString[expr, LeanForm]`: the same text as a string.

**Features**:
- `Protected`. An inert display wrapper, like `FortranForm`.
- Lean 4 precedences: `↔` 20, `→` 25, `∨` 30, `∧` 35, `¬` 40, relations 50, `+ -` 65, `* /` 70, `^` 75, application 1024. There is no line breaking.
- Arithmetic: `x ^ 2`, `a * b`, `x / (y * z)` for negative powers, `1 / x` for `x^-1`, `Real.sqrt x` for `Sqrt[x]`. A rational prints as `(3 / 4)`, so `3/4 x` is `(3 / 4) * x`. `Plus` prints negative terms as subtraction (`x - 3`, `a - 2 * b`), and leads with a non-negative term when the canonical order puts a negative one first. Unary negation is `-(…)` unless its operand is an atom (`-x`, `-(x ^ 2)`, `-(2 * x)`), and a negative literal as an operand is parenthesised (`f (-3)`).
- Relations and logic: `= ≠ < ≤ > ≥`; chains such as `a < b <= c` and `Inequality[...]` become conjunctions (`a < b ∧ b ≤ c`); `And Or Not Implies Equivalent` are `∧ ∨ ¬ → ↔`.
- Binders: `ForAll[{x, y}, body]` is `∀ x y, body`, `ForAll[x, cond, body]` is `∀ x, cond → body`, `Exists[x, cond, body]` is `∃ x, cond ∧ body`, and `Function[{x, y}, b]` is `fun x y => b`. A binder is parenthesised unless it is the whole expression.
- `List` is `[a, b]`. `f[x, y]` is the application `f x y`, with compound arguments parenthesised (`f (g x) (1 / 2)`).
- Mathlib names: `Log`, `Exp`, `Sin`, `Cos`, `Tan`, `ArcTan`, `Sinh`, `Cosh` are `Real.log` … `Real.cosh`, and `Log[b, x]` is `Real.logb b x`. `Pi` is `Real.pi`, `E` is `Real.exp 1`, `I` is `Complex.I`.
- Symbols: a plain identifier passes through. Anything else, such as a context mark, a `$`, or a Lean keyword like `fun`, is escaped as `«name»` (`` a`b `` is `` «a`b» ``).
- Strings are Lean string literals. Machine reals print as decimal literals.
- `Format[f[...], LeanForm] := rhs` display rules apply.

```mathematica
In[1]:= LeanForm[(x + y)^2 == x^2 + 2 x y + y^2]
Out[1]= (x + y) ^ 2 = x ^ 2 + 2 * x * y + y ^ 2

In[2]:= LeanForm[ForAll[{x, y}, x^2 + y^2 >= 2 x y]]
Out[2]= ∀ x y, x ^ 2 + y ^ 2 ≥ 2 * x * y

In[3]:= ToString[{Function[{x}, x + 1], 3/4 a, Log[b], f[-3]}, LeanForm]
Out[3]= "[fun x => x + 1, (3 / 4) * a, Real.log b, f (-3)]"
```

### `` Lean`Internal`ToLean ``
`` Lean`Internal`ToLean[expr, type, coreOnly, {heads...}] `` is the **strict** translator that the certificate generators use. It is internal to the package: the C pieces (`` Lean`Internal`ToLean ``, `` Lean`Internal`RunProcess ``, `` Lean`Internal`VerifySource ``) live in a context that is not on `$ContextPath`, so their short names never capture another package's symbols. `type` is `"Int"`, `"Rat"`, `"Real"` or `None`: it types the binders (`∀ x : Rat, …`), ascribes rational literals (`(3 / 4 : Rat)`), and ascribes the left side of a relation between numbers (`(2 + 2 : Int) = 4`) so that numerals do not default to `Nat`. `heads` lists the Lean names that may appear as applied functions. Unlike the display form, it **refuses** every node whose Lean text would not mean the same thing as the Mathilda expression, giving `$Failed` and `LeanForm::untrans`:

| Refused | Reason |
|---|---|
| `1.5`, any inexact number | `inexact number; Rationalize it first` |
| `x/2`, `x^-1` with `type` `"Int"` | `division over Int` (Int `/` is floor division) |
| `Log`, `Sqrt`, `Pi`, `I`, … with `coreOnly` `True` | needs Mathlib (Tier 1) |
| `x^(1/3)`, `x^n` | `non-integer exponent` |
| `Infinity`, `Indeterminate`, `ComplexInfinity` | `infinite or indeterminate quantity` |
| `f[x]` with `f` not in `heads` | `no Lean counterpart for this head` |
| `#^2 &` | `pure function with Slot` |

`Format` rules do not apply in the strict form, so a display rule can never change the statement handed to Lean.

## Certify
- `Certify[comp]`: evaluates the computation `comp` and certifies its result.
- `Certify[a == b]`: certifies the claim `a == b`.
- Supported in this release (`"Kind"`):
  - `"Identity"`: claims `a == b`, and `Expand`, `Collect`, `Together`, `Apart`, `Cancel`, `Simplify`, `FullSimplify`. For a computation, the claim is `input == result`, with the input taken from the held call, so it is the user's expression that is proved equal.
  - `"Factorization"`: `Factor`, `FactorTerms`. `"TrustedSteps" -> {"irreducibility of the factors (not certified)"}`.
- **Polynomial identities** are stated as `∀ x y : T, a = b` over the sorted `Variables`, with `T` being `Int` when every coefficient is an integer and there is no division, and `Rat` otherwise. The proof is `by intros; grind`.
- **Rational-function identities** are cleared: with `a = na/da` and `b = nb/db` (`Together`), the theorem is the polynomial identity `na * db = nb * da`. `"Hypotheses"` records `da ≠ 0` and `db ≠ 0` (informational in Tier 0), and `"TrustedSteps"` records the `Together` normal form.
- Both sides of a claim are evaluated by Mathilda first. A symbol that holds a claim (`c = a == b; Certify[c]`) is evaluated and then dispatched. Anything else (a transcendental claim, an unsupported head) issues `Certify::nocert` and returns `$Failed`.

**Options** (`Options[Certify]`; `SetOptions` works):

| Option | Default | Meaning |
|---|---|---|
| `Engine` | `"Lean"` | `"Lean"`, or `"Mathilda"`: a heuristic check through `PossibleZeroQ` / `Simplify` (`"Tier" -> "Heuristic"`), for contrast. Anything else issues `Certify::engine` and returns `$Failed`. |
| `Tier` | `Automatic` | `0` (core Lean), `1` (Mathlib; needs `$LeanProject`, otherwise `Certify::nomathlib`) or `Automatic`. |
| `Domain` | `Automatic` | `Integers`, `Rationals` or `Reals`. `Reals` needs Tier 1 (`Certify::tier` without a project). |
| `Method` | `Automatic` | The tactic: `"grind"`, `"decide"`, `"omega"`, or any Lean tactic string. It goes into the Lean source file only, never into a shell command. |
| `TimeConstraint` | `60` | Seconds, mapped to `set_option maxHeartbeats (40000 t) in` on the theorem. This is a deterministic resource bound, not a wall-clock kill. |
| `Assumptions` | `{}` | Extra hypotheses: `Assumptions -> {x == 1}` gives `∀ x : Int, x = 1 → …`. |
| `Verify` | `True` | `False` generates the source without running Lean (`"Verified" -> Missing["NotChecked"]`). |
| `LeanExecutable` | `Automatic` | Overrides `$LeanExecutable` for one call. |
| `LeanProject` | `Automatic` | Overrides `$LeanProject` for one call (`None` forces no project). |

**Messages**: `Certify::nocert` (no certificate for this input), `Certify::engine`, `Certify::tier`, `Certify::nomathlib`, and `Certify::rejected`, which is issued when Lean judged the certificate and rejected it.

## ProofObject
- `ProofObject[<|...|>]`: the result of `Certify`.
- `ProofObject[comp, opts]`, with `comp` not an `Association`: the same as `Certify[comp, opts]`, so `ProofObject[Factor[x^2 - 1], Engine -> "Lean"]` works.
- `obj["prop"]`: a property. `obj["Properties"]` gives the sorted list, and an unknown name gives `Missing["UnknownProperty", name]`.
- Attributes `HoldFirst`, `Protected` (and the same for `Certify`).
- It displays as `ProofObject[«Factorization» Verified: True | Lean Tier 0 | 0.51 s]`.

| Property | Content |
|---|---|
| `"Kind"` | `"Identity"` or `"Factorization"` (more kinds come with later releases) |
| `"Input"` | `HoldForm` of the computation or the claim |
| `"Result"` | The computation's value; for a claim, the claim with evaluated sides |
| `"Theorem"` | The certified statement as a held Mathilda expression (`ForAll[{x}, …]`) |
| `"LeanTheorem"` | The statement as Lean text |
| `"ProofScript"` | The tactic block, e.g. `by intros; grind` |
| `"LeanSource"` | The complete Lean file: header, prelude, `theorem mathilda_cert_1`, `#print axioms` |
| `"Certificate"` | Raw Mathilda data: `"Form"` (`"Polynomial"` or `"ClearedDenominators"`), `"Variables"`, and the numerators and denominators when cleared |
| `"Verified"` | `True`, `False`, `Missing["NotChecked"]` (`Verify -> False`) or `Missing["NotAvailable"]` (no Lean found) |
| `"Engine"`, `"Tier"`, `"Domain"` | As the options; `"Tier"` is `"Heuristic"` for `Engine -> "Mathilda"` |
| `"Hypotheses"` | Side conditions such as denominators `≠ 0` |
| `"TrustedSteps"` | What Lean did not check |
| `"Axioms"` | The theorem's axioms, from `#print axioms` |
| `"Logs"` | Lean's messages for the theorem: `<|"Severity", "Line", "Text"|>` |
| `"Timing"` | Lean's wall time in seconds |
| `"Cached"` | Whether the verdict came from the session cache (always `False` in this release) |
| `"LeanVersion"` | `lean --version` of the checker that ran |

**The `"Verified" -> True` contract.** All of these must hold: Lean's exit code is 0; no message has severity `error`; no message has kind `hasSorry`; and the theorem's `#print axioms` report (appended by the generator) names only `propext`, `Classical.choice` and `Quot.sound`. `sorryAx`, and the auxiliary axiom that `native_decide` adds, both fail the last test. An error in the theorem's own lines, or outside every theorem (the header or the prelude), makes it `False`.

```mathematica
In[1]:= Certify[Factor[x^4 - 1]]
Out[1]= ProofObject[«Factorization» Verified: True | Lean Tier 0 | 0.51 s]

In[2]:= %["LeanTheorem"]
Out[2]= "∀ x : Int, x ^ 4 - 1 = (x - 1) * (1 + x) * (1 + x ^ 2)"

In[3]:= Certify[(x + 1)^3 == x^3 + 3 x^2 + 3 x + 2]["Verified"]
Certify::rejected: Lean rejected the certificate for (x + 1)^3 == x^3 + 3 x^2 + 3 x + 2; see obj["Logs"]
Out[3]= False
```

**Lean catches the CAS.** `PossibleZeroQ[(1 + x)^20 - Expand[(1 + x)^20] + 10^-25 x]` is `True`, which is wrong: the expression is the polynomial `x/10^25`. The heuristic engine reproduces the mistake and Lean refutes it:

```mathematica
In[4]:= c = (1 + x)^20 == Expand[(1 + x)^20] - 10^-25 x;
In[5]:= {Certify[c, Engine -> "Mathilda"]["Verified"], Quiet[Certify[c]["Verified"]]}
Out[5]= {True, False}
```

## LeanVerify
- `LeanVerify[src]`: checks the Lean 4 source string `src` and gives `<|"Success", "ExitCode", "Messages", "Axioms", "Timing", "Command"|>`.

**Features**:
- `Protected`. Options `LeanExecutable`, `LeanProject` and `KeepFiles` (`True` keeps the temporary `.lean` file and returns its path as `"SourceFile"`).
- A `#print axioms` line is appended for every `theorem` the source declares and does not already report. `"Success"` applies the `Verified` contract to the whole file, and also requires an axiom report for every declared theorem. A `sorry` or a `native_decide` proof gives `False`. `example`s cannot be axiom-checked.
- `"Messages"` is a list of `<|"Severity", "Kind", "Line", "Text"|>` parsed from `lean --json`. `"Axioms"` maps each theorem name to its axioms.
- With a project (`LeanProject`), the command is `cd project && <lake> env <lean> --json FILE`, where `<lake>` is the `lake` binary beside the Lean executable.
- `LeanVerify::nolean` is issued when no Lean executable can be run. `LeanVerify::badpath` is issued when a path cannot be quoted safely for the shell (Windows `cmd.exe` only), and `LeanVerify::fail` when the process could not be started.

```mathematica
In[1]:= LeanVerify["theorem ok : 2 + 2 = 4 := by decide"]["Success"]
Out[1]= True

In[2]:= LeanVerify["theorem nd : 2 + 2 = 4 := by native_decide"]["Success"]
Out[2]= False
```

**Process execution.** Lean runs through ISO C `system()`. The command line holds only executable paths, flags written by Mathilda and temp-file paths, each single-quoted for the POSIX shell (with `set -C`, so the output redirections never follow a planted file). The Lean source, claims and tactic strings are always written to the `.lean` file. Temp files go in `$TemporaryDirectory` and are removed afterwards. `system()` blocks, and `TimeConstrained` measures this process's CPU time, so it does not interrupt a Lean run. The bound on a run is the deterministic `maxHeartbeats` of each theorem.

## Settings
- `$LeanExecutable`: `Automatic`, or a path. `Automatic` tries, in order, the environment variable `MATHILDA_LEAN`, the newest `~/.elan/toolchains/*/bin/lean` (called directly, so no `lean-toolchain` file is needed), and `lean` on `PATH`. The first hit is remembered for the session.
- `$LeanProject`: `None`, or the directory of a Lake project with Mathlib built, used for `Tier -> 1`. It defaults to the environment variable `MATHILDA_LEAN_PROJECT`.
- `$LeanVersion`: the `lean --version` string, or `Missing["NotAvailable"]`.
