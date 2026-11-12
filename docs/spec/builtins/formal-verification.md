# Formal verification (Lean 4)

Mathilda can hand its answers to the [Lean 4](https://lean-lang.org) theorem prover for checking. `Certify` computes a result, writes a self-contained Lean 4 file that states the result as a theorem with a proof script, and runs the Lean kernel on it. The answer is a `ProofObject`.

**Trust model.** Mathilda is an untrusted oracle and the Lean kernel is the trusted checker. `"Verified" -> True` means that Lean accepted the theorem with no `sorry` and with no axiom outside `{propext, Classical.choice, Quot.sound}`, and that a **kernel replay** confirmed it: the source is compiled to an `.olean` and a second, fresh `lean` process (`src/internal/lean/Replay.lean`, which runs none of the checked file's code) re-checks every declaration of the module in a new kernel environment and computes their axioms itself. A declaration the elaborator added without a kernel check (`debug.skipKernelTC`, `addDecl` from a user tactic) fails the replay. The replay reads the `.olean` that the elaborating process wrote, so it is only as trustworthy as that process: the foundation is that the checked file contains **no user-defined code** (the allowlist grammar of [LeanVerify](#leanverify), which every generated file passes too), and the theorem is the input **as written** (a computation inside a claim is never evaluated first; see `Certify::eval`). Each `ProofObject` also lists its `"TrustedSteps"`: the facts that Lean did **not** check for that kind of certificate, such as the `Together` normal form used to clear denominators, or the irreducibility of the factors of a `Factor` result. When the theorem is the whole claim, the list is empty.

**What "as written" leaves to Mathilda.** The claim is held, but its `Plus`, `Times`, `Power` and `Rational` nodes are put in Mathilda's canonical arithmetic form before translation (`(x + 1)^2 == x^2 + 2 x + 1` becomes `(1 + x) ^ 2 = 1 + 2 * x + x ^ 2`; a rational literal is reduced), and so is a claimed antiderivative `F`. Lean checks the canonicalised statement, not the characters the user typed. This arithmetic normalisation is therefore part of the trust base, alongside the Lean kernel, core Lean / Mathlib as built, and the per-kind `"TrustedSteps"`: Mathilda cannot make Lean accept a false theorem, and it can make the theorem differ from the claim only through the listed `"TrustedSteps"` and the arithmetic normalisation of the claim. `"Theorem"` and `"LeanTheorem"` show exactly what Lean proved. Only Lean can make `"Verified"` `True`: the heuristic engine (`Engine -> "Mathilda"`) reports `Missing["Heuristic", True | False]`, never `True`.

Lean is optional. With no Lean installed, `Certify[..., Verify -> False]` still generates the Lean source and the statement, and `"Verified"` is `Missing["NotChecked"]`. Tier 0 certificates need only core Lean 4 (v4.35 or newer, for `grind`) and no Mathlib. The package lives in `src/internal/lean/` (`Lean.m`, the prover-neutral core, which loads `LeanBackend.m`, the Lean backend, `ExactBackend.m`, the Exact backend, `Engines.m`, engine lists, the heuristic engine and consensus, `Refute.m`, the counterexample search and `Refute`, and `Conjecture.m`, the guessers and `Conjecture`; and `Prelude.lean`, which is inlined into every generated file). The C side is `src/print_lean.c` (LeanForm) and `src/lean_run.c` (the runner).

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
`` Lean`Internal`ToLean[expr, type, coreOnly, {heads...}, closed] `` (`closed` optional, default `False`) is the **strict** translator that the certificate generators use. It is internal to the package: the C pieces (`` Lean`Internal`ToLean ``, `` Lean`Internal`RunProcess ``, `` Lean`Internal`VerifySource ``) live in a context that is not on `$ContextPath`, so their short names never capture another package's symbols. `type` is `"Int"`, `"Rat"`, `"Real"` or `None`: it types the binders (`∀ x : Rat, …`), ascribes rational literals (`(3 / 4 : Rat)`), and ascribes the left side of a relation between numbers (`(2 + 2 : Int) = 4`) so that numerals do not default to `Nat`. `heads` lists the Lean names that may appear as applied functions. Unlike the display form, it **refuses** every node whose Lean text would not mean the same thing as the Mathilda expression, giving `$Failed` and `LeanForm::untrans`:

| Refused | Reason |
|---|---|
| `1.5`, any inexact number | `inexact number; Rationalize it first` |
| `x/2`, `x^-1` with `type` `"Int"` | `division over Int` (Int `/` is floor division) |
| `Log`, `Sqrt`, `Pi`, `I`, … with `coreOnly` `True` | needs Mathlib (Tier 1) |
| `x^(1/3)`, `x^n` | `non-integer exponent` |
| `Infinity`, `Indeterminate`, `ComplexInfinity` | `infinite or indeterminate quantity` |
| `f[x]` with `f` not in `heads` | `no Lean counterpart for this head` |
| `#^2 &` | `pure function with Slot` |
| a symbol whose name contains `«`, `»` or a control character | it could close the `«…»` escape early |
| two distinct symbols with the same Lean text | they would become one Lean variable |
| with `closed` `True`: a symbol not bound by an enclosing `ForAll` / `Exists` / `Function` and not in `heads` | `free variable (not bound by a quantifier)`; a free `default` or `id` would silently mean Lean's global |

Every theorem statement a generator emits is translated with `closed -> True`; open pieces (the `have` steps of a proof, after `intro`) are not.

`Format` rules do not apply in the strict form, so a display rule can never change the statement handed to Lean.

## Certify
- `Certify[comp]`: evaluates the computation `comp` and certifies its result.
- `Certify[a == b]`: certifies the claim `a == b`.
- `Certify[{c1, c2, …}]`: certifies a batch in one Lean run and gives a list.
- Supported in this release (`"Kind"`):
  - `"Identity"`: claims `a == b`, and `Expand`, `Collect`, `Together`, `Apart`, `Cancel`, `Simplify`, `FullSimplify`. For a computation, the claim is `input == result`, with the input taken from the held call, so it is the user's expression that is proved equal.
  - `"Factorization"`: `Factor`, `FactorTerms`. `"TrustedSteps" -> {"irreducibility of the factors (not certified)"}`.
- **Polynomial identities** are stated as `∀ x y : T, a = b` over the sorted `Variables`, with `T` being `Int` when every coefficient is an integer, there is no division **and there are no hypotheses**, and `Rat` otherwise. A hypothesis-free integer polynomial identity that holds over `Int` holds in every commutative ring, so `Int` loses nothing there.
- **Identities with hypotheses** (`Assumptions`) are stated over a **generic field**, matching Mathematica's semantics (variables are complex unless an inequality makes them real, and an inequality makes only **its own** variables real): `∀ {K : Type} [Lean.Grind.Field K] [Lean.Grind.IsCharP K 0], ∀ x y : K, hyps → a = b` — every field of characteristic 0, so ℚ, ℝ and ℂ alike (`"Domain" -> Complexes`). When the hypotheses contain an order relation the binders are those of every ordered field, `[Lean.Grind.Field K] [LE K] [LT K] [Std.LawfulOrderLT K] [Std.IsLinearOrder K] [Lean.Grind.OrderedRing K]` (`"Domain" -> Reals`) -- but only if **every** variable is real: it occurs in an order relation of the assumptions, or is declared with `Element[x, Reals]` or `Element[{x, y}, Reals]` in the assumptions. Otherwise there is no certificate: `Certify[x == 7, Assumptions -> x^2 + 4 == 0 && y > 0]` issues `Certify::realvars` (`x` is complex there, `x = 2 I`); `Assumptions -> x^2 + 4 == 0 && y > 0 && Element[x, Reals]` is the explicit real statement. With every variable declared real and no order relation the statement is over every ordered field too. Over `Int` or `Rat` a hypothesis with no solution there made the statement vacuous: `Certify[x == 7, Assumptions -> x^2 + 4 == 0]` would prove `∀ x : Rat, 4 + x^2 = 0 → x = 7` (false for `x = 2 I`); now it is `False`, and so is `Assumptions -> x^2 == 2`. The theorem takes the binders as parameters (`theorem mathilda_cert_1 {K : Type} [...] : …`); `"LeanTheorem"` and the footer show the full `∀ {K : Type} …` statement. `Domain -> Integers` / `Rationals` are honoured and shown in `"Domain"` (see the `Domain` table below). The proof is `by intros; grind`.
- **The `Domain` option** (v0.6522) has exactly one meaning per value, applied with or without hypotheses; any other value is refused with `Certify::domain` and `$Failed` (before v0.6522 an unknown value such as `Complexes` silently fell back to the automatic rules **without** their generic form, so `Certify[x == 7, Assumptions -> x^2 == -4, Domain -> Complexes, Method -> …]` proved `∀ x : Rat, x ^ 2 = -4 → x = 7` and reported `True`):

  | `Domain` | Statement | `"Domain"` shown |
  |---|---|---|
  | `Automatic` | Mathematica's default (variables complex unless an order relation or `Element[x, Reals]` makes them real). No hypotheses: `Int` for an integer polynomial identity, `Rat` otherwise -- equivalent, since a polynomial identity that holds on all of ℤ or ℚ holds in ℚ[x] and so in every field of characteristic 0. With hypotheses: as `Complexes`. | `Integers` / `Rationals`, or as `Complexes` |
  | `Complexes` | Every field of characteristic 0 (`[Lean.Grind.Field K] [Lean.Grind.IsCharP K 0]`), hypotheses or not; every ordered field when every variable is declared real; an order relation on a variable that is not real is `Certify::realvars`. | `Complexes` (or `Reals`) |
  | `Reals` | Every ordered field (all variables real by declaration), hypotheses or not. Tier 0, no Mathlib (before v0.6522: `Real`, Tier 1). | `Reals` |
  | `Rationals` | `∀ x : Rat` | `Rationals` |
  | `Integers` | `∀ x : Int` (the claim over ℤ) | `Integers` |

  `Element` assumptions are read the same way for every domain: `Element[x, Reals]` declares `x` real (over `Int` / `Rat` it is implied and dropped); `Element[x, Complexes]` is implied and dropped; `Element[x, Rationals]` is implied over `Rat` and `Int`, `Element[x, Integers]` over `Int`; any other `Element` (`Element[x, Integers]` with `Automatic`, for example) is `Certify::element`. A variable printed as `K` or `Rat` cannot be stated in the generic form of an identity (the type variable is `K`): `Certify::nocert`. Inequality, implication and refutation statements name the type variable `K1`, `K2`, ... instead when the claim uses `K` (v0.6528).
- **Rational-function identities** are cleared with the **unreduced** numerators and denominators of each side, read off the expression as written (`rawFraction`: sums over the product of their denominators, no cancellation): with `a = na/da` and `b = nb/db`, the theorem is the polynomial identity `na * db = nb * da`. `Certify[(x^2 - 1)/(x - 1) == x + 1]` proves `∀ x : Int, x ^ 2 - 1 = (x - 1) * (1 + x)` (with `Together`, which cancels, it was the empty `1 + x = 1 + x`). `"Hypotheses"` records `da ≠ 0` and `db ≠ 0` (informational in Tier 0), and `"TrustedSteps"` records the field algebra from the cross-multiplied identity to the rational one where the denominators are nonzero.
- `"Antiderivative"`: `Integrate[f, x]` with `f` a rational function of `x` over ℚ (see [Antiderivatives](#antiderivatives-tier-0)), and the claim `Integrate[f, x] == F` (either side), which certifies the **claimed** `F` -- not Mathilda's -- as an antiderivative of `f`. The claim as written is **not** what is certified (v0.6522): in Mathematica `Integrate[2 x, x] == x^2 + 5` is `False`, because `Integrate` returns one particular antiderivative, and Lean cannot decide which one without evaluating `Integrate`. What is certified is `D[F, x] == f`: `"Theorem"` shows `ForAll[x, D[F, x] == f]` (Tier 1: the `HasDerivAt` statement), `"Input"` stays the claim as written, `Certify::antideriv` says so, and `"TrustedSteps"` gains "the claim Integrate[f, x] == F is certified as D[F, x] == f …; Integrate[f, x] == F itself holds only up to an additive term with derivative 0, which is locally constant (constant on each interval where F and f are defined; it may differ across a point where they are not, as ArcTan[1/x] and -ArcTan[x] do at x = 0)". To certify Mathilda's own antiderivative, use `Certify[Integrate[f, x]]`; to state the derivative claim directly, `Certify[D[F, x] == f]` works when `F` is a polynomial.
  - `"Bezout"`: `PolynomialExtendedGCD[a, b, x]` → `{g, {s, t}}`. The theorem is `∀ x, s·a + t·b = g ∧ a = g·qa ∧ b = g·qb`, with the quotients from `PolynomialQuotient`. It proves that `g` is a common divisor and a combination of `a` and `b`, so every common divisor divides `g`: `g` is a gcd. `"TrustedSteps" -> {}`.
  - `"IdealMembership"`: `PolynomialReduce[q, {p1, …}, vars]` → `{{c1, …}, r}` is certified as the division identity `∀ vars, q = c1·p1 + … + r` (`"TrustedSteps"` notes that the remainder being fully reduced is not certified). The claim `Implies[p1 == q1 && …, p == q]` is stated over every field of characteristic 0 (as identities with hypotheses above; `"Domain" -> Complexes`) and certified with cofactors `c_i` such that `p - q = Σ c_i (p_i - q_i)`: from `PolynomialReduce` when its remainder is 0, otherwise from a degree-bounded ansatz (degree 0 to 4) whose coefficient equations are solved by `LinearSolve`. The proof is `have k_i : p_i - q_i = 0`, `have hc : p - q = Σ c_i * (p_i - q_i)`, `rw [k_i] at hc`, `grind`.
  - `"Roots"`: `Solve[eqs, vars]` with rational solutions. A single univariate polynomial equation that splits over ℚ is certified **complete**: `∀ x : Rat, lhs = rhs ↔ x = r1 ∨ … ∨ x = rk`, through the factorization `lhs - rhs = lc · Π (x - r_i)^m_i` (`"Certificate"` has `"Roots"` and `"Multiplicities"`). Otherwise (systems, inequality / `Unequal` constraints, a polynomial with irreducible non-linear factors, or division) only **soundness** is certified: every solution substituted into **every constraint** (`==`, `!=`, `<`, `<=`, `>`, `>=`, joined by `&&` or a list), checked by `decide +kernel`, with completeness listed in `"TrustedSteps"`. Lean's `x / 0 = 0` would otherwise accept a "solution" that zeroes a denominator (`Solve[(x-1)/(x^2-1) == (x-1)/(2x-2), x]` returns the spurious `x -> 1`), so every denominator of the system is also stated nonzero at every solution (`"Certificate"["Denominators"]`); the spurious root is then `False`. Completeness needs a single polynomial equation with no other constraint and no division. A system with `||`, `Element` or any other head, or a symbol other than the variables, has no certificate. Radical or complex roots have no certificate yet.
  - `"MatrixInverse"`, `"Determinant"`, `"LinearSolve"`, `"NullSpace"`: `Inverse[A]`, `Det[A]` (n ≤ 6), `LinearSolve[A, b]` and `NullSpace[A]` for matrices of rational numbers. The theorems use the prelude's `Mathilda.matMul`, `Mathilda.idMat`, `Mathilda.det` (Laplace expansion) and `Mathilda.matVec`, and the kernel evaluates them (`decide +kernel`): `matMul A B = idMat n ∧ matMul B A = idMat n`, `det A = d`, `matVec A x = b`, and `matVec A v = 0` for each basis vector (`NullSpace`'s completeness is a trusted step). These functions use `zipWith` / `getD`, which silently truncate or pad ragged input, so every statement also fixes the shapes with `Mathilda.HasShape A r c` (`r` rows of `c` entries) and `Mathilda.HasLength v n`: `HasShape A n n ∧ HasShape B n n ∧ …`, `HasShape A n n ∧ det A = d`, `HasShape A m n ∧ HasLength x n ∧ matVec A x = b`, and `HasLength v n ∧ matVec A v = 0` per basis vector.
  - `"Nonnegativity"`: claims `ForAll[vars, p >= q]` (or `<=`, or a bare `p >= q`), optionally with hypotheses `ForAll[vars, g1 >= 0 && …, p >= q]` or `Implies[g1 >= 0 && …, p >= q]` (see [Sums of squares and optimality](#sums-of-squares-and-optimality)).
  - `"Optimality"`: `Minimize[{f, cons…}, vars]` and `Maximize[…]` with a polynomial objective, linear equalities and polynomial (non-strict) inequalities.
  - `"Primality"`: `PrimeQ[n]` and `FactorInteger[n]` (see [Primality](#primality)).
  - `"ClosedFormSum"`: `Sum[t, {k, a, n}]` and claims `Sum[t, {k, a, n}] == closed` (see [Closed-form sums](#closed-form-sums)).
- **Batches.** `Certify[{c1, c2, …}]` gives a list of `ProofObject`s from **one** Lean run (per tier). The theorems are `mathilda_cert_1`, `mathilda_cert_2`, … in input order, every object's `"LeanSource"` is the shared file, and each theorem gets its own verdict from the messages in its own line span, so one bad member is `False` while its siblings stay `True`. An item with no certificate gives `$Failed` in its slot.
- **Session cache.** A Lean verdict (`True` or `False`) is remembered for the session, keyed by the theorem text, the tier, the executable and the project. Certifying the same thing again, alone or inside a batch, gives `"Cached" -> True` with the original `"Timing"` and does not run Lean.
- **The input is translated as written.** Mathilda does not evaluate a claim before translating it (since v0.6520; before, `Simplify[(x^2 - 1)/(x - 1)] == x + 1` became the tautology `1 + x = 1 + x`). The input is walked held:
  - numbers, symbols, `Plus`, `Times`, `Subtract`, `Divide`, `Minus`, integer powers, relations, `&&`, `Implies`, `ForAll` and lists are kept (Mathilda puts them in its canonical arithmetic form -- flattening, ordering, collecting numeric coefficients, reducing rationals -- before translation; Lean checks that canonical form, so this normalisation is part of the trust base, see the trust model above); numerators and denominators of an identity are read off the held sides, so `x/x == 1` is `∀ x, x = x` with the hypothesis `x ≠ 0`, not `1 = 1`;
  - `Expand`, `Together`, `Cancel`, `Apart`, `Factor`, `FactorTerms` and `Collect` are read through: they preserve the value of a rational function, so Lean checks the user's own expression (`Factor[x^2 - 1] == (x - 1)(x + 1)` proves `x ^ 2 - 1 = (x - 1) * (1 + x)`; a `Factor` / `FactorTerms` makes the kind `"Factorization"` with the irreducibility trusted step);
  - `D[F, x]` with `F` a polynomial is replaced by Mathilda's derivative, and `"TrustedSteps"` names it ("Mathilda's differentiation: D[F, x] = …");
  - a whole side `Integrate[f, x]` (indefinite) or `Sum[t, {k, a, n}]` goes to the antiderivative or closed-form-sum generator with the claimed other side;
  - the outermost head of a computation (`Certify[Simplify[e]]`) is the computation;
  - anything else inside the input (`Simplify`, `FullSimplify`, `Refine`, `Limit`, `Mod`, `Exp`, `Log`, `Sqrt`, `Sin`, a definite `Integrate`, `Det`, `N`, a user function, …) issues `Certify::eval` and returns `$Failed`: certify the computation itself, or state the claim with arithmetic only.

  `"Input"` is the input as written. `"Result"` and the heuristic engine use Mathilda's own evaluation of the two sides. A claim with no generator (a transcendental claim, an unsupported shape) issues `Certify::nocert` and returns `$Failed`.
- **Nothing in the claim is evaluated, and no value is substituted** (v0.6521). `Certify` and `ProofObject` are `HoldAll`, and the claim walker is `HoldAllComplete`, so an `Evaluate[...]` inside a claim is refused like any other computation (`ProofObject[Evaluate[Simplify[(x^2 - 1)/(x - 1)]] == x + 1]` issues `Certify::eval`; with `HoldAll` it used to be evaluated into the tautology `1 + x = 1 + x`). A symbol that is not a plain variable is refused with `Certify::value`, wherever it occurs (v0.6522: any `OwnValue`, `DownValue`, `UpValue`, `SubValue`, `NValue` or `DefaultValue`, or an attribute that changes evaluation such as `Constant`, `NumericFunction`, `Flat`, `Orderless`, `Listable` or a `Hold*`; `Protected`, `ReadProtected`, `Locked` and `Temporary` are harmless; the message names the kind, e.g. "t has UpValues ({HoldPattern[t^2] :> -1})"). Before v0.6522 `Protected` symbols were exempt and only `OwnValues` were looked at, so `e = Simplify[...]; Protect[e]; Certify[e == x + 1]` was a tautology and `t /: t^2 = -1; Certify[t^2 == -1]` proved `(-1 : Int) = -1`. ``System` `` symbols (the built-in heads and constants: `Plus`, `Power`, `Pi`, `I`, `Reals`, `True`, …) are accepted unless a user gave one an `OwnValue` or an `UpValue` (`I`'s own `I :> I` is not a value). Symbols in the contexts `` Lean`Private` `` and `` Lean`Internal` `` (Certify's generator symbols: surd variables, prelude placeholders, sum binders) are refused with `Certify::context`, so a claim cannot capture them (`Certify[Integrate[2, x] == Sqrt[2] Lean`Private`s1 x]` was `True`); the surd variables are also fresh symbols per call (`` Lean`Private`s1$n ``, printed `s1` in Lean). Examples: `x = 5; Certify[x^2 == 25]` and `e = Simplify[...]; Certify[e == x + 1]` give `$Failed` instead of certifying `(25 : Int) = 25` or a tautology, consistently with user functions (DownValues, `Certify::eval`). This includes a symbol that holds the whole claim (`c = a == b; Certify[c]`): its value was computed when it was assigned. `Clear` the symbol, or ask for the value-substituted claim explicitly with a top-level `Certify[Evaluate[c]]`. The `Assumptions` option is walked held with the same rules before the options are evaluated (`a = 0; Certify[y == 7, Assumptions -> y a == 1]` would otherwise make the hypothesis `False`): since v0.6522 it may contain arithmetic, relations, `&&`, `||`, `!`, `Implies`, `Element[vars, domain]`, lists and the value-preserving heads read through in claims; any other head (`Reduce`, `Simplify`, `D`, a user function, `Evaluate`) is `Certify::eval` -- `Assumptions -> Reduce[x^2 == 1 && x > 0, x]` used to give the hypothesis `x == 1`, computed by `Reduce`, with no trusted step. A default set with `SetOptions[Certify, Assumptions -> …]` was evaluated by `SetOptions` itself, so it cannot be checked as written: `Certify` then issues `Certify::optdefault` and returns `$Failed`, unless the call gives its own `Assumptions`.

### Antiderivatives (Tier 0)
`Certify[Integrate[f, x]]` computes `F = Integrate[f, x]` and proves, in core Lean, that `f` equals the formal derivative of `F`:

1. Each quadratic surd in `F` is written over primes, `Sqrt[d] = Π p^⌊k/2⌋ · Π s_p`, and `s_p` becomes a fresh variable (`s1, s2, …`; `"Certificate"["Surds"]` maps them back). Prime generators keep `Sqrt[2] Sqrt[3] = Sqrt[6]` visible.
2. `F` must be built from `x`, the surds and rational numbers by `Plus`, `Times`, integer powers and one level of `Log`, `ArcTan` or `ArcTanh`. It is split as `R + Σ c_i φ_i(u_i)` with each `c_i` free of `x`. Anything else (`RootSum`, `Root`, `I`, other functions) issues `Certify::tier`, because it needs a Tier 1 (Mathlib) certificate.
3. With `u_i = p_i/q_i`, `R = r/w` and `f = a/b` (`Together`), and derivatives of **polynomials only** taken by `D`, each term contributes `N_i/D_i`: `(p'q - pq')/(pq)` for `Log`, `(p'q - pq')/(q² + p²)` for `ArcTan`, `(p'q - pq')/(q² - p²)` for `ArcTanh`, and `(r'w - rw')/w²` for `R`.
4. The theorem is the cleared identity, with every piece written in its structured, unexpanded form so that it can be traced back to `F`:
   `∀ x s1 … : Rat, (Σ_i c_i N_i Π_{j≠i} D_j) b - a Π_j D_j = Σ_k (s_k^2 - d_k) K_k`.
   The cofactors `K_k` come from `PolynomialReduce` modulo `s_k^2 - d_k`. This is an **unconditional** polynomial identity in `ℚ[x, s…]`; it never assumes `s^2 = 2` over `ℚ` (which would be vacuous, since `√2 ∉ ℚ`). Mapping `s_k ↦ √d_k` is a ring homomorphism, so the identity holds for the real surds. The proof is `grind`.
5. Mathilda checks the identity itself first and records the result as `"Certificate"["PreCheck"]`. A failed pre-check does **not** suppress the theorem: it is still emitted, and Lean refutes it. A denominator `D_i` that vanishes identically would make the identity vacuous, so in that case there is no certificate.

What Lean proves is that the integrand equals the formal derivative of the antiderivative, as an identity of rational functions over `ℚ(√d…)`, cleared of denominators. Lean does **not** check the following, which `"TrustedSteps"` lists: the elementary derivative rules and the quotient rule (applied symbolically by Mathilda), the polynomial derivatives written into the statement (Mathilda's `D`), the `Together` normal forms `f = a/b`, `u_i = p_i/q_i` and of the rational part `R = r/w`, and that the cleared denominators are not identically zero (checked by Mathilda). `"Hypotheses"` lists the denominators `≠ 0`, with the surds written back as `Sqrt[d]`; they are informational in Tier 0. The Tier 1 certificate below proves the analytic statement instead, with nothing trusted.

```mathematica
In[1]:= c = Certify[Integrate[(x^4 - 3 x^2 + 6)/(x^6 - 5 x^4 + 5 x^2 + 4), x]];
In[2]:= {c["Result"], c["Verified"], c["Tier"]}
Out[2]= {ArcTan[x] + ArcTan[x^3] + ArcTan[(x - 3 x^3 + x^5)/2], True, 0}

In[3]:= Certify[Integrate[1/(x^2 + 2), x]]["LeanTheorem"]
Out[3]= "∀ x s1 : Rat, (1 / 2 : Rat) * s1 * (s1 * 2 - s1 * x * 0) * (2 + x ^ 2) - (2 ^ 2 + (s1 * x) ^ 2) = (s1 ^ 2 - 2) * 2"
```

`` Lean`Private`$antiderivOverride `` is a test hook: when it is not `None`, its value replaces the antiderivative that `Integrate` returned, so that a test can check that Lean rejects a wrong one.

```mathematica
In[4]:= Certify[Implies[x^2 + y^2 == 1 && x == y, 2 y^2 == 1]]["Certificate"]["Cofactors"]
Out[4]= {1, -x - y}

In[5]:= Certify[Solve[x^3 - 6 x^2 + 11 x - 6 == 0, x]]["LeanTheorem"]
Out[5]= "∀ x : Rat, 11 * x - 6 - 6 * x ^ 2 + x ^ 3 = 0 ↔ x = 1 ∨ x = 2 ∨ x = 3"

In[6]:= #["Verified"] & /@ Certify[{Factor[x^6 - 1], Inverse[{{2, 1}, {1, 1}}], (x + 1)^2 == x^2 + 1}]
Certify::rejected: Lean rejected the certificate for (x + 1)^2 == x^2 + 1; see obj["Logs"]
Out[6]= {True, True, False}
```

### Antiderivatives (Tier 1, Mathlib)
`Certify[Integrate[f, x], Tier -> 1]` proves `∀ x : Real, d1 ≠ 0 → … → HasDerivAt (fun x => F) f x` in Mathlib, with `"TrustedSteps" -> {}` when every hypothesis is a polynomial in `x` (see the witness rule below). It needs `$LeanProject` (or `LeanProject -> dir`, or the environment variable `MATHILDA_LEAN_PROJECT`): a Lake project with Mathlib built for the same toolchain as the Lean executable. Without one it issues `Certify::nomathlib`. Nothing in the default path depends on Mathlib. With `Tier -> Automatic`, an antiderivative outside the Tier 0 shape uses Tier 1 when a project is configured, and issues `Certify::tier` otherwise.

- The proof term is built from `F`'s tree, and the same walk writes the function body, so the two have the same shape: `hasDerivAt_id'`, `hasDerivAt_const`, `hasDerivAt_pow n x`, then `.add`, `.mul`, `.const_mul c`, `.div … h`, `.pow n`, `.log h` and `.arctan`. `ArcTanh[u]` is written as `(Real.log (1 + u) - Real.log (1 - u)) / 2`, because Mathlib has no `Real.artanh`.
- The hypotheses `d ≠ 0` are real binders. They cover every `Log` argument and denominator, the integrand's denominator, and each irreducible factor of these. Each factor hypothesis follows from its product, so adding them does not weaken the statement. `"Hypotheses"` lists them.
- **No vacuous statements** (v0.6523). A hypothesis that is zero at every real `x` under Lean's semantics would make the theorem true for any `F`: `Log[x^2] - 2 Log[x]` is nonzero over ℂ for `x < 0`, but `Real.log x = log |x|`, so it is identically 0 on ℝ, and `Integrate[x, x] == x^2/2 + 1/(Log[x^2] - 2 Log[x])` was `"Verified" -> True`. Mathilda now needs a witness: a rational point among `3, 7/2, 11, 1/3, …, -3, …` where every hypothesis is nonzero under `Real.log` semantics (evaluated at 40 digits, `|d| > 10^-20`). Without one the input is refused with `Certify::vacuous` and `$Failed`. The point is `"Certificate"["HypothesisWitness"]`; when some hypothesis is not a polynomial in `x` the numeric check is listed in `"TrustedSteps"`.
- The proof closes with `refine d.congr_deriv ?_`, a pointwise `simp only [Pi.add_apply, …]`, `eq_div_iff`, `field_simp`, then `ring1` (or `ring_nf` with `Real.sq_sqrt` for surd coefficients). A tactic that fails gives `"Verified" -> False` with the logs, never a wrong `True`.
- `import Mathlib` takes about 10–40 s per Lean run, so a batch of Tier 1 items shares one run.

```mathematica
In[1]:= c = Certify[Integrate[(1 + x)/(1 + x^2), x], Tier -> 1];
In[2]:= {c["Verified"], c["TrustedSteps"], c["Hypotheses"]}
Out[2]= {True, {}, {1 + x^2 != 0}}
```

### Sums of squares and optimality
`Certify[ForAll[vars, p >= q]]` looks for `p - q = Σ d_k ℓ_k² · m_k`, with every `d_k > 0` rational, every `ℓ_k` a polynomial and every `m_k` equal to `1` or to one of the hypotheses `g ≥ 0`. The squares come from Gram matrices through an exact rational LDLᵀ (`Q = L·diag(d)·Lᵀ`, `ℓ = Lᵀ·basis`). The Gram matrices are found in one of three ways (`"Certificate"["Method"]`):

- `"Quadratic"`: a quadratic form with no hypotheses has a unique Gram matrix over `{1, vars}`, so no numeric engine and no LAPACK is needed.
- `"SOSCertificate"`: [`` Reduce`SOSCertificate ``](solutions-of-equations.md#reducesoscertificate), the exact Putinar certificate of `Reduce`'s sum-of-squares stage.
- `"ExactGram"`: when that engine declines a problem with no hypotheses, the Gram matrix over the half-degree monomials is pinned by coefficient matching, by `Q·m(z) = 0` at every rational real zero `z`, and by zero rows where a diagonal entry is forced to 0. Free parameters are then set to 0. This handles problems whose SDP is not strictly feasible, such as `(x² - 1)² + y²` with its two zeros.

The proof states the hints `have : 0 ≤ (ℓ)^2 := Mathilda.sqnn` and `have : 0 ≤ (ℓ)^2 * (g) := Mathilda.mulnn Mathilda.sqnn hg`, then the explicit identity `have : p - q = Σ d (ℓ)^2 (g) := (by grind)`, and closes with `grind`.

`Certify[Minimize[{f, cons…}, vars]]`:
1. Linear equalities are solved for the last variable with a constant coefficient. In Lean this is `have : v = … := (by grind); subst`.
2. At the optimum `x*`, the multipliers `λ_j ≥ 0` of the active inequalities `g_j ≥ 0` solve `∇f = Σ λ_j ∇g_j` exactly (`LinearSolve`). Inactive constraints get 0, and a negative multiplier means no certificate.
3. `σ = f - m - Σ λ_j g_j` must be a sum of squares, as above.

The theorem is `(∀ x, cons → f ≥ m) ∧ cons(x*) ∧ f(x*) = m`, with `≤` for `Maximize`. `"Certificate"` holds `"Multipliers"` (one per inequality), `"Squares"`, `"Eliminated"` and `"LinearProgram"`. When `σ ≡ 0`, the problem is an LP and the multipliers are its dual certificate: weak duality plus attainment is optimality.

Both kinds are stated over **every ordered field** (`∀ {K : Type} [Lean.Grind.Field K] [LE K] [LT K] [Std.LawfulOrderLT K] [Std.IsLinearOrder K] [Lean.Grind.OrderedRing K], …`), which includes ℝ, so `"Domain" -> Reals` and `"TrustedSteps" -> {}` (an earlier version proved them over `Rat` and listed the unchecked transfer to ℝ). `Mathilda.sqnn` / `Mathilda.mulnn` are core's generic `Lean.Grind.OrderedRing` lemmas.

```mathematica
In[1]:= c = Certify[Maximize[{3 x + 2 y, x + y <= 4, x + 3 y <= 6, x >= 0, y >= 0}, {x, y}]];
In[2]:= {c["Verified"], c["Certificate"]["Multipliers"], c["Certificate"]["LinearProgram"]}
Out[2]= {True, {3, 0, 0, 1}, True}

In[3]:= Certify[ForAll[{x, y}, x^2 + y^2 >= 2 x y]]["ProofScript"]
Out[3]= "by
  intro x y
  have mth_s1 : 0 ≤ (x - y) ^ 2 := Mathilda.sqnn
  have mth_e : x ^ 2 - 2 * x * y + y ^ 2 = (1) * (x - y) ^ 2 := (by grind)
  grind"
```

### Primality
A flat Pratt certificate lists steps `(p, a, factors of p - 1)` in dependency order, with 2 as the base. The kernel evaluates the prelude's `Mathilda.prattCheck` with `decide +kernel`. The witness `a` is `PrimitiveRoot[p]`, or the first `a ≥ 2` that satisfies Lucas's conditions (checked with `PowerMod`). Each step checks `2 < p`, that the listed factors multiply out to exactly `p - 1` (`prodPow fs = p - 1`, so no prime factor can be left out), `a^(p-1) ≡ 1`, `a^((p-1)/q) ≢ 1` for every listed `q`, and that every `q` is an earlier step (or 2). `Mathilda.powMod` is exact for every exponent: its fuel is `e + 1`, which bounds the `log2 e + 1` halving steps (an earlier fixed fuel of 256 bits accepted the composite `2^256 + 1` with the certificate `[(2, 256)]`). The kernel only unfolds as many fuel steps as there are bits, so large exponents cost nothing extra: `PrimeQ[51 2^257 + 1]` checks in well under a second, and `p = NextPrime[2^256]; Certify[PrimeQ[p]]` is `True` after about 118 s (`Certify[PrimeQ[NextPrime[2^256]]]` is refused with `Certify::eval`: the input would have to be evaluated by Mathilda first), nearly all of it Mathilda finding the certificate (`PrimitiveRoot` for the 39-digit factor of `p - 1`); Lean takes 0.8 s.

- `Certify[FactorInteger[n]]` emits one lemma `theorem mathilda_cert_1_pratt_i : Mathilda.PrattCertifies p cert` per odd prime factor. The theorem is `n = Π p^k ∧ (∃ c, Mathilda.PrattCertifies p c) ∧ …`.
- `Certify[PrimeQ[n]]` proves `Mathilda.PrattCertifies n cert` for an odd prime `n`. A non-prime is stated with the prelude's definition `Mathilda.IsPrime n := 2 ≤ n ∧ ∀ d, d ∣ n → d = 1 ∨ d = n`: a composite `n` gets `¬ Mathilda.IsPrime n` from a factor `a` (`h.2 a ⟨n/a, _⟩` gives `a = 1 ∨ a = n`, refuted by `decide +kernel`), and `n < 2` gets `¬ Mathilda.IsPrime n` from `2 ≤ n`. `PrimeQ[2]` proves `Mathilda.IsPrime 2` outright. These have `"TrustedSteps" -> {}`.
- `"TrustedSteps"` of a Pratt certificate: Lucas's theorem. Lean checks the certificate, not the theorem that a valid certificate implies primality. Proving that through Mathlib is a follow-up.

```mathematica
In[1]:= c = Certify[FactorInteger[2^64 + 1]];
In[2]:= {c["Verified"], c["LeanTheorem"]}
Out[2]= {True, "18446744073709551617 = 274177 * 67280421310721 ∧ (∃ c, Mathilda.PrattCertifies 274177 c) ∧ (∃ c, Mathilda.PrattCertifies 67280421310721 c)"}
```

### Closed-form sums
`Certify[Sum[t, {k, a, n}]]` (or `Sum[t, {k, n}]`, which is read as `Sum[t, {k, 1, n}]`, as `Sum` defines it) takes a non-negative integer `a`, an upper bound `n` that is a variable (not `Infinity`, `Pi` or any other built-in constant: those have no certificate) and a term `t` that is a polynomial in `k`, or a polynomial times `c^k` for a rational `c`. A term with a denominator in `k` (or a closed form with a denominator in `n`) has no certificate: the denominator could vanish inside the range, where Lean's `x / 0 = 0` would give the sum a junk value. It certifies Mathilda's closed form. For a claim `Sum[t, {k, a, n}] == closed`, it certifies the claimed form instead. Any other iterator (`{k, 1, n + 1}`, `{k, 1, 2 n}`, `{k, n, 2 n}`, ...) has no certificate (`Certify::nocert`, and `Conjecture` reports `"Tested"`); before v0.6528 such a claim fell to the identity generator, which evaluated the `Sum` and certified Mathilda's own output against itself. Lean defines

```lean
def mathilda_cert_1_t (k1 : Nat) : Rat := …
def mathilda_cert_1_S : Nat → Rat
  | 0 => mathilda_cert_1_t a
  | m1 + 1 => mathilda_cert_1_S m1 + mathilda_cert_1_t (m1 + 1 + a)
```

and proves `∀ m1 : Nat, mathilda_cert_1_S m1 = closed(m1 + a)` by induction (`simp only`, `push_cast`, `grind`). Every `c^(α n + β)` becomes `(c^α)^m1 · c^(αa + β)`, so the exponent is always the `Nat` binder. Every other occurrence of the variable becomes `(m1 : Rat)`. `"Theorem"` is `ForAll[n, n >= a, Sum[t, {k, a, n}] == closed]`. `"TrustedSteps"` names the one gap between that and Lean's statement (v0.6522): "the closed form is proved for integer n >= a only (Lean's statement is over n = a + m with m : Nat)"; `n` is not checked for non-integer or smaller values.

```mathematica
In[1]:= Certify[Sum[k 2^k, {k, 0, n}]]["LeanTheorem"]
Out[1]= "∀ m1 : Nat, mathilda_cert_1_S m1 = 2 + 2 ^ m1 * 2 * ((m1 : Rat) - 1)"

In[2]:= Quiet[Certify[Sum[k^3, {k, 1, n}] == n^2 (n + 1)^2/4 + 1]["Verified"]]
Out[2]= False
```

**Fibonacci** (v0.6527). The term may contain `Fibonacci[k + j]` (an integer shift `j >= 0`) and the closed form `Fibonacci[n + j]` (`a + j >= 0`), e.g. `Sum[Fibonacci[k]^2, {k, 1, n}] == Fibonacci[n] Fibonacci[n + 1]`. The source then also defines

```lean
def mathilda_cert_1_F : Nat → Rat
  | 0 => 0
  | 1 => 1
  | n + 2 => mathilda_cert_1_F n + mathilda_cert_1_F (n + 1)
```

before `_t`, `Fibonacci[n + j]` becomes `mathilda_cert_1_F (m1 + a + j)` and both `simp only` lists gain `mathilda_cert_1_F`. `"TrustedSteps"` gains "the Lean function _F of the source (F 0 = 0, F 1 = 1, F (n + 2) = F n + F (n + 1)) is Mathilda's Fibonacci on the naturals (checked by Mathilda at n = 0..29, not by Lean)". Sums without `Fibonacci` are unchanged. Other arguments (`Fibonacci[2 k]`, a negative shift in the term) have no certificate; `Fibonacci` outside a sum claim is `Certify::eval` as before. The Exact engine checks the step with `Fibonacci[n + j] = F(j - 1) f0 + F(j) f1` (valid for every integer `j`), a polynomial identity in fresh `f0 = F(n)`, `f1 = F(n + 1)`.

**Options** (`Options[Certify]`; `SetOptions` works). Option names are matched by their text, as `OptionValue` does: `Engine -> "Mathilda"` (the user's own symbol `Engine`) and `"Engine" -> "Mathilda"` are the same option. `Options[Certify]` lists the package's names as strings (`"Engine"`, `"Tier"`, `"Domain"`, `"Verify"`, `"Counterexample"`, `"LeanExecutable"`, `"LeanProject"`) next to the System symbols `Method`, `TimeConstraint` and `Assumptions`. The package creates no option-name symbols (v0.6523; before, `Engine`, `Tier`, `Domain`, `Verify`, … lived in `` Lean` `` on `$ContextPath`, captured the user's symbols of those names and survived `ClearAll["Global`*"]`). An unknown name is `Certify::optx`. A name symbol with a value (`Engine = "zzz"`) would turn `Engine -> "Lean"` into `"zzz" -> "Lean"`, so `Certify` checks the held options first and issues `Certify::optvalue`; use `Clear[Engine]` or the string form.

| Option | Default | Meaning |
|---|---|---|
| `"Engine"` | `"Lean"` | `"Lean"`; `"Exact"`: Mathilda's exact re-check of the certificate, no Lean needed (`"Tier" -> "Exact"`, `"Verified" -> Missing["Exact", True \| False]` or `Missing["NotChecked"]`; see [Engines and consensus](#engines-and-consensus)); `"Mathilda"`: a heuristic check through `PossibleZeroQ` / `Simplify` (`"Tier" -> "Heuristic"`, `"Verified" -> Missing["Heuristic", True \| False]`; the summary box reads `Mathilda Heuristic`); or a non-empty list of these without duplicates (v0.6525), which runs every listed engine and combines the verdicts. Neither an exact re-check nor a heuristic answer is a verification: only a kernel (Lean) makes `"Verified"` `True`. Anything else issues `Certify::engine` and returns `$Failed`. `Method` and `"Tier"` are read by `"Lean"` (and by `"Mathilda"`, which shows the Lean encoding); given explicitly to engines that do not read them they issue `Certify::engineopt` (`"LeanExecutable"`, `"LeanProject"` and `TimeConstraint` are ignored silently). |
| `"Tier"` | `Automatic` | `0` (core Lean), `1` (Mathlib; needs `$LeanProject`, otherwise `Certify::nomathlib`) or `Automatic`. |
| `"Domain"` | `Automatic` | `Automatic`, `Complexes`, `Reals`, `Rationals` or `Integers` (table above; Tier 0). Anything else is `Certify::domain`. |
| `Method` | `Automatic` | The tactic: `"grind"`, `"decide"`, `"omega"`, or a tactic script string. It goes into the Lean source file only, never into a shell command. A script is refused (`Certify::method`, `$Failed`) when it contains a tab, a control character, `«`, `»`, a backquote, `$`, `#`, `@[`, `--`, `/-` or `"`, or uses a command keyword or meta-level escape as a word (`theorem`, `def`, `example`, `macro`, `syntax`, `elab`, `set_option`, `attribute`, `open`, `namespace`, `section`, `end`, `axiom`, `run_tac`, `native_decide`, `sorry`, `Lean`, `IO`, …). Every line is indented under `by`, and each item with a script runs in a Lean file of its own. |
| `TimeConstraint` | `60` | Seconds, mapped to `set_option maxHeartbeats (40000 t) in` on the theorem. This is a deterministic resource bound, not a wall-clock kill. |
| `Assumptions` | `{}` | Extra hypotheses (walked held like the claim; a `SetOptions` default other than `{}` is refused, `Certify::optdefault`): `Assumptions -> {x == 1}` gives `∀ {K : Type} [Lean.Grind.Field K] [Lean.Grind.IsCharP K 0], ∀ x : K, x = 1 → …` (hypotheses make the automatic domain every char-0 field, or every ordered field when every variable is real; variables that occur only in the assumptions are bound as well). `Element[x, Reals]` / `Element[{x, y}, Reals]` declare variables real (they are not translated as hypotheses); an order relation with a variable that is not real gives `Certify::realvars`. |
| `"Verify"` | `True` | `False` generates the source without running Lean (`"Verified" -> Missing["NotChecked"]`). |
| `"LeanExecutable"` | `Automatic` | Overrides `$LeanExecutable` for one call. |
| `"LeanProject"` | `Automatic` | Overrides `$LeanProject` for one call (`None` forces no project). |
| `"Counterexample"` | `Automatic` | What happens when a kernel engine **rejects** an item (v0.6526): `Automatic` searches a counterexample for at most 2 s and, if one is found, refutes the claim (see [Counterexamples](#counterexamples)); a positive number is the search budget in seconds; `False` turns the search off (the object then has the usual keys only). `Verify -> False`, a missing Lean and the non-kernel engines never trigger it. |

`Domain`, `Method` and `Assumptions` apply to identity claims (`a == b` and the `"Identity"` / `"Factorization"` computations). The other kinds choose their own proof script and number type (and issue `Certify::optignored` when one of these options is given), and take their hypotheses from the claim itself (`Implies[...]`, `ForAll[vars, cond, ...]`); a claim with hypotheses is stated over every field of characteristic 0 (equations) or every ordered field (inequalities, `Minimize` / `Maximize`), never over `Int` or `Rat` alone. (`Solve` soundness statements are closed rational facts and stay over `Rat`.)

**Messages**: `Certify::nocert` (no certificate for this input), `Certify::eval` (the input applies a computation Mathilda would have to evaluate before Lean sees the claim, `Evaluate` included), `Certify::value` (a symbol in the claim or its `Assumptions` has a value of some kind, or an evaluation-changing attribute), `Certify::context` (a symbol of `` Lean`Private` `` / `` Lean`Internal` ``), `Certify::domain` (an unsupported `Domain`), `Certify::element` (an `Element` assumption the domain cannot state), `Certify::optdefault` (an `Assumptions` default set with `SetOptions`), `Certify::antideriv` (a claim `Integrate[f, x] == F` is certified as `D[F, x] == f`), `Certify::fracpow` (an integrand with a fractional power of the variable: outside ℚ(x)), `Certify::vacuous` (a Tier 1 hypothesis with no witness point), `Certify::optx` / `Certify::optvalue` (an unknown option name, or a name symbol with a value), `Certify::realvars` (order relations, but a variable not declared real), `Certify::optignored` (`Domain` / `Assumptions` / `Method` given to a kind that ignores them), `Certify::method`, `Certify::engine`, `Certify::engineopt` (`Method` or `"Tier"` given to engines that do not read them), `Certify::exactfail` (the `"Exact"` re-check failed), `Certify::disagree` (the listed engines disagree), `Certify::tier`, `Certify::nomathlib`, `Certify::unsafe` (the generated file failed the source check -- an internal assertion; `"Verified" -> False`), and `Certify::rejected`, which is issued when Lean judged the certificate and rejected it.

### Engines and consensus
`Certify` runs its certificates through registered **backends** (engines). `"Lean"` is the default and the only one whose verdict can be `True`: it is a kernel (`"Trust" -> "Kernel"`). Since v0.6525 there is a second, always-available engine, `"Exact"`, and `"Engine"` may be a list.

**`Engine -> "Exact"`** re-checks the certificate with Mathilda's own exact rational arithmetic and polynomial `Expand`, without Lean. It is an exact re-check by the system that produced the certificate, **not** a kernel proof, so it never reports `True`:

| `"Verified"` | Meaning |
|---|---|
| `Missing["Exact", True]` | every check passed |
| `Missing["Exact", False]` | a check failed on a concrete exact computation; `Certify::exactfail` |
| `Missing["NotChecked"]` | no decisive exact check for this shape (an identity with an inequality or `Unequal` hypothesis, a check that timed out), or `Verify -> False` |

The checks are recomputed from the certified statement (the obligation behind `"Theorem"`) plus the certificate data, never from the generator's own `"PreCheck"`; `obj["Checks"]` lists each check's text with its result (`True`, `False`, `Missing["NotChecked"]`). Per kind: identities -- `Expand[lhs - rhs] === 0` (with equality hypotheses: cofactors with `Expand[lhs - rhs - c . p] === 0`; cleared denominators must not vanish identically); Bezout, `PolynomialReduce` and Tier 0 antiderivatives -- each equation of the statement expands to 0 (surd variables included); ideal membership -- the cofactor identity; `Solve` completeness -- `lhs - rhs` equals the factorization (a nonzero constant times linear factors) and its roots are the stated alternatives; soundness and the matrix kinds -- the closed statement evaluates to `True` with the prelude's `matMul`, `idMat`, `matVec`, `HasShape`, `HasLength` read as `Dot`, `IdentityMatrix`, `Length`, `Dimensions` and `det` as a cofactor expansion independent of `Det`; sums of squares -- the target recomputed from the statement equals `Σ d l² g` with every `d > 0` and every `g` 1 or a hypothesis polynomial; optimality -- feasibility and value at the optimum, multipliers `≥ 0`, each eliminated variable forced by an equality constraint, and the KKT identity; primality -- the Pratt certificates replayed (`PowerMod` Lucas checks, the factorization of `q - 1`, every odd prime factor certified earlier), the factorization's product, a divisor of a composite, `n < 2`; sums -- the base and the step `C(n + 1) - C(n) = t(n + 1)` as an exponential-polynomial identity (powers `c^(α n + β)` rewritten over prime bases), and when that fails the claim evaluated at `n = a, …, a + 10` (a mismatch is `False`, otherwise `Missing["NotChecked"]`); the Tier 1 shape -- `Together[D[F, x] - f] === 0` (no Mathlib needed). Each check runs under `TimeConstrained[…, Min[TimeConstraint, 10]]`. The object's properties are `Cached`, `Certificate`, `Checks`, `Domain`, `Engine`, `Hypotheses`, `Input`, `Kind`, `Result`, `Theorem`, `Tier`, `Timing`, `TrustedSteps`, `Verified`; `"TrustedSteps"` gains "Mathilda's exact rational arithmetic and polynomial Expand: an exact re-check by the system that produced the certificate, not a kernel proof" and, per kind, Lucas's theorem, the induction principle, Mathilda's `D` / `Together`, no zero divisors, or the meanings of the matrix operations. An `"Exact"`-only call issues no Lean message (no Lean encoding is made).

**A list of engines** runs every listed engine (in list order; Lean still runs when Exact fails, so both verdicts are visible) and combines them. With `K` the kernel engines (heuristic engines only add their entry):
1. a `False` or `Missing[_, False]`: `False`; `Certify::disagree` when a kernel says `True`, or a kernel says `False` while another engine says `Missing[_, True]` (a checker or a translation is wrong, or the kernel ran out of time);
2. `K` non-empty and every kernel `True`: `True`;
3. no kernel listed: `Missing["Consensus", True]` when every verdict is `Missing[_, True]`, else `Missing["NotChecked"]`;
4. otherwise the first kernel verdict that is not `True` (`Missing["NotAvailable"]` without Lean, `Missing["NotChecked"]` with `Verify -> False`).

The object has `"Engine" -> {…}` (the list), `"Engines" -> <|"Lean" -> True, "Exact" -> Missing["Exact", True]|>` and `"EngineDetails" -> <|engine -> <|"Verified", "Timing", "TrustedSteps", "Logs", …|>|>` (Exact's `"Checks"`, Lean's `"Axioms"` and `"LeanVersion"`), on top of the Lean object's properties when `"Lean"` is listed (21 keys), else of the Exact object's. `"Tier"` is Lean's tier when Lean is listed, else `"Exact"`; `"Timing"` is the sum; `"TrustedSteps"` is the certificate's plus the backend-level steps of the kernel engines (of every engine when none is listed). The summary box joins the engine labels. A one-element list gives this shape too.

```mathematica
In[1]:= ToString[Certify[Factor[x^10 - 1], Engine -> {"Lean", "Exact"}]]
Out[1]= "ProofObject[«Factorization» Verified: True | Lean Tier 0 + Exact | 2.41 s]"

In[2]:= Quiet[Certify[(1 + x)^20 == Expand[(1 + x)^20] - 10^-25 x, Engine -> {"Lean", "Exact"}]]["Engines"]
Out[2]= <|"Lean" -> False, "Exact" -> Missing["Exact", False]|>

In[3]:= Certify[FactorInteger[2^64 + 1], Engine -> "Exact"]["Verified"]
Out[3]= Missing["Exact", True]
```

### Counterexamples
Since v0.6526 a kernel **rejection** is followed by a bounded counterexample search (`"Counterexample" -> Automatic`, 2 s per item; `False` turns it off). When a witness is found and Mathilda's exact arithmetic confirms it, the claim is refuted: Lean's kernel checks the theorem `¬ S`, where `S` is the certified statement (what Certify tried to prove, generic binders included), proved from the witness by `fun h => absurd (h w…) (by decide +kernel)`. All refutations of one call run in one Lean file. The object keeps `"Verified" -> False` and the rejected `"LeanSource"`, and gains:

| Key | Content |
|---|---|
| `"Counterexample"` | The witness, `{x -> -1, …}` (`{}` for a closed statement such as a `Solve` soundness fact) |
| `"Refuted"` | `True` (the kernel accepted `¬ S`) or `Missing["Exact", True]` (the witness is checked by exact arithmetic only: an algebraic witness, a claimed antiderivative, or a refutation run that did not succeed) |
| `"RefutationTheorem"` | `HoldForm[!S]` |
| `"RefutationSource"` | The Lean file of the refutation (only when it ran) |

It issues `Certify::refuted` and displays as `ProofObject[«Identity» Refuted: {x -> -1} | Lean Tier 0 | …]`. The search and the faithfulness rules are those of [Refute](#refute). Nothing is searched for primality, Tier 1 (Mathlib) or `Minimize` / `Maximize` certificates (a rejection there is a certificate or translation bug, not a false claim).

```mathematica
In[1]:= ToString[Quiet[Certify[(x + 1)^3 == x^3 + 3 x^2 + 2 x + 1]]]
Out[1]= "ProofObject[«Identity» Refuted: {x -> -1} | Lean Tier 0 | 0.63 s]"

In[2]:= %%["RefutationTheorem"]   (* of the object *)
Out[2]= !ForAll[{x}, (1 + x)^3 == 1 + 2 x + 3 x^2 + x^3]
```

## ProofObject
- `ProofObject[<|...|>]`: the result of `Certify`.
- `ProofObject[comp, opts]`, with `comp` not an `Association`: the same as `Certify[comp, opts]`, so `ProofObject[Factor[x^2 - 1], Engine -> "Lean"]` works.
- `obj["prop"]`: a property. `obj["Properties"]` gives the sorted list, and an unknown name gives `Missing["UnknownProperty", name]`.
- Attributes `HoldAll`, `Protected` (and the same for `Certify`; `HoldAll` since v0.6521, so that the `Assumptions` option is checked before it evaluates).
- It displays as `ProofObject[«Factorization» Verified: True | Lean Tier 0 | 0.51 s]`; a refuted object (a kernel rejection with a counterexample, v0.6526) as `ProofObject[«Identity» Refuted: {x -> -1} | Lean Tier 0 | 0.63 s]`.
- **Authentication** (v0.6522). `ProofObject[assoc]` is an inert container, so anyone can write `ProofObject[<|"Verified" -> True, …|>]`. Every association `Certify` produces is registered by its `Hash` for the session, and `obj["Verified"]` (and the summary box) reports the stored verdict only for a registered object; a hand-made or edited one gives `Missing["Unauthenticated"]` (its other properties are still readable). Since v0.6528 the same holds for every verdict-bearing key: `"Status"` (Conjecture), `"Refuted"`, `"Counterexample"` and `"RefutationTheorem"`. This guards against minted objects in ordinary code, not against code that deliberately writes `` Lean`Private` `` internals, and an object from another session (saved and read back) is unauthenticated too: re-run `Certify`, or check `obj["LeanSource"]` with `LeanVerify`.

| Property | Content |
|---|---|
| `"Kind"` | `"Identity"`, `"Factorization"`, `"Antiderivative"`, `"Bezout"`, `"IdealMembership"`, `"Roots"`, `"MatrixInverse"`, `"Determinant"`, `"LinearSolve"`, `"NullSpace"`, `"Nonnegativity"`, `"Optimality"`, `"Primality"` or `"ClosedFormSum"`; `"Refutation"` for [Refute](#refute) and `"Conjecture"` for [Conjecture](#conjecture) (their own property lists are given there) |
| `"Input"` | `HoldForm` of the computation or the claim |
| `"Result"` | The computation's value; for a claim, the claim with evaluated sides |
| `"Theorem"` | The certified statement as a held Mathilda expression (`ForAll[{x}, …]`); for a claim `Integrate[f, x] == F`, `ForAll[x, D[F, x] == f]` |
| `"LeanTheorem"` | The statement as Lean text |
| `"ProofScript"` | The tactic block, e.g. `by intros; grind` |
| `"LeanSource"` | The complete Lean file: header, prelude, any auxiliary definitions or lemmas (`mathilda_cert_1_S`, `mathilda_cert_1_pratt_1`, …), `theorem mathilda_cert_1` (`mathilda_cert_k` in a batch), `#print axioms` |
| `"Certificate"` | Raw Mathilda data, per kind: e.g. `"Form"` and `"Variables"` for identities, `"Cofactors"` for Bezout and ideal membership, `"Terms"`, `"Surds"` and `"PreCheck"` for antiderivatives, `"Roots"` and `"Completeness"` for `Solve` |
| `"Verified"` | `True`, `False`, `Missing["NotChecked"]` (`Verify -> False`, or no decisive exact check), `Missing["NotAvailable"]` (no Lean found), `Missing["Exact", True \| False]` (`Engine -> "Exact"`), `Missing["Heuristic", True \| False]` (`Engine -> "Mathilda"`), `Missing["Consensus", True]` (a list of engines with no kernel) or `Missing["Unauthenticated"]` (not produced by `Certify` in this session) |
| `"Engine"`, `"Tier"`, `"Domain"` | As the options (`"Engine"` is the list for a list of engines); `"Tier"` is `"Heuristic"` for `Engine -> "Mathilda"` and `"Exact"` for `Engine -> "Exact"` |
| `"Checks"` | `Engine -> "Exact"` only: each exact check's text `->` its result |
| `"Engines"`, `"EngineDetails"` | A list of engines only: each engine's verdict, and its verdict, timing, backend-level trusted steps, logs and backend keys |
| `"Counterexample"`, `"Refuted"`, `"RefutationTheorem"`, `"RefutationSource"` | A kernel rejection with a counterexample only (see [Counterexamples](#counterexamples)); `"RefutationSource"` only when the Lean refutation ran |
| `"Hypotheses"` | Side conditions such as denominators `≠ 0` |
| `"TrustedSteps"` | What Lean did not check |
| `"Axioms"` | The theorem's axioms, from `#print axioms` |
| `"Logs"` | Lean's messages for the theorem: `<|"Severity", "Line", "Text"|>` |
| `"Timing"` | Lean's wall time in seconds |
| `"Cached"` | Whether the verdict came from the session cache |
| `"LeanVersion"` | `lean --version` of the checker that ran |

**The `"Verified" -> True` contract.** The generated file is: a header (`set_option autoImplicit false`), the prelude, the theorem chunks (`theorem mathilda_cert_k`, with any auxiliary declarations), and then a **trusted footer** after all chunks with, for each theorem, `example : <statement> := mathilda_cert_k` (the statement is re-elaborated and must be the theorem's type, so a chunk cannot pass by declaring something else under that name) and `#print axioms mathilda_cert_k`. A theorem is `True` only when all of these hold:
1. no message of severity `error` and no `hasSorry` message lies in its chunk or in its footer lines;
2. no such message lies outside every chunk and footer (header, prelude);
3. the `#print axioms` report printed **on its own footer line** (a report anywhere else is ignored) names only `propext`, `Classical.choice` and `Quot.sound` — `sorryAx` and the auxiliary axiom `native_decide` adds fail this;
4. Lean exited with code 0, or every error it reported lies in some theorem's lines (a non-zero exit is then explained by sibling theorems of the batch, which are `False`).

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
Out[5]= {Missing["Heuristic", True], False}
```

## Refute
`Refute[claim, opts]` (v0.6526) searches a counterexample to `claim` and proves the **negation** of its certified statement `S` from the witness. It gives a `ProofObject` of Kind `"Refutation"`, or `$Failed`.
- Claims: an identity `a == b` (with `"Domain"` and `Assumptions` as in `Certify`), an implication of equations `Implies[p1 == q1 && …, p == q]`, a `ForAll[vars, (cond,) p >= q]` or bare inequality, `Sum[t, {k, a, n}] == C`, and `Integrate[f, x] == F` (exact only). Anything else (`PrimeQ[7]`, a computation such as `Integrate[x, x]` or `Factor[…]`) is `Refute::nostmt`. The claim and the `Assumptions` are walked held exactly as in `Certify` (`Certify::eval`, `Certify::value`).
- `S` is the statement `Certify` would prove (`"RefutedStatement"`), `"Theorem"` is `HoldForm[!S]`, and `"LeanTheorem"` is `¬ (S)` in Lean text with `"ProofScript"` `fun h => absurd (h w…) (by decide +kernel)`: `(K := Rat)` first for the generic char-0 / ordered-field statements, then the witness values `(-1 : Int)` / `(1/2 : Rat)`, then one `(by decide +kernel)` per antecedent of `S`; a sum's statement `∀ m : Nat, S m = C (m + a)` is refuted by `h m₀` with the preamble's definitions; a closed statement by `decide +kernel`.
- **Faithfulness.** A witness must lie in `S`'s quantifier domain: integers for `Int`, rationals for `Rat` and for the generic forms (ℚ is a field of characteristic 0 and an ordered field, so a rational witness refutes "for every such field"; since ℚ ⊂ ℂ and ℚ ⊂ ℝ, the claim is refuted under Mathematica's semantics too). Every hypothesis of `S` holds there, its conclusion fails, and every denominator of the claim is nonzero there: a point where a rational claim is undefined is no counterexample, even when it falsifies the cleared statement (`Refute[1/x == 1]` gives `x -> -1`, not the pole `x -> 0`). For a claim whose statement omits the call's `"Domain"` / `Assumptions` (every kind but identities: a sum, an inequality, an implication), the witness must also satisfy them exactly -- `S` is then stronger than the claim, so a counterexample to `S` alone would not refute the claim as stated (`Refute[x^2 >= x, Assumptions -> Element[x, Integers]]` is `Refute::none`; with `Assumptions -> x > 0` the witness is `x -> 1/2`, and the condition is listed in `"TrustedSteps"`). The automatic search after a `Certify` rejection applies the same rule. Mathilda checks all of it by exact arithmetic before a refutation is built. `FindInstance`'s domain argument is not trusted: it can return `x -> -2 I` for `Rationals`.
- **Search** (deterministic): sums at `n = a, …, a + 30`; an exact grid by increasing height (integers `−5..5`, or `0, ±1, ±2, ±3, ±1/2, ±3/2, ±1/3, ±2/3` for the rational forms; at most 4000 points; points where a denominator vanishes are skipped); `FindInstance` over `Integers` / `Rationals`; then, exact only, `FindInstance` over ℂ (char-0 form) or the reals (ordered form). An algebraic witness (`-2 I`, `-Sqrt[2]`) has no kernel proof in core Lean: the result is the Exact object with `"Verified" -> Missing["Exact", True]` (nonzero tests by `RootReduce` and a 50-digit value, listed in `"TrustedSteps"`). A claimed antiderivative is refuted where `D[F, x]` and `f` are different rational numbers (never at a surd or transcendental value). No witness within the bounds: `Refute::none` and `$Failed` — not a proof of the claim.
- Options: `"Engine"` (`"Lean"`, default, or `"Exact"`: the witness re-checked exactly, `Missing["Exact", True]`; anything else `Refute::engine`), `"Domain"`, `Assumptions`, `TimeConstraint` (60), `"Verify"` (`False`: the source is built, nothing runs, `Missing["NotChecked"]`), `"SearchTime"` (5 s), `"LeanExecutable"`, `"LeanProject"`; an unknown name is `Refute::optx`. Attributes `HoldAll`, `Protected`.
- Properties (Lean): `Axioms, Cached, Certificate, Counterexample, Domain, Engine, Hypotheses, Input, Kind, LeanSource, LeanTheorem, LeanVersion, Logs, ProofScript, RefutedStatement, Result, Theorem, Tier, Timing, TrustedSteps, Verified`; the Exact variant drops the Lean keys and adds `"Checks"`. `"Certificate"` holds the witness, the form, the refuted kind, the domain form and the antecedent count.

```mathematica
In[1]:= r = Refute[(x + 1)^3 == x^3 + 3 x^2 + 2 x + 1];
In[2]:= {ToString[r], r["Counterexample"], r["LeanTheorem"]}
Out[2]= {"ProofObject[«Refutation» Verified: True | Lean Tier 0 | 2.05 s]", {x -> -1},
         "¬ (∀ x : Int, (1 + x) ^ 3 = 1 + 2 * x + 3 * x ^ 2 + x ^ 3)"}

In[3]:= Refute[ForAll[{x, y}, x^2 + y^2 >= 2 x y + 1/10]]["ProofScript"]
Out[3]= "fun h => absurd (h (K := Rat) (0 : Rat) (0 : Rat)) (by decide +kernel)"

In[4]:= Refute[Implies[x^2 + 4 == 0, x == 7]]["Verified"]   (* x -> -2 I: exact only *)
Out[4]= Missing["Exact", True]

In[5]:= Refute[ForAll[{x, y}, x^4 y^2 + x^2 y^4 - 3 x^2 y^2 + 1 >= 0]]
Refute::none: no counterexample to ... was found within the search bounds; this is not a proof of the claim
Out[5]= $Failed
```

## Conjecture
`Conjecture[expr, n, opts]` (v0.6527) runs an experiment, guesses a closed form in `n`, tests it and certifies it. It gives a `ProofObject` of Kind `"Conjecture"`, or `$Failed`. Attributes `HoldAll`, `Protected`.
1. **Experiment.** `expr` (held) is evaluated exactly at `n = s, s + 1, …` for `"FitPoints"` (12) plus `"TestPoints"` (40) values. `"Start" -> Automatic` is the lower limit `a` of an outer `Sum[t, {k, a, n}]` (1 for `{k, n}`), else 0. Every value must be an integer or a rational, else `Conjecture::nonexact`.
2. **Guess**, from the fit values, by the first form (in this order, `"Forms"` may restrict them) that reproduces every fit value exactly and then every test value: `"Polynomial"` (`InterpolatingPolynomial`, degree ≤ 10, `Factor`ed); `"Geometric"` (the smallest linear recurrence with constant coefficients, order ≤ 6, from the exact `NullSpace` of the Hankel matrix; a characteristic polynomial with nonzero rational roots `r` of multiplicity `μ` gives the ansatz `Σ p_r(n) r^n`, deg `p_r` < `μ`, fitted exactly); `"Fibonacci"` (only when `expr` contains `Fibonacci`: the smallest subset of `{1, n, (-1)^n, F(n), F(n+1), F(n)², F(n) F(n+1), F(n+1)²}` that fits); `"Rational"` (`p(n)/q(n)`, degrees ≤ 4, from an exact `NullSpace`); `"Hypergeometric"` (`Binomial[2 n, n]`, `Binomial[2 n, n]/(n + 1)`, `n!`, `(2 n)!` times `c^n`, else a term ratio that is a rational function with rational roots, written with `Pochhammer`). None: `Conjecture::noguess`, `$Failed`.
3. **Certify** the claim `expr == guess` as written (the user's iterator names are kept) with `Certify`'s pipeline and the options `"Engine"` (`{"Lean", "Exact"}`), `"Verify"`, `TimeConstraint`, `"LeanExecutable"`, `"LeanProject"`; every message of that run is quieted.

`"Status"`: `"Proved"` (a kernel accepted: `"Verified" -> True`), `"Checked"` (no kernel verdict, e.g. Lean absent or `"Engine" -> "Exact"`, and the exact re-check passed), `"Refuted"` (the certification found a counterexample), `"Tested"` (anything else: no certificate kind, `"Verify" -> False`). A guess is only ever reported as tested on the points in `"TestedOn"` unless it is certified. `"Verified"` is the certification's verdict, or `Missing["NoCertificate", "Certify::eval"]` (the name of the message that refused it) when there is no certificate; `"Proof"` is the nested `ProofObject` or `Missing["NoCertificate"]`; `"Theorem"` is the proof's, else `HoldForm[ForAll[n, n >= s, expr == guess]]`.

Properties: `Engine, FittedOn, Form, Guess, Input, Kind, Proof, Status, TestedOn, Theorem, Timing, Variable, Verified`. Other messages: `Conjecture::var` (`n` not a symbol, or a symbol with a value), `Conjecture::badopt` (an invalid option value), `Conjecture::optx` / `Conjecture::optvalue` (option names).

```mathematica
In[1]:= c = Conjecture[Sum[k^4, {k, 1, n}], n];
In[2]:= ToString[c]
Out[2]= "ProofObject[«Conjecture» Proved: Sum[k^4, {k, 1, n}] == 1/30 n (1 + n) (1 + 2 n) (-1 + 3 n + 3 n^2) | Lean Tier 0 + Exact | 2.27 s]"

In[3]:= {c["Form"], c["FittedOn"], c["TestedOn"], c["Proof"]["Engines"]}
Out[3]= {"Polynomial", {1, 12}, {1, 52}, <|"Lean" -> True, "Exact" -> Missing["Exact", True]|>}

In[4]:= Conjecture[Sum[Fibonacci[k]^2, {k, 1, n}], n]["Guess"]
Out[4]= Fibonacci[n] Fibonacci[1 + n]

In[5]:= c = Conjecture[Sum[Binomial[n, k]^2, {k, 0, n}], n]; {c["Guess"], c["Status"], c["Verified"]}
Out[5]= {Binomial[2 n, n], "Tested", Missing["NoCertificate", "Certify::eval"]}
```

## LeanVerify
- `LeanVerify[src]`: checks the Lean 4 source string `src` and gives `<|"Success", "ExitCode", "Messages", "Axioms", "AxiomReports", "Timing", "Command", "Replay", "ReplayAxioms", "ReplayTiming"|>`.

**Features**:
- `Protected`. Options `"LeanExecutable"`, `"LeanProject"` and `"KeepFiles"` (`True` keeps the temporary `.lean` file and returns its path as `"SourceFile"`), matched by name like `Certify`'s (`KeepFiles -> True` works too); an unknown name is `LeanVerify::optx`.
- **Accepted sources: a minimal surface, then an allowlist grammar** (v0.6521). LeanVerify's guarantee holds only when the checked source cannot run user-controlled code while Lean elaborates it -- code run at elaboration time can forge everything the elaborating process writes, including the `.olean` the replay reads. So, before Lean runs, `` Lean`Internal`CheckSource `` (`src/lean_check.c`) checks the source, and anything else is **not run** (`"Success" -> False`, `"ExitCode" -> Missing["NotRun"]`, `"Rejected" -> {what}`). The check does **not** try to reproduce Lean's lexer (the third review's lexer disagreed with Lean in three ways, each of which hid a column-0 `#eval` in what the checker took for a string, a character literal or a comment). Instead the accepted surface has no construct whose content Lean would skip, so every byte belongs to a token the checker sees:
  - **(A) Lexical surface** (`LeanVerify::lex`, with the construct and its line). Refused anywhere: `"` (no strings of any kind: plain, raw, interpolated, `s!` / `m!` / `f!`); `/-`, `-/` and `--` (no comments, no doc comments); `«` and `»`; the backquote; `#` (no hash commands at all -- not even `#print axioms`: the replay computes the axioms); `@[`; `$`; `%`; `\`; tabs, carriage returns and every other control character; a quote `'` unless it is an identifier prime (it follows an ASCII letter, digit, `_` or `'`, and the run of `[A-Za-z0-9_']` it ends starts with a letter or `_`; so `h'`, `h''` and `x1'` are fine, `1'`, `x!'`, `↦'` and every character literal are not); invalid or overlong UTF-8; and every non-ASCII character except `∀ ∃ → ← ↔ ∧ ∨ ¬ ≤ ≥ ≠ ⟨ ⟩ ↦ λ · ∘ ∈ ∉ ∣ × ⁻ ¹ ² ³ ℕ ℤ ℚ ℝ ℂ ⊢ ≡ ∑ ∏ ∞ ⌊ ⌋ ⌈ ⌉ ‖` and the Greek letters `α`-`ω`, `Α`-`Ω`.
  - **(B) Words** (`LeanVerify::unsafe`). A word is a maximal run of `[A-Za-z0-9_]`; every other character separates words. Lean's identifiers are coarser (they also absorb `'`, `!`, `?`, `.` and Unicode letters), so any keyword Lean sees starts a word here -- with two exceptions that are handled: after a numeric literal (`2theorem` is `2` then `theorem` for Lean: in a run that starts with a digit every suffix that starts with a letter is checked), and after a Mathlib notation token that ends in a letter (`→o`, `[MOD`, `≤i`, …: the word after it is checked). Refused anywhere: every command-leading word of Lean's parser -- the list in `src/lean_check_tokens.h` is **derived from Lean's parser tables** (`tools/lean_command_tokens.lean`) and tested against them for core Lean and Mathlib (ctest `lean_check_tokens_sync`, `lean_check_tokens_mathlib`) -- except an allowed command starter in its starting position; and the escape words `by_elab`, `run_tac`, `run_elab`, `run_meta`, `run_cmd`, `native_decide`, `native` (`decide +native`), `implemented_by`, `extern`, `unsafe` (and `unsafeBaseIO`, `unsafeIO`, …), `sorry`, `admit`, `ofReduceBool`, `ofReduceNat`, `trustCompiler`, `attribute`, `debug`, `elab`, `macro`, `eval`, `skipKernelTC`, `addDecl`, `IO`, `EIO`, `BaseIO`, `Elab`, `MetaM`, `TacticM`, `CoreM`, …, `plausible`, `slim_check`, `polyrith`, `bv_decide`. (So, for example, a variable named `end`, `open` or `local` is refused; with Mathlib so are `Polynomial.eval` and `alias`.)
  - **(C) Grammar.** **Imports**: none (core `Init` only), or, with a project (`LeanProject` / `MATHILDA_LEAN_PROJECT`), only `import Mathlib` and `import Mathlib.<Path>`, one module per line, before every other command; anything else (`import Lean`, `Lake`, `Mathilda`, `Qq`, `Aesop`, `Batteries`, …) is `LeanVerify::import`. **Commands** start at column 0 and their continuation lines are indented. A command must start with `theorem NAME`, `lemma NAME` (with a project only; `NAME` is `[A-Za-z_][A-Za-z0-9_']*` components joined by `.`), `example`, `section [A]`, `noncomputable section [A]`, `end [A]`, `variable …`, `universe …`, `set_option maxHeartbeats | maxRecDepth | synthInstance.maxHeartbeats <number> [in]` or `set_option autoImplicit | relaxedAutoImplicit | linter.<name> true | false [in]`, followed by a space or a newline (`variable?` is another command). `open` and `namespace` are refused: in Lean 4 a plain `open X` or `namespace X` activates X's scoped extensions, including scoped *command* syntax that the derived keyword list (computed with nothing open) cannot cover; write fully qualified names instead. Everything else (`def`, `abbrev`, `instance`, `axiom`, `macro`, `syntax`, `attribute`, `initialize`, `run_cmd`, `mutual`, `module`, …) is `LeanVerify::cmd` (or `LeanVerify::unsafe` when it is not the first word).
  - The trusted prelude (`src/internal/lean/Prelude.lean`, byte-identical, once, at the start of a line) is skipped, so `` LeanVerify[Lean`Private`prelude[] <> "theorem …"] `` works. The prelude itself contains no comment, string or character literal (its documentation is in `LeanBackend.m`), so Lean's lexer is at the top level where it ends; a prelude that does not meet this is refused.

  Without user definitions, without meta-level imports and without any span that the two lexers could read differently, there is no user-defined code to run. What remains trusted is stated in [the design note](../../design/lean_integration.md#the-guarantee-precisely): that Lean's built-in (and, with a project, Mathlib's) elaborators and tactics do not execute user-supplied code through ordinary term and tactic syntax. A development with its own definitions has to be checked by Lean / Lake directly.
- **Method scripts** of `Certify` satisfy the same rules (A) and (B): the script is checked by `` Lean`Internal`CheckSource `` as the body of a theorem, on top of `Certify`'s own list (`Lean`, `IO`, …).
- **Kernel replay** (defence in depth). The source is compiled (`lean --json -R <dir> -o <dir>/MathildaCert.olean`, in the private temp directory). After a clean run, the private directory must hold nothing but the source, Lean's output files and the `.olean` (a regular file); then `lean --run src/internal/lean/Replay.lean <dir>/MathildaCert.olean` (through `lake env` with a project) reads **only** the main module file (a module-system file, or a `.olean.server` / `.olean.private` beside it, is refused), replays every declaration into a fresh kernel environment on top of its imports (the per-module mode of the toolchain's `leanchecker`), and computes, from the replayed environment, the axioms each declaration depends on. `"Replay"` is `True` or the reason it failed (e.g. `(kernel) declaration type mismatch`, `unexpected file in the private directory after elaboration: …`); `"ReplayAxioms"` maps every declaration of the module to its axioms; `"ReplayTiming"` is its wall time (about 1.1 s on core Lean). Imported modules (core, Mathlib, …) are trusted as built.
- Every `example` is renamed to `noncomputable def _root_.mathilda_example_k` (same line), so its axioms are checked too. Every declared `theorem`, `lemma` and renamed example -- with its namespace prefix, tracked through `namespace` / `section` / `end` by the source checker -- gets a `#print axioms` line in a footer appended by LeanVerify **after** the check (the user's source has no `#`).
- `"Success"` is `True` only when the source passed the grammar, Lean's exit code is 0, no message has severity `error` or kind `hasSorry`, every reported axiom list is standard, every declared name has a standard axiom report **printed on its own footer line**, the replay accepted the module, **every** declaration of the module depends on standard axioms only according to the replay (an unused `sorry` fails too), and every theorem / example the **source text** declares is among the replayed declarations with standard axioms (a module that drops or renames one fails).
- `"Messages"` is a list of `<|"Severity", "Kind", "Line", "Text"|>` parsed from `lean --json`. `"Axioms"` maps each reported name to its axioms; `"AxiomReports"` lists every report as `<|"Theorem", "Line", "Axioms"|>`.
- With a project (`LeanProject`), the command is `cd project && <lake> env <lean> --json -R DIR -o DIR/MathildaCert.olean FILE`, where `<lake>` is the `lake` binary beside the Lean executable; the replay runs the same way. Nothing is written into the project.
- `LeanVerify::import`, `LeanVerify::cmd`, `LeanVerify::unsafe` and `LeanVerify::lex` (with the line) are the refusals above. `LeanVerify::nolean` is issued when no Lean executable can be run. `LeanVerify::badpath` is issued when a path cannot be quoted safely for the shell (Windows `cmd.exe` only), and `LeanVerify::fail` when the process could not be started.

```mathematica
In[1]:= LeanVerify["theorem ok : 2 + 2 = 4 := by decide"]["Success"]
Out[1]= True

In[2]:= LeanVerify["theorem nd : 2 + 2 = 4 := by native_decide"]["Success"]
LeanVerify::unsafe: the source uses native_decide, which can run code at elaboration time, bypass or forge the kernel check; not run
Out[2]= False

In[3]:= LeanVerify["import Lean\ntheorem t : True := trivial"]["Rejected"]
LeanVerify::import: the source has import Lean; LeanVerify accepts no imports (core Lean) or, with a LeanProject, only import Mathlib / import Mathlib.<Path>; not run
Out[3]= {"import Lean"}

In[4]:= LeanVerify["theorem t : True := trivial -- a comment"]["Rejected"]
LeanVerify::lex: line 1: the source contains a comment (--), which is outside LeanVerify's minimal lexical surface (no strings, comments, character literals, #, $, %, backquotes, backslashes, @[, French quotes or unlisted Unicode); not run
Out[4]= {"a comment (--)"}
```

**Process execution.** Lean runs through ISO C `system()`. The command line holds only executable paths, flags written by Mathilda and temp-file paths, each single-quoted for the POSIX shell (with `set -C`, so the output redirections never follow a planted file). The Lean source, claims and tactic strings are always written to the `.lean` file. Each run creates a private directory with `mkdtemp` (mode `0700`) under `$TemporaryDirectory`; the source is written there with `O_CREAT | O_EXCL | O_NOFOLLOW` (mode `0600`) and the shell's redirections create the output files there, so no other user can plant a file or a symlink, and the output parsed is the one this run wrote (a missing output file is a failure). The directory is removed afterwards with **every** file in it (the known names, then each remaining entry, unlinked -- a symlink is never followed, a subdirectory never recursed into); with `KeepFiles -> True`, the source and its directory are kept. Windows has no `mkdtemp` in its C runtime and keeps a probed fresh name in the temp directory. `system()` blocks, and `TimeConstrained` measures this process's CPU time, so it does not interrupt a Lean run. The bound on a run is the deterministic `maxHeartbeats` of each theorem.

## Settings
- `$LeanExecutable`: `Automatic`, or a path. `Automatic` tries, in order, the environment variable `MATHILDA_LEAN`, the newest `~/.elan/toolchains/*/bin/lean` (called directly, so no `lean-toolchain` file is needed), and `lean` on `PATH`. The first hit is remembered for the session.
- `$LeanProject`: `None`, or the directory of a Lake project with Mathlib built, used for `Tier -> 1`. It defaults to the environment variable `MATHILDA_LEAN_PROJECT`.
- `$LeanVersion`: the `lean --version` string, or `Missing["NotAvailable"]`.

## Using Mathilda from Lean

The other direction: the core-only Lake package in [`lean/`](../../../lean/README.md) gives Lean 4 tactics that ask Mathilda for a certificate and then check it in Lean. Mathilda is still the untrusted oracle. It only supplies data (cofactors, squares, a Pratt certificate), the tactic turns that data into ordinary proof steps, and the Lean kernel checks them. A wrong answer makes the tactic fail; it never proves a false goal. The package imports only `Lean` (no Mathlib), and its `Mathilda/Prelude.lean` is a byte-for-byte mirror of `src/internal/lean/Prelude.lean` (the `lean_prelude_sync` test enforces this).

| Tactic / command | Goal | Mathilda supplies | Lean checks |
|---|---|---|---|
| `mathilda_linear_combination` | `a = b`, from hypotheses `hᵢ : lᵢ = rᵢ` of the same type | cofactors `cᵢ` with `a - b = Σ cᵢ (lᵢ - rᵢ)` | `have kᵢ : lᵢ - rᵢ = 0`, the identity by `grind`, `rw [kᵢ] at hc`, `grind` |
| `mathilda_sos` | `0 ≤ p`, `p ≤ q` or `p ≥ q`, from hypotheses `0 ≤ gᵢ` | triples `(d, ℓ, i)` with `q - p = Σ d ℓ² gᵢ` (`g₀ = 1`) | hints `0 ≤ ℓ²` / `0 ≤ ℓ² gᵢ`, the identity by `grind`, `grind` |
| `mathilda_prime` | `∃ c, Mathilda.PrattCertifies p c` | a flat Pratt certificate | `decide +kernel` (no axioms at all) |
| `#mathilda "expr"` | — | the `LeanForm` of the evaluated `expr` | nothing: it is an information message and proves nothing |

```lean
import Mathilda

example (x y : Rat) (h1 : x^2 + y^2 = 1) (h2 : x = y) : 2*y^2 = 1 := by
  mathilda_linear_combination

example (x y : Rat) : 2*x*y ≤ x^2 + y^2 := by
  mathilda_sos

example : ∃ c, Mathilda.PrattCertifies 274177 c := by
  mathilda_prime

#mathilda "Factor[x^10-1]"
-- (x - 1) * (1 + x) * (1 - x + x ^ 2 - x ^ 3 + x ^ 4) * (1 + x + x ^ 2 + x ^ 3 + x ^ 4)
```

**Features**:
- Goals over any type with `Lean.Grind.CommRing` (and `Lean.Grind.OrderedRing` for `mathilda_sos`): `Int`, `Rat`, and `ℝ` when the user's project has Mathlib. The serializer accepts `+ - * /`, unary `-`, `^` with a natural-number literal exponent, numerals, `Nat.cast` / `Int.cast`, and local variables; anything else is `mathilda: unsupported term …`. Hypotheses that do not serialize are skipped.
- Local variables are sent as `v1, v2, …`. Mathilda answers with **data**: `[c1, c2]`, `[(d, ℓ, i), …]` or a Pratt list. The answer is never elaborated as Lean syntax. `Mathilda.Answer.parse` (`lean/Mathilda/Answer.lean`, a hand-written parser and a pure function on `String`) reads it into an AST of natural-number literals, variables `v<i>`, `+ - * /`, unary `-`, `^` with a literal exponent, lists and tuples; everything else (`by`, `fun`, `sorry`, any other identifier or character) is refused, and sizes are bounded (answer ≤ 10⁶ characters, numerals ≤ 2000 digits, exponents ≤ 10⁴, nesting ≤ 100). The tactic builds the terms from the AST with fixed quotations, numerals and the goal's own free variables (`v<i>` must name one), elaborates them at the goal's type, and then the kernel checks the certificate as before.
- When Mathilda finds no certificate, the tactic fails with `mathilda: Mathilda found no certificate: <reason>`. When the answer is not well-formed data, it fails with `<tactic>: Lean could not parse the certificate returned by Mathilda (<reason>); …` and nothing in the answer is run. When Lean rejects a certificate, it fails with `<tactic>: Lean rejected the certificate returned by Mathilda (…)`. Error recovery is switched off while the certificate is checked, so a failing step is an error, never a `sorry`.
- Each call writes a temporary script, runs `Mathilda -file <script>` with a null stdin through `IO.Process.output`, and reads the answer between the sentinels `@@MATHILDA-BEGIN@@` / `@@MATHILDA-END@@`. A tactic call costs about 0.05 s of Mathilda time.
- The Mathilda binary is the environment variable `MATHILDA_BIN`, else `set_option mathilda.path "…"`, else a `Mathilda` file in the working directory or its parent (the repository layout), else `Mathilda` on `PATH`. A binary that cannot be run gives an error that names these settings.

### `` Lean`Internal`Oracle ``
- `` Lean`Internal`Oracle[kind, args…] ``: the server side of the tactics, called by its full name (`` Lean`Internal` `` is not on `$ContextPath`, so the short name `Oracle` stays the user's; before v0.6523 it was `` Lean`Oracle ``). It returns plain data as a String (numerals, `v<i>`, `+ - * / ^`, lists, tuples), or `"MATHILDA-FAIL: <reason>"`. It needs no Lean installation.

| Query | Answer |
|---|---|
| `` Lean`Internal`Oracle["LinearCombination", {a, b}, {{l1, r1}, …}, {v1, …}] `` | `[c1, …]` |
| `` Lean`Internal`Oracle["SOS", {p, q}, {g1, …}, {v1, …}] `` | `[(d, ℓ, i), …]` for `p ≤ q` |
| `` Lean`Internal`Oracle["Pratt", p] `` | `[(q, a, [(r, k), …]), …]` |
| `` Lean`Internal`Oracle["Eval", "input"] `` | the `LeanForm` of the evaluated input |

The cofactors come from `PolynomialReduce`, else a degree-bounded ansatz (up to degree 4) solved with `LinearSolve`; the squares come from the same engine as `Certify[ForAll[vars, p >= q]]`; the Pratt certificate is the one `Certify[PrimeQ[p]]` uses. Coefficients carry no numeric type: Lean builds them at the goal's type.

```mathematica
In[1]:= Lean`Internal`Oracle["LinearCombination", {2 v2^2, 1}, {{v1^2 + v2^2, 1}, {v1, v2}}, {v1, v2}]
Out[1]= "[1, -v1 - v2]"

In[2]:= Lean`Internal`Oracle["Pratt", 274179]
Out[2]= "MATHILDA-FAIL: 274179 is not prime"
```

**Demo and design.** [`examples/lean/showcase.m`](../../../examples/lean/README.md) runs every certificate kind end to end with timings, including the "Lean catches the CAS" case, a tampered certificate and the tactics. [`docs/design/lean_integration.md`](../../design/lean_integration.md) explains the architecture, the shell-safety argument and the TrustedSteps of each kind. `tests/scripts/lean_leakcheck.sh` is the leak gate.

**Tests.** `lean_tactics_tests` (`tests/scripts/lean_tactics.sh`) builds the package and elaborates `lean/test/Tactics.lean` (every tactic, `#print axioms` of each, `#mathilda`, and false goals that must fail), `lean/test/Untrusted.lean` (run against `lean/test/bad_oracle.sh`, a fake Mathilda that returns wrong certificates, the reviewer's `by_elab` exploit with a file write, `sorry`, a foreign identifier, a huge exponent / numeral, malformed tokens and an unknown variable: Lean rejects or refuses every one, adds no theorem and writes no file) and `lean/test/Missing.lean` (a missing binary). It is Skipped (exit 77) when no `lake` or no built `Mathilda` is found.
