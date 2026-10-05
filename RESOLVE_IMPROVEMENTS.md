# Resolve improvements: simple decidable statements currently declined

`Resolve` / `Reduce` decide first-order formulas over the reals by cylindrical
algebraic decomposition (CAD, `src/solve/reduce_cad.c`) — a complete decision
procedure for the first-order theory of real-closed fields (Tarski). In
principle every **semialgebraic** sentence (built from polynomial/rational/
algebraic atoms and `ForAll`/`Exists` over `Reals`) is decidable. In practice a
small, structurally simple class is **declined** — returned unevaluated — when it
should return `True` or `False`.

This file tracks those gaps. It was motivated by the ε–δ limit-proof examples in
the book (§5.3.2 *Limits*), where the finite cases prove but the infinite ones do
not.

**Soundness first.** Every failure below is a *decline* (the input comes back
unevaluated). No wrong `True`/`False` has ever been observed for these. The gaps
are **completeness** gaps, not soundness bugs — consistent with the system's
"decline rather than guess" contract. Any fix must strictly *expand* what is
proved and must never turn a decline into a wrong verdict.

Observed on **v0.284** (2026-10-05). Reproduce with the block at the end.

> **RESOLVED in v0.289 (2026-10-05): the whole of failure class A below is fixed,
> and the §2 root-cause hypothesis was wrong.** The CAD handles the
> parameter×variable *product* perfectly. The real cause was a **missing
> denominator-clearing pass**: chained quantifier elimination emits a *rational*
> atom at an intermediate level (the inner `∀x` of the product case eliminates
> correctly to `del>0 && m <= 1/del`), and every multivariate real engine declined
> any atom with a non-constant denominator (`nonconst_denom`). Clearing
> `p/q REL 0` to the sign-exact polynomial form in the multivariate Reals
> preprocessing (`reduce_piecewise_preprocess`, plus `Abs[p/q] -> Abs[p]/Abs[q]`
> so the Abs sign-split stays sound at a pole) closes all three §1 gap rows and
> both §2 product cases. Class B (higher-arity *refutation*) and the transcendental
> declines (§4) are unaffected and remain open. The original analysis is kept below
> as the record; the frontier table's "declined" rows marked (A) are now decided.

---

## 1. The decidability frontier (observed)

Each row is the literal ε/M–δ (or –N) transcription of a limit definition,
`Resolve[..., Reals]`.

| Limit class | Example | Statement shape | Result |
|---|---|---|---|
| Finite value, finite point, prove | `lim_{x→2}(3x−1)=5` | `∀ε ∃δ ∀x (0<\|x−2\|<δ ⇒ \|3x−1−5\|<ε)` | **True** |
| Finite value, finite point, refute | `lim_{x→2}(3x−1)=6` | same, wrong `L` | **False** |
| Finite value, finite point (1-sided) | `lim_{x→0⁺} x = 0` | `∀ε ∃δ ∀x (0<x<δ ⇒ x<ε)` | **True** |
| Finite value, planar, prove | `lim_{(x,y)→0}(x+y)=0` | `∀ε ∃δ ∀{x,y} (0<x²+y²<δ² ⇒ \|x+y\|<ε)` | **True** |
| Infinite value, **at infinity**, prove | `lim_{x→∞} x² = ∞` | `∀M ∃N ∀x (x>N ⇒ x²>M)` | **True** |
| Infinite value, **at infinity**, refute | `lim_{x→∞}(−x²)=∞` | same, `−x²` | **False** |
| **Finite value, planar, refute** | `lim_{(x,y)→0}(x+y)=1` | `∀ε ∃δ ∀{x,y} (0<x²+y²<δ² ⇒ \|x+y−1\|<ε)` | **declined** (class B, open) |
| ~~Finite value, at infinity, prove~~ | `lim_{x→∞} 1/x = 0` | `∀ε ∃N ∀x (x>N ⇒ \|1/x\|<ε)` | **True (v0.289)** |
| ~~Infinite value, finite point, prove~~ | `lim_{x→0} 1/x² = ∞` | `∀M ∃δ ∀x (0<\|x\|<δ ⇒ 1/x²>M)` | **True (v0.289)** |

The last three rows *were* the gaps. Two (class A) are **fixed in v0.289** — see the
note at the top. The planar refute (class B) remains open. They split into two
independent causes, analysed below.

---

## 2. Failure class A — a parameter×variable product in the matrix

This is the main gap and the easiest to characterise. Hold the quantifier
pattern and the neighbourhood *fixed* and vary only the shape of the matrix
atom:

```
[1] finite pt, ADDITIVE   x < eps          → True
    Resolve[ForAll[eps, eps>0, Exists[del, del>0, ForAll[x, 0<x<del, x < eps]]], Reals]

[2] finite pt, PRODUCT    m x < 1           → declined
    Resolve[ForAll[m,   m>0,   Exists[del, del>0, ForAll[x, 0<x<del, m x < 1]]], Reals]

[3] at infinity, ADDITIVE x^2 > m           → True
    Resolve[ForAll[m,   m>0,   Exists[n,   n>0,   ForAll[x, x>n,     x^2 > m]]], Reals]

[4] at infinity, PRODUCT  eps x > 1         → declined
    Resolve[ForAll[eps, eps>0, Exists[n,   n>0,   ForAll[x, x>n,     eps x > 1]]], Reals]
```

[1] and [2] have the *same* shrinking neighbourhood `0<x<δ` and the *same*
`∀∃∀` prefix; [3] and [4] have the *same* ray `x>N`. The only difference is the
atom: `x < ε` / `x² > M` put the outer-∀ parameter **additively** against a
variable expression, whereas `M·x < 1` / `ε·x > 1` form a **product of the
outer-∀ parameter with the inner-∀ variable**. The additive ones decide; the
product ones decline. All four are true (or, with a wrong constant, false); none
is harder than first-year analysis by hand (`δ = 1/M`, `N = 1/ε`).

**Why this blocks real mathematics.** The product is not an artefact of a clumsy
transcription — it is intrinsic to two whole families of limit:

- **Infinite limits at a finite point** (vertical asymptotes), `lim_{x→a} f = ∞`.
  The definition is `f(x) > M`; for the archetype `f = 1/x²` this is `1/x² > M`,
  which clears to `M·x² < 1` — a parameter×variable product. No transcription
  avoids it (`1/x² > M` and `M x² < 1` and `x² < 1/M` are the same product), and
  it declines in every form tried: literal and cleared, one- and two-sided,
  `1/x`, `1/x²`, `1/|x|`, with `Abs` or with the squared neighbourhood
  `0<x²<δ²`, with `δ` bounded, and under `Reduce` as well as `Resolve`.
- **Limits at infinity with a finite value**, `lim_{x→∞} f = L`. The definition
  `|f − L| < ε` for a rational `f` clears to a product with `ε` (e.g.
  `1/x < ε ⟺ ε·x > 1`). Declines literal and cleared (`1/x → 0`,
  `(2x+1)/x → 2`).

So the single product gap removes **all** infinite limits at a point and **all**
finite limits at infinity — the bulk of a standard limits chapter — while the
additive cases (finite limit at a point, `x → ∞` of a function that *grows* to
∞) go through.

**Root-cause (CONFIRMED in v0.289 — the hypothesis below was WRONG).** The CAD does
*not* blow up on the product: `Reduce[Exists[del, del>0 && m del <= 1], {m}, Reals]`
→ `True`, `Reduce[m del <= 1 && del>0, {m, del}, Reals]` → the correct region, and
`Resolve[ForAll[m, m>0, Exists[del, del>0 && m del <= 1]], Reals]` → `True`. The
decline was a **front-end plumbing** gap: the `∃`-witness is non-polynomial in the
parameter, so the inner `∀x` *correctly* eliminates to a **rational** atom
(`del>0 && m <= 1/del`), and the next elimination level then declined it because
every multivariate real engine bails on an atom with a non-constant denominator
(`nonconst_denom` in `reduce_atom.c`; checked in `reduce_cad.c` / `reduce_fm.c` /
`reduce_sys.c` / `reduce_zerodim.c`). The Reals preprocessing cleared radicals and
split selectors but never cleared rational denominators. The fix (v0.289) adds that
clearing pass (`p/q REL 0 → p q REL 0 [&& q≠0]`, with `Abs[p/q]→Abs[p]/Abs[q]` to
keep the Abs sign-split sound at a pole) — a general preprocessing rewrite, not the
per-case `t = 1/x` substitution the superseded hypothesis below guessed at.

> *Original (unverified, now disproved) hypothesis, kept for the record:* "A product
> `p·v` … raises the total degree of the projection factors … CAD projection/lifting
> on these mixed product terms either blows up or hits a branch that gives up …
> the fix may be a normalising substitution in the front end rather than a change to
> the core CAD." — The CAD was never the problem; the missing denominator-clearing
> pass was.

---

## 3. Failure class B — refutation cost at higher arity

Separate from the product issue. Refuting a *wrong* limit is cheap in one
variable (the `=6` row in §1 returns `False`) but declines in the plane:

```
planar refute, lim_{(x,y)→0}(x+y) = 1  → declined
    Resolve[ForAll[eps, eps>0, Exists[del, del>0,
      ForAll[{x,y}, 0<x^2+y^2<del^2, Abs[x+y-1] < eps]]], Reals]
```

Proving a *true* planar limit is fast because an `∃δ` witness is found early;
refuting a wrong one must show that **no** `δ` works, which drives the full
decomposition in four variables (`x, y, δ, ε`), and CAD is doubly exponential in
the variable count. At present Mathilda runs it to its budget and declines. This
is a **scalability** gap (faster projection, better variable ordering, or a
specialised refutation search), not the matrix-shape gap of class A. It is
already noted in the book's §5.3.2 "Under the hood".

---

## 4. Out of scope (correct declines, listed to avoid confusion)

Transcendental targets such as `lim_{x→0} sin(x)/x = 1` are **not** semialgebraic
and have no CAD. Declining them is correct, not a bug; proving them would need a
different engine (and `Limit` already reaches them). These are not "simple
failures" and should not be conflated with classes A and B.

---

## 5. Suggested improvements, by payoff

1. **Class A (parameter×variable product). ✅ DONE in v0.289.** It was not a
   product/CAD problem at all but a missing **rational-denominator clearing** pass
   (the inner elimination emits `m <= 1/del`, which the engines declined). Clearing
   `p/q REL 0 → p q REL 0 [&& q≠0]` in the multivariate Reals preprocessing (with
   `Abs[p/q]→Abs[p]/Abs[q]`) unblocked both infinite limits at a point and finite
   limits at infinity. See the §2 root-cause note and the v0.289 changelog.
2. **Class B (higher-arity refutation).** Faster/space-bounded CAD or a dedicated
   "no witness exists" search so wrong planar (and higher) limits are *refuted*
   rather than declined.
3. **Longer term.** A non-CAD path (series/Gruntz-style or nullstellensatz
   certificates) for transcendental ε–δ statements — out of scope for the real
   QE engine, tracked here only for completeness.

---

## 6. How to reproduce

```
./Mathilda -file repro.m          # with Print[...] wrappers, since -file echoes only Print output
```

The verified reproducers are the four lines in §2 and the three in §1/§3. A
declined result prints the input expression back verbatim; a decided one prints
`True` or `False`.
