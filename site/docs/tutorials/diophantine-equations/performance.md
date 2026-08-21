# How Mathilda compares

How does Mathilda's integer solver stack up against the tools a working
number-theorist actually reaches for? Two of them are the natural yardsticks:
**sympy**, whose `diophantine` is the most widely used open-source Diophantine
solver, and **PARI/GP**, whose `thue()` is the gold-standard implementation of the
Thue equation. Mathilda ships a purpose-built benchmark against each
(`benchmarks/87-diophantine-integers/` and `benchmarks/88-thue-equations/`), and
the picture is genuinely two-sided: against sympy, Mathilda wins on *both*
coverage and speed; against PARI, it matches on *correctness* but trails on speed
and breadth. This page reports both honestly.

---

## 1. Versus sympy

The comparison covers 19 representative equations (`cases.py`, run against sympy
1.14.0). The headline:

- **Mathilda** returns the finite, constraint-satisfying set for **19 of 19**,
  all matching the known solution count.
- **sympy `diophantine`** answers **3** directly (all sums of squares), **3** only
  *parametrically* (it returns a symbolic family but not the bounded answer the
  question asked), rules 1 empty by a gcd argument, and raises
  **`NotImplementedError` on 11** — every mixed, cubic, and exponential form, and
  *every system*, because `diophantine` takes a single equation.

Timing is in-process on both sides (Mathilda `AbsoluteTiming`, sympy
`perf_counter`), excluding start-up.

### Where both succeed — Mathilda is 20–100× faster

| Case | Mathilda | sympy `diophantine` | count |
|---|--:|--:|--:|
| \(x^2 + y^2 = 25\), box \(\pm5\) | 0.090 ms | 10.4 ms | 12 |
| \(x^2 + y^2 = 10^6\)-box | 0.459 ms | 9.2 ms | 28 |
| \(x^2 + y^2 + z^2 = 29\) | 0.235 ms | 9.8 ms | 72 |

Even on its home turf — sums of squares, the one family `diophantine` returns
directly — sympy is one to two orders of magnitude slower.

### Where sympy answers only parametrically

| Case | Mathilda | sympy `diophantine` | count |
|---|--:|---|--:|
| Pell \(x^2 - 61y^2 = 1\), \(x < 10^{10}\) | 0.134 ms | parametric only (24.6 ms) | 1 |
| Negative Pell \(x^2 - 2y^2 = -1\) | 0.132 ms | parametric only | 4 |
| Linear \(3\)-var, positive box | 2.1 ms | parametric only | 703 |
| Pythagorean, perimeter 3000 | 0.510 ms | 252 raw, needs hand-elimination | 3 |

For \(x^2 - 61y^2 = 1\) the honest bounded answer needs the fundamental solution
\((1766319049, 226153980)\); sympy hands back a symbolic recurrence, and the
same-box Python search a sympy user would write instead **times out past 25
seconds**. Mathilda returns the pair in 0.134 ms.

### Where sympy has no solver at all (`NotImplementedError`)

| Case | Form | Mathilda | count |
|---|---|--:|--:|
| Mixed quadratic–cubic | \(x^2 + 2y^3\) | 0.123 ms | 3 |
| Two-quadratic *system* | coupled | 0.881 ms | 0 (proved) |
| Egyptian fractions | \(4/2027 = 1/x+1/y+1/z\) | 740.6 ms | 73 |
| Markov triples | \(\le 1000\) | 427.8 ms | 13 |
| Mordell | \(y^2 = x^3 + k\) | 47.5 ms | 1 |
| Taxicab | two cubes, two ways | 22.0 ms | 10 |
| **Lander–Parkin** | \(x^5+y^5+z^5+w^5 = r^5\) | 201.2 ms | 1 |
| Catalan / exponential | variable exponents | 0.126 ms | 1 |

Every cubic, exponential, mixed, and multi-equation form is a `NotImplementedError`
for sympy — including the two-quadratic system, which `diophantine` cannot even
express. On the Lander–Parkin quintic, the brute-force Python fallback takes
**15 seconds**; Mathilda takes 201 ms.

!!! note "Soundness, checked on held-out equations"
    Speed means nothing without correctness. A separate **held-out gate**
    (`make check-diophantine-heldout`) runs Mathilda on 40 equations drawn from
    standard references that the engine was *not* developed against, cross-checking
    each against an independent brute-force oracle. The result: **35 OK, 5 DECLINE,
    0 WRONG** — never a wrong or silently-incomplete answer, only honest declines
    on out-of-reach inputs.

---

## 2. Versus PARI/GP

For **Thue equations**, PARI/GP's `thue()` is the reference implementation, and
here the story reverses: PARI is faster and covers more. Mathilda's job is to
match it *where it answers* and to decline safely elsewhere. Over an adversarial
corpus of 113 Thue equations (`benchmarks/88-thue-equations/`), Mathilda scores
**111 CORRECT, 2 DECLINE, 0 WRONG**, and a randomised 400-case grid (degrees 3–6,
mixed \(m\)) adds **281 CORRECT with 0 Mathilda bugs**.

### PARI is faster on the hard forms

| Form \(= m\) | Mathilda | PARI `thue()` | ratio |
|---|--:|--:|--:|
| \(131y^3 - 211xy^2 + 97x^2y + x^3 = 1\) | 410.5 ms | 11 ms | 37× |
| \(x^6 - 3y^6 = 1\) | 252.4 ms | 10 ms | 25× |
| \(x^5 - 5y^5 = 1\) | 238.9 ms | 10 ms | 24× |
| \(x^3 - 41y^3 = 1\) | 146.0 ms | 8 ms | 18× |

Simple cubics run in a few milliseconds either way, but on the hardest forms PARI
is several to \(\sim\)40× quicker. It is also **unconditionally broader**:
`thue()` solves *any* Thue equation regardless of degree, field, or \(m\), whereas
Mathilda still declines a handful — very large regulators, some non-monogenic
fields, some \(|m| \ne 1\). Those declines are safe non-answers, and the roadmap
to close them lives in `docs/design/thue_completion_plan.md`.

### But Mathilda is sound — and caught PARI in an error

The randomised grid does not blindly trust the oracle: on any disagreement it
verifies the disputed points directly against \(F(x,y) = m\). That caught a
genuine **incompleteness in PARI**. On
\[
y^4 - 3xy^3 + 4x^2y^2 - 2x^3y + x^4 = 5 \quad\text{over } \mathbb{Q}(\zeta_5),
\]
PARI `thue()` returned the empty set, but \((1, 2)\) and \((-1, -2)\) are genuine
solutions — which Mathilda found and brute-force-verified. So while PARI is the
faster and more general tool, "faster" is not "infallible": Mathilda's insistence
that every reported set be provably complete let it flag a case the reference
implementation missed.

---

## 3. The bottom line

- **Versus sympy** — Mathilda wins on **coverage and speed**. It solves every
  case sympy's `diophantine` cannot (systems, cubics, exponentials, mixed forms)
  and is 20–100× faster even where sympy succeeds.
- **Versus PARI/GP** — **correctness parity plus one oracle catch**, with a real
  **speed and coverage deficit** on Thue equations that is an active roadmap, not
  a hidden limitation. Where Mathilda cannot yet match PARI's generality, it
  declines rather than guess.

Both benchmarks are reproducible: `benchmarks/87-diophantine-integers/run.py`
(needs a virtualenv with sympy) and `benchmarks/88-thue-equations/run.py` plus
`grid.py` (need `gp` on the `PATH`). The development-driving document
`SOLVE_INTEGERS.md` tracks the open tiers, and `docs/design/thue_completion_plan.md`
the Thue-specific completion plan.

---

You have now toured the whole integer solver — from the extended Euclidean
algorithm to the Tzanakis–de Weger method, and from Pell's equation to
Lander–Parkin. Return to the **[Diophantine index](index.md)** for the family
map, revisit the general solver in **[Solutions of
equations](../07-solutions-of-equations.md)**, or consult the
[`Solve` reference](../../documentation/solutions-of-equations/Solve.md) for
every option in full.
