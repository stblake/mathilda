# `benchmarks/osp_famous/` — the math.osp.no "famous" problems

The [math.osp.no](https://math.osp.no/problems?famous=any) CAS benchmark curates
222 "famous" problems (Basel, Ahmed's integral, Apéry's recurrence, Painlevé I,
Burgers' equation, …). Each one has a statement and a reference answer. This
directory grades Mathilda on all of them **locally**. The site grades only
runs that are uploaded to it.

```bash
make bench-osp-famous                                  # all 222: results/math.{md,jsonl}
python3 benchmarks/osp_famous/osp_famous.py --only basel_problem,ahmed_integral -v
python3 benchmarks/osp_famous/osp_famous.py --binary /tmp/Mathilda_snapshot   # keep building meanwhile
python3 benchmarks/osp_famous/osp_famous.py --runner mathilda-run.jsonl       # grade a math-run file
```

| file | what it is |
|---|---|
| `famous.jsonl` | the 222 problems with reference answers, from the site's `quick.jsonl.gz` (suite 2026.1, CC BY 4.0, © math.osp.no editors) |
| `baseline_v0.6222.tsv` | the site's official grades for Mathilda v0.6222 (run `r-hkums2rrgh`) |
| `osp_famous.py` | translator (samara neutral syntax → Mathilda input), driver, report |
| `check_famous.m` | the grader, run by `./Mathilda -file` (Mathilda grades itself) |
| `results/` | `math.md` / `math.jsonl` (our translation), `runner.*` (a math-run file) |

## Two scores

The site runs Mathilda through its runner, **math-run 0.1.0**. Its Mathilda
adapter shares its translation with another system's adapter, and that translation has no form
for four of the suite's own functions inside an expression: `harmonic`,
`integrate`, `sum` and `product`. It passes them on as unknown lower-case
heads (`harmonic[k]`, `integrate[f, y, 0, 1]`).

As a result, nine famous problems cannot be answered in the form the site
poses them: three Euler sums, three Watson integrals, two Beukers integrals and
the Rogers–Ramanujan continued fraction. That is a defect of the runner, not of
any system.

`osp_famous.py` therefore translates every problem itself into the Mathilda
input the runner intends. It uses the same heads, the same
`Assumptions -> …`, and the same `mr` prefix on the single capitals C D E I N O
K. It grades that translation as the **math** score. `--runner` grades a math-run
upload file as the **runner** score, which is what the site would report.

## Grading

Every check is done by Mathilda at 60 digits, with a relative tolerance of 1e-30.
Free parameters are sampled inside the problem's assumptions.

- **Values** (integrate, sum, product, limit, laplace, determinant, …): the answer
  is free of the task's head and agrees with the reference. Indefinite integrals
  are compared by derivative.
- **simplify / evaluate**: the value check, plus the site's "worked out" rule.
  No head of the input that the reference eliminated (Cos, ArcTan, PolyLog,
  EllipticK, Root) may survive. Radical nesting may be no deeper than the
  reference's. The answer may be at most 1.5 × the reference's LeafCount + 4.
- **ODE/PDE**: the answer is substituted into the equation and the conditions,
  with the constants `C[i]` sampled. A condition at a removable singularity is
  approached from inside. An ODE without conditions must carry `order` constants.
  An implicit first-order answer is checked by implicit differentiation.
- **Recurrences**: the answer satisfies the recurrence for n = n₀ … n₀+10, and
  it meets the conditions.
- **Solve**: every answer satisfies the equations. A finite solution set must
  match the reference in size and numerically. A parametric family must contain
  the reference's k = 0 member.
- **Open or impossible problems** (no reference answer: Painlevé, Blasius,
  `MinimalPolynomial[Pi]`, …) count as solved when left unevaluated, and wrong
  when an answer is claimed.

Grades are `solved`, `wrong`, `unsolved` (left unevaluated, or equal but not
worked out), `unclear` (the checker could not decide), `timeout`, `crash` and
`error`.
