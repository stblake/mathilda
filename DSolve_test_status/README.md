# DSolve test status

Cross-session dashboard for `DSolve` development against Nasser Abbasi's
*Solving ODEs* corpora (<https://12000.org/my_notes/solving_ODE/current_version/>).
Each section there lists ODEs that Maple **and** Mathematica solve; we track how
many `DSolve` currently solves and drive that number to parity, one general
method-wave at a time (see `../DSOLVE_PLAN.md`, milestones M15+).

The point of this folder: **pick up DSolve development between sessions with a
clear, self-contained context** — the corpora, the self-verifying harness, and a
living scoreboard all in one place.

## Contents

| File | Role |
|---|---|
| `DE_examples_2.m` | Section **2.1.2** corpus — 1204 ODE records (1000 scalar + 204 systems). |
| `DE_examples_3.m` | Section **2.1.3** corpus. |
| `test_dsolve_corpus.c` | Fork-per-case runner (compiled via `tests/CMakeLists.txt`). |
| `dsolve_corpus_prelude.m` | Self-verifier: runs `DSolve` under `TimeConstrained` and numerically back-substitutes each branch. |
| `STATUS.md` | **The scoreboard** — per-section, per-bucket solve counts + wave history. Update after every wave. |
| `reports/` | Per-section bucketed gap reports (regenerated from the run TSV). |

A corpus record is `{"label", equation(s), function(s), indVar, "MapleClassif", sympy?}`.
Systems (`classif` contains `system_of_ODEs`, or a `List` function) are **skipped** —
this is the scalar-first campaign.

Verdict codes: `0 PASS` (verified closed form), `1 FAIL` (a branch is
demonstrably nonzero — a wrong answer), `2 UNEVAL` (declined / timed out / `{}`),
`3 SKIP` (system).

## Run a section (progress dashboard)

```bash
cd tests/build && cmake .. && make dsolve_corpus_tests -j
ctest -R dsolve_corpus_2_1_2_tests --output-on-failure     # gated at STATUS.md baseline

# Manual full run, capturing the per-case TSV for the bucket report:
DSOLVE_CORPUS_BASELINE=100000 ./dsolve_corpus_tests \
    ../../DSolve_test_status/DE_examples_2.m \
    ../../DSolve_test_status/dsolve_corpus_prelude.m \
    > /tmp/corpus.tsv 2>/tmp/corpus.err
python3 ../../tools/dsolve_corpus_report.py /tmp/corpus.tsv \
    > ../../DSolve_test_status/reports/2.1.2.md
```

The `ctest` gate fails only when non-PASS scalar cases exceed the per-section
baseline (argv[3] in `tests/CMakeLists.txt`). Each landed method-wave must
**lower** that baseline; a newly-failing case trips the test.

## Regenerate / add a corpus

The site 403s automated fetch, so save the section HTML manually, then convert:

```bash
curl -sL -A "Mozilla/5.0" \
  'https://12000.org/my_notes/solving_ODE/current_version/indexsubsectionN.htm' \
  -o /tmp/sectionN.html
python3 tools/latex_ode_to_mathilda.py /tmp/sectionN.html \
  DSolve_test_status/DE_examples_K.m --label 2.1.K
```

Then add a `dsolve_corpus_<section>_tests` entry in `tests/CMakeLists.txt` and a
section block in `STATUS.md`. The converter is section-agnostic (auto-selects the
problems table and reads its column layout from the header row, which varies between
sections).

**Initial-value problems.** When the source rows carry initial conditions (`y(0)=3`,
`y'(0)=10`, symbolic `y(a)=b`), the converter emits the record's equation slot as the
DSolve-native list `{ode, ic1, …}` (the function slot stays a bare symbol, so the scalar
harness still handles it). The prelude solves `DSolve[{ode, ics}, y, x]` and verifies the
ODE residual **and every initial condition**; an IVP whose solved branch still carries a
generated constant `C[k]` (the general solution with the IC unfitted) scores UNEVAL, not
PASS. §2.2.1 is the first IVP-carrying section; §2.1.2 has no ICs (bare equations).

## After a method-wave

1. Re-run the affected section(s), regenerate `reports/<section>.md`.
2. Update the counts and the wave-history line in `STATUS.md`.
3. Lower the section baseline (argv[3]) in `tests/CMakeLists.txt` to the new
   non-PASS count.
