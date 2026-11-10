# Gradshteyn & Ryzhik: Integrate corpus, table stage, triage — campaign status

Status as of **v0.6493 (2026-11-10)**. This is the hand-off document for the G&R
campaign: what exists, how to rebuild every artefact, the measured numbers, and the
open work in priority order. Read it before touching `tools/gr/`,
`src/internal/GRTable.m` or `src/internal/GRTables/`.

Source: I. S. Gradshteyn & I. M. Ryzhik, *Table of Integrals, Series, and Products*,
7th ed. (2007), chapters 2–7 and 17 (the integrals; chapters 8+ are special-function
identities). Transcription contract: `tools/gr/TRANSCRIBE.md`.

## 1. What exists

| Piece | Where | State |
|---|---|---|
| Transcriptions (104 page batches) | `tools/gr/raw/p*.json` | done; gap audit clean |
| Merged entries (9761 formulas) | `merge.py` → `all_entries.json` (snapshot) | done |
| Numeric oracle | `tools/gr/verify.m` via `verify.py` | done |
| Python cross-check (mpmath/scipy) | `tools/gr/pycheck.py` | done |
| String repair (pre-v0.6488 printer bug) | `tools/gr/restring.{m,py}` | done; 34 records repaired |
| Integrate corpora (7157 records) | `tools/gr/build_corpus.py` | built in scratch, **not yet in `tests/`** |
| Book check (book value in place of Integrate) | `build_corpus.py --book` | done: no formula wrong as stated |
| Condition audit | `tools/gr/condaudit.{m,py}` | first pass done (old script); **re-run pending** |
| Table engine | `src/internal/GRTable.m` | done |
| Table records (3745) | `src/internal/GRTables/*.m` via `build_table.py` | committed |
| Integrate wiring | `src/calculus/integrate.c` (`try_grtable`, `try_grtable_definite`, `METHOD_GRTABLE`) | done; **last stage of both cascades** |
| Unit test | `tests/test_grtable.c` | passes |
| Wrong Integrate answers | GitHub issue #132 | filed: 46 wrong results across 42 records |

The pipeline state that is not derivable from the repo is snapshotted gzipped in
`tools/gr/corpus/state-2026-11-10/` (`gunzip -k` before use):

| File | What |
|---|---|
| `all_entries.json` | every transcribed formula (`merge.py` output) |
| `{ch2,ch34,ch5to17}_entries.json` | the per-batch entries `build_table.py` pairs with |
| `ch2_py2r.json`, `ch34_py3r.json`, `ch5to17_py3r.json` | oracle + Python verdicts, string-repaired: THE verified inputs |
| `names.txt` | `Names["*"]` of the binary (head-existence check) |
| `run2.jsonl` (ch2), `runB.jsonl` (ch3–7, 17), `run_logabs.jsonl` | **baseline** Integrate runs with an EMPTY table (the CAN-DO filter) |
| `runT.jsonl` | a preview pass WITH the table (early slot, pre-fix engine) |
| `book.jsonl` | the book check |

The earlier snapshot `tools/gr/corpus/{ch2,ch34,ch5to17}/` (2026-11-09) holds the
pre-repair verifier outputs; `state-2026-11-10/` supersedes it as input.

## 2. Rebuilding

Paths below use `S=tools/gr/corpus/state-2026-11-10` after `gunzip -k $S/*.gz`.

```sh
# corpora (Integrate questions; S = symbolic + Assumptions, N = numeric instance)
python3 -I tools/gr/build_corpus.py $S/all_entries.json \
    $S/ch2_py2r.json,$S/ch34_py3r.json,$S/ch5to17_py3r.json OUT $S/names.txt \
    --book-wrong tools/gr/audit/condaudit_wrong.tsv --cond-fix tools/gr/cond_fixes.tsv
# book check: same, plus --book (writes GRBook<ch>Corpus.m)

# table
python3 -I tools/gr/build_table.py \
    --batch $S/ch2_entries.json:$S/ch2_py2r.json,$S/ch34_entries.json:$S/ch34_py3r.json,$S/ch5to17_entries.json:$S/ch5to17_py3r.json \
    --runs $S/run2.jsonl,$S/runB.jsonl,$S/run_logabs.jsonl --names $S/names.txt \
    --book-wrong tools/gr/audit/condaudit_wrong.tsv --cond-fix tools/gr/cond_fixes.tsv \
    --out src/internal/GRTables

# condition audit of a built table (definite records); --ids for a subset
python3 -I tools/gr/condaudit.py src/internal/GRTables OUT.tsv --jobs 8
#   -> OUT.tsv, OUT.tsv.fix (automatic corrections), OUT.tsv.wrong (removals)

# a corpus pass (runner: tests/test_mathar_corpus.c, target mathar_corpus_tests)
MATHAR_ONLY=id,id,... MATHAR_JOBS=10 MATHAR_TIMEOUT=60 MATHAR_CORPUS_JSON=run.jsonl \
    tests/build/mathar_corpus_tests OUT/GR3IntegralsCorpus.m
# a baseline (empty table): point MATHILDA_HOME at a copy of src/internal whose
# GRTables/*.m keep only the header/footer
```

`--runs` order matters: a later file's S status overrides an earlier one's
(`run_logabs.jsonl` re-baselines the 7 records whose integrand the log convention
changed).

## 3. Selection and soundness rules (all enforced in code)

A formula becomes a table record only when (`build_table.py` docstring):
1. The oracle verified it, either strictly or loosely with the Python check also
   verifying, and no Python probe disagrees.
2. The baseline `Integrate` with an empty table did not return CORRECT in S
   mode. The table never replaces an algorithm.
3. Every head exists in Mathilda.
4. The value has no `Sum`/`Product`/`Integrate`/`Limit`, and no test of its own
   parameters (`If`, `Boole`, `OddQ`, …). Such a test decides on a bare symbol.
5. The formula is not on a removal list (`--book-wrong`), and its corrected
   conditions (`--cond-fix`) are used.

Rules the engine enforces:
- **Log convention.** `Integrate[1/x, x]` is `Log[x]`. Indefinite records have
  `Log[Abs[u]]` rewritten to `Log[u]` in both the integrand and the value
  (`tools/gr/logconv.py`, also used by `build_corpus.py`). Other `Abs` forms are
  left as the book prints them, for example 2.449.2a's `ArcTan[Abs[..]]` and the
  2.813/2.833 `Sign[a] ArcSin[x/Abs[a]]` family. The derivative self-check
  decides whether those are served (§6, item 6).
- **Principal root.** When a power slot has several solutions, the one that is
  real and positive for positive query symbols is tried first. Across all
  records in the bucket, the unconditional answer with the most principal
  solution wins.
- **Conditional answers.** A `ConditionalExpression` is returned only under
  `Method -> "GRTable"`. In the automatic cascade it is a decline.
- **Contexts.** Engine names live in ``GRTable`Private` ``. Record symbols live in
  ``GRTable`R` ``. Never generate records into the engine context: a `Do`
  iterator is not renamed the way a `Module` local is. 17.13.30's parameter `s`
  hit `$RecursionLimit` that way.
- **Self-check.** A definite answer with no free symbols is compared with
  `NIntegrate` (3 s). An antiderivative is differentiated at generic parameter
  values, with `Abs[u]` read as `Sqrt[u^2]` for the derivative. A symbolic
  definite answer has no numeric check, so its soundness rests on the book's
  conditions. That is why the condition audit exists.

## 4. Book errors: corrected from first principles, or removed

User rule (2026-11-10): an incorrect G&R formula is corrected from first
principles or removed from the corpus and the table. It is never kept on the
book's authority.

1. **Book check.** Each formula's value replaced Integrate in the runner and was
   checked at every instantiation and x point.
   - 14155 of 14314 S/N checks passed. None of the failures was a book error.
   - The WRONGs were three Mathilda bugs (fixed in v0.6487–v0.6489) plus a
     substitution-order artefact (fixed in `build_corpus.py`).
   - The NONNUMERICs are the runner unable to differentiate `Abs`.
   - 6.666 is an NSum bug (§5).
   - Every unconfirmed antiderivative in the table passed a 60-digit
     central-difference check.
2. **Condition audit** (`condaudit.py`). At up to four sampled points that the
   stated conditions admit, the book's value is compared with `NIntegrate`. A
   quadrature that fails to converge with an estimate more than 10⁶ times the
   value counts as divergence.
   - First pass: 3580 definite records, **268 disagreed**, 554 had a point where
     the quadrature failed, and 291 had no checkable point.
   - Of the 268, 266 are removed (`tools/gr/audit/condaudit_wrong.tsv`) and 2
     are corrected by hand (`tools/gr/cond_fixes.tsv`, each with its
     derivation):
     - 3.197.7: the endpoint powers x^(μ−1/2) at 0 and x^(−μ−1/2) at ∞ give
       μ > 1/2, where the book prints Re μ > 0.
     - 3.165.2: the integrand is even in a and the value is odd, so a > 0.
   - The first pass ran the old script, which had an 8 s quadrature budget and
     no correction step.
   - The current script adds the automatic correction: when the integrand and
     the limits are even in a parameter p, the integral is even in p. It tries
     p > 0, one parameter at a time and then all together. It also uses a 20 s
     budget, so a divergent quadrature reports instead of timing out.

## 5. Measured

The baseline is `Integrate` with an empty table (S correct, N correct):

| Chapter | S | N |
|---|---|---|
| 2 | 1052/1651 | 1364/1651 |
| 3 | 483/2527 | 1128/2527 |
| 4 | 293/1074 | 414/1074 |
| 5 | 18/102 | 22/102 |
| 6 | 268/1077 | 344/1077 |
| 7 | 9/485 | 149/485 |
| 17 | 161/241 | 195/241 |

There is a preview pass with the table, `runT.jsonl`. It ran with the table in its
early slot, before the root, conditional, context and audit fixes, on the 4017
records of that first table. It closed 3223 S and 3464 N of them, all integrals
the baseline had not closed.
- WRONG in N: 1. That one is 7.244.2, already WRONG in the baseline.
- WRONG in S: 52. After the fixes, re-running every record that was WRONG in
  that pass leaves only baseline Integrate bugs, all in issue #132:
  3.354.5, 3.724.3, 3.731.1–4, 3.741.3a, 3.828.2, 4.234.6/8, 6.693.1a/1b/2a,
  17.23.14/15, 17.34.21a.
- **The final pass with the v0.6493 table, now the last cascade stage, has not
  been run** (§6, item 2).

Mathilda bugs found by the campaign and not yet fixed:
- **Issue #132:** 46 wrong Integrate results. The families are a 0 returned
  for a nonzero integral, `Piecewise` falling through at a = b, a parameter
  ignored (17.34.21), `ExpIntegralEi[±I…]`, and a spurious `I` in elliptic
  antiderivatives.
- **NSum** on 6.666's series Σ (−1)^(n−1) n^(ν+1) sin(nα) K_ν(nβ) at
  β = 3/8 returns 0.113461 with `NSum::accgl`. The correct value is 0.113176:
  `Method -> "AlternatingSigns"` and direct summation get it. The
  elementary analogues are fine.
- **A partial `Inactive[Integrate][0, 0]`** leaks from some definite stage on
  the 3.162–3.172 elliptic family in N mode (NONNUMERIC). It became visible
  once the table moved last.
- **The corpus runner** cannot check an answer containing `Abs` by
  differentiation, because `Abs'` has no numeric value. It should
  differentiate the `Sqrt[u^2]` form, as the table's self-check does.

## 6. Next steps, in order

1. **Finish the condition audit.**
   - Re-run `condaudit.py --ids` with the current script on the 266 removed
     records plus the 554 with failed quadratures, against a table built
     WITHOUT `condaudit_wrong.tsv`.
   - Merge its `.fix` into `tools/gr/cond_fixes.tsv`, marked as automatic.
   - Replace `condaudit_wrong.tsv` with its `.wrong`.
   - Rebuild, then re-audit the corrected ids. They must come back clean.
   - Hand-derive corrections where the cause is clear, as for 3.197.7;
     otherwise leave the record removed.
2. **Final corpus pass.**
   - Regenerate the corpora with `--book-wrong` and `--cond-fix`. The rule
     applies to the corpus too: a formula wrong as stated is not a valid
     question until corrected.
   - Run all seven chapters with `MATHAR_ONLY` set to the table ids. The gate
     is zero WRONG coming from table answers; a WRONG already WRONG in the
     baseline belongs in #132.
   - Measure the gain against the baseline, and the TIMEOUT cost of running
     last: an integrand the earlier stages grind on to the 25 s budget may
     never reach the table.
3. **Register the corpora.**
   - Copy them to `tests/GR{2,3,4,5,6,7,17}IntegralsCorpus.m`.
   - Add `add_test(NAME gr<ch>_integrals_corpus_tests COMMAND mathar_corpus_tests <file> <S floor> <N floor>)`
     in `tests/CMakeLists.txt` next to the Mathar corpora. A `*corpus*` name
     makes it nightly.
   - Set the floors from step 2 on a quiet machine.
4. **Triage chapters 3–7 and 17.** Bucket the remaining non-CORRECT records
   as for ch2 (CAN-DO / INSTANCE-ONLY / elliptic / risch / symbolic-order /
   BUG; `tools/gr/triage.py`), and add the table to this document.
5. **Fix the bugs in §5:** #132, NSum, the Inactive leak, and the runner's
   `Abs` check.
6. **Convention review of the other `Abs` forms** (§3). For example,
   2.813.1's `Sign[a](x ArcSin[x/Abs[a]] + Sqrt[a^2 - x^2])` is
   `x ArcSin[x/a] + a Sqrt[1 - x^2/a^2]` analytically. This needs a derived
   rewrite per family, not a regex.
7. **Later (user):** G&R integrals as a `benchmarks/` gap experiment, after
   all dev campaigns.

Lessons are in `tasks/lessons.md` ("G&R table" entries). Memory notes:
`feedback_table_entry_wrong_correct_or_remove`, `feedback_antiderivative_log_not_log_abs`.
