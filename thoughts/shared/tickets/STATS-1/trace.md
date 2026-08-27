---
type: trace
ticket: STATS-1
date: 2026-08-27
flow: qrispy
traced: 0
partial: 0
gaps: 0
basis: not-applicable
---

# Spec-check (QRISPY S) — STATS-1

## trace.py report (verbatim)

```
NOT APPLICABLE — no prd.md under the tickets tree. There are no product requirements to
trace; this says nothing about test coverage (the verification ladder owns that).
```

Command: `python3 skills/spec-as-test/scripts/trace.py --ticket STATS-1` (exit 0).

This is the documented outcome for an engineering-only ticket (commands/rpi-loop.md:
"S on an engineering-only ticket ... `trace.py` reports NOT APPLICABLE, and that is the
honest outcome, not a failed step"). `gaps: 0` here means "nothing to trace", NOT
"everything traced" — recorded that way deliberately.

## The join that DOES exist: plan Acceptance Criteria -> test functions

The plan's AC table is the requirement side; tests/test_quantile_family.c is the proof
side, one assert per row. Manual join (verified by reading both):

| AC | Test function | Status |
|---|---|---|
| AC-1 | test_quantile_default_half | covered |
| AC-2 | test_quantile_default_quarter | covered |
| AC-3 | test_quantile_q_list | covered |
| AC-4a/4b | test_quantile_edges | covered |
| AC-5 | test_quantile_explicit_params | covered |
| AC-6 | test_quantile_reals | covered |
| AC-7 | test_quantile_unsorted | covered |
| AC-8 | test_iqr_basic | covered |
| AC-9 | test_mean_deviation | covered |
| AC-10 | test_median_deviation | covered |
| AC-11a/11b | test_quantile_symbolic_declines | covered |
| AC-12 | test_quartiles_regression | covered |
| AC-13 | test_quantile_matrix_columnwise | covered |
| AC-14 | test_mean_deviation (second assert) | covered |
| AC-15 | test_quantile_empty_declines | covered |
| AC-16 | test_quantile_out_of_range_q | covered |
| AC-17 | test_iqr_matrix_columnwise | covered |
| AC-18 | test_quantile_ndarray_exact | covered |
| AC-19 | test_meandeviation_ndarray | covered |

19/19 acceptance criteria have a named test. Extra tests beyond the AC table:
test_quantile_singleton, test_deviation_declines, test_quartiles_param_integer_h (pins
the deliberate semantic change).

**gaps: 0** on this join as well — asserted from reading both artifacts, not from a
tool, since no tool joins a plan AC table to C test functions in this repo.
