# SIAM 100-digit holdout A (experiments 196-225)

Ten generalised problems per SIAM problem (C1-C10, D1-D20), each with a
`machine` and a `wp50` row: 600 cases. Every truth was established to at least
55 digits by two independent methods; the reference scripts are in
`siam100_refs/hc*_ref.*` and `hd*_ref.*`, and each folder's `truth.tsv` names
its provenance. Methodology, scoring and the target are those of
[`SIAM100.md`](SIAM100.md).

## The pre-fix run (out of sample)

Run once, on v0.6377, **before any fix of the campaign**, all five systems,
scored by `siam100_table.py`. This is the honest out-of-sample number; after it
holdout A is training data for the fixes. Sealed holdout B (226-255, branch
`siam100-holdout-b`) plays that role from here on.

| system | correct | inexact | wrong | none |
|---|---|---|---|---|
| Mathilda | 252 | 117 | 30 | 201 |
| WL | 257 | 179 | 65 | 99 |
| scipy / mpmath | 330 | 58 | 44 | 168 |
| Maxima | 194 | 23 | 5 | 378 |
| FriCAS | 113 | 14 | 12 | 461 |

Mathilda on target (digits >= best other system and time < 1.5x the fastest
correct one): **262 / 600 (44%)**. Weakest families: D10, D11, D12, D20 (1/20),
C8, D18 (2/20), D7 (3/20), D19 (4/20).

Caveat: the machine was loaded during the run (holdout-B reference jobs and
fix-branch builds, load average peaking above 30 on 16 cores), so the DIGITS
are exact but the SLOWER verdicts are approximate. To re-time the pre-fix
binary on a quiet machine: `git worktree add <dir> v0.6377`, build there and
rerun with `--systems mathilda --from` the archived JSONs.

Files: `siam100_results/holdoutA_v0.6377.md` (full per-case table),
`siam100_results/raw_v0.6377.tar.xz` (the raw per-system JSONs of the
baseline 166-195 and of this run, for `siam100_table.py --from`).

The per-family findings (the fix queue) are in `tasks/todo.md`, section
"SIAM 100-Digit Challenge".
