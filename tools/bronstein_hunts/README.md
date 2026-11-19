# Bronstein / Risch false-certificate hunts

Seeded hunts for the C transcendental Risch decision (`Risch`ElementaryIntegralQ`,
`Risch`RischDECertified`). They were scratch files until 2026-11-19, when a machine
crash lost them and they had to be recovered from a session transcript.

Each hunt builds `F` at random in a tower and asks about `D[F]`, whose integral is
elementary by construction. A `False` (or `"NoSolution"`) is then a **false
certificate**, the one outcome that must never happen.

| Hunt | Field |
|---|---|
| `hunt.m`, `hunt2.m` (second seed) | Q, Q(i), Q(a), exponential, nested logarithms |
| `hunt_exp.m` | exponential towers |
| `f8hunt.m` | one hypertangent level, plus random integrands (ground truth unknown) |
| `f9hunt.m`, `f9hunt2.m` | numeric algebraic constants (`Sqrt[2]`, `2^(1/3)`) |
| `f12hunt.m` | composite exponential classes |
| `dephunt.m` | dependent logarithms (`Log[-u]` beside `Log[u]`) |
| `rdehunt.m` | `Risch`RischDECertified` on `g = D[y] + f y` |
| `primhunt.m`, `primhunt2.m` | Liouvillian primitives (Erf, Erfi, Erfc, Ei, li) |

The primhunt oracle: a `False` is wrong only when every Liouvillian kernel of `F`
also occurs in `D[F]`. `D[2 ExpIntegralEi[u]] = 2 E^u/u`, so `False` is right there.

`run.sh [binary] [outdir]` runs them all in parallel (about 10 minutes), with each
call wrapped to print `CASE` / `VERDICT` lines around the `rde_cert_mark` trace. For
each hunt it prints the verdict counts. `solo.py` groups the undecided cases by their
exact set of marks, and `blocking.py` ranks the marks. A group shown as
`(NO MARK: silent decline)` is undecided with no mark at all, which
`tools/check_bronstein_completeness.py` cannot see.

Run the A/B for a campaign step against a scratch build of the previous tag
(`git worktree add`), and compare the two `run.sh` outputs.
