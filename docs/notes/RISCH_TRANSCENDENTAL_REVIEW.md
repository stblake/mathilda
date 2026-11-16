# Adversarial review: the C transcendental Risch integrator

Date: 2026-11-16. Binary: `./Mathilda` as built at 20:13 on 2026-11-16, commit `892fb40e8` (v0.6587), working tree as found. The review was read-only. Every probe was run against the existing binary, and every reproducer below was run at least twice.

Oracles used:
- every antiderivative was checked with `PossibleZeroQ[D[F,u]-f]` and numerically at u = 7/10, 13/10, 23/10 (a = 3/7, b = 5/3);
- every verdict was cross-checked against FriCAS 1.3.x (`fricas -nosman`, `unparse(r::InputForm)`): an `integral(...)` in FriCAS's output, or a special function in it (Ei, li, erf, dilog, polylog), means the integral is not elementary;
- `Integrate`ParallelMixedTower` was used as a second certificate.

Terms: **RT** is `Integrate`RischTranscendental` and **EIQ** is `Risch`ElementaryIntegralQ`.

## 1. Summary

| id | severity | one line | minimal reproducer | source |
|----|----------|----------|--------------------|--------|
| F1 | **WRONG** | SolveAlways in the single-extension ansätze also solves for the user's symbolic parameters. It silently assumes `a -> 0` (or `a -> -1`), and the answer is shipped without a differentiation check, out of top-level `Integrate` too | `Integrate[(1 + a + a E^u + u)/((1 + E^u)(1 + u)), u]` gives `u - Log[1 + E^u]`; `Integrate`RischTranscendental[a/(1 + E^u), u]` gives `0` | `risch_singleext.c:405-409,453` (rt_frac_try); same pattern at `:787`, `:1039`; `integrate.c:974-990` |
| F2 | **WRONG** (decision) | EIQ returns `True` for non-elementary integrands. It takes F1's wrong antiderivative as an existence proof | `Risch`ElementaryIntegralQ[1/(1 + E^u) + a E^u/u, u]` gives `True` (FriCAS: needs Ei) | `integrate_risch_transcendental.c:690-699` |
| F3 | **WRONG** (decision) | The same parameter capture inside the log fold of the field decision (`rt_arctan_fold_coeff`) | `Risch`ElementaryIntegralQ[Log[u]/(u + a), u]` gives `True`; also `Log[a+u]/u`, `Log[u]^2/(u+a)`, `Log[u]/(u+1+a)` | `risch_field_integrate.c:1392-1466` (SolveAlways at `:1433`) |
| F4 | **WRONG** (decision) | The field decision accepts the field integrator's Q without verifying it. Limited integration also lets a *higher* tower variable into a lower-field coefficient | `Risch`ElementaryIntegralQ[Log[u] Log[1 + u], u]` gives `True`; also `u^k Log[u] Log[1+u]`, `Log[u] Log[1-u]`, `Log[u] Log[2+u]` | `integrate_risch_transcendental.c:659-660`; `risch_field_integrate.c:1490,1531-1545` |
| F5 | UNDECIDED-IN-SCOPE | Limited integration with a new logarithm (Bronstein §5.8 / §7.2) is a silent decline. This is the root of the known `ArcTan[u]^k` issue, and it also covers `Log[u]^2/(1+u)` with no ArcTan at all | `Risch`ElementaryIntegralQ[Log[u]^2/(1 + u), u]` is undecided (PMT certifies it in 6 ms) | `risch_field_integrate.c:1531-1565` (no branch for a foreign Log); fold restriction at `:1400-1401` |
| F6 | PERF / HANG | The flat SolveAlways ansatz fallback (`rt_log_tower_case`) has no cap on its size. It runs after the deterministic path declines, twice per Integrate call, and becomes a Gröbner-basis computation when a parameter is present | `Integrate`RischTranscendental[TrigToExp[u^5 ArcTan[u]^3], u]`: 16–18 s, undecided (about 175 unknowns). `R21_21` (parametric): more than 30 s | `risch_field_integrate.c:163-213` ("No cap"), SolveAlways `:267` |
| F7 | UNDECIDED-IN-SCOPE + PERF | The tower log part declines when the residue polynomial is irreducible of degree ≥ 3, or quadratic over Q(i): LogToReal declines and there is no RootSum fallback. The decline then costs 17–68 s | `Integrate`RischTranscendental[1/(u (1 + Log[u] + Log[u]^3)), u]`: 44–68 s, unevaluated (`Integrate` answers it in 4 ms) | `intrat.c:2283` via `risch_field_integrate.c:803-824` |
| F8 | UNDECIDED-IN-SCOPE | No Risch DE over a hypertangent lower field ("Gap 3"). The decline raises no `rde_cert_mark`, so the completeness gate does not see it | `Integrate`RischTranscendental[E^Tan[u] Sec[u]^2, u]` is unevaluated | `risch_field_integrate.c:1712-1720` |
| F9 | UNDECIDED-IN-SCOPE | The structure theorem fails on algebraic constants: `Risch`RationalSpan[2, {Sqrt[2]}, {u}]` gives `$Failed`, so the tower is marked dependent and the verdict is UNKNOWN | `Risch`ElementaryIntegralQ[E^(Sqrt[2] u) + E^(2 u)/u, u]` is undecided | `risch_structure.c` (RationalSpan); `risch_tower.c:727-728`; `integrate_risch_transcendental.c:671-672` |
| F10 | UNSOUND-SHORTCUT (PLAUSIBLE) | `g_rt_decision` is a process global and is not scoped around probe integrations inside the Risch DE. A probe's non-elementarity flag can survive an uncertain decline and become a certificate | not reproduced | `integrate_risch_rde.c:695-711` (rc_integrate_in_K), `:887`, `:915`; `integrate_risch_transcendental.c:653-676` |
| F11 | UNSOUND-SHORTCUT | Stage 2 of `Integrate` and the Bronstein `.m` pipeline use RT answers unverified, on the strength of a "correct by construction" claim that F1 falsifies | see F1 | `integrate.c:974-990`; `src/internal/bronstein/LimitedIntegrate.m:34-38`, `Residues.m:201-205` |
| F12 | UNDECIDED-IN-SCOPE (PLAUSIBLE cause) | The answer depends on the input's form. A merged exponent `E^(2u + u^2)` defeats the tower, while the split form integrates | `Integrate`RischTranscendental[(6 E^(2u) + 6 E^(2u + u^2) - 2u)/(1 + E^(u^2)), u]` is unevaluated; `6 E^(2u) - 2u/(1 + E^(u^2))` integrates | tower build (`risch_tower.c`, rt_expand_exp_sums / dependence check) |
| F13 | DEFICIENCY (tooling) | `make check-bronstein-completeness` counts only `rde_cert_mark` sites. The silent declines behind F5, F7 and F8 are invisible to it, so the "13 uncertain paths" figure undercounts the gaps | — | `tools/check_bronstein_completeness.py:35,168-176` |

**Not found:** no FALSE_CERT. Every `False` that EIQ returned in about 400 probes was confirmed non-elementary, by FriCAS or by a textbook family. No crash, and no leak (`leaks --atExit`: 0 leaks, both on a baseline script and on a 30-call RT/EIQ script).

## 2. Findings

### F1: WRONG. SolveAlways solves for the user's parameters, and the answer is shipped unverified

**Reproducers** (CONFIRMED, run twice; `scratchpad/risch_review/confirm.m`):

```
Integrate[(1 + a + a E^u + u)/((1 + E^u) (1 + u)), u]          -> u - Log[1 + E^u]        (residual -a/(1+u))
Integrate[(a + E^u + a E^u + E^u u^2)/((1 + E^u) (1 + u^2)), u]  -> Log[1 + E^u]            (drops a ArcTan[u])
Integrate[(a E^u + a E^(2 u) + u)/((1 + E^u) u), u]             -> u - Log[1 + E^u]        (drops a ExpIntegralEi[u]!)
Integrate[(1 + a + a E^(2u) + u)/((1 + E^(2u)) (1 + u)), u]      -> u - Log[1 + E^(2u)]/2
Integrate`RischTranscendental[a/(1 + E^u), u]                   -> 0
Integrate`RischTranscendental[a + 1/(1 + E^u), u]               -> -Log[1 + E^u]           (a -> -1 assumed)
Integrate`RischTranscendental[a/(u (1 + Log[u]^2)), u]          -> 0
Integrate`RischTranscendental[1/(u Log[u]) + a/(1 + u), u]      -> Log[Log[u]]
Integrate`RischTranscendental[a/(1 + E^u)^2 + b/(1 + u), u]     -> a u + (a + a u)/(...) - a Log[1 + E^u]   (b -> 0)
Integrate`RischTranscendental[a/Log[u] + 1/(u (1 + Log[u])), u] -> Log[1 + Log[u]]         (drops a li(u))
```

Each residual is non-zero, numerically and by PossibleZeroQ. With a numeric value in place of `a` (e.g. `a -> 3`), every one of them is correct.

**Root cause.** `rt_frac_try` (`risch_singleext.c:308`) calls `SolveAlways[residual == 0, {rmT, u}]` (`:405-409`). SolveAlways solves for every symbol that is not in the variable list. That covers the residue unknowns `rmK_i`, and also the user's constants `a`, `b`, `c`. It returns `{{a -> 0, rmK0 -> -1}}`, or for `a + 1/(1+E^u)` the solution `a -> -1`.

The code then:
- checks only that the rules' right-hand sides are free of x and t (`:418-431`);
- substitutes them into `Σ rmK_i Log[g_i]` (`:453`), where `a` does not occur, so the constraint on `a` disappears;
- returns the result with no differentiation check.

The repeated-pole and coupled ansätze (`rt_hermite_try`, SolveAlways at `:787`; `rt_hyperexp_case`, `:1039`) have the same shape. The `b -> 0` reproducer comes from one of these two; I did not isolate which.

`rt_log_tower_case` and `rt_exp_tower_case` do diff-back (`risch_field_integrate.c:322,633`), so there F1 shows only as F6's hang.

**Bronstein.** §5.6 (residue criterion) and §2.5 require the constant field to contain the parameters, C = Q(a). The ansatz instead treats `a` as an unknown to solve for. That is a different problem: integrability for *some* value of a.

**Fix direction.**
- Pass the ansatz unknowns explicitly. Match coefficients in `{t, x}` (CoefficientList) and `Solve` for the `rmK`/`rmW`/`rmHh` symbols only, so every other symbol stays a constant.
- Alternatively, reject any solution that binds a symbol which is not an ansatz unknown.
- Independently, diff-back verify these three cases, as the flat-tower cases already do.

### F2: WRONG (decision). EIQ takes an unverified antiderivative as an existence proof

`rt_decide` (`integrate_risch_transcendental.c:690-699`) returns ELEMENTARY whenever `rt_integrate` returns something with elementary heads. With F1, this gives false "True" verdicts (CONFIRMED; FriCAS: `-log(exp(x)+1)+(a*Ei(x)+x)` and `log(log(x)+1)+a*li(x)`):

```
Risch`ElementaryIntegralQ[1/(1 + E^u) + a E^u/u, u]       -> True   (needs Ei)
Risch`ElementaryIntegralQ[a/Log[u] + 1/(u (1 + Log[u])), u] -> True   (needs li)
```

**Fix.** Gate the existence proof on `PossibleZeroQ[D[anti,x] - f]`. The cost is negligible next to the decision itself.

### F3: WRONG (decision). Parameter capture in the field decision's log fold

```
Risch`ElementaryIntegralQ[Log[u]/(u + a), u]      -> True     FriCAS: integral(log(x)/(x+a),x)
Risch`ElementaryIntegralQ[Log[a + u]/u, u]        -> True     FriCAS: integral(log(x+a)/x,x)
Risch`ElementaryIntegralQ[Log[u]^2/(u + a), u]    -> True     FriCAS: integral(log(x)^2/(x+a),x)
Risch`ElementaryIntegralQ[Log[u]/(u + 1 + a), u]  -> True     FriCAS: integral(log(x)/(x+(a+1)),x)
Risch`ElementaryIntegralQ[Log[u]/(1 + u), u]      -> undecided (the same family without a parameter)
```

RT itself leaves all four unevaluated, so these verdicts come from `rt_decide_field`.

**Root cause.** For `Log[u]/(u+a)`, the tower is {t0 = Log u}. The t0^1 coefficient 1/(u+a) integrates to the foreign `Log[a+u]`. `rt_arctan_fold_coeff` (`risch_field_integrate.c:1392`) then solves `SolveAlways[1/(a+u) - c0/u == 0, {u, t0}]` (`:1433`). SolveAlways solves for `a` as well and returns `{c0 -> 1, a -> 0}`. The fold "succeeds" with c = 1, which is exactly the elementary `Log[u]^2/2` case at a = 0, so the field integrator returns a Q and the verdict is ELEMENTARY (see F4: Q is not verified). `Log[u]/(u+1+a)` gets `a -> -1` the same way.

**Fix.** As for F1: solve for `rtatc$j` only. Also require every returned rule's left-hand side to be one of the unknowns.

### F4: WRONG (decision). An unverified Q in the field decision, plus higher tower variables leaking into lower coefficients

```
Risch`ElementaryIntegralQ[Log[u] Log[1 + u], u]     -> True   FriCAS: integral(log(x)*log(x+1),x)
Risch`ElementaryIntegralQ[u Log[u] Log[1 + u], u]   -> True   FriCAS: integral(x*log(x)*log(x+1),x)
Risch`ElementaryIntegralQ[u^2 Log[u] Log[1 + u], u] -> True   FriCAS: integral(...)
Risch`ElementaryIntegralQ[Log[u] Log[1 - u], u]     -> True   FriCAS: integral(log(-x+1)*log(x),x)
Risch`ElementaryIntegralQ[Log[u] Log[2 + u], u]     -> True
```

All are CONFIRMED. RT is unevaluated on all of them, because the recursive tower case verifies D_tower Q == F at `risch_field_integrate.c:1817-1840` and rejects Q. PMT on `Log[u] Log[1+u]` says "failed: not independently confirmed", not "elementary".

**Root cause.** There are two parts.

1. **Q is accepted unverified.** `rt_decide_field_canon` (`integrate_risch_transcendental.c:659-660`) sets `verdict = RT_DEC_ELEMENTARY` whenever `rt_field_integrate` returns a Q with no remainder. It never checks Q, although the recursive integrate case runs the same integrator and verifies its result.

2. **A higher tower variable leaks into a lower coefficient** (PLAUSIBLE: read from the code; consistent with the probes). In `rt_limited_field_integrate`, the lower integral R is rewritten with `T->subrules` (`:1490`), and those rules alias *every* tower kernel, including t_{L+1..n-1}.

   For F = t0·t1 (t0 = Log u, t1 = Log(1+u)), the bottom coefficient is integrated over K_0. Its t0^1 coefficient -u/(1+u) integrates to `-u + Log[1+u]`. That `Log[1+u]` becomes t1, a variable *above* the current level, so the result no longer looks "foreign" (`:1502-1503`). The branch at `:1531-1545` then accepts `s = -u + t1` as an element of K_0. Only the coefficient c is checked for tower variables (`:1534-1536`); s is never checked. The deeper recursion then treats t1 as a constant.

   Bronstein §5.8 (IntegratePrimitivePolynomial) needs s ∈ k (= K_{L-1}). Here a new logarithm that is not t_L sits at an i ≥ 1 coefficient, which proves the integral is not elementary.

   A single-kernel control behaves as this predicts: `Log[u]^2/(1+u)` has no higher kernel to alias to, and it stays undecided (F5), not True.

**Fix.**
- Verify Q in the decision path, reusing the exact tower check at `risch_field_integrate.c:1817-1838`.
- In `rt_limited_field_integrate`, substitute only kernels t_0..t_L, or reject `s` and `c` that depend on t_{>L}. Treated as foreign logs, these then fall into F5's (correct) limited-integration branch.

### F5: UNDECIDED-IN-SCOPE. Limited integration with a new logarithm is a silent decline (the known ArcTan issue, generalised)

These are confirmed non-elementary by PMT or FriCAS, and EIQ leaves all of them undecided:

```
Log[u]^2/(1+u), Log[u]^3/(1+u), Log[u]/(1+u), Log[1+u]/u, Log[1+u^2]/u, u Log[u]/(1+u^2),
Log[u+I]/u, Log[u]/(u-Sqrt[2]), ArcTan[u]/u, ArcTan[u]^2, ArcTan[u]^3, u^2 ArcTan[u]^3, u^5 ArcTan[u]^3,
TrigToExp[ArcTan[u]^2] (two logs over Q(i))
```

PMT certifies `Log[u]^2/(1+u)` in 6 ms, `Log[1+u]/u` in 5 ms, and `ArcTan[u]^2` and `u^5 ArcTan[u]^3` in under 0.05 s.

These elementary inputs are also undecided:

```
-2 Log[Cos[u]] Tan[u]   (antiderivative Log[Cos[u]]^2)
D[ArcTan[E^u]^2, u]
```

**Decline site.** This was verified with `MATHILDA_RDE_CERT_DEBUG=1`: no `rde_cert_mark` fires for `ArcTan[u]^2` or for `Log[u]^2/(1+u)`. The decline is the plain `rc = -1` of `rt_limited_field_integrate` (`risch_field_integrate.c:1469`).

For a θ^i coefficient with i ≥ 1 whose integral carries a new Log, the logic runs:
- the fold is tried first (`:1510-1529`);
- then `if (!foreign && …)` (`:1531`);
- then `else if (rt_has_arctan_of_x)` (`:1549`).

There is no branch for "foreign Log, fold failed". The comment at `:1553-1555` says it is left as a plain decline on purpose: it "can be a content / combination artifact".

Two more restrictions:
- **The fold is restricted to logs in x alone.** `rt_arctan_fold_coeff` returns NULL as soon as the new-log part A contains any tower variable (`:1400-1401`). That makes `Log[1 - I E^u]`-type new logs undecidable; this is the cause for `ArcTan[E^u]^2` and `Log[Cos[u]]^2`.
- **Only the syntactic spelling is folded.** The new-log rewrite at `:1484-1488` folds only `Log[u_L in t-vars]`. Any other spelling (`Log[u + I]` against the kernel `Log[1 - I u]`) has to go through the fold.

For `ArcTan[u]^k` the decision runs over the complex tower `{Log[1 - I u], Log[1 + I u]}` from TrigToExp. There the "new" logarithms are spelled as `Log[1 + u^2]` / `Log[u ± I]`, and some of them are also caught by F4's aliasing. Exactly which coefficient dies is PLAUSIBLE, not traced (lldb/sample cannot unwind through the evaluator's stack; see F6). The class is the same: `Log[u]^2/(1+u)`, with no ArcTan and no Q(i), fails identically.

**Bronstein.** This deviates from §5.8 IntegratePrimitivePolynomial and §7.2 LimitedIntegrate. When `∫p_i ∉ c·t + k`, the result is "not elementary". That needs the parametric version of the residue test (are the new logs' residues in the span of t's?), or equivalently `Risch`LogReducible` on the new-log arguments over the tower.

**Fix direction.** Replace the fold with LimitedIntegrate (Phase 2 item 4/5 of the plan). As an interim, decidable step:
1. write R = r_K + Σ c_k Log w_k (w_k over the tower);
2. call `Risch`LogReducible` on each w_k;
3. if any w_k is independent of {t_0..t_L}, prove non-elementary;
4. otherwise fold.

Allow the fold to accept A with tower variables, by computing `D_tower[A]` instead of `D[A,x]`.

### F6: PERF / HANG. The flat ansatz fallback is uncapped and runs after the exact path has declined

**Profile.** lldb backtraces, 16 of them, 1 s apart, on `TrigToExp[u^5 ArcTan[u]^3]` (`perf1.m`): every sample was in

```
rt_solve_always_quiet <- rt_log_tower_case <- rt_transcendental_case <- rt_integrate <- builtin_rischtranscendental
```

(and later the same under `rt_decide <- builtin_elementaryintegralq`), with Together/Cancel inside SolveAlways at the leaves.

**Size of the system.** `rt_log_tower_case` (`risch_field_integrate.c:83`) builds a dense ansatz with `nunk = (Ntop+1)·Π(bd_j+1) + ng` unknowns (`:163-213`), where each bound is `bd_j = deg num + deg den + 1` and the code says "no cap". For `u^5 t0^3`-type inputs that is 5·(7·5) = 175 unknowns. Timings:
- 16.3 s in RT (18 s on a rerun), undecided;
- 1.2 s for `TrigToExp[u^2 ArcTan[u]^3]`;
- about 0.01 s for the `ArcTan` spelling, which this case does not see.

**Wasted work.** The case runs only after the exact recursion (`rt_recursive_tower_case`) has declined. Its decline proves nothing, so on a non-elementary integrand it always spends the maximum time. `builtin_rischtranscendental` and `rt_decide` each call `rt_integrate`, so `Integrate` + `Risch`ElementaryIntegralQ` pay the cost twice.

**With a parameter it becomes a Gröbner basis.** On `R21_21` (`u^2/(1 + Log[1 + E^E^u]^2) + Log[1 + Log[u]^2]^2/(-a + u)`, non-elementary per FriCAS), RT and EIQ each hit their 30 s TimeConstrained. Backtraces show:

```
builtin_groebner_basis <- builtin_solve <- builtin_solvealways <- rt_log_tower_case
```

SolveAlways treats `a` as an unknown (F1), which makes the system nonlinear.

**Fix.**
- Cap `nunk`, or skip the flat ansatz whenever the recursive case declined for a reason the ansatz cannot overcome (a foreign-log decline, F5).
- Pass the unknowns explicitly so the system stays linear (F1).
- Cache the `rt_integrate` result between RT and the decision.

### F7: UNDECIDED-IN-SCOPE + PERF. Tower log part with an irreducible residue polynomial of degree ≥ 3, or quadratic over Q(i)

```
Integrate`RischTranscendental[1/(u (1 + Log[u] + Log[u]^3)), u]          -> unevaluated after 44-68 s
Integrate`RischTranscendental[(1 + Log[u]^2)/(u (1 - 3 Log[u] + Log[u]^3)), u] -> undecided, 17 s
Risch`ElementaryIntegralQ[1/(u (I + Log[u]^2)), u]                       -> undecided
```

FriCAS integrates all three (rootOf / `%i^(1/2)` logs). `Integrate` (Automatic) answers the first in 3.6 ms with a RootSum. The exponential analog `E^u/(1 + E^u + E^(3u))` takes 2 ms, because `rt_exp_ratreduce_case` kernelizes it to a rational integral; there is no log analog of that case.

**Root cause.** `Integrate`TranscendentalLogPart[1, 1 + t + t^3, t, z, 1 + 3t^2, {u}]` returns unevaluated in about 1 ms. `intrat_transcendental_log_part` declines when LogToReal cannot convert a residue factor (`intrat.c:2283`, `/* LogToReal declined */`), and there is no RootSum fallback. `rt_field_lrt_logpart` (`risch_field_integrate.c:817-824`) turns that into a plain NULL.

Bronstein §5.6 allows residues in the algebraic closure: a `RootSum[Q(z), z Log[S(z,t)]]` is a valid elementary log part. The 44–68 s are spent in later fallbacks, which the stack switch hid from the profiler (lldb unwinding stops at `builtin_together_compute`).

**Fix.** Return the RootSum form when LogToReal declines, with a real form attempted afterwards. Add a log analog of `rt_exp_ratreduce_case`: when F = g(t)·Dt with g ∈ C(t), substitute and integrate the rational g.

### F8: UNDECIDED-IN-SCOPE. No Risch DE over a hypertangent lower field

`Integrate`RischTranscendental[E^Tan[u] Sec[u]^2, u]` is unevaluated (antiderivative `E^Tan[u]`), and EIQ is undecided.

The tower is {Tan u (RT_TAN), E^Tan u}. The RDE `q' + Dcoef q = p` over C(u)(tan u) reaches `rt_field_rde`'s Gap-3 return NULL (`risch_field_integrate.c:1712-1720`). That is the correct thing for soundness (see the A4 note in the comment), but it raises **no `rde_cert_mark`**, so the completeness gate's 13 does not count it (F13).

**Bronstein.** The missing pieces are the §6.2 special denominator for a hypertangent t (`RdeSpecialDenomTan`, the `t^2+1` part) and the §6.3 bound for a nonlinear t (δ(t) = 2).

### F9: UNDECIDED-IN-SCOPE. Algebraic constants in the structure theorem

`Risch`RationalSpan[2, {Sqrt[2]}, {u}]` and `Risch`RationalSpan[Sqrt[2], {2}, {u}]` both return `$Failed`; the correct answer is `False`, because √2 ∉ Q. The tower builder turns that into `T->dependent = true` (`risch_tower.c:727-728`), and the decision demotes the verdict to UNKNOWN (`integrate_risch_transcendental.c:671-672`).

Effects:
- `E^(Sqrt[2] u) + E^(2u)/u` and `E^(Sqrt[2] u)(2 + u) + E^(2u)/(u - a)` are undecided. Both are non-elementary (FriCAS: Ei). With 3 in place of √2, both are False.
- `E^u u^Sqrt[2]/(1 + u)` is undecided.

**Fix.** Decide Q-linear independence over Q(√2…), by writing coefficients in a Q-basis of the number field (the Phase 4 CF record).

### F10: UNSOUND-SHORTCUT (PLAUSIBLE, not reproduced). The decision flag is not scoped around probe integrations

`rt_dec_nonelem()` writes the global `g_rt_decision` (`risch_field_integrate.c:696-698`) whenever `g_rt_decide_mode` is on. The Risch DE calls the field integrator as a *probe*: `rc_integrate_in_K` (`integrate_risch_rde.c:695-711`), used by `rc_limited_integrate_prim` in `rde_bound_degree_prim` (`:887`, `:915`) and by `rde_cancel_logderiv` (`:1098`).

That call is wrapped in an `rde_cert` scope, but `g_rt_decision` is not saved and restored. A probe whose own integrand is non-elementary therefore sets the flag. The bound step then just marks "limited integration undecided", and `rt_field_rde` correctly declines to raise the certificate (`:1703-1708`), but the flag is already NONELEMENTARY. If the outer integration later fails for any reason, `rt_decide_field_canon` (`integrate_risch_transcendental.c:662`) returns NONELEMENTARY: a false certificate.

Reaching this needs a primitive-over-primitive RDE with deg b = deg a − 1 and an alpha whose integral over K_{m-1} is not elementary. I did not construct one; no false certificate appeared in about 400 probes.

**Fix.** Save and restore `g_rt_decision` around `rc_integrate_in_K`, or make every probe run with `g_rt_decide_mode = false`.

### F11: UNSOUND-SHORTCUT. Consumers trust RT answers as "correct by construction"

- **`integrate.c`.** The stage-2 comment says RT "fires only behind an exact certificate, so it emits no wrong closed forms and needs no differentiation check" (`integrate.c:974-981`). `try_rischtranscendental` (`:982-990`) returns the result unchecked. F1's four top-level `Integrate` wrong answers come through this path. (Inputs that `Integrate` splits by linearity before stage 2 come out right, e.g. `Integrate[a/(1 + E^u), u]`.)
- **Builtin comment.** `builtin_rischtranscendental`'s own comment says "no differentiation check is applied" (`integrate_risch_transcendental.c:737-739`).
- **Bronstein `.m` pipeline.** `brIntegrateInK` and `brIntegrateInKLogLike` (`LimitedIntegrate.m:34-38`, `Residues.m:201-205`) use `G` unverified, so a wrong RT answer on a parametric integrand propagates into limited integration and residues of the algebraic pipeline. That impact is PLAUSIBLE; I did not construct an algebraic reproducer.

**Fix.** Verify at the RT builtin boundary (one PossibleZeroQ), or fix F1 and re-audit every SolveAlways site for parameters.

### F12: UNDECIDED-IN-SCOPE. Form dependence on a merged exponent (cause PLAUSIBLE)

`Together[D[3 E^(2u) + Log[1 + E^(-u^2)], u]]` is `(6 E^(2u) + 6 E^(2u + u^2) - 2u)/(1 + E^(u^2))`. RT leaves it unevaluated and EIQ is undecided. The split spelling `6 E^(2u) - 2u/(1 + E^(u^2))` integrates in 0.0 s.

The likely cause is a third kernel `E^(2u + u^2)` beside E^(2u) and E^(u^2): either `rt_expand_exp_sums` does not split it inside the quotient, or the dependence check (`risch_tower.c:713-736`) then marks the tower dependent. Not traced further.

**Fix.** Rewrite dependent towers (plan item "Rewrite dependent towers") rather than declining.

### F13: Completeness-gate blind spot

`tools/check_bronstein_completeness.py` scans `integrate_risch_rde.c` and `risch_field_integrate.c` for `rde_cert_mark("…")` only (`:35`, `:168-176`). The declines behind F5 (`rc = -1`), F7 (LogToReal NULL) and F8 (Gap 3) are silent NULL returns, so the "13 uncertain Risch DE paths" figure does not measure them.

**Fix.** Give each of these declines a mark, even outside a cert scope (the mark is a no-op there), so the ratchet counts them.

## 3. Coverage: cases probed and found clean

These areas showed no wrong answer and no false certificate. A case listed as "decided" was decided correctly, by the oracles above.

- **Rational base case** (`Integrate`BronsteinRational`), and Hermite reduction over the tower for log tops, exp tops and nested levels: `H1-H16` in `batB.m`. All 16 integrated and verified, including repeated poles in Log, E^u, E^(u^2), E^(E^u), Log Log, and mixed `E^u + Log[u]` denominators.
- **Rothstein-Trager / LRT with algebraic residues**:
  - quadratic and quartic over Q, with ArcTan real forms: `R1-R12`, `R16-R18` in `batB.m`, including `Log[u]^4+1`, `E^(4u)+1` and `E^(3u)-2`;
  - parametric residues (`Log[u]^2 - a`, `a^2 + Log[u]^2`, `E^(2u) - a`);
  - Q(i) residues in the exponential case (`E^u/(I + E^(2u))`).
- **Residue criterion with non-constant residues** (Thm 5.6.1(ii)): `1/Log[u]`, `1/(1+Log[u])`, `1/(u+Log[u])`, `E^u/(u Log[u])`, and `u/(1+E^u)^2` are all decided False.
- **Primitive polynomial (§5.8) without new logs**: `Log[u]^k/u`, `Log[Log[u]]/u`, depth-3 `Log[Log[Log[u]]]/(u Log[u])`, and the Dc ≠ 0 certificate (`Log[Log[u]]`, `Log[u]^2 Log[1+u]`, `u^3 Log[u]^2 Log[1+u]`) are all correct.
- **Real ArcTan obstruction** (dilogarithm, real tower): `Log[u]/(1+u^2)` and `Log[1+u^2]/(1+u^2)` are False; `ArcTan[u]/(1+u^2)` is True.
- **Hyperexponential polynomial and RDE (Ch. 6)**: `E^u/u`, `E^(u^2)`, `E^(1/u)`, `E^(E^u)`, `E^(E^(E^u))`, `u^u`, `E^(Log[u]^2)`, `E^(u^3)`, `E^(u - u^2)`, and the elementary `E^u(1/u - 1/u^2)`, `E^(u^2)(2u^2+1)`, `E^(u/(1+u))/(1+u)^2`.
- **Parametric RDE with a in exponents and powers**: `u^a E^u`, `E^(a u)/(u+b)`, `E^(a u^2)`, `u^a Log[u]`, `E^(a/u)/(1+u)`, and the elementary derivatives `Q1-Q20` in `p12.m`.
- **Hypertangent (§5.10) single monomial**: `Tan[u]^2`, `Tan[u]/(1+Tan[u])`, `u Tan[u]` (False), `Tan[u]/u` (False), `Sec[u]`, `a/(2 + Cos[u])`.
- **Structure theorem / dependent kernels**: `Log[u^2]/u`, `Log[2u]/u`, `E^(u+u^2)(1+2u)`, `E^(Log[u]^2)` against `Log[u]`, `Log[u] Log[2u]/u`.
- **Transcendental constants** (E, Pi, EulerGamma, Log[2]): all elementary cases in `p13.m` (`C1-C11`) integrated. `rt_has_transcendental_constant` correctly blocks certificates there.
- **Q(i) logs**: `Log[u+I] Log[u-I]`, `Log[1+I u] Log[1-I u]`, `Log[(u-I)/(u+I)]^2`, `E^(I u)/(u+I)`.
- **Randomised D[F] for elementary F** (seeds 1, 2, 11; 190 cases): one wrong answer (F1, `G1_29`); otherwise every answer verified, and no false certificate.
- **Randomised integrands against FriCAS** (seed 21, 80 cases): every False (56) and every True (4) agree with FriCAS; 7 undecided and 1 timeout, all non-elementary per FriCAS.
- **Memory**: `leaks --atExit` reports 0 leaks on a 30-call RT+EIQ script (`leak1.m`).

Probed only lightly:
- special-function outputs (`RT_SF_*`, Cherry), seen only as Ei/Erf/li answers;
- the Knowles kernel-radical path;
- the trig front end with non-tangent kernels;
- tower depth ≥ 4.

## 4. Appendix: scripts and tallies

All scripts are in `/private/tmp/claude-501/-Users-user-Mathilda-dev-Mathilda/e3830d4b-e34d-4b5f-a4ea-f39cb06ea566/scratchpad/risch_review/`.

The harness, `harness.m`, defines `chk[id, f, expect]`. It runs RT and EIQ under `TimeConstrained[…, 30]` with `AbsoluteTiming`, then applies PossibleZeroQ plus the 3-point numeric check, and tags each case WRONG? / FALSE_CERT / FALSE_TRUE / UNDEC / TIMEOUT / ok. An "INCONSISTENT" tag only means that RT returned a special-function answer (Ei, Erf) beside EIQ False; it is not a bug. Early "nonnum" tags came from the unbound parameter, and were fixed by substituting a = 3/7 and b = 5/3.

| script | content | tally |
|---|---|---|
| `batA.m` | 40 known non-elementary + 40 elementary hand cases | FALSE_TRUE 2 (`Log[a+u]/u`, `Log[u] Log[1+u]`); undecided non-elementary 9; undecided elementary 1 (`E^Tan[u] Sec[u]^2`) |
| `gen.py` → `genE1.m`, `genE2.m` (seeds 1, 2) | 120 × D[F], F random over 20 kernels | 98 verified (21 of them after the harness fix), 1 WRONG (F1), 7 undecided |
| `gen2.py E 11` → `genH11.m` | 70 × D[F] with params / Q(i) / √2 / deeper kernels | 70 ok |
| `gen2.py R 21` → `genR21.m`, `casesR21.tsv`, `fricas_run.py` → `fricasR21.tsv` | 80 random integrands, FriCAS ground truth | 56 False and 4 True, all agree; 7 undecided + 1 timeout, all non-elementary |
| `batB.m` | 18 LRT/algebraic-residue cases + 16 Hermite D[F] cases | 31 ok; `R13`, `R15` undecided; `R14` timeout (F7) |
| `p5.m`, `p6.m`, `p7.m`, `p8.m`, `p10.m`, `p11.m`, `p16.m`, `p17.m` | parameter-capture family (F1–F2) | 22 wrong RT answers, 4 wrong top-level `Integrate` answers, 2 false True |
| `p2.m`, `p20.m` | two-log and parametric decision family (F3–F4) | 11 false True |
| `p12.m`, `p13.m`, `p14.m`, `p15.m`, `p21.m` | parametric RDE, constants, Q(i), √2 span, merged exponent | F8, F9, F12 |
| `perf1.m` … `perf4.m`, `p18.m`, `p19.m` | profiling (lldb `bt` sampling) and timing | F6, F7 |
| `casesX.tsv`, `casesY.tsv` → `fricasX.tsv`, `fricasY.tsv` | FriCAS confirmation of every reported verdict | as quoted above |
| `confirm.m` | one-shot rerun of every reproducer in §2 | all reproduced |

The two generators are seeded with Python's `random.seed(<seed>)`, so rerunning them with the same arguments regenerates the same cases.
