# Proposal: making scoping substitution cheap (Module / With / rule application / Function)

Status: **implemented in v0.6653–v0.6657**; the design as built, and three corrections to the prototype (the epoch wrap, the capture-rename skip, root frozenness), are in `docs/design/scoping_substitution.md`. Original status: proposal, for decision. Nothing in the repository was changed. All measurements
were made on 2026-11-21 against `main` at `78f8ba58a` (v0.6650), on the shared development
Mac (16 logical CPUs, x86-64, `gcc-16 -O3`). Every number below is labelled **measured**,
**simulated** (computed exactly by instrumentation over a real run, but for a mechanism that
was not built), or **estimate**. Section 12 lists what was not measured.

The artefacts referred to throughout (scripts, raw profiles, the prototype diff) are in
`docs/design/eval_change/`, abbreviated `$S` below: `$S/bench/` holds the scripts and
raw results, and `$S/prototype.diff` is the prototype's full diff against `78f8ba58a`. The
instrumented scratch build tree (`$S/inst` in the text) was not kept; rebuild it by
applying the prototype diff.

---

## 1. Executive summary

- **What the profile shows.** On the slow integration cases that run the Trager and Bronstein
  `.m` packages, scoping substitution is real but bounded. I instrumented every scoping walk
  and timed it exactly. Scoping-related work is **17% of the Trager case** (`Method -> "Trager"`
  on `((1+x^2)(2+x^4)^(1/6))/x`) and **18–20% of the `Log`-tower case** (`cM-1-292`). It is
  **6–11% on lccontent, logxy, texpp and the M-class cases**, and **0–1% on the rest of the slow
  tail** (`trager-ratradicand-007`, `cM-3-21`, `mperturb-2-10`). Those last three spend their
  time in algebra or are capped by wall-clock budgets. So even a zero-cost scoping mechanism
  would make the worst `.m`-bound cases at most about **1.2–1.25× faster**. The 43-case tail
  will not become 2× faster from this change. `substitute_scoping` tops the self-time list
  only because the rest of the time is spread over hundreds of functions.
- **Why the walks are slow.** The cost is full walks, not allocation. Each `Module`/`With`
  call walks its whole body. Each rule application walks its whole right-hand side.
  Unchanged subtrees are already *shared*: v0.6057 (`src/modular.c:196-210`,
  `src/match.c:3090-3112`). But they are still *visited*. In the Trager case, **97% of the
  101 million nodes `substitute_scoping` visits come back unchanged**. **59% of all those
  visits come from one function**, `goodPrime` (`src/internal/trager/TorsionModP.m:57`). Its
  pattern variable `ctx` is bound to a large context `Association`. `rb_rec` injects that
  value into the `Module`, and `substitute_scoping` then walks it again on every call: 29,449
  nodes per call, 2,024 calls. After capture avoidance an injected value cannot contain the
  construct's locals, so all of that walking finds nothing.
- **The prototype.** I built it in a scratch copy, not the repository. It is a 64-bit
  per-node free-symbol summary that is filled lazily on shared ("frozen") nodes. Bits are
  assigned only to *binder names*, meaning the names that scoping constructs and patterns
  bind, with an epoch for late registrations. Results:
  - The Trager case's `substitute_scoping` visits drop from **101.1 M to 10.1 M (−90%)**.
  - Its instructions retired drop from **78.13 G to 65.11 G (−16.7%)**, and wall time on a
    quiet machine from **15.0 s to 12.8 s**.
  - `cM-1-292` drops **26.44 G → 22.47 G (−15%)**.
  - The micro-benchmarks gain 10–200×. A rule whose right-hand side holds a 10,000-leaf
    branch that is never taken drops from 232 µs to 1.1 µs per call.
  - Twenty-seven scoping-semantics probes give identical output with the fast paths on and off,
    and every corpus answer is unchanged.
- **Recommendation**, in this order:
  1. **(D)** Stop bumping the evaluation clock for `$ModuleNumber` and for a fresh local's
     first binding. This is small and independent, and worth −3.5% instructions on Trager
     and −10% on `cM-1-292`.
  2. **(A′)** A per-node free-symbol summary that is lazily filled, binder-indexed and
     epoch-validated. Store it **in the existing `data.function.index` slot as a tagged
     word**. That slot has seven uses in two files. Storing it there keeps `sizeof(Expr)`
     at 56 bytes. Then use it in `substitute_scoping`, `rb_rec`, the capture-avoidance
     "danger set" and the pure-function walkers.
  3. Option (b), a pointer-keyed side cache, is fragile: the node pool recycles addresses
     at once. Reject it. Moving `.m` stages to C, option (c), does not address this cost:
     the cost is spread across every `Module` in two packages. Keep (c) for algebra
     hotspots, as project policy already says.
  4. Reject lazy or dynamic renaming: it breaks escaping closures and recursion.
  5. Defer fused rule-and-`Module` substitution and per-rule substitution plans. Once (A′)
     exists they add little, and plans cannot see injected values, which are the dominant
     cost.

---

## 2. Problem statement and measurements

### 2.1 The reported symptom

Your `sample` of the Trager case put `substitute_scoping.part.0` first by self time (82),
ahead of `_nanov2_free` (53), `evaluate_step_core` (47) and `nanov2_malloc_type` (38). An
earlier sample of `2/((1+2x) Log[1+2x] (1+Log[1+2x])^(1/6))` showed the same. The
2026-11-16 changelog records it as a known gap: `docs/spec/changelog/2026-11-16.md:194-198`
says "Skipping subtrees that hold no local needs a per-node symbol summary in `Expr`, a core
layout change left for a decision." I reproduced the profile: Appendix A.1, 10,759 samples
over 20 s. Self time was `substitute_scoping` 1045, `evaluate_step_core` 632, `_nanov2_free`
594, `nanov2_malloc_type` 453, `expr_free` 337, `rb_rec` 294, `node_intern_ex` 164,
`rb_collect_symbol_names` 125. That puts `substitute_scoping` at 9.7% self and 11.8%
inclusive.

### 2.2 Measurement method, and why wall time was not used

The machine was shared, with load averages between 12 and 43 on 16 CPUs while this work ran.
Under that load:
- the same Trager run took between 14.7 s and 62.8 s of wall time;
- three cases (`loglog`, `cM-3-152`, `lccontent-*`) **gave different answers on different
  runs**, because the Integrate cascade's stages are wall-clock budgeted
  (`TimeConstrained`).

The run-to-run variation in instruction counts makes this visible. `cM-3-152` retired 23 G
instructions under load and 36 G on a quiet machine, with the same answer hash. So I used
four instruments:

1. **Instructions retired** (`/usr/bin/time -l`). This is deterministic for an unchanged
   evaluation path; the answer hash (`Hash[r]`) shows whether the path was unchanged. It is
   the primary A/B metric. Startup (`Print[1]`) costs 2.96 G instructions; subtract it when
   quoting relative gains on the integration part.
2. **An instrumented scratch build** that counts every node each walker visits and times
   each scoping phase (Appendix C). This attributes time exactly rather than by sampling.
   Its counters cost about 3% (78.1 G instrumented versus 75.5 G for the stock binary on the
   Trager case).
3. **CPU time** (`TimeUsed[]`), best of three, for the micro-benchmarks.
4. **`sample`**, for the qualitative picture. Under load, and with deep recursion, its call
   graphs are truncated at the root, so I use only leaf-near inclusive numbers.

### 2.3 How much of each slow case is scoping (measured, instrumented build)

The scoping share is the sum of five timed phases:
- **subst**: the `substitute_scoping` walks (Module, With, rename);
- **rb**: the `replace_bindings` walks (all rule applications, including the danger set and
  capture renames);
- **prologue**: Module setup, i.e. `$ModuleNumber` and locals, without initialiser
  evaluation;
- **epilogue**: env free and local reclamation;
- **free**: freeing the substituted body.

| case | CPU (s) | Module calls | subst nodes | subst | rb | prologue | epilogue | free | **scoping share** |
|---|---|---|---|---|---|---|---|---|---|
| Trager `((1+x^2)(2+x^4)^(1/6))/x` | 14.66 | 117,823 | 101.1 M | 1.392 | 0.706 | 0.200 | 0.140 | 0.068 | **17.1%** |
| `loglog` (Automatic) | 4.42 | 96,175 | 18.9 M | 0.250 | 0.317 | 0.140 | 0.076 | 0.037 | **18.5%** |
| `cM-1-292` (Bronstein) ‡ | 9.53 | 96,175 | 18.9 M | 0.966 | 1.276 | 0.658 | 0.246 | 0.174 | 34.8% ‡ |
| `cM-3-152` (Bronstein) ‡ | 7.49 | 13,839 | 4.1 M | 0.251 | 0.250 | 0.084 | 0.235 | 0.018 | 11.2% |
| `logxy-3-6` (Bronstein) ‡ | 6.38 | 9,622 | 3.4 M | 0.184 | 0.180 | 0.064 | 0.107 | 0.035 | 8.9% |
| `lccontent-2-3` (Bronstein) ‡ | 4.40 | 6,447 | 2.2 M | 0.114 | 0.135 | 0.027 | 0.070 | 0.034 | 8.6% |
| `texpp-3-2` (Automatic) ‡ | 6.72 | 9,448 | 2.3 M | 0.235 | 0.176 | 0.042 | 0.069 | 0.009 | 7.9% |
| expcube `(2E^x^2 x)/(6(1+E^x^2)^(2/3)+3(1+E^x^2))` | 1.23 | 6,412 | 1.9 M | 0.026 | 0.027 | 0.011 | 0.007 | 0.003 | **6.0%** |
| `mperturb-2-10` ‡ | 11.40 | 2,956 | 1.0 M | 0.037 | 0.035 | 0.013 | 0.016 | 0.005 | 0.9% |
| `cM-3-21` ‡ | 7.59 | 98 | 0.1 M | — | — | — | — | — | 0.1% |
| `trager-ratradicand-007` ‡ | 10.80 | 0 | 0 | — | — | — | — | — | 0.0% |

‡ These rows ran as eight concurrent processes. Section timers measure wall time inside a
section, so preemption inflates them. `cM-1-292`'s 34.8% is an overestimate; the same
integrand under Automatic measured 18.5% alone. Use the instruction A/B of §2.5 for those
rows.

**Ideal speedup if scoping cost nothing (estimate, computed as 1/(1 − share)).**
- Trager: 1.21×
- `loglog` / `cM-1-292`: about 1.2–1.25×
- `cM-3-152`, `logxy`, `lccontent`, `texpp`: 1.08–1.12×
- expcube: 1.06×
- the algebra-bound or budget-bound tail: none

This is a ceiling. It is not a forecast: some of the work, such as building changed paths
and evaluating the substituted body, is inherent.

**Not every case in the tail behaves like this.** The 256/257 tail lists many cases at
about 10, 12, 20 or 30 s that end `UNEVAL` or `TIMEOUT`. They are capped by budgets. For
them a faster evaluator changes *what finishes inside the budget*, not the elapsed time
(§2.2).

### 2.4 Where the visits go: one function dominates the Trager case (measured)

Instrumented counts for the Trager case:
- `substitute_scoping` ran 119,265 top-level walks and visited 101,148,312 nodes, about 848
  per call.
- 98,178,134 of those visits (97.1%) returned the node unchanged (shared).
- Only 1,956,203 nodes were rebuilt.
- There were 42.9 M symbol visits, and each scanned the `ScopingEnv` linked list 9.06
  entries deep on average.
- `rb_rec` visited 22.7 M nodes in 273,802 rule applications.

Grouping the `substitute_scoping` visits by `Module` (keyed by its local names) gives this
(Appendix C.3):

| nodes walked | calls | nodes/call | Module (locals) |
|---|---|---|---|
| 59,605,466 | 2,024 | 29,449 | `goodPrime[ctx_, p_] := Module[{B = ctx["B"], x, m, G, Gs, Bp = ctx["Bp"], Ep = ctx["Ep"], F = ctx["F"], …}, …]`, TorsionModP.m:57 |
| 10,193,885 | 30,194 | 338 | `Module[{r}, {kc[<|"Type" -> "Q"|>, Numerator[r]], Denominator[r]}]` |
| 9,395,028 | 57,994 | 162 | `Module[{e1}, If[tragerParams =!= {} && …]]` |
| 4,582,726 | 2,580 | 1,776 | `Module[{t, m, dd, ords}, …]` |
| 3,088,680 | 430 | 7,183 | `Module[{t, R, n, cols, ords, clas, …}, …]` |

`goodPrime` alone is 59% of all scoping visits. Its body is not large. What is large is the
value of its pattern variable `ctx`: the component context, an `Association` holding the
curve basis, residues and packed matrices. `rb_rec` substitutes it into the right-hand side
in four places (`ctx["B"]`, `ctx["Bp"]`, `ctx["Ep"]`, `ctx["F"]`). `builtin_module` then
walks the result, and so walks `ctx` four times, on each of the 2,024 candidate primes.
This **"walking an injected value"** pattern is the main cost, and it is invisible in the
`.m` source.

### 2.5 Instruction-count A/B of the prototype (measured)

Variants are written `[symsum rbsum nobump binder]`:
- `0000` is the current algorithm (instrumented build);
- `0010` is clock-bump removal (D) only;
- `1201` is summaries in `substitute_scoping` and `rb_rec`, binder-indexed;
- `1211` is both.

The final round ran on a quiet machine, where wall time ≈ CPU time. Values are G instructions
retired for the whole process, including 2.96 G of startup. The answer hash was identical
across all variants in every row.

| case | 0000 | 0010 (D) | 1201 (A′) | 1211 (A′+D) | Δ total | Δ excl. startup | wall 0000 → 1211 |
|---|---|---|---|---|---|---|---|
| trager-torsion-009 (Trager) | 78.134 | 75.390 | 67.873 | 65.110 | −16.7% | **−17.3%** | 15.02 → 12.83 s |
| cM-1-292 (Bronstein) | 26.436 | 23.771 | 25.117 | 22.472 | −15.0% | **−16.9%** | 6.72 → 5.64 s |
| lccontent-1-12 | 31.494 | 30.931 | 30.970 | 29.730 | −5.6% | −6.2% | 6.23 → 5.97 s |
| logxy-3-6 | 35.741 | 35.701 | 37.281 | 35.296 | −1.2% | −1.3% | ≈ |
| cM-3-152 | 36.211 | 36.940 | 35.372 | 35.677 | −1.5% | −1.6% | ≈ |
| expcube | 10.375 | 10.362 | 10.238 | 10.230 | −1.4% | −2.0% | 1.69 → 1.64 s |
| texpp-3-2 | 47.594 | 47.568 | 47.384 | 47.343 | −0.5% | −0.6% | ≈ |
| lccontent-2-3 | 31.325 | 30.282 | 31.435 | 31.708 | +1.2% | +1.3% | ≈ |
| mperturb-2-10 | 55.419 | 57.143 | 55.495 | 56.027 | +1.1% | | ≈ |
| cM-3-21 | 54.242 | 54.791 | 55.453 | 54.026 | −0.4% | | ≈ |
| trager-ratradicand-007 | 43.362 | 43.368 | 43.364 | 43.368 | 0.0% | | ≈ |

Rows inside ±2% cannot be resolved. The budget-dependent cases vary by several percent
between runs of the *same* variant: compare `lccontent-2-3` and `mperturb` across variants.

A second datum is the stock repository binary against the full prototype (instrumented),
on the Trager case: 75.52 G → 65.13 G (−13.8%). Peak RSS was 53.2 MB → 55.1 MB, of which about
1 MB is the instrumentation's own tables.

### 2.6 Micro-benchmarks (measured; CPU µs per call, best of 3 × 2,000 calls)

The script is `$S/bench/micro.m`, reproduced in Appendix B. Each function's body holds a
dead branch, `If[cond, live, <m statements>]`. `If` is `HoldRest`, so the dead branch is
never evaluated, but every substitution walker still visits it. This is the usual shape of
`.m` library code: large bodies with many `If`/`Which` branches, of which one is taken.

The stock binary (`./Mathilda`, quiet run):

| dead statements (≈ leaves) | rule only `f[x_] := If[…]` | Module | With | Block | Function `{a}` |
|---|---|---|---|---|---|
| 0 (2) | 1.0 | 7.2 | 1.5 | 2.6 | 1.8 |
| 10 (103) | 2.6 | 6.8 | 4.4 | 6.2 | 3.7 |
| 100 (1,003) | 19.1 | 41.3 | 39.5 | 22.4 | 18.8 |
| 1,000 (10,003) | 190.8 | 332.2 | 373.2 | 195.7 | 171.7 |

Reading this table:
- **rule only.** This is the `rb_rec` walk, about 19 ns per leaf even though the pattern
  variable occurs only in the condition.
- **Module.** It adds a second walk, about 14 ns per leaf.
- **With.** It is the same as Module.
- **Block.** It does no substitution, so it costs the same as rule only.
- **Function.** This is the `substitute_names` walk.
- **Nesting.** Wrapping the 1,003-leaf body in 1, 2, 4 and 8 nested `Module`s gives 40, 57,
  84 and 160 µs. Each level walks everything below it, so the cost is O(depth × size).
- **Locals.** With 1, 4, 16 and 32 locals the costs are 41, 46, 48 and 82 µs: the
  linear `ScopingEnv` scan at every symbol.

The prototype, all on (`1211`), against the same instrumented binary all off (`0000`):

| row | 0000 | 1211 | ratio |
|---|---|---|---|
| rule only, 10k leaves | 231.7 | 1.1 | 210× |
| Module, 10k | 386.8 | 4.2 | 92× |
| With, 10k | 422.9 | 2.0 | 210× |
| Block, 10k | 252.6 | 2.6 | 97× |
| Function, 10k | 183.1 | 191.5 | not prototyped |
| nested depth 8 | 209.6 | 22.9 | 9.2× |
| 32 locals | 100.6 | 39.3 | 2.6× |

The instrumented build runs 15–20% slower than stock on these rows, because of the counters.

---

## 3. How the current mechanism works

### 3.1 The `Expr` node (src/expr.h:138-235)

`sizeof(Expr)` is **56 bytes** (measured). macOS's allocator hands out 64 for it, also
measured with `malloc_size`. The fields are:
- `type`: 4 bytes, an enum with 11 values.
- `refcount`: an `unsigned`.
- `last_evaluated_at` (`uint64_t`): the evaluation stamp. Bit 63 is GROUND and bit 62 is
  INERT (`src/eval.c:266-299`).
- `hash_cache` (`uint64_t`): the lazily memoised structural hash. "MUST be reset to 0
  (`expr_invalidate_hash`) anywhere a live node's head/args/arg_count/payload is rewritten
  in place" (expr.h:172-184).
- A 32-byte union. The FUNCTION arm is `head`, `args`, `arg_count` (`size_t`) and
  `struct AssocIndex* index`. `index` is non-NULL only on canonical `Association` nodes
  (expr.h:202-214), and the union has no slack left.

How nodes behave:
- **Sharing.** `expr_copy` is a refcount increment (expr.c:878-882). A node with
  `refcount > 1` is immutable: mutators call `expr_unshare` (expr.c:660), which builds a
  one-level private copy with stamp and hash reset.
- **The node pool.** Nodes come from a bounded free-list of 8,192 nodes (expr.c:46-143).
  The link to the next free node lives in `data.function.head`. Freed addresses are
  **reused immediately**, which matters for option (b).
- **Constructors.** All of them go through `expr_alloc_node` (expr.c:113), which zeroes
  `hash_cache`.
- **Symbol nodes.** Symbols have one canonical shared node per `SymbolDef`
  (`expr_new_symbol_def`, expr.c:276-300), except Module locals. Each node naming a local
  is a fresh counted node (`temp_nodes`), which is how locals are reclaimed.
- **Small integers.** Integers in [−256, 1024] are shared too.
- **In-place mutation.** I found 121 direct `data.function.args[i] = …` stores and 11 head
  stores under `src/`, but only 10 `expr_invalidate_hash` calls. Most stores act on freshly
  built nodes, as the M3 audit comment at expr.c:868-877 says. The discipline is "mutate
  only what you own exclusively". It is not "call the invalidator", and a per-node summary
  must not assume otherwise (§5.1).

### 3.2 Rule application: `replace_bindings` / `rb_rec` (src/match.c:2948-3150)

Every DownValue, OwnValue (non-trivial), SubValue, UpValue and `ReplaceAll` application
ends in `replace_bindings(rule->replacement, env)`. The call sites are symtab.c:2422, 2675,
2761 and 2859, and replace.c:174, 647 and 1166.

- **No bindings.** The right-hand side is returned shared (match.c:3011-3021, the issue #112
  fix).
- **Otherwise, a full walk.** `rb_rec` (match.c:3023) visits **every** node of the
  right-hand side:
  - A symbol bound in `env` is replaced by `expr_copy(bound)`, a refcount bump: the value is
    injected by pointer.
  - A `Sequence` value is spliced, but only where a pattern variable stood (match.c:3054-3066).
  - A node whose head and arguments came back pointer-identical is shared
    (match.c:3090-3112). This shared-unchanged rule is what keeps evaluation stamps and
    hashes alive on the unchanged parts.
- **Capture avoidance** (`docs/notes/MATHILDA_DIVERGENCES.md` A11; changelog 2026-09-21,
  v0.212).
  - **When it acts.** On meeting a *binding scope*, `rb_rec` builds the "danger set" once
    per application, lazily (match.c:2967-3007). A binding scope is
    `Module`/`With`/`Function`-with-names; Block is excluded because it is dynamic,
    modular.c:352-381. The danger set is every distinct symbol name in every binding
    value, collected by `rb_collect_symbol_names` (match.c:2952, a full walk of each value
    with a linear de-dup of up to 256 names), plus the pattern variables' own names.
  - **The rename.** If a local of the construct is in the danger set,
    `scoping_capture_avoid` (modular.c:481-517) mints `L$n` for it and alpha-renames the
    construct with `scoping_apply_rename` (modular.c:434-476). Bodies are renamed with
    `substitute_scoping`; binding occurrences are renamed directly.
  - **The invariant this establishes.** After `rb_rec`, no injected value contains a local
    of any binding scope in the right-hand side, so the later `Module` walk over an injected
    value can never find anything. §4 uses this.

### 3.3 `Module` (src/modular.c:519-631)

1. **`$ModuleNumber`.** It is read by `evaluate($ModuleNumber)`, so a user's reassignment is
   honoured. The counter is incremented and written back as an **OwnValue**
   (modular.c:532-541, `symtab_add_own_value`). That write goes through `add_rule_in`
   (symtab.c:1878), which bumps the evaluation clock (`eval_clock_bump_binding`,
   eval.c:183) and runs `pattern_canonicalize` and duplicate detection, **on every Module
   call**.
2. **Locals.** For each local the code formats `name$N` with `snprintf` and interns it with
   `symtab_intern_def` (symtab.c:300). For a fresh name this costs a hash, a `calloc` of a
   `SymbolDef` and a `malloc` of the name: `node_intern_ex` at symtab.c:182 is 2% of the
   Trager profile. The local is marked `Temporary`, and when it is fresh,
   `symtab_temp_begin` counts it for reclamation (modular.c:573-584). The env node holds
   `expr_new_symbol_def(ld)`. An initialiser is evaluated **before** the body is substituted,
   and stored with `symtab_add_own_value_def`, which is another clock bump (modular.c:589-591).
3. **Body.** `substitute_scoping(body, env)` runs (modular.c:599), then `evaluate`, then
   `expr_free` of the substituted body.
4. **`Return`.** It is trapped (modular.c:604-616).
5. **Epilogue.** The env and init values are freed, then `symtab_temp_release` runs for each
   counted local (modular.c:631-632). Releasing a local walks its own rule lists with
   `temp_walk` (symtab.c:394-417) to decide whether anything outside still refers to it. The
   reclamation design and its invariants are in changelog 2026-10-26 (v0.6056-0.6059) and in
   the `project_module_local_reclamation` notes.

### 3.4 `substitute_scoping` (src/modular.c:74-224)

This is a recursive walk over the body with a `ScopingEnv` linked list of
(interned old name → replacement) pairs (modular.c:29-33).
- **Symbol nodes.** A symbol is looked up by scanning the env list with a pointer compare
  (modular.c:76-86).
- **Atoms.** They are returned shared.
- **Function nodes.**
  - **Shadowing.** If the node is a scoping construct that binds in argument 0
    (modular.c:45-72), the names it binds are removed from the env pushed into its body.
    The constructs are `Module`, `Block`, `With`, and `Function` with two or more arguments.
    A filtered copy of the env is `malloc`ed node by node (modular.c:91-127). Note that
    `Block` *does* shadow here, even though it is excluded from capture avoidance.
  - **The binding list.** Only the right-hand sides of `Set`/`SetDelayed` entries are
    substituted, and against the *outer* env (simultaneous-binding semantics). The
    left-hand-side names are left alone (modular.c:129-181).
  - **Rebuild or share.** Children are substituted. If the head and all arguments come back
    pointer-identical, the node itself is shared. Otherwise a new node is built
    (modular.c:196-210, v0.6057).

`Table`/`Do`/`Sum`/`Product` are deliberately **not** scoping constructs. They localise
dynamically, as `Block` does (modular.c:35-43; changelog 2026-11-09, v0.6525). The walk
descends into held expressions: `Module[{x}, Hold[x]]` gives `Hold[x$n]`, as in the
reference. It allocates only along changed paths, which the 1.96 M rebuilt nodes against
98 M shared visits confirms.

### 3.5 `With` (src/modular.c:1098-1171) and `Block` (modular.c:1072-1096)

`With` evaluates each `Set` initialiser; a `SetDelayed` entry binds its right-hand side
unevaluated. A `Throw` from an initialiser propagates at once (modular.c:1118-1131). `With`
then does the same `substitute_scoping` walk with values in place of fresh symbols
(modular.c:1149).

`Block` does no substitution. `mth_block_enter`/`mth_block_leave` (modular.c:974-1070) save
and restore OwnValues, DownValues, UpValues, attributes and builtin pointers, and keep a
frame stack for `TimeConstrained` unwinds.

### 3.6 Pure functions (src/purefunc.c)

- **Named arguments.** `substitute_names` (purefunc.c:284-350) is a full walk with a
  per-call name array. A nested `Function` that rebinds a name drops it from the set
  (`function_binds_name`, purefunc.c:237).
- **Slots.** `substitute_slots` (purefunc.c:165-229) is a full walk that stops at a nested
  slot-form `Function`.
- **Sharing.** Both share unchanged subtrees through `subst_children` (purefunc.c:129-156).
- **Call site.** They are invoked from `apply_pure_function` (purefunc.c:379).

### 3.7 The evaluation clock and Module

`last_evaluated_at` stamps a fixed point. A stamp is reusable if any one of these holds
(eval.c:296-298):
- it equals the live clock `g_eval_clock`, an exact hit;
- it is GROUND and newer than the rule epoch;
- it is INERT and newer than `g_inert_epoch` (v0.6053, v0.6119).

Every OwnValue write bumps `g_eval_clock` (symtab.c:1879-1887). So **every Module call
invalidates every exact-clock stamp in the system at least once**, through the
`$ModuleNumber` write, and once more per initialised local. This also fails every open
dependency frame ("the clock moving", eval.c:155-200). The INERT/GROUND mechanisms limit
the damage for inert data, but not for the rest.

### 3.8 Invariants any change must keep

- Lexical renaming for `Module`, `With` and named `Function`, including inside held
  expressions. Dynamic scope for `Block` and the iterator heads.
- Shadowing by inner constructs, including a bare `Function[x, …]` parameter
  (modular.c:96-100). A slot-form `Function[body]` binds no names (modular.c:51-66).
- Simultaneous binding: initialisers see the outer scope (modular.c:135-146).
- Capture avoidance in rule application (A11). It applies to Module/With/Function, not Block.
- `x$nnn` naming from `$ModuleNumber`, including a user's reassignment of it. `Unique`
  shares the counter (modular.c:237-286).
- Local reclamation: every symbol node naming a local is counted. Never store a raw local
  name without an `Expr` that names it (`project_module_local_reclamation`).
- Sharing: return unchanged nodes as themselves, so their stamps and hashes survive (issue #112).
- Message, `Throw` and `Return` propagation from initialisers and bodies (modular.c:563-568,
  604-616).

---

## 4. Root-cause analysis

1. **Every walk is O(body) regardless of where the locals are.** The v0.6057 sharing removed
   the *allocation* from unchanged paths but not the *visit*. In the Trager case 97% of
   visits are to unchanged nodes, and a perfect oracle would visit about 6% of today's nodes
   (simulated: rebuilt 1.96 M + symbol hits 1.01 M + 3.11 M unchanged roots ≈ 6.1 M of 101 M).
2. **Injected values are walked again.** `rb_rec` puts binding values into the right-hand
   side by pointer, and the `Module`/`With` walk then traverses them, though capture
   avoidance guarantees they hold no local. This is 59% of the Trager case's visits through
   `goodPrime` alone. The same pattern appears wherever a `.m` function takes a big record
   (`Association`, matrix, polynomial) as an argument and opens a `Module`, which is the
   house style of the Trager and Bronstein packages.
3. **Two walks per rule-plus-Module call, and O(depth × size) for nesting.**
   `f[x_] := Module[…]` costs an `rb_rec` walk plus a `substitute_scoping` walk, and every
   inner `Module` walks again what its parent just walked.
4. **The danger set walks every binding value** (`rb_collect_symbol_names`, with a linear
   256-entry de-dup) on 119,197 rule applications in the Trager case. That is 0.14–0.21 s
   just to decide that nothing needs renaming, which happened only 30 times.
5. **Linear env scans at every symbol.** The average scan depth is 9 in the Trager case,
   because many of its Modules have 10–20 locals. The `locals = 32` micro row shows the
   cost.
6. **Per-call fixed costs in `Module`.**
   - The `$ModuleNumber` OwnValue rewrite: `evaluate` + `add_rule_in` +
     `pattern_canonicalize` + a clock bump.
   - Local interning: `snprintf` + hash + `calloc` + `malloc`.
   - A clock bump per initialised local.
   - `temp_walk` at exit.

   Together these were 0.34 s of the Trager case, about 2.9 µs per Module. The clock bumps
   also have a **non-local** cost: they invalidate exact stamps everywhere and fail open
   dependency frames. Removing just two of them saved 2.7 G instructions on Trager and 2.7 G
   on `cM-1-292` (variant `0010`, §2.5). That is more than the measured prologue time, which
   is consistent with the evaluator re-deriving fixed points it would otherwise have reused.
7. **Allocation churn is *not* mainly scoping.** I attributed the `malloc`/`free` leaf
   samples of the Trager profile to their nearest caller (`$S/tools/callers.py`). The
   largest are FLINT and GMP through the `TimeConstrained` accounting allocators:
   `tc_flint_deallocate` 220, `tc_gmp_deallocate` 163, `tc_gmp_allocate` 120,
   `tc_flint_allocate` 105, then `to_mpoly` 79. `substitute_scoping` (22), `builtin_module`
   (17) and `node_intern_ex` (12) are small. Allocation-focused fixes, such as pools or
   inline arguments, would not move this problem.

---

## 5. Design options

The common goal is to visit only the parts of a tree that can contain a name being
substituted, and to stop paying non-local costs per `Module`.

### 5.1 Option A — per-node free-symbol summary (recommended in the A′ form)

**Idea.** Each FUNCTION node can carry a 64-bit summary: a one-hash Bloom filter of the
symbols in its subtree, heads included. A walker computes a mask from the names it
substitutes. At any node whose summary is valid and disjoint from the mask, it returns the
node shared without descending. A Bloom filter only ever answers "definitely absent" or
"maybe present", so a skip is exact. False positives only cost a descent.

**Three sub-choices decide whether it pays off. All three were measured or simulated.**

*(i) Which symbols get bits: all symbols, or binder names only.* The obvious version gives
every symbol a bit by hashing. It **saturates**. Code and values contain dozens of distinct
symbols, and every head (`Plus`, `If`, `List` …) is a symbol. I simulated it exactly over
the Trager run (Appendix D), as the fraction of today's visits that remain:
- 16 bits: 79%
- 32 bits: 52%
- 63 bits: 14%, if summaries are trusted on every node
- 63 bits, only on shared nodes (see (ii)): 68%

Built (hashing, excluding builtins and `Protected` symbols), it cut visits only to 62%
(101 M → 63 M). The many-local `goodPrime` mask collided with the many-symbol `ctx` value.
Hash bits also depend on pointer values, which vary from run to run with ASLR: the same
micro row varied 50 → 115 µs between runs.

The fix is to **assign bits only to binder names**: names that occur as Module/With/Function
locals or as pattern variables. All other symbols contribute nothing. A value like `ctx`
contains almost no binder names, so its summary is nearly empty and it is skipped at the
top. Built this way, visits fell to **10%** (101 M → 10.1 M) and `goodPrime` from 29,449 to
616 nodes per call. Bits are handed out round-robin over 55 positions in first-seen order,
which is deterministic for a given definition order.

Soundness needs one rule. A symbol that was summarised *before* it became a binder is
missing from older summaries. Each summary therefore carries an **epoch**, and registering
a new binder advances the epoch, which invalidates every older summary. There were 294
registrations in the whole Trager run, most of them early. Two refinements reduce the
bumps further:
- Register binders eagerly when a rule is stored, by scanning the stored right-hand side
  for binding lists and the left-hand side for pattern names, so registrations happen at
  package-load time.
- Bump only if the symbol was previously *seen* by a summary, using a per-`SymbolDef`
  flag.

*(ii) Which nodes may carry a summary: every node, or frozen nodes only.* The risk is an
in-place mutation that *adds* a symbol to a node whose summary is already filled. §3.1
shows that mutation is controlled by ownership, not by invalidation calls. The prototype
therefore fills a summary only on **frozen** nodes:
- nodes reached during a walk through a node with `refcount ≥ 2`, which the M3 rule makes
  immutable together with everything below it;
- or nodes of a stored rule, which the prototype reaches by treating `replace_bindings`'
  argument as frozen.

Fresh, rebuilt nodes get no summary. Simulated, frozen-only costs little: 16.7% of visits
remain, against 11.4% when every node is trusted, with a comparable exemption rule. The
summary is computed bottom-up from the children's summaries *during the walk that is
happening anyway*, with no extra pass, and a node whose children are not all valid is left
unfilled.

Invalidation:
- `expr_alloc_node` and `expr_unshare` start at "invalid";
- `expr_invalidate_hash` also clears the summary, so every place that already resets the
  hash resets this too;
- a verify build (§8.3) checks every skip with an exact occurs test.

*(iii) Where the 64 bits live.*
- **A1, a new field.** `sizeof(Expr)` goes 56 → 64. On macOS this is free in practice:
  `malloc(56)` already occupies 64 bytes (measured). The prototype used A1 and peak RSS rose
  from 53.2 MB to 55.1 MB, about 1 MB of it instrumentation tables. On glibc,
  `malloc(56)` uses a 64-byte chunk and `malloc(64)` an 80-byte chunk, so the cost there
  would be +25% per node (estimate from glibc's chunk rule, not measured). Every TU must be
  rebuilt, and `bench_assoc` is sensitive to node locality (`project_evaluator_allocation_perf`).
- **A2, steal bits from existing fields.** Narrowing `type` to 8 bits frees 24 bits, which
  per (i) is too few. The clock in `last_evaluated_at` has no safe spare bits.
- **A3, a tagged `data.function.index` (recommended).** That 8-byte slot is used only by
  canonical Associations, with **seven uses in two files** (expr.c:372, 741 and 823;
  assoc.c:307, 325, 341 and 343). Read it as a word:
  - `0` means none;
  - an even value is an `AssocIndex*` (malloc'd, so aligned);
  - an odd value is a summary. Bit 0 is the tag, bits 1–8 are the epoch and bits 9–63 are
    55 binder bits.

  `sizeof(Expr)` is unchanged, so there is no memory cost on any platform, no locality
  change and no full rebuild beyond what `expr.h` already forces. An Association that has an
  index carries no summary, so walks descend into its keys, each of which is summarised. The
  `goodPrime` `ctx` is such an Association, so A3 will skip one level lower than the A1
  prototype did. That costs O(number of keys) per visit, about 10–20 here (estimate; A3
  was not built). The cost is a slightly subtle accessor pair,
  `expr_fn_assoc_index()` / `expr_fn_symsum()`, plus a static assert that pointers are 64
  bits.

**Consumers, all exact.**
1. **`substitute_scoping`.** Mask = the bits of the env's old names. The prototype also
   showed that the per-symbol env scan can be pre-filtered by the symbol's own bit.
2. **`rb_rec`.** Mask = the bits of the pattern-variable names. Skipping a subtree with no
   pattern variable is exact for every `rb_rec` duty: splicing happens only at pattern
   symbols, and capture avoidance only matters where something is injected.
3. **The danger set.** Replace "collect all names in all values" by "for this binding
   scope, does any local's bit intersect the OR of the value summaries?". Collect names
   exactly only if it does. Values are shared and therefore frozen.
4. **`substitute_names` / `substitute_slots`.** Register `Slot` and `SlotSequence` as
   pseudo-binders, so a subtree with no `#` is skipped. Not prototyped.
5. **Later.** `temp_walk` and `FreeQ`-style queries are possible consumers; out of scope here.

**Correctness risks and how each is handled.**
- **A stale summary after an in-place mutation** is the one real risk. It is handled three
  ways: summaries are filled on frozen nodes only; `expr_invalidate_hash` clears them; and
  the verify and differential modes (§8) check every skip. A mutation that only *removes*
  symbols, such as flattening, is harmless, because a superset is conservative.
- **Variable capture.** Capture logic is unchanged, because skipping only drops visits that
  would have returned the node unchanged.
- **Held expressions** are walked exactly as today, since holds are not special to the
  walker.
- **Nested scopes and shadowing.** Skipping with the unfiltered outer mask is conservative:
  the mask is a superset.
- **Symbols created at runtime**, such as `Unique`, `Module` locals and `ToExpression`
  output, are not binders unless they are bound, so they contribute no bits. A fresh local
  name appears in a summary only if something later binds that name, which goes through the
  epoch.
- **Module and Unique counters** are untouched.
- **Block** is unaffected: it has no walk.
- **OwnValues.** A value read through an OwnValue is shared, hence frozen and summarisable.
- **Reclaimed `SymbolDef`s** leave stale bits behind, which only produces false positives.
- **Epoch wrap.** With 8 epoch bits, the 256th bump must disable the fast path for the
  session; a 16-bit epoch with 47 binder bits also fits. The trigger is one global flag.
- **Threads** (`MATHILDA_THREADS`). Writing a summary into a shared node is the same benign
  race the hash cache already has: idempotent, word-sized. The binder registry must be
  guarded, or run on the main thread only.

**Expected speedup.** Estimate: what the prototype measured, −17% instructions on the
Trager case and −17% on `cM-1-292` with (D) included, and 0–6% on the other slow cases.
Pure-function consumers and the danger set were not prototyped and should add a little.

**Effort.**
- New `src/symsum.{c,h}`: summary compute/validate, the binder registry and epoch, about
  250 lines.
- Edits to `expr.h`/`expr.c` (accessors, invalidation), `assoc.c` (4 lines),
  `modular.c`, `match.c`, `purefunc.c` and `symtab.c` (the binder-bit field and eager
  registration in `add_rule_in`).
- A unit suite.

In total, about 2–4 days with verification.

**Blast radius.** Every walker that substitutes, plus the Association index accessors.
There is no change to evaluation order, printing or hashing.

### 5.2 Option B — a side cache keyed by node pointer

**Design.** An open-addressing table mapping `Expr*` to a summary, filled like A.

**Problems.**
- Every lookup is a hash probe at each visited node, which costs about as much as the visit
  it is meant to save (about 14 ns per node, §2.6).
- Invalidation must happen in `expr_free` for every FUNCTION node, the hottest free path,
  or the table holds dead keys.
- The node pool reuses a freed address **at once** (expr.c:113-143), so a stale entry is
  silently wrong, not merely missing.
- Restricting it to stored-rule nodes, with entries owned by the `Rule`, is safer but misses
  injected values, which are the dominant cost (§2.4).

**Verdict.** It is more fragile and slower than A3, which costs zero bytes. Reject it.

### 5.3 Option C — move the hottest `.m` stages to C

**What it would mean.** The scoping cost is spread across about 440 `Module`s in the two
packages (Bronstein 213, Trager 227). The Trager case's top function, `goodPrime`, is 59% of
*its* scoping visits but irrelevant to the Bronstein cases.

**Why it does not fit.** Porting removes the overhead only for the functions ported, and it
conflicts with project policy: algorithms in `.m`, the algebra substrate in C
(`feedback_new_integrate_method_in_dotm`, `feedback_algebra_substrate_in_c`), no Wolfram
Language templates in C. A cheap `.m`-level workaround also exists for the single hot spot:
bind `ctx` parts once in a caller and pass only the parts.

**Verdict.** Use C ports for algebra hotspots, which the same profile shows are FLINT
determinant/adjugate and `nf_elem_mul` work, not as the fix for evaluator overhead.
Effort grows with every stage ported, and each port carries its own correctness risk.

### 5.4 Option D — no evaluation-clock bump per Module call (recommended, independent)

**D1, a fresh local's first binding.** A local created from a fresh name has had no node
evaluated anywhere, so no stamp can depend on its binding. The prototype skipped
`eval_clock_bump_binding` for that one `add_rule_in` call only (`symtab_add_own_value_def`
when `fresh`). Later `Set`s to the local in the body still bump.

**D2, `$ModuleNumber`.** The C counter (`module_number`, modular.c:13) becomes
authoritative, so `Module` no longer rewrites an OwnValue. User reads and writes must still
work, with two possible shapes:
- install `$ModuleNumber := Internal`ModuleNumberValue[]` as a delayed OwnValue that never
  changes, and have `Module`/`Unique` check whether it has been replaced by a user `Set`;
- or add a narrow read hook in `apply_own_values` for that one symbol.

**A trap found while measuring.** The first prototype of D2 stopped updating the OwnValue,
but `Unique` still read it (modular.c:309-315) and so moved the counter **backwards**.
Later `Module`s then reused existing `x$n` names. That changed answers (`cM-3-152`) and made
two cases look 5× faster. Every reader and writer of the counter must change together, and
§8.1 has the regression probe.

**Measured.** Variant `0010`: Trager −3.5% instructions, `cM-1-292` −10.1%, others within
noise.

**Effort and risk.** About 60 lines plus tests. The risk is the dependency-frame and INERT
machinery (eval.c:155-200): the change must keep `eval_timestamps_tests` green, and the
argument that a fresh symbol cannot have been observed needs a test.

**Blast radius.** `modular.c` and `symtab.c`.

### 5.5 Option E — fuse rule substitution and Module renaming into one walk

When a rule's right-hand side is `Module[…]`, mint the locals during rule application and do
one walk with env = bindings ∪ renames. This saves one of two walks. But it moves
`$ModuleNumber` side effects and initialiser evaluation into rule application, changes
evaluation order and needs a hidden "pre-renamed Module" head. After A′ both walks cost
little. **Defer.**

### 5.6 Option F — per-rule substitution plans

At rule store time, record the paths to pattern-variable and local occurrences, then rebuild
only along those paths. A plan is exact and cheap for the rule's own right-hand side, but it
says nothing about **injected values**, which are the dominant cost, or about code built at
runtime. `Sequence` splicing also shifts positions. A′ covers all of these. **Reject.**

### 5.7 Option G — lazy renaming, or Module without rewriting the body

Evaluate the body under a dynamic binding of the original name, as `Block` does, and rename
only on escape. This breaks:
- recursion: two activations would share one name;
- escaping closures: `Module[{x}, Function[y, x + y]]` must return a function mentioning
  `x$n`;
- held results: `Module[{x}, Hold[x]]`;
- `Names[]`;
- the A11 capture semantics.

Detecting escape soundly is as expensive as the walk it saves. **Reject.**

### 5.8 Option H — cheaper env lookup

Replace the linear `ScopingEnv` scan (modular.c:76-86) with a pre-filter on the symbol's
binder bit, which is free once A′ exists, falling back to the list. Alternatively, put a
per-walk scratch index in `SymbolDef` for envs with more than 8 entries. After A′, symbol
visits fell 91% in the Trager case (42.9 M → 4.6 M) and the average scan from 9 to 3, so H's
remaining value is small. **Take the pre-filter as part of A′.**

### 5.9 Option I — summary-driven danger set

This is consumer 3 of A′, listed separately because it is a self-contained change to
match.c:2967-3007 and removes the 0.14–0.21 s of `rb_collect_symbol_names` in the Trager
case. **Part of A′, phase 3.**

---

## 6. Recommendation

Adopt **D**, then **A′ with A3 storage** (binder-indexed, epoch-validated, frozen-only,
tagged in the `index` slot), with consumers in this order:
1. `substitute_scoping`, with the H pre-filter;
2. `rb_rec`, with eager summaries at rule store;
3. the danger set (I);
4. pure functions.

Justification:
- **Evidence.** The prototype of exactly this design, A′ + D on A1 storage, cut the
  Trager case by 17% of instructions and the `Log`-tower case by 17%. Those are the two
  cases where scoping is a large share. The other cases are unaffected, as their measured
  shares predict, and answers are unchanged.
- **Cost.** A3 costs no memory and changes no node size. Only two files touch the slot.
- **Generality.** Every substituting walker benefits, including `.m` code not yet written.
  The rule-only micro row (−99.5%) shows that plain DownValue application with large
  right-hand sides benefits, not just `Module`.
- **Safety.** Exactness is structural: a Bloom miss proves absence. The only new failure
  mode, a stale summary, is confined by frozen-only filling, can be checked by a verify
  build, and can be switched off at run time.
- **Expectations.** Do not expect the 43-case tail to halve. Most of it is algebra-bound or
  budget-bound (§2.3). The right gate for this work is the instruction count on the
  `.m`-bound cases, not the tail's wall time.

If A3's tagged slot is judged too clever, A1, a new field, is the fallback. It is free on
macOS and costs about +25% node memory on glibc (estimate). Measure glibc RSS before
choosing it.

---

## 7. Implementation plan

Each phase ends at a checkpoint that must pass before the next starts. Version bumps and
tags follow CLAUDE.md: each substantive commit gets +0.0001 and a tag. Docs go in
`docs/spec/changelog/2026-11-16.md` for this week, or the week the work lands.

**Phase 0 — measurement harness (no behaviour change).**
- Land the counters as a compile-time option (`-DMTH_SCOPE_STATS`), with a report at exit:
  per-walker nodes, rebuilt, unchanged, unchanged roots, timers, and the hottest-Module
  table. Port it from `$S/inst/src/scopestat.{c,h}`, using `mth_mono_seconds`
  (`src/mtime.c`) rather than `mach_absolute_time`.
- Add `tools/scope_bench/`: `micro.m`, the semantic probes `semantics.m`, the 11-case corpus
  and an instruction-count A/B driver. On macOS it uses `/usr/bin/time -l`; on Linux,
  `perf stat -e instructions` or `valgrind --tool=callgrind` counts.

*Checkpoint:* the harness reproduces §2.3–2.6 within noise on a quiet machine.

**Phase 1 — Option D.**
- D1: in `builtin_module`, bind a *fresh* local's initial value without the clock bump. Add a
  `bool bump` parameter (or a dedicated `symtab_bind_fresh_local`) rather than a global flag.
- D2: make `module_number` authoritative, change `builtin_module`, `builtin_unique` and
  `unique_sync_module_number` together, and keep `$ModuleNumber = n` (a user `Set`) working.

*Checkpoint:*
- `eval_timestamps_tests`, `ctest -LE nightly -E dsolve` and the new probes (§8.1), including
  `$ModuleNumber = 1000; Module[{x}, x]` and the `Unique`-then-`Module` sequence.
- Measure `0010`: Trager −3 to −4% instructions.

**Phase 2 — summary substrate (dormant).**
- `src/symsum.{c,h}`:
  - the word layout (A3: tag, epoch, 55 bits);
  - `symsum_get(e)` / `symsum_set(e, s)` / `symsum_valid(s)`;
  - `symsum_symbit(SymbolDef*)`;
  - `symsum_register_binder(SymbolDef*)`, which bumps the epoch if the symbol was already
    seen;
  - the epoch-wrap kill switch.
- `SymbolDef` gains `uint8_t ss_bit` (0 = none) and `uint8_t ss_seen`. Check that it fits
  in existing padding next to `is_system`/`has_def`, symtab.h:137-138.
- `expr.c`/`assoc.c`: the tagged slot accessors; `expr_free` frees the index only for an
  even non-zero value; `expr_unshare` and `expr_new_function` clear the slot;
  `expr_invalidate_hash` clears a summary.
- Switches: `MATHILDA_NO_SYMSUM=1` (runtime) and `-DMTH_NO_SYMSUM` (compile time).
- `tests/test_symsum.c`: layout, tag, epoch, invalidation, an Association with an index,
  wrap.

*Checkpoint:* the suites are green with the substrate present and unused; `bench_assoc` is
unchanged.

**Phase 3 — consumers, one commit each, each with its own A/B and probes.**
- 3a. `substitute_scoping`. A wrapper computes the mask and registers the env names
  (modular.c:74), and the walk is threaded with a `frozen` flag. Fill on unchanged frozen
  nodes, using the H pre-filter for symbol lookups. *Checkpoint:* the Module/With micro rows
  are ≥10× faster, Trager instructions −8% or better, the probes are identical, and the
  verify build is clean on all suites.
- 3b. `rb_rec`. The mask comes from `env->symbols`, which must be **interned**: `env_set`
  may hold a literal such as `"$OptionsPattern$"`, so canonicalise with `intern_symbol`.
  Register pattern variables when the rule is stored, in `add_rule_in`, symtab.c:1878, and
  compute the right-hand side's summaries eagerly then. *Checkpoint:* the rule-only micro row
  is ≥50× faster, and the match and pattern suites are clean in verify mode.
- 3c. The danger set (match.c:2967-3007), using summaries. *Checkpoint:* the A11 capture
  probes pass, and `rb_collect_symbol_names` is gone from the Trager profile.
- 3d. `substitute_names` / `substitute_slots` (purefunc.c), with `Slot`/`SlotSequence` as
  pseudo-binders. *Checkpoint:* the Function micro row is ≥10× faster, and the purefunc
  suite is clean.

**Phase 4 — verification sweep and docs.** Everything in §8. Update:
- `docs/design/performance.md`, or a new `docs/design/scoping_substitution.md` holding this
  proposal's design section;
- `SPEC.md` §3.1, to document the slot's two meanings;
- the changelog.

**Phase 5 — optional refinements.**
- Eager binder registration for every stored definition, at package load.
- `ss_seen`-gated epoch bumps.
- A `temp_walk` consumer.
- A16 over A3 if more bits are needed (16-bit epoch, 47 binder bits).

---

## 8. Test and verification plan

### 8.1 Semantic probes

These are the tricky cases. Twenty-seven probes are in `$S/bench/semantics.m`. Their output is
identical between the stock binary and all prototype variants (`sem.base.txt`,
`sem.0.txt`, `sem.1.txt`, `sem.r2.txt`, `sem.b1.txt`). Turn them into a ctest. The cases:
- Capture (A11): `g[v_] := Module[{e = 1}, v + e]; g[e + 1]` must be `2 + e`.
  Shadowing: `sh[e_] := Module[{e = 1}, e + 1]; sh[99]` must be `2`.
- Block is dynamic: `gb[v_] := Block[{e = 1}, v + e]; gb[e + 1]` must be `3`.
- Simultaneous binding: `With[{q = 12}, With[{k = q}, k]]` must be `12`.
- Iterators are not scoping constructs:
  - `With[{k = 2}, Hold[Table[k, {k, 3}]]]` must be `Hold[Table[2, {2, 3}]]`;
  - `gs[f_] := Sum[f, {k, 1, 3}]; gs[1/k]` must be `11/6`.
- Slot functions inside named functions: `Function[x, (x + #) &][a]` must be `a + #1 &`.
- Closures over locals: `Module[{c = 0}, Scan[(c = c + 1) &, …]; c]` must be `3`.
- Renaming inside holds: `Module[{x}, Hold[x]]`, and the `x$n` form.
- Inner rebinding: `Hold[With[{a = 2}, a], a]` and `Hold[Module[{a}, a], a]` under an outer
  `With[{a = 1}, …]`.
- Nested shadowing: `Function[{x}, Function[{x}, x][2] + x][10]` must be `12`.
- A value carrying the local's name: `cap[v_] := Module[{t}, t = v; {t, Hold[v]}]; cap[t]`.
- Redefinition after summaries exist: define `k2`, call it, redefine it with different
  locals, call again.
- System-named locals: `With[{Plus = Times}, Plus[2, 3]]` must be `6`;
  `Module[{List = 7}, List]` must be `7`.
- Large dead branches, with repeated application after summaries are filled.

Add:
- Epoch tests: a name first used as a free symbol inside a summarised value, *then* used as
  a local:

  ```
  v = Hold[{u, w}]; v;
  f[] := Module[{u = 1}, {u, v}]; f[]
  ```

  The `v` inside must not be renamed, and the inner `u` must be.
- `ToExpression`-built `Module`s with names never seen before.
- `Unique` / `$ModuleNumber` interplay: `$ModuleNumber = 5000; Module[{x}, x]`;
  `Unique[x]` followed by `Module[{x}, …]` gives distinct names.
- `TimeConstrained` aborting mid-walk; `Throw` out of initialisers.

### 8.2 Differential mode

Behind `MATHILDA_SYMSUM_DIFF=1`, every fast walk also runs the slow walk and compares the
two results with `expr_eq`. It aborts with the expression and the env on any mismatch. Run:
- the full `ctest -LE nightly -E dsolve`;
- the DSolve suite separately, since it is slow;
- the integration suites;
- the shootout `--quick` subset.

A mismatch is a release blocker.

### 8.3 Verify mode

Behind `-DMATHILDA_SYMSUM_VERIFY`, modelled on `MATHILDA_HASH_VERIFY` (expr.c:1265), every
skip runs an exact occurs check of the mask's names in the skipped subtree, and every stored
summary is recomputed on read. This catches the one real failure mode: a summary made stale
by an in-place mutation. Run all ctest suites under it.

### 8.4 Memory safety

- An ASan + UBSan build (`MTH_EXPR_NO_POOL` is automatic under ASan) of the scoping, match,
  purefunc, eval and integration suites, with `MATHILDA_TEMP_SWEEP_EAGER=1`, as required
  after touching scoping (`project_module_local_reclamation`).
- `valgrind --leak-check=full ./Mathilda -file tools/scope_bench/semantics.m` and a short
  Trager case, compared against the macOS startup baseline
  (`feedback_mathilda_startup_fixed_init_leak_baseline`).
- `leaks` on macOS.
- The Windows/MinGW leg, because of the pointer-tag assumption.

### 8.5 Regression suites

- `cd tests/build && ctest -LE nightly -E dsolve -j4`, then the dsolve suite.
- `make check-c99`, `check-messages` and `check-packed-aware`. These are not expected to
  change, but cheap.
- `bench_assoc`, several runs because it is noisy. It is sensitive to node locality, and
  this matters only for the A1 fallback.

### 8.6 Performance gates

- `tools/scope_bench` micro rows. Thresholds: Module, With and rule-only 10k-leaf rows ≥ 20×
  faster than today; nested depth 8 ≥ 5×.
- Instruction counts on the corpus: Trager ≤ −12% and `cM-1-292` ≤ −12% against the same
  build with `MATHILDA_NO_SYMSUM=1`, every other case within ±3%, and every answer hash
  identical.
- `python3 tools/integration_shootout/shootout.py --classes A,M,B --fresh` and
  `--bench benchmarks/257-compact-integrals --fresh` on a **quiet** machine: no new
  failures, and the time distribution of the A/M classes shifts down. Report the
  budget-bound cases separately (§2.2).
- `benchmarks/run_all.py --only 40,41,42,43,44,45,46` (the evaluator suite): no regression.
  The pure-dispatch rows should improve slightly.

---

## 9. Rollout and fallback

- **Runtime kill switch.** `MATHILDA_NO_SYMSUM=1` disables all consumers and is read once
  at startup. It is the A/B baseline for every gate and the fallback if a bug surfaces in
  the field.
- **Compile-time switch.** `-DMTH_NO_SYMSUM` compiles the consumers out, for bisecting.
- **Differential mode.** `MATHILDA_SYMSUM_DIFF=1` runs in the nightly CI job (see
  `.github/workflows/tests.yml`) for at least a week after each consumer lands.
- **Order.** Land D first: it has its own kill switch, `MATHILDA_MODULE_CLOCK_BUMP=1`. Then
  land the substrate dormant. Then each consumer separately, so a regression can be bisected
  to one walker.
- **Epoch exhaustion** turns the fast path off for the session, with a one-time
  `Internal`SymsumDisabled` flag visible to diagnostics.

---

## 10. Risks and open questions

**Risks.**
- **A stale summary from an in-place mutation of a node filled earlier.** This is the only
  way to get a wrong answer. It is mitigated by frozen-only filling, invalidation in
  `expr_invalidate_hash`, and the verify and differential modes.
- **Binder churn from generated code.** `ToExpression`-built `Module`s with new names each
  time would bump the epoch constantly. They are rare in the current packages (294 bumps
  per Trager case), but should be monitored. The counter goes into the stats report.
- **Bit collisions among binder names** (55 bits for about 1,000 local names). These cost
  performance, never correctness. Assignment is deterministic in definition order, so the
  benchmarks are reproducible, but package load order changes which names collide.
- **The tagged slot.** It assumes 64-bit pointers and malloc alignment ≥ 2. Add a
  `static_assert`-style check (a C99 negative-array trick) and guard the MinGW build.
- **Benchmark noise.** Wall-clock budgets in the Integrate cascade make timing and even
  answers load-dependent (§2.2). The gates must use instruction counts and answer hashes, on
  a quiet machine.

**Open questions.**
1. Should binder registration happen eagerly for *all* stored definitions at load time
   (fewest bumps), or lazily (simplest)? The measured lazy count was 294 per run.
2. Should `Block` locals be registered? Block does not walk, but its names shadow in
   `substitute_scoping` (modular.c:44-49). Today they are harmless either way.
3. Should summaries also be filled at **construction**, as in the original option (a)? That
   would cover fresh rebuilt nodes, which are about 2 M per Trager run, but it widens the
   mutation hazard to every node and costs on every `expr_new_function`. I recommend against
   it unless the verify build proves the mutation discipline airtight.
4. A3 against A1 on Linux: is A1's +25% node memory on glibc real? I did not measure it. If
   it is not, A1 is simpler.
5. `goodPrime`-style code: is it worth a `.m` refactor anyway (bind `ctx` parts in the
   caller)? After A′ it should not matter. Re-measure.
6. Can more of `Module`'s per-call cost go: lazy interning of `x$n`, a cheaper `temp_walk`?
   The earlier attempt "fully lazy names would be materialised by ordering, hashing and
   printing anyway" (changelog 2026-10-26, item 2) suggests not much.

---

## 11. Prior art in the repository

- **v0.6057** (changelog 2026-10-26, "Symbol-table follow-ups"). `substitute_scoping` shares
  subtrees with no local and compares names by pointer. `rb_rec`/`substitute_scoping` use
  stack buffers for arity ≤ 8. A rule-plus-Module call went 2.60 → 1.81 µs. Integer symbol
  ids were declined on measurements, and "fully lazy" local names were judged not worth it.
  `docs/notes/SYMBOL_TABLE_IMPROVEMENTS.md` items 2–5 record the outcomes.
- **v0.6056/0.6058/0.6059.** Module locals are reclaimed, names are freed, and symbol nodes
  are canonical and shared, with invariants for any code that creates or renames symbol
  nodes.
- **v0.6053/v0.6119.** INERT fixed points and dependency-recorded fixed points. Any OwnValue
  binding bumps the clock. `eval_clock_bump_binding` limits INERT invalidation to symbols
  relied on, which Option D builds on.
- **v0.212** (changelog 2026-09-21). Capture-avoiding `replace_bindings`; the lazy one-time
  danger set came from a Charlwood regression when it was first built per node.
- **v0.6525** (changelog 2026-11-09). Table/Do/Sum/Product removed from the scoping list:
  they are dynamic. Getting this wrong gave `g[1/k] = 3/k`.
- **2026-11-16 changelog.** The present gap, with the per-node summary named as the needed
  decision.
- **`project_evaluator_allocation_perf` / `project_evaluator_vs_reference_perf (memory note)`.** The
  node pool is bounded at 8,192 because an unbounded pool broke `bench_assoc` scaling; a
  slab allocator was rejected; adopting a per-step args array was rejected as noise. The
  lesson recorded there is to profile and decompose before optimising the evaluator, which
  is why §2 decomposes by phase.
- **`tasks/lessons.md:1085`** (2026-05-22): an old note that Module did not reach into
  HoldAll bodies. That is superseded; the current walker does (probe `Hold[x$n]`).

---

## 12. What was not measured

- A3 storage itself. The prototype used A1, a new field. A3's extra descent into indexed
  Associations is an estimate.
- glibc memory for A1, and anything on Linux or Windows.
- The pure-function consumers (`substitute_names`/`substitute_slots`) and the summary-driven
  danger set (I). Both are in the plan but were not built.
- ASan, valgrind and the full ctest suites on the prototype. Only the 27 probes and the 11
  corpus answers were compared.
- The full 256/257 shootout; only 11 cases were A/B'd.
- Wall time on a reliably quiet machine for more than the final round.
- `loglog` under Automatic in the A/B: it hit a 10 s budget under load and its answer varied.
- The cost of the H pre-filter by itself, and that of the eager rule-store summaries.
- Multi-threaded evaluation.

---

## Appendix A — raw profiles

### A.1 Trager case, `sample` 20 s, stock binary (`$S/bench/trager.sample.txt`)

Top of stack (self samples), 10,759 samples in total:

```
substitute_scoping.part.0        1045
evaluate_step_core                632
_nanov2_free                      594
nanov2_malloc_type                453
expr_free.part.0                  337
rb_rec                            294
expr_new_function                 200
node_intern_ex.part.0             164
expr_copy                         162
_platform_bzero                   160
expr_free                         157
evaluate_core                     153
rb_collect_symbol_names           125
_free                             111
builtin_times                     110
expr_compare.part.0               108
_platform_strcmp                   93
nanov2_calloc_type                 91
match_internal_impl                82
fmpz_add                           78
```

Inclusive (`$S/tools/incl.py`). Roots are truncated by stack depth, so only frames near the
leaves are reliable:

```
substitute_scoping 11.8%   rb_rec 5.6%   node_intern_ex 2.0%   rb_collect_symbol_names 1.2%
temp_walk 0.5%   builtin_flint_polynomialdetadjugate 7.3%   gr_mat_det_berkowitz 7.1%
_nf_elem_mul_red 5.0%   flint_expand_polynomial 4.3%   to_mpoly 6.8%
```

`malloc`/`free` leaves by nearest Mathilda caller (`$S/tools/callers.py`):

```
tc_flint_deallocate 220   tc_gmp_deallocate 163   tc_gmp_allocate 120   tc_flint_allocate 105
to_mpoly 79   tc_flint_calloc 51   get_coeff_expanded 40   rat_canon_reduce 30 ...
substitute_scoping 22   builtin_module 17   node_intern_ex 12
```

### A.2 `loglog` case, `sample` 8 s (`$S/bench/loglog.sample.txt`)

```
self: evaluate_step_core 135  rb_rec 86  expr_free 85  substitute_scoping 58  _nanov2_free 55
      expr_new_function 40  evaluate_core 37  nanov2_malloc_type 35  is_packed_list 34 ...
inclusive: rb_rec 11.8%  builtin_module 11.0%  expr_free 8.2%  substitute_scoping 6.2%
           node_intern_ex 1.7%
```

### A.3 A rule application with a large right-hand side

The script is `$S/bench/ruleonly.m`: `fP[x_] := If[x > 0, x + 1, <10k-leaf dead branch>]`,
30,000 calls. In the profile, `rb_rec` is 97.7% inclusive, with `expr_free` 11%,
`expr_copy` 8.3% and `expr_is_binding_scope` 7.7%. The walk itself is the cost.

## Appendix B — scripts

Everything is under `$S/bench/`.

**`micro.m`**, the micro-benchmarks of §2.6:

```
dead[m_] := StringJoin[Table["q" <> ToString[k] <> "[u, v, {1, 2, 3}, w + " <> ToString[k] <> "]; ", {k, m}]];
time[expr_, n_] := Module[{best = Infinity, t0}, Do[t0 = TimeUsed[]; Do[expr, {n}];
    best = Min[best, TimeUsed[] - t0], {3}]; Round[10^6 best/n, 0.1]];
SetAttributes[time, HoldFirst];
n = 2000;
Do[
  ToExpression["fM[x_] := Module[{a = x, b}, If[x > 0, b = a + 1, " <> dead[m] <> "Null]; b]"];
  ToExpression["fW[x_] := With[{a = x}, If[x > 0, a + 1, " <> dead[m] <> "Null]]"];
  ToExpression["fB[x_] := Block[{a = x, b}, If[x > 0, b = a + 1, " <> dead[m] <> "Null]; b]"];
  ToExpression["fF = Function[{a}, If[a > 0, a + 1, " <> dead[m] <> "Null]]"];
  ToExpression["fP[x_] := If[x > 0, x + 1, " <> dead[m] <> "Null]"];
  Print[...timings of fP[1], fM[1], fW[1], fB[1], fF[1]...], {m, {0, 10, 100, 1000}}];
(* nesting: d nested Modules around a 100-statement dead branch; locals: L locals *)
```

The remaining scripts:
- `semantics.m`: the 27 semantic probes of §8.1.
- `corpus/*.m`: the 11 A/B cases. Each prints CPU time, wall time, `Hash[r]` and
  `LeafCount[r]`.
- `corpusv/*.m`: the same cases with a numeric `D[r, x] - f` check at 3/7, 5/4 and 7/3.
- `corpus_ab.sh VARIANT`: an instruction-count A/B driver, where the variant digits are
  `symsum rbsum nobump binder`.
- `$S/tools/incl.py`, `$S/tools/callers.py`: inclusive-time and allocation attribution from
  `sample` output.

## Appendix C — instrumented build and its output

### C.1 What the instrumentation adds

The diff is `$S/prototype.diff`. The instrumentation adds:
- `src/scopestat.{c,h}`: counters, timers, the hot-body table and the report at exit, gated
  by `MTH_SCOPESTAT=1`;
- counters in `substitute_scoping`, `rb_rec`, `rb_danger_build` and `substitute_names`;
- phase timers in `builtin_module` and `builtin_with`.

### C.2 Trager case, baseline (all fast paths off) — `$S/bench/trager_stat_baseline.txt`

```
Module calls 117823 (locals 190943)  With calls 1442 (bindings 2668)
substitute_scoping: top 119265 nodes 101148312 (fn 31028432 sym 42867436 sym_hits 1013975)
  rebuilt 1956203 unchanged 98178134 unchanged_roots 3114082 avg_env_len 9.06 filters 4665
replace_bindings: top 273802 nodes 22731251 rebuilt 4422743 unchanged 3732335
  unchanged_roots 3306967 scopes_seen 129314 danger_builds 119197 renames 30
pure Function named-arg subst: calls 232 nodes 261778
(timers from the 14.66 s quiet run: ss(Module) 1.392 s, prologue 0.200, epilogue 0.140,
 free 0.068, replace_bindings 0.706 [danger 0.144, capture rename 0.028])
```

The same run with the prototype, hash bits (`1211`, binder=0), `$S/bench/trager_hot_proto.txt`:

```
substitute_scoping: nodes 63040573 ... avg_env_len 10.87 ; replace_bindings: nodes 13729713
```

The same run with the prototype, binder bits (`1211`, binder=1), `$S/bench/trager_binder.txt`:

```
substitute_scoping: nodes 10078824 (sym 4593105) ... avg_env_len 3.01
replace_bindings: nodes 14086885 ; epoch_bumps 294 ; ss_skip 969934 ; rb_skip 1115585
```

### C.3 Hottest Module bodies, Trager case, baseline (`$S/bench/trager_hot_bodies.txt`)

Columns are nodes walked, calls, nodes per call, then the local names and the start of the
body.

```
59605466   2024  29449  {B x m G Gs Bp Ep F rho lin rhos alphas a} :: x = B["x"]; m = B["m"]; ...   (goodPrime)
10193885  30194    338  {r} :: {kc[<|"Type" -> "Q"|>, Numerator[r]], Denominator[r]}
 9395028  57994    162  {e1} :: If[tragerParams =!= {} && !FreeQ[e1, Alternatives @@ tragerParams], ...
 4582726   2580   1776  {t m dd ords} :: dd = kDeg[<|"Type" -> "GF", "p" -> 17|>, ...
 3088680    430   7183  {t R n cols ords clas ...} :: If[nativeQ[...], nat = FLINT`Polynomial...
 2879430    430   6696  {t rows ords best h s} :: rows = reducedRows[...]
 2656160    434   6120  {x M det adj1 nrm cof nat} :: nat = If[nativeQ[...], FLINT`PolynomialDetAdjugate[...
```

With binder summaries, `goodPrime` drops to 616 nodes per call (1,246,522 nodes in total).

## Appendix D — Bloom-filter simulation (`$S/bench/*.sim.txt`)

For every Module/With call in a real run, I computed the visits a summary-guided walk would
make, exactly. Each node's filter was precomputed bottom-up per call. Shadowing was ignored,
which is conservative. Values are the percentage of today's visits that remain.

| case (calls, baseline visits) | trust | 16-bit | 32-bit | 63-bit | 63-bit, `is_system` excluded † |
|---|---|---|---|---|---|
| Trager (119,265; 101.2 M) | all nodes | 79.3 | 51.9 | 14.1 | 11.4 |
| | refcount ≥ 2 only | 94.0 | 78.3 | 68.4 | 16.7 |
| loglog run (52,169; 8.0 M) | all nodes | 55.7 | 33.4 | 32.0 | 31.8 |
| | refcount ≥ 2 only | 71.9 | 47.3 | 45.4 | 44.9 |
| expcube (6,539; 1.9 M) | all nodes | 52.1 | 33.1 | 42.9 | 20.7 |
| | refcount ≥ 2 only | 77.2 | 53.8 | 61.1 | 28.0 |

† `is_system` turned out to be too broad an exemption. `intern_mark_all_system`
(symtab.c:520) marks every name interned at startup, including `a`, `x` and `t`. The built
prototype therefore exempted builtins and `Protected` symbols instead (101 M → 63 M
measured), and then moved to binder-indexed bits (→ 10.1 M measured). That supersedes the
hash-based columns.
