# UpValues for Mathilda — implementation plan

## Context

Mathilda today has **OwnValues** (`x = 5`) and **DownValues** (`f[x_] := …`) but **no
UpValues** — the Mathematica mechanism that lets a symbol `g` carry a rule for how it
behaves when it appears *inside* another head's call (`f[…, g[…], …]`). This blocks the
idiomatic "associate a rule with the operand, not the operator" pattern (modular
arithmetic `mod/: mod[a_,p_] + mod[b_,p_] := …`, symbolic properties
`area[square] ^= s^2`, per-symbol databases `rate[chf] ^= 0.812`). The subsystem is
confirmed **fully greenfield** — no stubs anywhere in the tree (only prose in
`plans/EVAL_IMPROVEMENTS_PLAN.md` marking it "intentionally not implemented").

This is a **trial**: every change is gated behind `#if UP_VALUES` so the feature can be
compiled out to exactly today's behaviour. The goal is a faithful, highly-efficient,
leak-free (valgrind-clean) recreation of UpValues plus its five assignment/inspection
heads, with extensive unit tests.

**Scope decisions confirmed with the user:**
- **Build default: ON** — `-DUP_VALUES=1` in the makefile and test build by default;
  `make UP_VALUES=0` compiles it out.
- **Inspection: full** — add a `Definition[f]` builtin and extend `?name`/`Information`
  to display own/down/up values (today they show only the docstring).
- **List setters: all three** — `UpValues[f]=list`, **and** `DownValues[f]=list` /
  `OwnValues[f]=list` (both previously unsupported), via shared new infrastructure.

**Heads to implement:** `UpValues`, `UpSet` (`^=`), `UpSetDelayed` (`^:=`),
`TagSet` (`f /: lhs = rhs`), `TagSetDelayed` (`f /: lhs := rhs`),
`TagUnset` (`f /: lhs =.`), `Definition`.

> **First implementation step:** copy this plan to `UPVALUES_PLAN.md` at the repo root
> (the user asked for it there; plan-mode only permits editing this plan file). Keep the
> two in sync.

---

## Semantics to reproduce (from the Mathematica reference the user supplied)

1. **Storage & dispatch.** An upvalue associated with symbol `s` fires when `s` appears
   at **level one** of an enclosing call — i.e. `s` is a direct argument, or `s` is the
   **head of a direct argument**. Matching upvalues are tried **before** the enclosing
   head's downvalues/builtin.
2. **`UpSet` (`lhs ^= rhs`) / `UpSetDelayed` (`lhs ^:= rhs`)** install the rule on
   **every distinct symbol** at level one of `lhs` (both plain-symbol args and
   arg-heads). `prop[a, b[c]] ^= v` defines upvalues on **both** `a` and `b`.
3. **`TagSet` (`f /: lhs = rhs`) / `TagSetDelayed` (`f /: lhs := rhs`)** install on the
   **single tag `f`**, which must occur in `lhs`. Classification ("defines upvalues,
   downvalues, or subvalues as appropriate"):
   - `f` *is* `lhs` (bare symbol) → **OwnValue** (tag redundant).
   - `f` is the **head of `lhs`** → **DownValue** (tag redundant).
   - `f` occurs as an **element of `lhs`, or the head of an element** → **UpValue**.
   - `f` is the **head of the head** (`f[x][y]`) → SubValue: **out of scope** (Mathilda
     has no SubValues). Emit a message and leave unevaluated; documented limitation.
   - `f` does not occur in `lhs` → message `TagSet::tagnf` (verify exact tag text), no
     change, still yield `rhs`/`Null` per WMA.
4. **`TagUnset` (`f /: lhs =.`)** removes the single matching rule from the appropriate
   list (own/down/up, classified as above) of tag `f`. Plain `lhs =.` (Unset, already
   implemented) does **not** reach upvalues — matching WMA's `Unset::norep` for
   `f[h[x_]] =.`.
5. **`UpValues[f]`** returns `{HoldPattern[lhs] :> rhs, …}` (same shape/readback as
   `DownValues[f]`). Accepts a **symbol** (HoldAll) **or a string** (`UpValues/@Names["x*"]`
   must work); a string naming a non-existent symbol emits `UpValues::sym`.
6. **`UpValues[f] = list`** (and the two siblings) replace the whole rule list.
7. **HoldAllComplete suppresses upvalues; HoldAll does NOT.** (`Hold[g[x]]` still fires
   `g`'s all-heads upvalue; `HoldComplete[g[x]]` does not.)

**Documented divergences (pre-existing Mathilda conventions, not regressions):**
- Rule order follows Mathilda's **descending-specificity sort** (insertion order as the
  tie-break), identical to DownValues. The spec's reorder example
  (`UpValues[x]=Reverse[…]`) relies only on equal-specificity tie-breaking, which this
  preserves correctly.
- **Immediate vs delayed is not recorded** on a `Rule` (a long-standing Mathilda fact:
  `DownValues[f]` already renders everything as `RuleDelayed`). So `Definition`/`?name`
  render all values in **delayed** form (`:=`, `^:=`, `f /: lhs := rhs`). Adding a
  `delayed` bit to `Rule` to make `=` vs `:=` faithful is noted as an optional
  follow-on (touches core `add_rule` and all readback — deliberately out of scope for
  the trial).

---

## Architecture (where each piece lands)

All new code is wrapped in `#if UP_VALUES … #endif`. `#if UP_VALUES` evaluates an
undefined macro as 0, so an `UP_VALUES=0` build is byte-for-byte today's behaviour.

### Efficiency guard (critical)
`apply_up_values` runs before downvalues on **every** function call, so the no-upvalues
case must be ~free. Maintain a global counter `symtab_up_value_count` (incremented on
install, decremented on remove/clear/set-replace). The eval hook early-outs on
`symtab_up_value_count == 0` **before** even calling `apply_up_values`, and the function
re-checks at entry. Per-candidate we skip any def with `up_values == NULL`, dedupe
candidate defs (small fixed stack array), and reuse the existing arity/first-arg-head
dispatch pre-filter from `apply_down_values_def`.

---

## Phase 0 — Build flag, symbols, struct field, teardown

**Files:** `makefile`, `tests/CMakeLists.txt`, `src/symtab.h`, `src/symtab.c`,
`src/sym_names.h`, `src/sym_names.c`.

- `makefile` (after CFLAGS at line 70): `UP_VALUES ?= 1` then
  `ifeq ($(UP_VALUES),1)` → `CFLAGS += -DUP_VALUES=1`. Consistent flag across all TUs.
- `tests/CMakeLists.txt`: `option(UP_VALUES … ON)` + `add_compile_definitions(UP_VALUES=1)`
  (mirrors the `USE_MPFR` block near line 40) so struct layout matches the lib.
- `src/symtab.h` `SymbolDef` (lines 48–137): add `#if UP_VALUES Rule* up_values; #endif`.
  `calloc`'d node (`node_intern`, symtab.c:68) zero-inits it. Declare `symtab_add_up_value`,
  `symtab_get_up_values`, `apply_up_values`, `symtab_set_up_values/_down_values/_own_values`,
  and extend the removal API with a 3-way selector enum `{OWN,DOWN,UP}` (replace the
  `bool own_value` of `symtab_remove_matching_rule`, symtab.c:960).
- `src/symtab.c` teardown — add an `up_values` walk/free everywhere own/down are handled:
  `reset_node_payload` (:41), `symtab_clear_symbol` (:842). Add
  `extern size_t symtab_up_value_count;` (definition + reset in `symtab_init`).
- `src/sym_names.{h,c}`: intern `SYM_UpValues`, `SYM_UpSet`, `SYM_UpSetDelayed`,
  `SYM_TagSet`, `SYM_TagSetDelayed`, `SYM_TagUnset`, `SYM_Definition`, and
  `SYM_DownValues`, `SYM_OwnValues` (last two not currently cached). Three edits each
  (extern decl, `= NULL` def, `= intern_symbol("…")` in `sym_names_init`), guarded.

**Gate:** `make UP_VALUES=1` and `make UP_VALUES=0` both build clean; no behaviour yet.

## Phase 1 — Symtab layer (storage, application, setters, removal)

**File:** `src/symtab.c` (+ `.h`).

- `symtab_add_up_value(target, pattern, replacement)` — mirror `symtab_add_down_value`
  (:833): `add_rule(&def->up_values, …)` (reuses canonicalize/dedup/specificity-sort),
  then `eval_rule_epoch_mark()` and `symtab_up_value_count++`. `add_rule` already does
  `eval_clock_bump()` — mandatory so a memoized `g[f[x]]` re-evaluates (eval.c fixed-point
  cache at :2360/:2435).
- `symtab_get_up_values(name)` — mirror `symtab_get_down_values` (:1022).
- `apply_up_values(Expr* res)` — borrows `res`, returns a fresh owned rewrite or NULL
  (same contract as `apply_down_values_def`, :1050). Steps: early-out on
  `symtab_up_value_count==0`; collect level-one candidate defs from `res` args (arg itself
  when `EXPR_SYMBOL`; arg's head when `EXPR_FUNCTION` with symbol head — a dedicated
  helper, *not* `input_arg_head_canon` which maps bare symbols to the literal "Symbol");
  dedupe; for each def with `up_values`, scan rules with the arity/first-head pre-filter,
  `match(res, rule->pattern, env)`, on hit `replace_bindings(rule->replacement, env)` and
  return. First match across candidates (left-to-right) wins.
- `symtab_remove_matching_rule` (:960): widen `bool own_value` → selector enum; pick
  `&def->{own,down,up}_values`; decrement `symtab_up_value_count` on an up removal.
- `symtab_set_{up,down,own}_values(name, list)` — new: free the existing target list,
  then for each element of `list` (`Rule`/`RuleDelayed[lhs,rhs]`, stripping an outer
  `HoldPattern` on `lhs`) call the matching `add_rule`/`symtab_add_*`. Adjust
  `symtab_up_value_count` by (new − old) for the up list. `eval_clock_bump` +
  `eval_rule_epoch_mark`.

## Phase 2 — Parser & printer

**Files:** `src/parse.c`, `src/print.c`.

- `src/parse.c`: add `OP_UPSET`, `OP_UPSETDELAYED`, `OP_TAGSET` to the `OperatorType`
  enum (~:915). In `get_operator` (:984): add `^:=` (3-char) **and** `^=` (2-char)
  **before** the bare `^` test at :1117 (else `^:=` mis-lexes as `^` then `:=`); add `/:`
  **before** the bare `/` Divide at :1109. `UpSet`/`UpSetDelayed` take `prec=500,
  right_assoc=1, head_name="UpSet"/"UpSetDelayed"` and ride the **generic infix else**
  (:1799) exactly like Set — no new led branch. `TagSet` (`/:`) needs a **dedicated led
  branch** among the special forms before :1671 (model on `OP_SPAN`, :1638): `left` is the
  tag; parse the middle `lhs` with `parse_expression_prec(s, 501)` so it stops at the
  trailing `=`/`:=`/`=.`; then inspect `s->pos` → build `TagUnset[f,lhs]` (on `=.`, with
  the `!isdigit(pos[2])` guard), `TagSetDelayed[f,lhs,rhs]` (on `:=`), or
  `TagSet[f,lhs,rhs]` (on `=`). All via `expr_new_function(expr_new_symbol("…"), …)`.
- `src/print.c`: `get_expr_prec` (:84) → `SYM_UpSet`/`SYM_UpSetDelayed` = 500, TagSet
  family ≈ 500 (for parenthesisation). `print_standard` (:673): add `SYM_UpSet` → `" ^= "`
  and `SYM_UpSetDelayed` → `" ^:= "` to the infix op-chain (2-arg, joins fine). Add a
  **dedicated branch** for the 3-arg TagSet/TagSetDelayed (`f /: lhs = rhs` / `… := …`)
  and 2-arg `TagUnset` (`f /: lhs =.`), modelled on the `SYM_MessageName` special case
  (:667). Mirror in `print_tex` (:1818) if TeX round-trip is wanted. (Note: `Unset` has
  no postfix printer today — the TagUnset branch is new code, no precedent to copy.)

## Phase 3 — Evaluator wiring

**File:** `src/eval.c`.

- **UpValues application hook** — insert immediately before `apply_down_values_def(hdef,
  res)` at :1865 (after the Orderless block closing at :1855, inside the
  `head->type==EXPR_SYMBOL` block):
  ```c
  #if UP_VALUES
      if (!hold_all_complete && symtab_up_value_count &&
          !is_assignment_primitive(head_name)) {   /* skip Set/Up*/Tag*/Unset wrappers */
          Expr* up = apply_up_values(res);
          if (up) { expr_free(res); *changed = true; return up; }
      }
  #endif
  ```
  `hold_all_complete` is the bool at :1379 (HoldAllComplete ⇒ no upvalues; HoldAll ⇒
  allowed, exactly as required). `is_assignment_primitive` tests `head_name` against the
  assignment SYM_ set so upvalues never fire on an assignment wrapper's held LHS.
- **Assignment dispatch** — extend the inline Set/SetDelayed block (:1967). Factor the
  LHS-target construction at :1977–:2035 into a reusable helper
  `build_assignment_target(Expr* lhs, bool* free_out)` (used by Set **and** the new
  heads — keeps it DRY). Add branches:
  - `UpSet`/`UpSetDelayed` (2-arg): build target from arg0, call new
    `apply_up_assignment(target_lhs, rhs, is_delayed)`; return evaluated `rhs` (UpSet) /
    `Null` (UpSetDelayed).
  - `TagSet`/`TagSetDelayed` (3-arg): arg0 held tag, build target from arg1, evaluate
    arg2 for the immediate case; call `apply_tag_assignment(tag, target_lhs, rhs,
    is_delayed)`.
  - `TagUnset` (2-arg): `apply_tag_unset(tag, target_lhs)` → `Null`.
- **New handlers in `apply_assignment`'s neighbourhood** (eval.c, near :968):
  - `apply_up_assignment`: collect distinct level-one symbols of `lhs`; per symbol check
    `ATTR_PROTECTED` (message `UpSet`/`UpSetDelayed` `::write`/`wrsym`) then
    `symtab_add_up_value(sym, lhs, rhs)`.
  - `apply_tag_assignment`: validate `tag` is a symbol occurring in `lhs`; classify
    own/down/up (per Semantics §3); protection-check the tag; route to
    `symtab_add_own_value`/`_down_value`/`_up_value`. SubValue case → message + no-op.
  - `apply_tag_unset`: classify, then `symtab_remove_matching_rule(tag, lhs, SELECTOR)`;
    `Unset::norep`-style message when nothing matched.
- **`{Up,Down,Own}Values[sym] = list` intercept** — in `apply_assignment` add a branch
  **before** the Protected guard (:1014), mirroring the `Options[sym]=…` intercept
  (:975–:1005): when `lhs` is `UpValues[sym]`/`DownValues[sym]`/`OwnValues[sym]` with a
  symbol arg and `rhs` a List of rules → `symtab_set_{up,down,own}_values`; `eval_clock_bump`.

## Phase 4 — Builtins, attributes, docstrings, scoping, clear/remove

**Files:** `src/core.c`, `src/attr.c`, `src/info.c`, `src/modular.c`.

- **`builtin_up_values`** (template: `rules_to_list` + `builtin_down_values`,
  core.c:2230/:2278): accept a symbol or a string; string naming a symbol with no
  definition → `mth_message("UpValues","sym", …)`. Return
  `{RuleDelayed[HoldPattern[lhs], rhs], …}`.
- Register all new heads in `core_init` next to Clear/Unset/DownValues (core.c:273–:326).
  `UpSet`/`UpSetDelayed`/`TagSet`/`TagSetDelayed`/`TagUnset` are recognised inline in the
  evaluator (like Set), so their registration is a stub builtin returning NULL **or**
  just a symtab/attr entry — but they still need docstrings + attributes.
- **Attributes** in `src/attr.c` `builtin_attrs[]` (format `{"Name", ATTR_…}`, :22):
  `UpValues` = `HoldAll|Protected`; `UpSet` = `HoldFirst|Protected|SequenceHold`;
  `UpSetDelayed` = `HoldAll|Protected|SequenceHold`; `TagSet` = `HoldAll|Protected|SequenceHold`;
  `TagSetDelayed` = `HoldAll|Protected|SequenceHold`; `TagUnset` = `HoldAll|Protected|SequenceHold`;
  `Definition` = `HoldAll|Protected`. (TagSet/TagSetDelayed get HoldAll — the evaluator
  evaluates the rhs explicitly for the immediate case.)
- **Docstrings** (terse, via `symtab_set_docstring` — likely in `src/info.c` near the
  DownValues/OwnValues strings at :3598) for every new head.
- **Clear / ClearAll / Remove / Unset** must also drop up_values:
  `core_clear_all_one` (core.c:1217), `builtin_clear` (:1127), `builtin_unset` (:1145 —
  note plain Unset stays own/down only, per WMA). `symtab_clear_symbol`/`reset_node_payload`
  already extended in Phase 0 cover Clear/Remove's symtab side.
- **Block / Module** save-restore in `src/modular.c`: add `Rule* old_up;` to
  `BlockSavedVar` (:682), save+null in `builtin_block` (:772), free+restore in
  `blk_restore_frame` (:711). Otherwise local upvalues leak.

## Phase 5 — Inspection display (Definition + ?name/Information)

**Files:** `src/core.c`, `src/print.c`, `src/info.c`.

- **`Definition[f]`** is an inert object: `builtin_definition` validates the arg is a
  symbol and returns NULL (stays `Definition[f]`; FullForm shows `Definition[f]`).
  Special-case `head == SYM_Definition` in `print.c` (standard/OutputForm only, **not**
  FullForm) to render `f`'s values in assignment syntax via a shared helper:
  - OwnValues → `f := rhs`
  - DownValues → `lhs := rhs`
  - UpValues → `f /: lhs := rhs` (one canonical form; we don't replicate WMA's `^:=`-vs-`/:`
    display heuristic — documented).
  All in **delayed** form (immediate/delayed not tracked — see divergences).
- **`?name` / `Information`** (`builtin_information`, core.c:2885): after the docstring,
  append the same rendered definitions (reuse the helper). Keep output minimal
  (docstring + values); attributes line optional.

## Phase 6 — Tests & valgrind

**Files:** new `tests/test_upvalues.c`, registered in `tests/CMakeLists.txt`.

Cover **every** example from the user's reference, grouped:
- UpSet/UpSetDelayed: single + multi-symbol install (`prop[a,b[c]]^=v` on both `a`,`b`);
  immediate vs delayed RHS evaluation; `area[sq[s_]]^:=s^2`; all-heads `_[g[x_]]^=1`.
- TagSet/TagSetDelayed: upvalue, redundant-tag downvalue, redundant-tag ownvalue;
  tag-not-in-lhs message; modular-arithmetic end-to-end
  (`mod[2,5]+3 mod[3,5]-mod[1,5] → mod[0,5]`).
- TagUnset / `f/:lhs=.`: removal; `Unset::norep` for plain `f[h[x_]]=.`.
- UpValues reader: symbol arg, `UpValues/@Names["x*"]`, empty list, `UpValues["x"]` message.
- Setters: `UpValues[h]={…}` then evaluate; `UpValues[x]=Reverse[UpValues[x]]` reorder;
  `UpValues[h]=UpValues[g]/.g->h` copy; plus `DownValues[f]=list` / `OwnValues[f]=list`.
- Precedence: upvalue beats outer downvalue (`f[g[2]]→h[2]` with `f[x_]:=1` present).
- Hold interaction: `Hold[g[x]]` fires all-heads upvalue; `HoldComplete[g[x]]` does not.
- Scoping: Block/Module restore upvalues; `Clear`/`Remove` drop them.
- **Efficiency:** assert `symtab_up_value_count==0` fast-path (sanity) and that a hot
  loop with zero upvalues is unaffected.

Run the suite, then `valgrind --leak-check=full` over a script exercising install →
fire → clear → Block, confirming zero leaks / no double-free. Also run existing
`make check-messages` (new diagnostics route through `mth_message`), `make check-c99`,
and the full existing test suite (no regressions with UP_VALUES on **and** off).

## Phase 7 — Docs, version, changelog

- `docs/spec/builtins/assignment-and-rules.md`: document UpSet/UpSetDelayed/TagSet/
  TagSetDelayed/TagUnset. `docs/spec/builtins/expression-information.md`: UpValues,
  Definition, extended Information.
- `docs/spec/changelog/2026-10-05.md` (exists): add a "UpValues subsystem" section.
- `src/version.h`: bump `MATHILDA_VERSION_NUMBER`/`_STRING` `0.294 → 0.295` (both lines).
- Commit message notes `; v0.295`; **commit/tag `v0.295`/push only on the user's
  go-ahead** (not done automatically).
- **Not required for the trial:** NDArray/packed/Compile surfaces — UpValues/UpSet/TagSet
  are structural assignment/inspection heads operating on expressions, not numeric
  kernels, so they are genuinely exempt (state this explicitly; no audit entries needed).
  Book index regen deferred (no book change in this trial).

---

## Critical files (summary)

| File | Change |
|------|--------|
| `makefile`, `tests/CMakeLists.txt` | `-DUP_VALUES=1` default (toggle `UP_VALUES=0`) |
| `src/symtab.{h,c}` | `up_values` field; add/get/apply/set/remove; `symtab_up_value_count`; teardown |
| `src/sym_names.{h,c}` | `SYM_UpValues/UpSet/UpSetDelayed/TagSet/TagSetDelayed/TagUnset/Definition/DownValues/OwnValues` |
| `src/parse.c` | `OP_UPSET/UPSETDELAYED/TAGSET`; lex `^=`,`^:=`,`/:`; TagSet led branch |
| `src/print.c` | prec + `print_standard`/`print_tex` for UpSet/UpSetDelayed/TagSet/TagUnset; Definition render |
| `src/eval.c` | upvalue hook (:1865); assignment dispatch (:1967); handlers + list-setter intercepts |
| `src/core.c` | builtins (UpValues, Definition); register; Clear/ClearAll/Remove/Unset up_values |
| `src/attr.c`, `src/info.c` | attributes + docstrings |
| `src/modular.c` | Block/Module save-restore of up_values |
| `tests/test_upvalues.c` | full spec coverage |
| `docs/spec/...`, `src/version.h` | docs, changelog, version bump |

## Verification
- `make UP_VALUES=1 -j` and `make UP_VALUES=0 -j` both build clean (C99, `-Wall -Wextra`).
- `cd tests/build && cmake .. && make -j && ./test_upvalues` (+ full suite) green.
- Manual REPL check of the reference transcript (modular arithmetic, `area[square]^=s^2`,
  `UpValues[g]`, `Definition[g]`, `?mod`).
- `valgrind --leak-check=full` clean on an install→fire→Block→clear script.
- `make check-messages`, `make check-c99` pass; no regression in the existing suite
  with UP_VALUES both on and off.
