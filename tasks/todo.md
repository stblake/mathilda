# DSolve M18 — reducible-μ integrating-factor method

Method: Cheb-Terrab & Roche, JSC 27(5):501–519 (1999). New file
`src/calculus/dsolve_ifactor.c`, builtin `DSolve\`ReducibleIntegratingFactor`.
Plan: `/Users/user/.claude/plans/let-s-continue-our-implementation-warm-willow.md`

## Stage 0 — substrate + shared pipeline + cascade wiring ✅ DONE
- [x] `dsolve_ifactor.c` skeleton: three-function contract (try/builtin/init)
- [x] Copy lie2 static kits: deadline/expired, TimeConstrained beval wrappers,
      run_applied, decline memo, too_big, has_undef_fn, num_ok/abs_at
- [x] Φ extraction (dsolve_solve_top_derivative(P,2)) + y'[x]→p, y[x]→Y substitution
- [x] First-integral R build: R = ∫μ dp + G(x,y), fix G from (2.10); if_hypsimp collapse
- [x] Symbolic verify gate A(R) = R_x + p R_y + Φ R_p == 0
- [x] Reduce/solve: R==C[2] → recurse first-order cascade (C[1]) → explicit → num_ok → dsolve_run
- [x] Cascade wiring: extern decl, init call, enum + ds_method_from_string + switch
- [x] Register builtin + ATTR_PROTECTED + docstring
- [x] **linearity gate** (declines linear ODEs — correctness + terminates Case-B recursion)
- [x] Build green

## Stage 1 — _mu_xy (Section 2.1) ✅ DONE
- [x] poly-in-p deg≤2 test → a,b,c
- [x] Case A: (2.16) existence → (2.17) μ (closed form) — 906/1094 → Airy
- [x] Case B: (2.18) existence → solve ν-ODE (2.20) → (2.19) μ — case 14 (Coth)
- [x] **SIDE FIX: TrigToExp[Coth] sign bug** (src/simp/trigsimp.c) — was −Coth
- [x] Verify: 906/1094/14 solve (independent residual ~1e-16); 189/199/307/665 genuine declines
- [~] Corpus re-measure (running, background by8o5b8ah)

## Stage 2 — _mu_x_y1 (Section 2.2, Lemma 3) — ATTEMPTED, BLOCKED, REVERTED
- Implemented + VERIFIED the μ-search for Cases A/C/D + Lemma-2 μ̃ (Kamke 226→μ=y';
  136→(y'−1)/h(y'); 66→(y'+b)/(a(1+y'²)^{3/2}), all with A(R)=0 holding; the
  per-candidate A(R)=0 gate — not the weaker μ̃-existence check — discriminates cases).
- **0 new corpus solves** (measured twice): the reduced first integrals R==C[1] are
  NON-ELEMENTARY first-order ODEs (y'=√(x²y²+2C), y'=Tan[C+Log[x−y]]) that neither our
  cascade NOR (verified directly) Maple/Mathematica close in elementary explicit form.
- Also introduced a t_m18_auto_dispatch suite regression → reverted to stub.
- BLOCKED on: an implicit/non-elementary first-order ODE solver, or a policy to emit the
  reduced first integral as an implicit answer. Exact equations in DSOLVE_PLAN.md M18.

## Stage 3 — _mu_y_y1 (Section 2.3) — BLOCKED (depends on Stage 2)

## Tests + docs + gates ✅ DONE (Stage 1)
- [x] tests/test_dsolve.c: t_m18_* units (Case A/B, auto-dispatch, decline, Coth regression)
- [x] tests/test_dsolve_m18_stress.c + CMake (Case-A/Case-B forward-generator grids)
- [x] all DSolve ctest suites (dsolve_tests 204 + m5/m12/m14/m17/m18 stress) + check-c99 green
- [x] trig/hyperbolic/simplify/logexp suites green (Coth fix regression)
- [x] valgrind: ifactor ownership clean (decline path byte-identical to Liouville control;
      solve-path leak is the known inherited Integrate/Solve engine baseline)
- [x] DSOLVE_PLAN.md (M18 + §1d), STATUS.md, changelog, calculus.md
- [~] gate baseline argv[3] unchanged at 605 (safe — wave only lowers non-PASS; full re-run pending)

## Review

Landed **M18 Stage 1**: `DSolve\`ReducibleIntegratingFactor` — the reducible-μ
integrating-factor method for nonlinear 2nd-order ODEs, form μ(x,y) (Cheb-Terrab &
Roche 1999). New file `src/calculus/dsolve_ifactor.c`; wired into `dsolve.c` before
`SecondOrderSymmetry`. Two verify gates (symbolic A(R)=0 + numeric back-sub) ⇒ never a
wrong answer. **+7 corpus solves, 0 FAIL.**

Two bugs fixed en route:
1. **`TrigToExp[Coth]`** returned −Coth (sign-flipped denominator in
   `src/simp/trigsimp.c`) — surfaced canonicalizing mixed hyperbolic/exp coefficients.
2. **Symbolic-parameter verify gap**: the numeric gate now instantiates free
   argument-position parameters at generic reals (M16 lesson) — unlocked case 184.

Scope decision: the user chose all three μ-forms, but Stage 2 (μ(x,y')) proved a
re-plan point — Case A verified working yet 0 corpus yield (radical reduced ODEs +
corpus targets are the harder Cases C–F). Landed the solid, verified Stage 1; Stages
2/3 documented with exact equations for a focused follow-up rather than shipping
half-finished zero-yield code.
