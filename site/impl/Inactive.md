---
source: src/core.c
---
**What it is.** `Inactive[h]` is an inert wrapper for a head. It has no C evaluation
handler; it is registered in `core_init` (`src/core.c`) only so the symbol is official and
`ATTR_PROTECTED`. Inertness is automatic: `Inactive[h][args]` has the *compound* head
`Inactive[h]`, which carries no DownValues and no builtin, so the evaluator evaluates the
arguments normally but never fires `h`'s rules — the expression stays a symbolic fixed
point (the same mechanism as `Derivative[n][f][x]`).

**Interaction.** Only the head is held inert; the arguments still evaluate, so
`Inactive[Plus][1 + 1, 3]` first reduces to `Inactive[Plus][2, 3]`. `D` recognises inert
heads — it applies the fundamental theorem of calculus to `Inactive[Integrate]` (see
`src/calculus/deriv.c`) — and `Activate` reverses the wrapper by rewriting `Inactive[h] ->
h` throughout and re-evaluating. The interned pointer `SYM_Inactive` is the handle the
DSolve and derivative machinery use to build and detect inert integrals.

**Complexity / limits.** No cost of its own — it is a structural marker, matched by a
pointer comparison on the head. Its single subtlety is the head/argument asymmetry above:
the head's rules are suppressed, the arguments' are not.
