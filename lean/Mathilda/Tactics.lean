/-
  Tactics that ask Mathilda for a certificate and check it in Lean.

  Mathilda is an untrusted oracle: it only ever supplies *data* (cofactors,
  squares, a Pratt certificate).  Each tactic turns that data into ordinary
  Lean proof steps (`have … := by grind`, `decide +kernel`) and the kernel
  checks them.  A wrong or malicious answer makes the tactic fail; it can never
  produce a proof of a false goal.

  * `mathilda_linear_combination` — goal `a = b`; hypotheses `hᵢ : lᵢ = rᵢ` of
    the same type in the local context.  Mathilda returns cofactors `cᵢ` with
    `a - b = Σ cᵢ (lᵢ - rᵢ)` (ideal membership); Lean checks the identity with
    `grind` and closes the goal from it.
  * `mathilda_sos` — goal `0 ≤ p`, `p ≤ q` or `p ≥ q`; hypotheses `hᵢ : 0 ≤ gᵢ`.
    Mathilda returns `(d, ℓ, i)` with `q - p = Σ d ℓ² gᵢ` (`g₀ = 1`, `d ≥ 0`);
    Lean adds the hints `0 ≤ ℓ²` / `0 ≤ ℓ² gᵢ` and checks with `grind`.
  * `mathilda_prime` — goal `∃ c, Mathilda.PrattCertifies p c`.  Mathilda
    returns a flat Pratt certificate; `decide +kernel` evaluates
    `Mathilda.prattCheck` on it.
  * `#mathilda "expr"` — evaluates a Mathilda expression and logs the result in
    `LeanForm`.  It proves nothing.
-/
import Lean
import Mathilda.Prelude
import Mathilda.Serialize
import Mathilda.Oracle

namespace Mathilda

open Lean Elab Tactic Meta

/-- Ask Mathilda (with the current options, for `mathilda.path`). -/
def oracle (query : String) : TacticM String := do
  callMathilda query (← getOptions)

/-- Parse a term returned by Mathilda. -/
def parseAnswer (ans : String) : TacticM Syntax := do
  match Parser.runParserCategory (← getEnv) `term ans with
  | .ok stx => return stx
  | .error e => throwError "mathilda: cannot parse Mathilda's answer `{ans}`: {e}"

/-- The elements of a list literal `[a, b, …]` (after `whnfR` at each cons cell). -/
partial def listElems (e : Expr) : MetaM (List Expr) := do
  let e ← whnfR (← instantiateMVars e)
  match e.getAppFnArgs with
  | (``List.nil, #[_]) => return []
  | (``List.cons, #[_, a, rest]) => return a :: (← listElems rest)
  | _ => throwError "mathilda: expected a list literal from Mathilda, got {e}"

/-- Elaborate Mathilda's answer `fun (v1 … vn : _) => (list : List τ)` against
    `∀ vᵢ : type(fvᵢ), List elemTy`, apply it to the fvars and return the list
    elements (beta-reduced, so they mention the real fvars). -/
def elabAnswer (ans : String) (fvs : Array Expr) (elemTy : Expr) : TacticM (List Expr) := do
  let stx ← parseAnswer ans
  let listTy ← mkAppM ``List #[elemTy]
  let ty ← fvs.foldrM (fun fv acc => do mkArrow (← inferType fv) acc) listTy
  let f ← Tactic.elabTermEnsuringType stx (some ty)
  Term.synthesizeSyntheticMVarsNoPostponing
  let f ← instantiateMVars f
  listElems (mkAppN f fvs).headBeta

/-- The components of a pair `(a, b)` (`Prod.mk a b`, after `whnfR`). -/
def pairOf? (e : Expr) : MetaM (Option (Expr × Expr)) := do
  match (← whnfR (← instantiateMVars e)).getAppFnArgs with
  | (``Prod.mk, #[_, _, a, b]) => return some (a, b)
  | _ => return none

/-- Mathilda text of the serialized fvars: `{v1, v2, …}`. -/
def varList (n : Nat) : String :=
  "{" ++ ", ".intercalate ((List.range n).map (fun i => s!"v{i + 1}")) ++ "}"

/-- A Mathilda list literal. -/
def mlist (xs : List String) : String := "{" ++ ", ".intercalate xs ++ "}"

/-- Serialize `e` with the variable state `st`; `none` (state unchanged) if unsupported. -/
def trySerialize (e : Expr) (st : Array Expr) : MetaM (Option (String × Array Expr)) := do
  try
    let (s, st') ← (toMathilda e).run st
    return some (s, st')
  catch _ => return none

/-- `e = rhs` / `l ≤ r` facts of the local context, serialized with the shared
    variable state.  `pick` returns the two sides of a usable hypothesis type. -/
def collectHyps (ty : Expr) (pick : Expr → MetaM (Option (Expr × Expr)))
    (st : Array Expr) : MetaM (Array (Expr × Expr × Expr × String × String) × Array Expr) := do
  let mut out := #[]
  let mut st := st
  for d in (← getLCtx) do
    if d.isImplementationDetail then continue
    let t ← instantiateMVars d.type
    let some (l, r) ← pick t | continue
    unless ← isDefEq (← inferType l) ty do continue
    let some (ls, st1) ← trySerialize l st | continue
    let some (rs, st2) ← trySerialize r st1 | continue
    st := st2
    out := out.push (d.toExpr, l, r, ls, rs)
  return (out, st)

/-- Run a list of proof steps; on failure restore the state and report that
    Lean rejected Mathilda's certificate. -/
def checkCertificate (tac : String) (what : MessageData) (steps : Array (TSyntax `tactic)) :
    TacticM Unit := do
  let s ← saveState
  try
    -- No error recovery: a failing step (also a nested `by grind`) must throw,
    -- not be logged and replaced by `sorry`.
    Term.withoutErrToSorry <| withoutRecover do
      for t in steps do evalTactic t
      unless (← getUnsolvedGoals).isEmpty do throwError "goals remain"
  catch _ =>
    s.restore
    throwError "{tac}: Lean rejected the certificate returned by Mathilda ({what}); \
      Mathilda is an untrusted oracle and its answer does not check"

/-- An identifier `mth_<base><i>`. -/
def hname (base : String) (i : Nat) : Ident := mkIdent (Name.mkSimple s!"mth_{base}{i}")

/-- `mathilda_linear_combination` closes `a = b` from equations `lᵢ = rᵢ` in the
    local context, using cofactors computed by Mathilda and checked by Lean. -/
elab "mathilda_linear_combination" : tactic => withMainContext do
  let goal ← instantiateMVars (← getMainTarget)
  let some (ty, a, b) := goal.eq?
    | throwError "mathilda_linear_combination: the goal must be an equation a = b, got{indentExpr goal}"
  let (as, st) ← (toMathilda a).run #[]
  let (bs, st) ← (toMathilda b).run st
  let (hyps, st) ← collectHyps ty (fun t => return (t.eq?.map fun (_, l, r) => (l, r))) st
  let query := "Lean`Oracle[\"LinearCombination\", " ++ mlist [as, bs] ++ ", " ++
    mlist (hyps.toList.map fun (_, _, _, ls, rs) => mlist [ls, rs]) ++ ", " ++
    varList st.size ++ "]"
  let cs ← elabAnswer (← oracle query) st ty
  unless cs.length == hyps.size do
    throwError "mathilda_linear_combination: Mathilda returned {cs.length} cofactors \
      for {hyps.size} hypotheses"
  let mut steps : Array (TSyntax `tactic) := #[]
  let mut terms : Array Term := #[]
  let mut ks : Array Ident := #[]
  for i in [0:hyps.size] do
    let (_, l, r, _, _) := hyps[i]!
    let ls ← Term.exprToSyntax l
    let rs ← Term.exprToSyntax r
    let cstx ← Term.exprToSyntax cs[i]!
    let k := hname "k" (i + 1)
    ks := ks.push k
    steps := steps.push (← `(tactic| have $k : $ls - $rs = 0 := (by grind)))
    terms := terms.push (← `($cstx * ($ls - $rs)))
  let asx ← Term.exprToSyntax a
  let bsx ← Term.exprToSyntax b
  let sum ← if h : terms.size > 0 then
      terms[1:].foldlM (fun acc t => `($acc + $t)) terms[0]
    else `(0)
  let hc := mkIdent `mth_hc
  steps := steps.push (← `(tactic| have $hc : $asx - $bsx = $sum := (by grind)))
  for k in ks do
    steps := steps.push (← `(tactic| try rw [$k:ident] at $hc:ident))
  steps := steps.push (← `(tactic| grind))
  checkCertificate "mathilda_linear_combination" m!"cofactors {cs}" steps

/-- The two sides of `0 ≤ g` (`g ≥ 0`), if `t` is such a fact. -/
def nonnegFact? (t : Expr) : Option (Expr × Expr) :=
  match t.getAppFnArgs with
  | (``LE.le, #[_, _, z, g]) => if natLit? z == some 0 then some (z, g) else none
  | (``GE.ge, #[_, _, g, z]) => if natLit? z == some 0 then some (z, g) else none
  | _ => none

/-- `mathilda_sos` proves `0 ≤ p`, `p ≤ q` or `p ≥ q` from hypotheses `0 ≤ gᵢ`,
    with a sum-of-squares certificate computed by Mathilda and checked by Lean. -/
elab "mathilda_sos" : tactic => withMainContext do
  let goal ← instantiateMVars (← getMainTarget)
  let (ty, lo, hi) ← match goal.getAppFnArgs with
    | (``LE.le, #[ty, _, lo, hi]) => pure (ty, lo, hi)
    | (``GE.ge, #[ty, _, hi, lo]) => pure (ty, lo, hi)
    | _ => throwError "mathilda_sos: the goal must be `p ≤ q`, `p ≥ q` or `0 ≤ p`, got{indentExpr goal}"
  let (los, st) ← (toMathilda lo).run #[]
  let (his, st) ← (toMathilda hi).run st
  let (hyps, st) ← collectHyps ty (fun t => return nonnegFact? t) st
  let query := "Lean`Oracle[\"SOS\", " ++ mlist [los, his] ++ ", " ++
    mlist (hyps.toList.map fun (_, _, _, _, gs) => gs) ++ ", " ++ varList st.size ++ "]"
  let triple ← mkAppM ``Prod #[ty, ← mkAppM ``Prod #[ty, mkConst ``Nat]]
  let ts ← elabAnswer (← oracle query) st triple
  let mut steps : Array (TSyntax `tactic) := #[]
  let mut terms : Array Term := #[]
  let mut k := 0
  for t in ts do
    k := k + 1
    let some (d, rest) ← pairOf? t | throwError "mathilda_sos: malformed certificate entry {t}"
    let some (l, idx) ← pairOf? rest | throwError "mathilda_sos: malformed certificate entry {t}"
    let some i ← (evalNat idx).run | throwError "mathilda_sos: non-literal index {idx}"
    let ds ← Term.exprToSyntax d
    let lsx ← Term.exprToSyntax l
    let s := hname "s" k
    if i == 0 then
      steps := steps.push (← `(tactic| have $s : 0 ≤ ($lsx) ^ 2 := Mathilda.sqnn))
      terms := terms.push (← `($ds * ($lsx) ^ 2))
    else
      let some (hfv, _, g, _, _) := hyps[i - 1]?
        | throwError "mathilda_sos: certificate index {i} names no hypothesis"
      let gs ← Term.exprToSyntax g
      let hs ← Term.exprToSyntax hfv
      steps := steps.push (← `(tactic| have $s : 0 ≤ ($lsx) ^ 2 * ($gs) :=
        Mathilda.mulnn Mathilda.sqnn $hs))
      terms := terms.push (← `($ds * (($lsx) ^ 2 * ($gs))))
  let los' ← Term.exprToSyntax lo
  let his' ← Term.exprToSyntax hi
  let sum ← if h : terms.size > 0 then
      terms[1:].foldlM (fun acc t => `($acc + $t)) terms[0]
    else `(0)
  let e := mkIdent `mth_e
  steps := steps.push (← `(tactic| have $e : $his' - $los' = $sum := (by grind)))
  steps := steps.push (← `(tactic| grind))
  checkCertificate "mathilda_sos" m!"squares {ts}" steps

/-- The prime `p` of a goal `∃ c, Mathilda.PrattCertifies p c`. -/
def prattGoal? (goal : Expr) : MetaM (Option Expr) := do
  match goal.getAppFnArgs with
  | (``Exists, #[_, .lam _ _ body _]) =>
    match body.getAppFnArgs with
    | (``Mathilda.PrattCertifies, #[p, .bvar 0]) =>
      return if p.hasLooseBVars then none else some p
    | _ => return none
  | _ => return none

/-- `mathilda_prime` proves `∃ c, Mathilda.PrattCertifies p c` (p prime) with a
    Pratt certificate computed by Mathilda and checked by `decide +kernel`. -/
elab "mathilda_prime" : tactic => withMainContext do
  let goal ← instantiateMVars (← getMainTarget)
  let some p ← prattGoal? goal
    | throwError "mathilda_prime: the goal must be `∃ c, Mathilda.PrattCertifies p c`, got{indentExpr goal}"
  let some n ← (evalNat p).run
    | throwError "mathilda_prime: `{p}` is not a natural-number literal"
  let ans ← oracle s!"Lean`Oracle[\"Pratt\", {n}]"
  let cert : Term := ⟨← parseAnswer ans⟩
  checkCertificate "mathilda_prime" m!"Pratt certificate {ans}"
    #[← `(tactic| exact ⟨$cert, by decide +kernel⟩)]

/-- `#mathilda "expr"` evaluates the Mathilda expression `expr` and logs the
    result in `LeanForm` as an information message.  It proves nothing. -/
elab "#mathilda " s:str : command => do
  let ans ← callMathilda ("Lean`Oracle[\"Eval\", " ++ mathildaString s.getString ++ "]")
    (← getOptions)
  logInfo ans

end Mathilda
