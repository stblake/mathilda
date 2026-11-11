/-
  Tactics that ask Mathilda for a certificate and check it in Lean.

  Mathilda is an untrusted oracle: it only ever supplies *data* (cofactors,
  squares, a Pratt certificate).  Its answer is never elaborated as Lean
  syntax: `Mathilda.Answer.parse` (a hand-written parser, a pure function on
  `String`) reads it into a data AST of numerals, variables `v<i>`, `+ - * / ^`,
  lists and tuples, and refuses everything else (`by`, `fun`, `sorry`, any other
  identifier, ...).  Each tactic builds the certificate's terms from that AST
  (fixed quotations, numerals and the goal's own free variables), turns them
  into ordinary Lean proof steps (`have … := by grind`, `decide +kernel`) and
  the kernel checks them.  A wrong, malformed or malicious answer makes the
  tactic fail with an error; it can never produce a proof of a false goal, and
  nothing in it is run.

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
import Mathilda.Answer

namespace Mathilda

open Lean Elab Tactic Meta

/-- Ask Mathilda (with the current options, for `mathilda.path`). -/
def oracle (query : String) : TacticM String := do
  callMathilda query (← getOptions)

/-- The error for an answer that is not a well-formed certificate. -/
def badAnswer {α} (tac why : String) : TacticM α :=
  throwError "{tac}: Lean could not parse the certificate returned by Mathilda ({why}); \
    Mathilda is an untrusted oracle and its answer is not used"

/-- Parse Mathilda's answer as data (never as Lean syntax). -/
def parseAns (tac ans : String) : TacticM Answer.Ans :=
  match Answer.parse ans with
  | .ok a => pure a
  | .error e => badAnswer tac e

/-- The syntax of a parsed term: fixed quotations over numerals and the
    goal's variables `vars` (`v<i>` is `vars[i-1]`). -/
partial def polySyntax (tac : String) (vars : Array Term) : Answer.Poly → TacticM Term
  | .num n => pure (Syntax.mkNumLit (toString n))
  | .var i =>
    if h : 0 < i ∧ i ≤ vars.size then pure vars[i - 1]
    else badAnswer tac s!"v{i} names no variable of the goal"
  | .add a b => do `(($(← polySyntax tac vars a) + $(← polySyntax tac vars b)))
  | .sub a b => do `(($(← polySyntax tac vars a) - $(← polySyntax tac vars b)))
  | .mul a b => do `(($(← polySyntax tac vars a) * $(← polySyntax tac vars b)))
  | .div a b => do `(($(← polySyntax tac vars a) / $(← polySyntax tac vars b)))
  | .neg a => do `((-$(← polySyntax tac vars a)))
  | .pow a k => do `(($(← polySyntax tac vars a) ^ $(Syntax.mkNumLit (toString k))))

/-- Elaborate a parsed term at type `ty`. -/
def elabPoly (tac : String) (vars : Array Term) (ty : Expr) (p : Answer.Poly) :
    TacticM Expr := do
  let stx ← polySyntax tac vars p
  let e ← Tactic.elabTermEnsuringType stx (some ty)
  Term.synthesizeSyntheticMVarsNoPostponing
  instantiateMVars e

/-- The goal's variables as terms, for `polySyntax`. -/
def varTerms (fvs : Array Expr) : TacticM (Array Term) := fvs.mapM fun e => Term.exprToSyntax e

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
  let tac := "mathilda_linear_combination"
  let .list xs ← parseAns tac (← oracle query)
    | badAnswer tac "expected a list of cofactors"
  let vars ← varTerms st
  let cs ← xs.toList.mapM fun
    | .poly p => elabPoly tac vars ty p
    | _ => badAnswer tac "a cofactor must be a term"
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
  let tac := "mathilda_sos"
  let .list xs ← parseAns tac (← oracle query)
    | badAnswer tac "expected a list of (d, l, i) triples"
  let vars ← varTerms st
  let ts ← xs.toList.mapM fun
    | .tuple #[.poly d, .poly l, .poly (.num i)] => do
      return (← elabPoly tac vars ty d, ← elabPoly tac vars ty l, i)
    | _ => badAnswer tac "an entry must be a triple (d, l, i) with a numeral i"
  let mut steps : Array (TSyntax `tactic) := #[]
  let mut terms : Array Term := #[]
  let mut k := 0
  for (d, l, i) in ts do
    k := k + 1
    let ds ← Term.exprToSyntax d
    let lsx ← Term.exprToSyntax l
    let s := hname "s" k
    if i == 0 then
      steps := steps.push (← `(tactic| have $s : 0 ≤ ($lsx) ^ 2 := Mathilda.sqnn))
      terms := terms.push (← `($ds * ($lsx) ^ 2))
    else
      let some (hfv, _, g, _, _) := hyps[i - 1]?
        | badAnswer tac s!"index {i} names no hypothesis"
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
  let shown : List MessageData := ts.map fun (d, l, i) => m!"({d}, {l}, {i})"
  checkCertificate "mathilda_sos" m!"squares {shown}" steps

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
  let tac := "mathilda_prime"
  let num (a : Answer.Ans) : TacticM Term := match a with
    | .poly (.num k) => pure (Syntax.mkNumLit (toString k))
    | _ => badAnswer tac "expected a natural-number literal"
  let .list steps ← parseAns tac ans
    | badAnswer tac "expected a list of Pratt steps"
  let steps ← steps.mapM fun
    | .tuple #[q, a, .list fs] => do
      let fs ← fs.mapM fun
        | .tuple #[r, k] => do `(($(← num r), $(← num k)))
        | _ => badAnswer tac "a factor must be a pair (q, k)"
      `(($(← num q), $(← num a), [$fs,*]))
    | _ => badAnswer tac "a Pratt step must be a triple (p, a, [(q, k), ...])"
  let cert ← `([$steps,*])
  checkCertificate "mathilda_prime" m!"Pratt certificate {ans}"
    #[← `(tactic| exact ⟨$cert, by decide +kernel⟩)]

/-- `#mathilda "expr"` evaluates the Mathilda expression `expr` and logs the
    result in `LeanForm` as an information message.  It proves nothing. -/
elab "#mathilda " s:str : command => do
  let ans ← callMathilda ("Lean`Oracle[\"Eval\", " ++ mathildaString s.getString ++ "]")
    (← getOptions)
  logInfo ans

end Mathilda
