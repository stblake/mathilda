/-
  Lean `Expr` → Mathilda input syntax.

  Only ring/field arithmetic is serialized: `+ - * /`, unary `-`, `^` with a
  natural-number literal exponent, numerals, and local variables (fvars).
  `Nat.cast` / `Int.cast` are transparent.  Every other term is rejected with
  `mathilda: unsupported term …`.

  Local variables are sent as `v1, v2, …` (never by their Lean names, which may
  be `x₁`, shadowed, or otherwise not valid Mathilda symbols).  The fvars are
  returned in the same order, so a closed lambda `fun v1 v2 => …` sent back by
  Mathilda can be applied to them directly (see `Mathilda.Oracle`).
-/
import Lean

namespace Mathilda

open Lean Meta

/-- The fvars met so far, in the order of their `v<i>` names. -/
abbrev SerM := StateT (Array Expr) MetaM

/-- The Mathilda name of the fvar `e`, allocating `v<n+1>` on first use. -/
def varName (e : Expr) : SerM String := do
  let vs ← get
  match vs.findIdx? (· == e) with
  | some i => return s!"v{i + 1}"
  | none =>
    set (vs.push e)
    return s!"v{vs.size + 1}"

/-- A natural-number literal (`2`, `OfNat.ofNat 2`, `instOfNatNat 2`), if `e` is one. -/
partial def natLit? (e : Expr) : Option Nat :=
  match e with
  | .lit (.natVal n) => some n
  | .mdata _ b => natLit? b
  | _ =>
    match e.getAppFnArgs with
    | (``OfNat.ofNat, #[_, n, _]) => natLit? n
    | (``instOfNatNat, #[n]) => natLit? n
    | _ => none

/-- Serialize a ring/field term to Mathilda syntax (fully parenthesised). -/
partial def toMathilda (e : Expr) : SerM String := do
  let e ← instantiateMVars e
  if e.isFVar then return (← varName e)
  if let some n := natLit? e then return toString n
  match e.getAppFnArgs with
  | (``HAdd.hAdd, #[_, _, _, _, a, b]) => bin "+" a b
  | (``HSub.hSub, #[_, _, _, _, a, b]) => bin "-" a b
  | (``HMul.hMul, #[_, _, _, _, a, b]) => bin "*" a b
  | (``HDiv.hDiv, #[_, _, _, _, a, b]) => bin "/" a b
  | (``Neg.neg, #[_, _, a]) => return s!"(-{← toMathilda a})"
  | (``HPow.hPow, #[_, _, _, _, a, n]) =>
    match natLit? n with
    | some k => return s!"({← toMathilda a}^{k})"
    | none => throwError "mathilda: unsupported exponent {n} (only natural-number literals)"
  | (``Nat.cast, #[_, _, a]) => toMathilda a
  | (``NatCast.natCast, #[_, _, a]) => toMathilda a
  | (``Int.cast, #[_, _, a]) => toMathilda a
  | (``IntCast.intCast, #[_, _, a]) => toMathilda a
  | _ => throwError "mathilda: unsupported term {e}"
where
  bin (op : String) (a b : Expr) : SerM String := do
    return s!"({← toMathilda a} {op} {← toMathilda b})"

end Mathilda
