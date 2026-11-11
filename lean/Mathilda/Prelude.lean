namespace Mathilda

def powModAux : Nat → Nat → Nat → Nat → Nat → Nat
  | 0, _, _, _, acc => acc
  | fuel+1, b, e, m, acc =>
    if e = 0 then acc
    else powModAux fuel (b * b % m) (e / 2) m (if e % 2 = 1 then acc * b % m else acc)

def powMod (b e m : Nat) : Nat := powModAux (e + 1) (b % m) e m (1 % m)

def prodPow : List (Nat × Nat) → Nat
  | [] => 1
  | (q, k) :: r => q ^ k * prodPow r

def lucasStep (known : List Nat) (p a : Nat) (fs : List (Nat × Nat)) : Bool :=
  decide (2 < p) && prodPow fs == p - 1 && powMod a (p - 1) p == 1 &&
  fs.all (fun qk => powMod a ((p - 1) / qk.1) p != 1 && known.contains qk.1)

def prattCheck : List Nat → List (Nat × Nat × List (Nat × Nat)) → Bool
  | _, [] => true
  | known, (p, a, fs) :: rest => lucasStep known p a fs && prattCheck (p :: known) rest

def IsPrime (n : Nat) : Prop := 2 ≤ n ∧ ∀ d, d ∣ n → d = 1 ∨ d = n

abbrev PrattCertifies (p : Nat) (cert : List (Nat × Nat × List (Nat × Nat))) : Prop :=
  prattCheck [2] cert = true ∧ (cert.getLast?.map (·.1)) = some p

abbrev Mat := List (List Rat)

def matMul (A B : Mat) : Mat :=
  A.map fun r => (List.range ((B.headD []).length)).map fun j =>
    (List.zipWith (· * ·) r (B.map (·.getD j 0))).sum

def matVec (A : Mat) (v : List Rat) : List Rat :=
  A.map fun r => (List.zipWith (· * ·) r v).sum

abbrev HasShape (A : Mat) (r c : Nat) : Prop := A.length = r ∧ ∀ row ∈ A, row.length = c

abbrev HasLength (v : List Rat) (n : Nat) : Prop := v.length = n

def idMat (n : Nat) : Mat :=
  (List.range n).map fun i => (List.range n).map fun j => if i = j then 1 else 0

def dropCol (j : Nat) (A : Mat) : Mat :=
  A.map fun r => r.eraseIdx j

def detAux : Nat → Mat → Rat
  | 0, _ => 1
  | _, [] => 1
  | fuel+1, r :: rest =>
    ((List.range r.length).map fun j =>
      (if j % 2 = 0 then (1 : Rat) else -1) * r.getD j 0 * detAux fuel (dropCol j rest)).sum

def det (A : Mat) : Rat := detAux A.length A

abbrev sqnn := @Lean.Grind.OrderedRing.sq_nonneg

abbrev mulnn := @Lean.Grind.OrderedRing.mul_nonneg

end Mathilda
