/-
  Mathilda core-only prelude (no imports, no Mathlib).

  Inlined (by src/internal/lean/Lean.m) into every generated certificate
  file, and mirrored byte-for-byte as lean/Mathilda/Prelude.lean for the
  Lean-side package. Everything here is a
  computable definition that the kernel can evaluate under `decide +kernel`:
  recursion is structural (an explicit fuel argument where the natural measure
  is not structural), never well-founded, so the kernel reduces it.
-/
namespace Mathilda

/-- Square-and-multiply `b ^ e % m`; `fuel` bounds the number of loop steps. -/
def powModAux : Nat → Nat → Nat → Nat → Nat → Nat
  | 0, _, _, _, acc => acc
  | fuel+1, b, e, m, acc =>
    if e = 0 then acc
    else powModAux fuel (b * b % m) (e / 2) m (if e % 2 = 1 then acc * b % m else acc)

/-- `b ^ e % m` for every exponent `e`.  The loop halves `e` and stops at 0, so it
    runs `log2 e + 1 ≤ e + 1` steps: fuel `e + 1` never runs out (a fixed fuel would
    silently truncate large exponents).  The kernel only peels as many successors off
    the fuel literal as there are steps, so a huge `e` costs nothing extra. -/
def powMod (b e m : Nat) : Nat := powModAux (e + 1) (b % m) e m (1 % m)

/-- `Π q ^ k` over a factor list `[(q, k), ...]`. -/
def prodPow : List (Nat × Nat) → Nat
  | [] => 1
  | (q, k) :: r => q ^ k * prodPow r

/-- One Lucas step: `a` has order `p - 1` modulo `p`, the factor list multiplies out to
    `p - 1`, and every prime factor `q` is already known to be prime. -/
def lucasStep (known : List Nat) (p a : Nat) (fs : List (Nat × Nat)) : Bool :=
  decide (2 < p) && prodPow fs == p - 1 && powMod a (p - 1) p == 1 &&
  fs.all (fun qk => powMod a ((p - 1) / qk.1) p != 1 && known.contains qk.1)

/-- Flat Pratt certificate: steps `(p, a, factors of p - 1)` in dependency order;
    2 is known prime. -/
def prattCheck : List Nat → List (Nat × Nat × List (Nat × Nat)) → Bool
  | _, [] => true
  | known, (p, a, fs) :: rest => lucasStep known p a fs && prattCheck (p :: known) rest

/-- `cert` is a checked Pratt certificate whose last step is `p`. -/
abbrev PrattCertifies (p : Nat) (cert : List (Nat × Nat × List (Nat × Nat))) : Prop :=
  prattCheck [2] cert = true ∧ (cert.getLast?.map (·.1)) = some p

/-- Dense matrices over `Rat` as row lists. -/
abbrev Mat := List (List Rat)

/-- Matrix product (the column count is taken from the first row of `B`). -/
def matMul (A B : Mat) : Mat :=
  A.map fun r => (List.range ((B.headD []).length)).map fun j =>
    (List.zipWith (· * ·) r (B.map (·.getD j 0))).sum

/-- Matrix-vector product. -/
def matVec (A : Mat) (v : List Rat) : List Rat :=
  A.map fun r => (List.zipWith (· * ·) r v).sum

/-- `A` has exactly `r` rows, each of exactly `c` entries (`matMul` / `matVec` /
    `det` use `zipWith` and `getD`, which silently truncate or pad ragged input,
    so every matrix statement also states the shapes). -/
abbrev HasShape (A : Mat) (r c : Nat) : Prop := A.length = r ∧ ∀ row ∈ A, row.length = c

/-- The vector `v` has exactly `n` entries. -/
abbrev HasLength (v : List Rat) (n : Nat) : Prop := v.length = n

/-- The `n × n` identity matrix. -/
def idMat (n : Nat) : Mat :=
  (List.range n).map fun i => (List.range n).map fun j => if i = j then 1 else 0

/-- Drop column `j` of every row. -/
def dropCol (j : Nat) (A : Mat) : Mat :=
  A.map fun r => r.eraseIdx j

/-- Laplace expansion along row 0; `fuel` (the dimension) makes the recursion structural. -/
def detAux : Nat → Mat → Rat
  | 0, _ => 1
  | _, [] => 1
  | fuel+1, r :: rest =>
    ((List.range r.length).map fun j =>
      (if j % 2 = 0 then (1 : Rat) else -1) * r.getD j 0 * detAux fuel (dropCol j rest)).sum

/-- Determinant of a square matrix. -/
def det (A : Mat) : Rat := detAux A.length A

/-- `0 ≤ t ^ 2` in any ordered ring (`Lean.Grind.OrderedRing.sq_nonneg`). -/
abbrev sqnn := @Lean.Grind.OrderedRing.sq_nonneg

/-- `0 ≤ a → 0 ≤ b → 0 ≤ a * b` (`Lean.Grind.OrderedRing.mul_nonneg`). -/
abbrev mulnn := @Lean.Grind.OrderedRing.mul_nonneg

end Mathilda
