#!/bin/sh
# A stand-in for the Mathilda binary that returns WRONG, malformed or hostile
# certificates, used by test/Untrusted.lean to show that the tactics never
# trust the oracle.  Called as: bad_oracle.sh -file <query script>.
# Goals carrying a marker constant (1001, 1002, ...) get a hostile answer;
# the others get a well-formed certificate for a different (true) statement.
q=$(cat "$2")
echo "@@MATHILDA-BEGIN@@"
case "$q" in
  *1001*)
    # Lean code (the reviewer's by_elab exploit plus a file write): must never run.
    cat <<'ANS'
(by_elab do
   IO.FS.writeFile "/tmp/mathilda_untrusted_pwned" "x"
   let ty ← Lean.Elab.Term.elabTerm (← `(∀ a b : Rat, a - b = 0)) none
   Lean.withOptions (fun o => o.setBool `debug.skipKernelTC true) <|
     Lean.addDecl (.thmDecl { name := `bogus.evil, levelParams := [], type := ty, value := Lean.mkConst ``True.intro })
   Lean.Elab.Term.elabTerm (← `(fun (v1 : Rat) => ([] : List (Rat × Rat × Nat)))) none)
ANS
    ;;
  *1002*) echo '[(1, sorry, 0)]' ;;
  *1003*) echo '[(1, x, 0)]' ;;
  *1004*) echo '[(1, v1 ^ 100000, 0)]' ;;
  *1005*) echo '[(1, v1 +* $, 0)]' ;;
  *1006*) echo "[(1, $(printf '9%.0s' $(seq 1 2500)), 0)]" ;;
  *1007*) echo '[(1, v7, 0)]' ;;
  *1009*) echo '[(by decide, 2, [])]' ;;
  *LinearCombination*) echo '[1, -v1 - v2]' ;;
  *SOS*) echo '[(1, v1 - v2, 0)]' ;;
  *Pratt*) echo '[(3, 2, [(2, 1)]), (7, 3, [(2, 1), (3, 1)]), (17, 3, [(2, 4)]), (274177, 5, [(2, 8), (3, 2), (7, 1), (17, 1)])]' ;;
  *) echo 'garbage' ;;
esac
echo "@@MATHILDA-END@@"
