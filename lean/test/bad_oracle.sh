#!/bin/sh
# A stand-in for the Mathilda binary that returns WRONG certificates, used by
# test/Untrusted.lean to show that the tactics never trust the oracle.
# Called as: bad_oracle.sh -file <query script>.
q=$(cat "$2")
echo "@@MATHILDA-BEGIN@@"
case "$q" in
  *LinearCombination*) echo 'fun (v1 v2 : _) => ([1, -v1 - v2] : List _)' ;;
  *SOS*) echo 'fun (v1 v2 : _) => ([(1, v1 - v2, 0)] : List (_ × _ × Nat))' ;;
  *Pratt*) echo '[(3, 2, [(2, 1)]), (7, 3, [(2, 1), (3, 1)]), (17, 3, [(2, 4)]), (274177, 5, [(2, 8), (3, 2), (7, 1), (17, 1)])]' ;;
  *) echo 'garbage' ;;
esac
echo "@@MATHILDA-END@@"
