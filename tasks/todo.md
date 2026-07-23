# Fix: Limit[2 ArcTan[Sqrt[(1+x)/(1-x)]], x -> 1] returns -Pi (should be Pi)

## Root-cause chain
1. `ArcTan[ComplexInfinity]` folds to `-Pi/2` because `exact_arctan` spuriously
   matches `Tan[-Pi/2] == ComplexInfinity` (n=-1 is tried before n=+1). WL gives
   `Indeterminate`.
2. Limit's layer1 numeric substitution does `ArcTan[Sqrt[2/0]]` =
   `ArcTan[ComplexInfinity]` = `-Pi/2`, accepts the (finite, wrong) value, and
   never reaches the series/compose layers. `2 * -Pi/2 = -Pi`.
3. `layer_compose_at_infinity` only accepts *real* +/-Infinity inner limits, so
   even after deferral it cannot fold `ArcTan[DirectedInfinity[I Sqrt[2]]]`.
4. No constant-factor linearity, so `2 ArcTan[...]` never reduces to
   `2 Limit[ArcTan[...]]` (`Limit[3 ArcTan[t], t->Infinity]` is also unevaluated).

Both one-sided limits genuinely -> Pi, so the two-sided answer is Pi.

## Fixes
- [x] trig.c `builtin_arctan`: handle infinite args before exact folds.
      `ArcTan[ComplexInfinity] -> Indeterminate`;
      `ArcTan[DirectedInfinity[dir]] -> +/-Pi/2` by the quadrant of `dir`
      (`arctan_direction_sign` via N + get_approx, so `I Sqrt[2]` resolves).
- [x] limit.c `layer_compose_at_infinity`: also accept ComplexInfinity /
      DirectedInfinity inner limits, applying the head; accept only an
      unambiguous folded value (reject Indeterminate / residual infinity).
- [x] limit.c new `layer_constant_factor` (late in cascade): pull x-free
      factors out of a Times and recurse.

## Verify
- [x] Target -> Pi; FromBelow/FromAbove -> Pi.
- [x] `Limit[3 ArcTan[t], t->Infinity] -> 3Pi/2`, `ArcTan[ComplexInfinity] ->
      Indeterminate`, `ArcTan[DirectedInfinity[I]] / [I Sqrt[2]] -> Pi/2`,
      `ArcTan[DirectedInfinity[-1]] / [-I] -> -Pi/2`.
- [x] Unchanged real-infinity paths (`ArcTan[x]`, `Tanh[x]`, `Erf[x]`,
      `ArcTan[x^2-x^4]` at ±Infinity) and constant-factor finite limits
      (`2x@3`, `5 Sin[x]/x @0`) still correct.
- [x] limit/trig/hyperbolic suites: identical soft-fail counts to pristine
      main (2 pre-existing, unrelated; 0 new). numeric_tests pass.
- [x] valgrind: leak profile identical to trivial-input baseline
      (13,440B/420blk = macOS dyld/libobjc noise); no stack in my code.

## Review
Root defect was a wrong-branch fold: `exact_arctan` matched the `Tan[-Pi/2]`
pole for `ArcTan[ComplexInfinity]`, and `Limit`'s numeric-substitution fast
path trusted the resulting `-Pi/2`. Fixed at the source (correct `ArcTan`
values at all infinities) rather than by guarding the substitution, then
taught `Limit` to (a) fold heads over directed/complex inner infinities and
(b) factor out constants. The composite now resolves because both one-sided
limits legitimately reach `Pi`. Minimal blast radius: the constant-factor
layer runs last (only rescues previously-unresolved shapes), and the compose
change preserves the existing real-infinity path byte-for-byte.
