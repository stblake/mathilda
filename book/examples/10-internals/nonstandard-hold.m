# Hold attributes let a head decide its arguments are not evaluated.
Hold[1 + 1]
Attributes[Hold]
# Evaluate punches through a Hold; the other argument stays held:
Hold[Evaluate[1 + 1], 2 + 2]
# Unevaluated passes an argument through unevaluated, then vanishes:
Length[Unevaluated[1 + 2 + 3]]
# ReleaseHold strips the wrapper so the contents evaluate:
ReleaseHold[Hold[1 + 1]]
# Set evaluates its right-hand side now; SetDelayed defers it:
a = 1 + 1
b := 1 + 1
a
b
# If has HoldRest, so the branch not taken is never evaluated (no error here):
If[True, ok, boom = 1/0]
boom
