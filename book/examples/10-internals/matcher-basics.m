# MatchQ drives the pattern matcher; these show what it accepts.
MatchQ[5, _Integer]
MatchQ[5.0, _Integer]
MatchQ[f[a, b], f[_, _]]
# A non-linear pattern: a repeated name must bind to equal expressions.
MatchQ[f[a, a], f[x_, x_]]
MatchQ[f[a, b], f[x_, x_]]
# Sequence blanks match runs of arguments, found by backtracking:
MatchQ[f[a, b, c], f[x__, y__]]
Cases[{1, a, 2, b, 3}, _Integer]
Count[{1, a, 2, b, 3}, _Integer]
# Orderless matching tries argument orderings:
MatchQ[a + b, x_ + y_]
# PatternTest and Condition constrain a match, calling the evaluator back:
MatchQ[4, _?EvenQ]
MatchQ[5, _?EvenQ]
MatchQ[6, x_ /; x > 5]
# Optional supplies a default when the argument is absent:
g[x_, y_: 10] := {x, y}
g[1, 2]
g[1]
