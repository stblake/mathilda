# The rule engine transforms an expression without defining a symbol.
x^2 + y^2 /. x -> a
{1, 2, 3, 4} /. n_Integer :> n^2
# ReplaceAll rewrites top-down and does not re-descend into a match:
f[f[x]] /. f[a_] -> a
# ReplaceRepeated iterates the same rule to a fixed point:
f[f[x]] //. f[a_] -> a
# Replace acts only at the top level by default:
Replace[f[x], f[a_] -> a]
Replace[{f[x], f[y]}, f[a_] -> a]
# ReplaceList returns EVERY way a pattern can match:
ReplaceList[a + b + c, x_ + y_ -> {x, y}]
