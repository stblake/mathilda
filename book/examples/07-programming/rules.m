# Ch.7 -- Rule/RuleDelayed, Replace, ReplaceAll, ReplaceRepeated
Replace[x^2, x^2 -> a + b]
1 + x^2 /. x^2 -> a + b
{x, x^2, y} /. x -> 2
a + b c /. {a -> 1, c -> 10}
Replace[x, {{x -> 1}, {x -> 2}}]
f[f[f[f[x]]]] /. f[a_] :> a
f[f[f[f[x]]]] //. f[a_] :> a
i = 0; {Replace[{p, q, r}, z_ -> (i = i + 1), {1}], i}
i = 0; {Replace[{p, q, r}, z_ :> (i = i + 1), {1}], i}
