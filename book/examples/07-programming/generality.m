# Ch.7 -- DownValues are ordered by pattern generality (v0.6022), not a specificity score
k[x_Integer] := "integer"
k[x_?EvenQ] := "even"
k[4]
p[a_, b_] := "any pair"
p[a_, a_] := "equal pair"
{p[1, 2], p[3, 3]}
DownValues[p]
s[x_ /; x > 0] := "positive"
s[x_] := "anything"
s[x_?NumberQ] := "number"
{s[5], s[-5], s[a]}
