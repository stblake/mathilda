# Ch.7 -- more kinds of rule and definition (U2)
# Dispatch indexes a large rule set by left-hand side
sq = Table[i -> i^2, {i, 4000}];
d = Dispatch[sq];
data = RandomInteger[{1, 4000}, 2000];
{First[AbsoluteTiming[r1 = data /. sq;]], First[AbsoluteTiming[r2 = data /. d;]]}
r1 === r2
# ReplaceRepeated stops after MaxIterations
ReplaceRepeated[x, x -> x + 1, MaxIterations -> 5]
# Replace matches the whole expression; only /. and //. rewrite part of a Flat head
Replace[a + b + c, a + b -> z]
a + b + c /. a + b -> z
# SubValues: a definition on a compound head f[a][x]
scale[a_][x_] := a x
scale[3][5]
SubValues[scale]
# NValues: a numerical value for a symbol
N[g0] = 9.80665;
{g0, N[2 g0]}
NValues[g0]
# FormatValues: change how an expression prints, not what it is
Format[cx[a_, b_]] := a + b I
cx[1, 2]
FullForm[cx[1, 2]]
