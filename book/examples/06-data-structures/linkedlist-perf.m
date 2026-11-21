# n appends: a List copies on every AppendTo, a LinkedList does not
n = 5*10^4;
$AutoCompilation = False;
tList = First[AbsoluteTiming[s = {}; Do[AppendTo[s, N[i]], {i, n}]]]
tLL = First[AbsoluteTiming[ll = LinkedList[]; Do[AppendTo[ll, N[i]], {i, n}]]]
Normal[ll] === s
# The same loop at 10^6, interpreted and then compiled behind the scenes
ll = LinkedList[]; tI = First[AbsoluteTiming[Do[AppendTo[ll, N[i]], {i, 10^6}]]]
$AutoCompilation = True;
ll = LinkedList[]; tC = First[AbsoluteTiming[Do[AppendTo[ll, N[i]], {i, 10^6}]]]
Round[tI/tC]
# A compiled function mutates its caller's list in place
sq = Compile[{{b, _LinkedList, _Real}, {k, _Integer}}, Do[Push[b, N[i]^2], {i, k}]; Length[b]];
q = LinkedList[]; sq[q, 5]
Normal[q]
