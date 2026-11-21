# Push and Pop make a stack (last in, first out); Push and PopFirst a queue
st = LinkedList[]; Push[st, a]; Push[st, b]; Push[st, c];
{Pop[st], Pop[st], Normal[st]}
q = LinkedList[{1, 2}]; Push[q, 3];
{PopFirst[q], PopFirst[q], Normal[q]}
e = LinkedList[]; Pop[e]
# A stack checks bracket matching: an opener is pushed, a closer must match the popped opener
closer = <|")" -> "(", "]" -> "[", "}" -> "{"|>;
step[st_, c_] := If[KeyExistsQ[closer, c], Length[st] > 0 && Pop[st] === closer[c], If[MemberQ[Values[closer], c], Push[st, c]]; True]
balancedQ[str_String] := Module[{st = LinkedList[]}, AllTrue[Characters[str], step[st, #] &] && Length[st] == 0]
{balancedQ["{[()()]}"], balancedQ["([)]"], balancedQ["(("]}
# A queue drives a breadth-first search: distance of every vertex from vertex 1
g = <|1 -> {2, 3}, 2 -> {4}, 3 -> {4, 5}, 4 -> {6}, 5 -> {6}, 6 -> {}|>;
bfs[g_, s_] := Module[{q = LinkedList[{s}], dist = <|s -> 0|>, u}, While[Length[q] > 0, u = PopFirst[q]; Scan[If[! KeyExistsQ[dist, #], dist[#] = dist[u] + 1; Push[q, #]] &, g[u]]]; dist]
bfs[g, 1]
