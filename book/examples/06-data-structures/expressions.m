# 6.1 Everything is an expression
# FullForm strips away the surface syntax and shows the raw head-and-arguments tree.
FullForm[a + b c]
FullForm[{1, 2, 3}]
FullForm[1/2]
# Head names the kind of every value, atom or compound.
Head[a + b c]
Head[b c]
Head[42]
Head[3.5]
Head[1/2]
Head[2 + 3 I]
Head["text"]
Head[{1, 2, 3}]
Head[f[x]]
# AtomQ separates the indivisible leaves from the compound nodes.
AtomQ[42]
AtomQ[x]
AtomQ[1/2]
AtomQ[f[x]]
AtomQ[{1, 2}]
# Depth and Level read the shape of the tree itself.
Depth[f[g[h[x]]]]
Depth[{{1, 2}, {3, 4}}]
Level[a + f[x, y^n], {-1}]
Level[a + f[x, y^n], 2]
