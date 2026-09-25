# 6.2 Lists: construction, parts, spans, shape
v = {2, 3, 5, 7, 11, 13};
v[[3]]
v[[-1]]
v[[2 ;; 4]]
v[[1 ;; -1 ;; 2]]
m = {{1, 2, 3}, {4, 5, 6}};
m[[2, 3]]
m[[All, 2]]
Length[v]
Length[m]
Dimensions[m]
Dimensions[{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}}]
Length[f[a, b, c, d]]
