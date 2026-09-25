# Rule chains compose by themselves: the loop just runs until nothing changes.
g[x_] := h[x]
h[x_] := x + 1
g[10]
Trace[g[10]]
