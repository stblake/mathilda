# A rank-1 array body: the element-wise chain is FUSED into one strip-mined
# pass over 64-wide tiles (VLOAD_R/VPOWI_R/VADD_R/VSTORE_R) driven by a single
# parallel loop (APAR). This is where array speed comes from -- not from
# removing interpretation, but from removing per-element dispatch.
Compile[{{u, _Real, 1}}, u^2 + 1.]
