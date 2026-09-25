# An integer accumulation loop: the counted Do compiles to a single LOOP
# opcode over integer opcodes (POWI_I, ADD_I) -- no all-Real fast path here.
Compile[{{n, _Integer}}, Module[{s = 0}, Do[s = s + i^2, {i, 1, n}]; s]]
