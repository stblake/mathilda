# A conditional: control flow lowers to conditional/unconditional jumps
# (JZ / JMP), with the '>' markers on the instructions that are jump targets.
Compile[{{x, _Real}}, If[x > 0., Sqrt[x], -Sqrt[-x]]]
