# A scalar real polynomial in Horner form (Chapter 8, section 8.2).
# CompilePrint dumps the bytecode: monomorphic typed opcodes, constants
# folded into ADD_RK/MUL_R immediates, three registers, RET.
Compile[{{x, _Real}}, (((x - 3.) x + 2.) x - 1.) x + 0.5]
