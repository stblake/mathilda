# Module renames its locals to fresh unique names (alpha-renaming):
Module[{x}, Hold[x]]
Module[{x}, Hold[x]]
# Block is dynamic scoping: it saves the global value and restores it after:
n = 10;
Block[{n = 2}, n^2]
n
# With substitutes constants lexically into the body:
With[{a = 5}, a + a]
# Pure functions bind slots by position:
(#^2 &)[7]
(#1 + #2 &)[3, 4]
