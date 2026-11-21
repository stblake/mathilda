# Every symbol lives in a context; $Context is where new ones are created.
$Context
Context[Sin]
y = 5;
Context[y]
# Begin switches the current context; a symbol created here is namespaced:
Begin["MyPkg`"];
Context[here]
End[]
$Context
# Kernel names that have no reference counterpart live in Mathilda` (v0.565):
Context[NDArray]
