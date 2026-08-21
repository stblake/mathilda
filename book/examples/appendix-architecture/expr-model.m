# Everything is an expression: FullForm reveals the underlying tree.
FullForm[a + b c]
FullForm[{a, b, c}]
FullForm[1/2]
FullForm[a/b]
FullForm[x_]
FullForm[a -> b]
# Even atoms have a head.
Head[2]
Head[2.5]
Head[{1, 2}]
Head["hi"]
