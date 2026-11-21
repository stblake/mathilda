# 4.4.2 More products: KroneckerProduct, TensorProduct, Projection, Orthogonalize
KroneckerProduct[{{1, 2}, {3, 4}}, {{0, 1}, {1, 0}}]
KroneckerProduct[{1, 2}, {3, 4}]
TensorProduct[{a, b}, {x, y}]
Dimensions[TensorProduct[{1, 2}, {3, 4, 5}, {6, 7}]]
Projection[{1, 1}, {3, 4}]
{1, 1} - Projection[{1, 1}, {3, 4}]
Orthogonalize[{{1, 1, 0}, {1, 0, 1}, {0, 1, 1}}]
Orthogonalize[{{1, 2}, {2, 4}}]
