# Computational geometry: Area of a polygon via the shoelace formula.
# Coordinates given exactly are kept exact -- the area is an exact number.
square = Polygon[{{0, 0}, {4, 0}, {4, 3}, {0, 3}}];
Area[square]
Area[Polygon[{{0, 0}, {1, 0}, {1/2, 1/2}}]]
Area[Polygon[{{0, 0}, {4, 0}, {4, 4}, {2, 1}, {0, 4}}]]
