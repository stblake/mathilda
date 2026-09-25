# Convex hull of a 2D point set (Andrew's monotone chain).
# Duplicate, interior, and collinear-middle points are discarded.
pts = {{0, 0}, {2, 0}, {1, 0}, {2, 2}, {0, 2}, {1, 1}};
ConvexHullRegion[pts]
Area[ConvexHullRegion[pts]]
ConvexHullRegion[{{0, 0}, {1, 0}, {1, 1}, {0, 1}, {1/2, 1/2}}]
ConvexHullRegion[{{0, 0}, {1, 1}, {2, 2}, {3, 3}}]
ConvexHullRegion[{{1, 2}}]
