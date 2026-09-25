# Point-in-polygon test. True inside OR on the boundary; exact input decides exactly.
sq = Polygon[{{0, 0}, {2, 0}, {2, 2}, {0, 2}}];
RegionMember[sq, {1, 1}]
RegionMember[sq, {2, 1}]
RegionMember[sq, {3, 1}]
RegionMember[Polygon[{{0, 0}, {2, 0}, {0, 2}}], {1, 1}/2]
RegionMember[Polygon[{{0, 0}, {4, 0}, {4, 4}, {2, 1}, {0, 4}}], {2, 3}]
