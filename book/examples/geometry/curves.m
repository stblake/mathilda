# Arc length, the measure of disks and annuli, and rotations
ArcLength[Line[{{0, 0}, {3, 4}, {3, 0}}]]
ArcLength[{Cos[t], Sin[t], t}, {t, 0, 2 Pi}]
ArcLength[{t, t^2}, {t, 0, 1}]
ArcLength[{t, Sin[t]}, {t, 0, Pi}]
RegionMeasure[Disk[{0, 0}, r]]
Area[Annulus[{0, 0}, {r1, r2}]]
Perimeter[Annulus[{0, 0}, {1, 2}]]
{RegionMember[Annulus[{0, 0}, {1, 2}], {3/2, 0}], RegionMember[Annulus[{0, 0}, {1, 2}], {1/2, 0}]}
RegionCentroid[Annulus[{1, 2}, {1, 2}]]
RotationMatrix[t]
RotationMatrix[Pi/3] . {1, 0}
RotationMatrix[t, {0, 0, 1}]
m = RotationMatrix[{{1, 2, 3}, {3, 2, 1}}]
{m . {1, 2, 3}, Det[m]}
