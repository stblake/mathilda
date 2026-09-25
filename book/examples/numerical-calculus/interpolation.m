# 4.5 Interpolation from data
h = Interpolation[{{1, 1}, {2, 8}, {3, 27}, {4, 64}}]
h[2.5]
f = ListInterpolation[{1, 4, 9, 16, 25}]
f[2.5]
ListInterpolation[{1, 4, 9, 16, 25}, {{0, 2}}][1]
