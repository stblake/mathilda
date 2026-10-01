# Storage: pixels live in a machine buffer, and ImageData hands back an ordinary List.
g = GaussianFilter[Image[{{0.1, 0.9}, {0.4, 0.6}}], 1];
Head[First[g]]
Head[ImageData[g]]
Head[First[Image[{{1, 2}, {3, 4}}]]]
CompileDiagnostics[{{v, _Real, 2}}, ImageData[GaussianFilter[Image[v], 1]]]
