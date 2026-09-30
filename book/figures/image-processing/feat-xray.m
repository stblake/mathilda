# format: png
# Featured: outline a bone in a noisy radiograph. A hollow ellipse (cortex around marrow) in noise.
SeedRandom[6]; bone = Image[Table[0.2 + 0.5 Boole[((x - 384)/280.)^2 + ((y - 192)/56.)^2 < 1] - 0.25 Boole[((x - 384)/240.)^2 + ((y - 192)/24.)^2 < 1] + RandomReal[{-0.1, 0.1}], {y, 384}, {x, 768}]];
fig = ImageAssemble[{{bone}, {EdgeDetect[bone, 6]}}]
