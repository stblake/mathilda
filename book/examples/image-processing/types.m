# Type inference reads the values; malformed data does not make an image.
ImageType /@ {Image[{{0, 1}, {1, 0}}], Image[{{0, 7}, {200, 255}}], Image[{{0., 1.}, {1., 0.}}]}
ImageData[Image[{{0, 1}, {1, 0}}]]
ImageQ[Image[{{1, 2}, {3}}]]
ImageQ[Image["hello"]]
ImageData[Image[{{-0.5, 1.5}}]]
