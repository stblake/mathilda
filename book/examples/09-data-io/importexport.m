# 9.7 Import and Export: raster images and graphics only, not general data.
# An image round-trips through a PNG file (see the Graphics chapter for depth).
img = Image[Table[{N[i/8], N[j/8], 0.5}, {i, 1, 8}, {j, 1, 8}], "Real"];
Export["/tmp/mathilda_book_swatch.png", img]
ImageDimensions[Import["/tmp/mathilda_book_swatch.png"]]
# A plain text/data file is not a format Import claims, so the call is simply
# left unevaluated -- Import is not a general data reader.
d = OpenWrite["/tmp/mathilda_book_notimg.txt"];
WriteString[d, "1 2 3\n"];
Close[d];
Head[Import["/tmp/mathilda_book_notimg.txt"]]
# Likewise Export declines a bare list -- there is no CSV/JSON writer yet.
Head[Export["/tmp/mathilda_book_out.csv", {{1, 2}, {3, 4}}]]
