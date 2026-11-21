# 9.7 Import and Export: raster images, graphics, and plain text.
# An image round-trips through a PNG file (see the Graphics chapter for depth).
img = Image[Table[{N[i/8], N[j/8], 0.5}, {i, 1, 8}, {j, 1, 8}], "Real"];
Export["/tmp/mathilda_book_swatch.png", img]
ImageDimensions[Import["/tmp/mathilda_book_swatch.png"]]
# Plain text: write a small file, then Import it three ways.
d = OpenWrite["/tmp/mathilda_book_notes.txt"];
WriteString[d, "alpha 1\nbeta 2\ngamma 3\n"];
Close[d];
Import["/tmp/mathilda_book_notes.txt"]
Import["/tmp/mathilda_book_notes.txt", "Lines"]
StringLength[Import["/tmp/mathilda_book_notes.txt"]]
StringLength[Import["/tmp/mathilda_book_notes.txt", "String"]]
ToExpression[Last /@ StringSplit /@ Import["/tmp/mathilda_book_notes.txt", "Lines"]]
# A structured data format is not one Import claims: the call stays unevaluated.
Head[Import["/tmp/mathilda_book_notes.txt", "CSV"]]
# Nor does Export write one -- there is no CSV/JSON writer yet.
Head[Export["/tmp/mathilda_book_out.csv", {{1, 2}, {3, 4}}]]
# A missing file is a failure, not a decline.
Import["/tmp/mathilda_book_absent_98765.txt"]
