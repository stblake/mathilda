# 9.7 Import and Export: raster images, graphics, and plain text.
# An image round-trips through a PNG file (see the Graphics chapter for depth).
img = Image[Table[{N[i/8], N[j/8], 0.5}, {i, 1, 8}, {j, 1, 8}], "Real"];
Export["/tmp/mathilda_book_swatch.png", img]
ImageDimensions[Import["/tmp/mathilda_book_swatch.png"]]
# Plain text: Export a small file, then Import it three ways.
Export["/tmp/mathilda_book_notes.txt", "alpha 1\nbeta 2\ngamma 3"]
Import["/tmp/mathilda_book_notes.txt"]
Import["/tmp/mathilda_book_notes.txt", "Lines"]
StringLength[Import["/tmp/mathilda_book_notes.txt"]]
StringLength[Import["/tmp/mathilda_book_notes.txt", "String"]]
ToExpression[Last /@ StringSplit /@ Import["/tmp/mathilda_book_notes.txt", "Lines"]]
# Export writes "Text", "Lines" and "String"; Text round-trips exactly.
s = "first line\nsecond line\n";
Import[Export["/tmp/mathilda_book_round.txt", s], "Text"] === s
Export["/tmp/mathilda_book_list.dat", {"a", "b", "c"}, "Lines"]
Import["/tmp/mathilda_book_list.dat", "String"]
Export["/tmp/mathilda_book_expr.txt", 1 + x^2];
Import["/tmp/mathilda_book_expr.txt"]
Export["/tmp/mathilda_book_no_such_dir/out.txt", "x"]
# A structured data format is not one Import claims: the call stays unevaluated.
Head[Import["/tmp/mathilda_book_notes.txt", "CSV"]]
# Nor does Export write one -- there is no CSV/JSON writer yet.
Head[Export["/tmp/mathilda_book_out.csv", {{1, 2}, {3, 4}}]]
# A missing file is a failure, not a decline.
Import["/tmp/mathilda_book_absent_98765.txt"]
