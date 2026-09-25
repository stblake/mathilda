# 9.4 The stream as a cursor: Read one object at a time, and move the point.
s = OpenWrite["/tmp/mathilda_book_seq.dat"];
WriteString[s, "1 2 3\n4.5 6\n"];
Close[s];
rd = OpenRead["/tmp/mathilda_book_seq.dat"]
Read[rd, Number]
Read[rd, Number]
Read[rd, {Number, Number}]
StreamPosition[rd]
Streams[]
SetStreamPosition[rd, 0]
Read[rd, Number]
Close[rd]
# A nested type structure reads a whole matrix in one call.
m = OpenWrite["/tmp/mathilda_book_matrix.dat"];
WriteString[m, "1 2\n3 4\n"];
Close[m];
mat = OpenRead["/tmp/mathilda_book_matrix.dat"];
Read[mat, {{Number, Number}, {Number, Number}}]
Close[mat];
# Reading past the end returns EndOfFile.
e = OpenRead["/tmp/mathilda_book_matrix.dat"];
Read[e, Number]; Read[e, Number]; Read[e, Number]; Read[e, Number];
Read[e, Number]
Close[e];
