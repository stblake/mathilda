# 9.2 Type-directed reading with ReadList.
# A small numeric table, including scientific-notation tokens.
w = OpenWrite["/tmp/mathilda_book_types.dat"];
WriteString[w, "1 2 3\n4.5 6\n2e5 1.5e-3\n"];
Close[w];
ReadList["/tmp/mathilda_book_types.dat", Number]
ReadList["/tmp/mathilda_book_types.dat", Real]
ReadList["/tmp/mathilda_book_types.dat", Word]
ReadList["/tmp/mathilda_book_types.dat", String]
ReadList["/tmp/mathilda_book_types.dat", Byte, 5]
ReadList["/tmp/mathilda_book_types.dat", Character, 5]
ReadList["/tmp/mathilda_book_types.dat", Number, 3]
# A type LIST reads one object of each type per pass, grouping each pass.
p = OpenWrite["/tmp/mathilda_book_people.dat"];
WriteString[p, "alice 30\nbob 25\ncarol 41\n"];
Close[p];
ReadList["/tmp/mathilda_book_people.dat", {Word, Number}]
# When end of file arrives mid-pass, the unread slots become EndOfFile.
o = OpenWrite["/tmp/mathilda_book_odd.dat"];
WriteString[o, "1 2 3\n"];
Close[o];
ReadList["/tmp/mathilda_book_odd.dat", {Number, Number}]
