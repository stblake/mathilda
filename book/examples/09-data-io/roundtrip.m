# Chapter 9 -- Data I/O
# 9.1 The round trip: open a file, write to it, close it, read it back.
# Every example in this chapter is self-contained -- it writes to a scratch
# path under /tmp and then reads that same file, so the whole book rebuilds
# from nothing.
s = OpenWrite["/tmp/mathilda_book_grades.dat"]
WriteString[s, "88 91 79\n95 82 100\n"];
Close[s]
FileExistsQ["/tmp/mathilda_book_grades.dat"]
ReadList["/tmp/mathilda_book_grades.dat", Number]
