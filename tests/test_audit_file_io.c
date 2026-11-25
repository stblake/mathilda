/* Builtin-audit regressions: file I/O (review/issues/<Head>.md; tests/test_audit.h).
 * Every file lives under $TemporaryDirectory and is deleted by the check that made it. */
#include "test_audit.h"

/* Export#1 (critical, crash): Export of a Graphics to PNG/JPEG segfaulted inside raylib's
 * InitWindow whenever the window or its GL context was unusable -- InitPlatform with a NULL
 * window under a sandbox that denies IOKit, or rlglInit -> rlLoadTexture jumping to a never-loaded
 * GL entry point (address 0) in a tool process that has a GUI session.  Nothing in-process predicts
 * either, so the first raster export now does the whole render once in a spawned child
 * (gl_probe.c) and declines with Export::nogl and $Failed if the child does not succeed.  In a
 * usable session the file is written.  Either outcome passes; a crash kills the suite.  (Checked by
 * hand under `sandbox-exec (deny iokit-open)`: SIGSEGV before, Export::nogl + $Failed after.) */
static void test_export_graphics_raster_no_crash(void) {
    audit_true("Export#1",
               "Module[{p = FileNameJoin[{$TemporaryDirectory, \"mth_audit_export1.png\"}], r, c},"
               " c = Quiet[Check[r = Export[p, Graphics[Disk[]]]; \"ok\", \"nogl\", Export::nogl]];"
               " If[FileExistsQ[p], DeleteFile[p]];"
               " (r === p && c === \"ok\") || (r === $Failed && c === \"nogl\")]");
    audit_true("Export#1",
               "Module[{p = FileNameJoin[{$TemporaryDirectory, \"mth_audit_export1.jpg\"}], r, c},"
               " c = Quiet[Check[r = Export[p, Graphics[Disk[]]]; \"ok\", \"nogl\", Export::nogl]];"
               " If[FileExistsQ[p], DeleteFile[p]];"
               " (r === p && c === \"ok\") || (r === $Failed && c === \"nogl\")]");
    /* sanity: a raster Image export (no renderer involved) still writes and round-trips */
    audit_true("Export#1",
               "Module[{p = FileNameJoin[{$TemporaryDirectory, \"mth_audit_export1b.png\"}], r, d},"
               " r = Export[p, Image[{{0., 1.}, {1., 0.}}]]; d = ImageData[Import[p]];"
               " DeleteFile[p]; r === p && d == {{0., 1.}, {1., 0.}}]");
}

/* Read#1 (high): Read[s, Number] parsed each token as Mathilda code and EVALUATED it --
 * Print["PWN"] in the file ran, 1/2 became a Rational and 2^3 became 8.  A Number is a
 * literal (sign, digits, point, e/E or *^ exponent); anything else is Read::readn and
 * $Failed, and is never evaluated. */
static void test_read_number_literal_only(void) {
    audit_run("audit$f = FileNameJoin[{$TemporaryDirectory, \"mth_audit_read1.txt\"}];"
              " audit$o = OpenWrite[audit$f];"
              " WriteString[audit$o, \"1/2 audit$side[] 2^3 -4 2.5 1.5e3 -2.*^-2 .5 7\\n\"];"
              " Close[audit$o]; audit$hit = 0; audit$side[] := (audit$hit = 1; 99)");
    audit_chk("Read#1",
              "Quiet[Read[audit$f, {Number, Number, Number, Number, Number, Number, Number, Number, Number}]]",
              "{$Failed, $Failed, $Failed, -4, 2.5, 1500.0, -0.02, 0.5, 7}");
    audit_chk("Read#1", "audit$hit", "0");                 /* the token was not run */
    audit_run("Close[audit$f]");
    audit_chk("Read#1", "Quiet[Check[Read[audit$f, Number], fired, Read::readn]]", "fired");
    audit_run("Close[audit$f]");
    audit_chk("Read#1", "Quiet[Read[audit$f, Real]]", "$Failed");
    audit_run("Close[audit$f]; DeleteFile[audit$f]");
    /* sanity: an ordinary numeric file reads as before */
    audit_chk("Read#1",
              "Module[{f = FileNameJoin[{$TemporaryDirectory, \"mth_audit_read1b.txt\"}], o, r},"
              " o = OpenWrite[f]; WriteString[o, \"10 -3.25\\n+6 1E2\\n\"]; Close[o];"
              " r = ReadList[f, Number]; DeleteFile[f]; r]",
              "{10, -3.25, 6, 100.0}");
}

/* DirectoryName#2 (low): the OperatingSystem option. */
static void test_s14_directoryname_os(void) {
    audit_chk("DirectoryName#2", "DirectoryName[\"a\\\\b\", OperatingSystem -> \"Windows\"] === \"a\\\\\"",
              "True");
    audit_chk("DirectoryName#2", "DirectoryName[\"a\\\\b\", OperatingSystem -> \"Unix\"] === \"\"",
              "True");
    audit_chk("DirectoryName#2", "DirectoryName[\"a/b/c\", 2, OperatingSystem -> \"Unix\"] === \"a/\"",
              "True");
}

/* Wave 3 S12/S13: declared argument counts and option checks
 * (src/arity.h, src/arity_table.c).  The evaluator checks a call
 * against its head's declared arity before the C code runs and
 * prints argx/argr/argrx/argt(u)/argb(u)/argm(u), nonopt or optx;
 * the call stays unevaluated. */
static void test_wave3_arity_options(void) {
    audit_chk("Environment#2", "Quiet[Check[Environment[Sequence[]], fired, Environment::argx]]", "fired");
    audit_chk("Environment#2", "Quiet[Environment[Sequence[]]]", "Environment[]");
    audit_chk("OpenRead#2", "Quiet[Check[OpenRead[], fired, OpenRead::argx]]", "fired");
    audit_chk("OpenRead#2", "Quiet[OpenRead[]]", "OpenRead[]");
    audit_chk("FileNameSplit#3", "Quiet[Check[FileNameSplit[\"a\", \"b\"], fired, FileNameSplit::nonopt]]", "fired");
    audit_chk("FileNameSplit#3", "Quiet[FileNameSplit[\"a\", \"b\"]]", "FileNameSplit[\"a\", \"b\"]");
    audit_chk("ReadList#4", "Quiet[Check[ReadList[StringToStream[\"a b\"], Word, Bar -> 1], fired, ReadList::optx]]", "fired");
}

/* Wave 4 G8: file names, streams and the file reader.
 *  - File["name"] is accepted wherever a file name is (AbsoluteFileName, DeleteFile,
 *    SetDirectory, Get, Put, PutAppend);
 *  - a stream object resolves only when its id, direction AND name match an open stream,
 *    and an id beyond the int range names none (it was truncated onto stream 1);
 *  - Put / PutAppend's noopen goes through the message funnel (Quiet, Check);
 *  - a file is read with line-ended statements, as the reference does: "5" and "-3" on
 *    two lines are two expressions, not 5 - 3;
 *  - a syntax error in Read / ReadList is a Syntax::sntx message and $Failed, not
 *    EndOfFile (or a silently dropped object);
 *  - Get and LoadModule refuse a directory; LoadModule's at-most-once memo is keyed by
 *    the normalised path;
 *  - Character reads a UTF-8 character, not a byte; export format names are
 *    case-insensitive and CSV / TSV / .m are written. */
static void test_wave4_g8_files_streams(void) {
    audit_run("w4$d = $TemporaryDirectory;"
              " w4$f[n_] := FileNameJoin[{w4$d, \"mth_audit_w4g8_\" <> n}]");
    /* File[...] arguments */
    audit_run("WriteString[w4$f[\"d.txt\"], \"x\"]; Close[w4$f[\"d.txt\"]]");
    audit_true("AbsoluteFileName#1", "StringQ[Quiet[AbsoluteFileName[File[w4$f[\"d.txt\"]]]]]");
    audit_true("File#1", "Quiet[Check[AbsoluteFileName[File[w4$f[\"d.txt\"]]], fired]] =!= fired");
    audit_chk("File#1", "{DeleteFile[File[w4$f[\"d.txt\"]]], FileExistsQ[w4$f[\"d.txt\"]]}", "{Null, False}");
    audit_run("Put[1, w4$f[\"a2.txt\"]]");
    audit_chk("DeleteFile#1", "{DeleteFile[File[w4$f[\"a2.txt\"]]], FileExistsQ[w4$f[\"a2.txt\"]]}", "{Null, False}");
    audit_true("SetDirectory#1", "Module[{r = SetDirectory[File[w4$d]]}, ResetDirectory[];"
               " StringQ[r] && r === AbsoluteFileName[w4$d]]");
    audit_run("Put[{1, 2, 3}, w4$f[\"f.m\"]]");
    audit_chk("File#2", "Get[File[w4$f[\"f.m\"]]]", "{1, 2, 3}");
    audit_chk("File#2", "{Put[7, File[w4$f[\"p.m\"]]], FileExistsQ[w4$f[\"p.m\"]], Get[w4$f[\"p.m\"]]}", "{Null, True, 7}");
    audit_chk("File#2", "{PutAppend[9, File[w4$f[\"f.m\"]]], ReadList[w4$f[\"f.m\"]]}", "{Null, {{1, 2, 3}, 9}}");
    /* Put / PutAppend noopen through the funnel */
    audit_chk("Put#1", "Quiet[Put[1, \"/nonexistent/dir/x.m\"]]", "$Failed");
    audit_chk("Put#1", "Quiet[Check[Put[1, \"/nonexistent/dir/x.m\"], failed]]", "failed");
    audit_chk("PutAppend#1", "Quiet[PutAppend[x, \"/nonexistent/dir/f\"]]", "$Failed");
    audit_chk("PutAppend#1", "Quiet[Check[PutAppend[x, \"/nonexistent/dir/f\"], caught, PutAppend::noopen]]", "caught");
    /* line-ended statements */
    audit_run("PutAppend[5, -3, w4$f[\"pb1\"]]");
    audit_chk("PutAppend#2", "Get[w4$f[\"pb1\"]]", "-3");
    audit_chk("PutAppend#2", "ReadList[w4$f[\"pb1\"]]", "{5, -3}");
    /* syntax errors in Read / ReadList */
    audit_run("w4$o = OpenWrite[w4$f[\"r2.txt\"]]; WriteString[w4$o, \"a + )\\n5\\n\"]; Close[w4$o];"
              " w4$s = OpenRead[w4$f[\"r2.txt\"]]");
    audit_chk("Read#2", "Quiet[Check[Read[w4$s], fired, Syntax::sntx]]", "fired");
    audit_chk("Read#2", "Read[w4$s, Number]", "5");
    audit_chk("Read#2", "Read[w4$s]", "EndOfFile");
    audit_run("Close[w4$s]");
    audit_chk("ReadList#2", "Quiet[ReadList[StringToStream[\"x+1\\n1+\"]]]", "{1 + x, $Failed}");
    audit_chk("ReadList#2", "Quiet[Check[ReadList[StringToStream[\"{1,2\"]], fired, Syntax::sntx]]", "fired");
    /* Get / LoadModule on a directory */
    audit_chk("Get#1", "Quiet[Get[w4$d]]", "$Failed");
    audit_chk("Get#1", "Quiet[Check[Get[w4$d], fired, Get::noopen]]", "fired");
    audit_chk("LoadModule#1", "Quiet[LoadModule[\"simp\"]]", "False");
    audit_chk("LoadModule#2", "{LoadModule[\"simp/FullSimplify.m\"], Check[LoadModule[\"simp//FullSimplify.m\"], fired]}",
              "{True, True}");
    /* stream objects resolve by id, direction and name */
    audit_run("w4$s4 = OpenWrite[w4$f[\"c4.txt\"]]");
    audit_chk("Close#1", "Quiet[Close[OutputStream[\"bogus\", w4$s4[[2]]]]]", "$Failed");
    audit_chk("Close#1", "Length[Streams[w4$f[\"c4.txt\"]]]", "1");
    audit_chk("OutputStream#1", "Quiet[Close[InputStream[w4$s4[[1]], w4$s4[[2]]]]]", "$Failed");
    audit_chk("OutputStream#1", "Close[w4$s4] === w4$f[\"c4.txt\"]", "True");
    audit_run("w4$ss = StringToStream[\"1 2 3\"]");
    audit_chk("InputStream#1", "Quiet[Close[InputStream[\"String\", 4294967297 + w4$ss[[2]] - 1]]]", "$Failed");
    audit_chk("InputStream#1", "Read[w4$ss, Number]", "1");
    audit_run("Close[w4$ss]; w4$ss = StringToStream[\"abcdef\"]; Read[w4$ss, Character]");
    audit_chk("StreamPosition#1", "Quiet[StreamPosition[OutputStream[\"bogus\", w4$ss[[2]]]]]", "$Failed");
    audit_chk("StreamPosition#1", "StreamPosition[w4$ss]", "1");
    audit_run("Close[w4$ss]");
    /* characters, not bytes */
    audit_chk("StringToStream#1", "Module[{s = StringToStream[\"\\[Alpha]\\[Beta]\\[Gamma]\"], r},"
              " r = ReadList[s, Character]; Close[s]; r === {\"\\[Alpha]\", \"\\[Beta]\", \"\\[Gamma]\"}]", "True");
    /* export formats */
    audit_chk("Export#2", "Export[w4$f[\"b\"], Image[{{0., 0.5}, {1., 0.25}}], \"png\"] === w4$f[\"b\"]", "True");
    audit_chk("Export#2", "ImageDimensions[Import[w4$f[\"b\"], \"PNG\"]]", "{2, 2}");
    audit_chk("Export#3", "Export[w4$f[\"t.csv\"], {{1, 2}, {3, \"a,b\"}}] === w4$f[\"t.csv\"]", "True");
    audit_chk("Export#3", "ReadList[w4$f[\"t.csv\"], String]", "{\"1,2\", \"3,\\\"a,b\\\"\"}");
    audit_chk("Export#3", "Export[w4$f[\"t.m\"], x^2 + 1] === w4$f[\"t.m\"]", "True");
    audit_chk("Export#3", "Get[w4$f[\"t.m\"]]", "1 + x^2");
    audit_chk("Export#3", "Export[w4$f[\"t.txt\"], \"hello\"] === w4$f[\"t.txt\"]", "True");
    audit_run("Quiet[DeleteFile[Map[w4$f, {\"f.m\", \"p.m\", \"pb1\", \"r2.txt\", \"c4.txt\", \"b\","
              " \"t.csv\", \"t.m\", \"t.txt\"}]]]");
}

/* FileNameSplit#1: the OperatingSystem default set by SetOptions is honoured. */
static void test_wave4_g8_filenamesplit_options(void) {
    audit_run("SetOptions[FileNameSplit, OperatingSystem -> \"Windows\"]");
    audit_chk("FileNameSplit#1", "FileNameSplit[\"a\\\\b\"]", "{\"a\", \"b\"}");
    audit_chk("FileNameSplit#1", "FileNameSplit[\"a\\\\b\", OperatingSystem -> \"Unix\"]", "{\"a\\\\b\"}");
    audit_run("SetOptions[FileNameSplit, OperatingSystem -> $OperatingSystem]");
    audit_chk("FileNameSplit#1", "FileNameSplit[\"a/b\"]", "{\"a\", \"b\"}");
}

/* Wave 5 G8: file-io findings -- File[] names, argument messages, bignum counts,
 * the standard streams, stream options and positions, Import/Export messages. */
static void test_wave5_g8_files(void) {
    audit_run("w5$d = $TemporaryDirectory;"
              " w5$f[n_] := FileNameJoin[{w5$d, \"mth_audit_w5g8_\" <> n}]");
    /* File["name"] in the name heads */
    audit_chk("FileBaseName#1", "FileBaseName[File[\"a/b.txt\"]]", "\"b\"");
    audit_chk("FileExtension#1", "FileExtension[File[\"a/b.txt\"]]", "\"txt\"");
    audit_chk("FileNameTake#1", "FileNameTake[File[\"/a/b/c\"]]", "\"c\"");
    audit_chk("DirectoryName#3", "DirectoryName[File[\"/a/b/c\"]]", "\"/a/b/\"");
    audit_chk("ParentDirectory#1", "ParentDirectory[File[\"/usr/bin\"]] === ParentDirectory[\"/usr/bin\"]", "True");
    audit_chk("DirectoryName#1", "Quiet[Check[DirectoryName[], fired, DirectoryName::argt]]", "fired");
    audit_chk("DirectoryName#1", "Quiet[Check[DirectoryName[\"a\", \"b\", \"c\"], fired, DirectoryName::argt]]", "fired");
    audit_chk("DirectoryName#1", "Quiet[DirectoryName[\"a\", \"b\", \"c\"]]", "DirectoryName[\"a\", \"b\", \"c\"]");
    audit_run("Put[1, w5$f[\"fs.txt\"]]");
    audit_chk("FileSize#1", "FileSize[File[w5$f[\"fs.txt\"]]] === FileSize[w5$f[\"fs.txt\"]] === 2", "True");
    audit_run("Put[{1, 2, 3}, w5$f[\"f.m\"]]");
    audit_chk("File#3", "Import[File[w5$f[\"f.m\"]], \"Text\"]", "\"{1, 2, 3}\"");
    audit_chk("Import#3", "Export[w5$f[\"a.png\"], Image[{{0., 0.5}, {1., 0.25}}, \"Real\"]] === w5$f[\"a.png\"]", "True");
    audit_chk("Import#3", "ImageDimensions[Import[File[w5$f[\"a.png\"]]]]", "{2, 2}");
    audit_chk("Import#2", "Quiet[Import[w5$f[\"nope.png\"]]]", "$Failed");
    audit_chk("Import#2", "Quiet[Check[Import[w5$f[\"nope.png\"]], fired, Import::nffil]]", "fired");
    audit_chk("Export#4", "Quiet[Export[\"/nonexistent/dir/a.png\", Image[{{0.5}}]]]", "$Failed");
    audit_chk("Export#4", "Quiet[Check[Export[\"/nonexistent/dir/a.png\", Image[{{0.5}}]], fired, Export::noopen]]", "fired");
    audit_chk("Export#5", "Export[w5$f[\"n.png\"], Image[{{0., 0.5}, {1., 0.25}}], ImageSize -> 10] === w5$f[\"n.png\"]", "True");
    audit_chk("Export#5", "ImageDimensions[Import[w5$f[\"n.png\"]]]", "{10, 10}");
    /* FilePrint: the funnel, bignum counts, File[] */
    audit_chk("FilePrint#1", "Quiet[FilePrint[w5$f[\"nofile.txt\"]]]", "$Failed");
    audit_chk("FilePrint#1", "Quiet[Check[FilePrint[w5$f[\"nofile.txt\"]], caught]]", "caught");
    audit_chk("FilePrint#1", "Quiet[Check[FilePrint[w5$f[\"nofile.txt\"]], caught, FilePrint::noopen]]", "caught");
    audit_chk("FilePrint#2", "FilePrint[w5$f[\"fs.txt\"], 10^30]", "Null");
    audit_chk("FilePrint#3", "FilePrint[File[w5$f[\"fs.txt\"]]]", "Null");
    /* names */
    audit_chk("Environment#1", "Environment[\"HOME=\"]", "$Failed");
    audit_chk("Environment#1", "Environment[\"\"]", "$Failed");
    audit_chk("FileNameJoin#1", "FileNameJoin[{\"a\", \"b\"}, OperatingSystem -> \"Windows\", OperatingSystem -> \"Unix\"]", "\"a\\\\b\"");
    audit_chk("FileNameSplit#2", "FileNameSplit[\"a/b\", {OperatingSystem -> \"Windows\"}]", "{\"a\", \"b\"}");
    audit_chk("FileNameSplit#2", "FileNameSplit[\"a\\\\b\", {OperatingSystem -> \"Windows\"}]", "{\"a\", \"b\"}");
    audit_chk("ToFileName#1", "ToFileName[\"a\", \"/c\"]", "\"a/c\"");
    audit_chk("ToFileName#1", "ToFileName[{\"/a/\", \"/b/\"}, \"c\"]", "\"/a/b/c\"");
    audit_chk("ToFileName#1", "ToFileName[{\"a/\", \"b\"}]", "\"a/b/\"");
    /* FileNames */
    audit_chk("FileNames#1", "FileNames[\"mth_audit_w5g8_fs*\", w5$d, 2^63] === FileNames[\"mth_audit_w5g8_fs*\", w5$d, 2^62]", "True");
    audit_chk("FileNames#1", "Length[FileNames[\"mth_audit_w5g8_fs*\", w5$d, 2^63]]", "1");
    audit_chk("FileNames#2", "Block[{$MessageList = {}}, Quiet[FileNames[RegularExpression[\"[\"], w5$d]]; 0]; "
              "Block[{$MessageList = {}}, FileNames[RegularExpression[\"[\"], w5$d]; Length[$MessageList]]", "1");
    /* streams */
    audit_chk("OpenAppend#2", "Quiet[OpenAppend[5]]", "$Failed");
    audit_chk("OpenAppend#2", "Quiet[Check[OpenAppend[5], fired, OpenAppend::fstr]]", "fired");
    audit_chk("OpenRead#3", "Quiet[Check[OpenRead[123], fired, OpenRead::fstr]]", "fired");
    audit_chk("OpenRead#1", "Quiet[OpenRead[w5$d]]", "$Failed");
    audit_chk("OpenRead#1", "Quiet[Check[OpenRead[w5$d], fired, OpenRead::noopen]]", "fired");
    audit_chk("OpenWrite#1", "Module[{s = OpenWrite[], r}, r = {Head[s], StringStartsQ[s[[1]], w5$d]}; DeleteFile[Close[s]]; r]",
              "{OutputStream, True}");
    audit_chk("OutputStream#2", "Take[Streams[], 2]", "{OutputStream[\"stdout\", 1], OutputStream[\"stderr\", 2]}");
    audit_chk("OutputStream#2", "Close[OutputStream[\"stdout\", 1]]; MemberQ[Streams[], OutputStream[\"stdout\", 1]]", "True");
    audit_chk("Write#1", "Quiet[Check[Write[w5x, 2], fired, General::stream]]", "fired");
    audit_chk("Write#1", "Quiet[Write[w5x, 2]]", "Write[w5x, 2]");
    audit_run("w5$o = OpenWrite[w5$f[\"ow3.txt\"]]; WriteString[w5$o, \"hello\"]");
    audit_chk("SetStreamPosition#1", "SetStreamPosition[w5$o, 100]", "5");
    audit_run("Close[w5$o]");
    audit_chk("SetStreamPosition#1", "FileSize[w5$f[\"ow3.txt\"]]", "5");
    audit_run("w5$o = OpenWrite[w5$f[\"w3.txt\"], FormatType -> OutputForm]; Write[w5$o, x^2, \"q\"]; Close[w5$o]");
    audit_chk("Write#2", "ReadString[w5$f[\"w3.txt\"]]", "\"x^2q\\n\"");
    audit_run("w5$o = OpenAppend[w5$f[\"oa12.txt\"], FormatType -> FortranForm, PageWidth -> 20];"
              " Write[w5$o, Expand[(a + b)^6]]; Close[w5$o]");
    audit_chk("OpenAppend#1", "Max[StringLength /@ StringSplit[ReadString[w5$f[\"oa12.txt\"]], \"\\n\"]] <= 20", "True");
    audit_run("w5$o = OpenWrite[w5$f[\"p22.m\"]]");
    audit_chk("Put#2", "Put[1 + x, w5$o]", "Null");
    audit_run("Close[w5$o]");
    audit_chk("Put#2", "ReadString[w5$f[\"p22.m\"]]", "\"1 + x\\n\"");
    audit_chk("Put#2", "Quiet[Put[1, w5$o]]", "$Failed");
    /* Read / ReadList / Skip */
    audit_chk("Read#3", "Quiet[Check[Read[], fired, Read::argt]]", "fired");
    audit_chk("Read#3", "Quiet[Check[Read[123], fired, General::stream]]", "fired");
    audit_chk("Read#3", "Quiet[Read[123]]", "Read[123]");
    audit_chk("Read#3", "Quiet[Check[Read[w5$f[\"w3.txt\"], Foo], fired, Read::readf]]", "fired");
    audit_run("Close[w5$f[\"w3.txt\"]]");
    audit_chk("ReadList#3", "Module[{s = StringToStream[\"a b\"], r}, r = ReadList[s, Word, 2^70]; Close[s]; r]", "{\"a\", \"b\"}");
    audit_chk("Skip#1", "Module[{s = StringToStream[\"1 2 3\"], r}, r = Skip[s, Number, 2^70]; Close[s]; r]", "EndOfFile");
    audit_chk("ReadString#1", "StringContainsQ[Information[ReadString], \"CRLF line ends stay\"]", "True");
    audit_run("Quiet[DeleteFile[Map[w5$f, {\"fs.txt\", \"f.m\", \"a.png\", \"n.png\", \"ow3.txt\", \"w3.txt\","
              " \"oa12.txt\", \"p22.m\"}]]]");
}

int main(void) {
    audit_init();
    TEST(test_export_graphics_raster_no_crash);
    TEST(test_read_number_literal_only);
    TEST(test_s14_directoryname_os);
    TEST(test_wave3_arity_options);
    TEST(test_wave4_g8_files_streams);
    TEST(test_wave4_g8_filenamesplit_options);
    TEST(test_wave5_g8_files);
    return audit_finish("audit_file_io");
}
