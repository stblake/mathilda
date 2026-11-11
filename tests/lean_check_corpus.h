/*
 * lean_check_corpus.h -- shared inputs of the LeanVerify source-check tests
 * (test_lean_certify.c: Lean`Internal`CheckSource only; test_lean_verify.c:
 * the same sources judged by a real Lean kernel).
 *
 *   LEAN_CORPUS     legitimate LeanVerify sources (core Lean 4.35, no
 *                   imports): forall / exists / anonymous constructors /
 *                   fun / mapsto / have / show / calc / match / decide /
 *                   omega / grind / simp / rw / cases / namespaces / sections
 *                   / variables / universes / open / set_option.  Every one
 *                   is true, accepted by the checker and proved by Lean.
 *   LEAN_REFUSED    constructs outside LeanVerify's surface or word rules;
 *                   inserting any of them anywhere (at a character boundary)
 *                   into a corpus source must make the checker refuse it.
 *   REVIEW4_X1..X4  the fourth review's reproducers (each hides a column-0
 *                   #eval inside what an earlier checker took for a string,
 *                   character literal or comment).  The lexical structure is
 *                   the reviewer's, byte for byte; only the payloads' file
 *                   paths are neutralised (/dev/null, /usr/bin/true), so a
 *                   regression could not touch anything.  Never give these to
 *                   a plain `lean`: they must be refused before Lean runs.
 */
#ifndef MATHILDA_LEAN_CHECK_CORPUS_H
#define MATHILDA_LEAN_CHECK_CORPUS_H

#include <stdlib.h>
#include <string.h>

static const char* const LEAN_CORPUS[] = {
    "theorem c1 : 1 + 1 = 2 := rfl\n",
    "theorem c2 : \xE2\x88\x80 n : Nat, n + 0 = n := fun n => rfl\n",
    "theorem c3 (a b : Nat) : a + b = b + a := Nat.add_comm a b\n",
    "theorem c4 : \xE2\x88\x83 n : Nat, n * n = 16 := \xE2\x9F\xA8" "4, rfl\xE2\x9F\xA9\n",
    "theorem c5 (p q : Prop) (hp : p) (hq : q) : p \xE2\x88\xA7 q := \xE2\x9F\xA8hp, hq\xE2\x9F\xA9\n",
    "theorem c6 (p q : Prop) : p \xE2\x88\xA7 q \xE2\x86\x92 q \xE2\x88\xA7 p := fun \xE2\x9F\xA8hp, hq\xE2\x9F\xA9 \xE2\x86\xA6 \xE2\x9F\xA8hq, hp\xE2\x9F\xA9\n",
    "theorem c7 : (fun x : Nat \xE2\x86\xA6 x + 1) 2 = 3 := rfl\n",
    "example : 2 ^ 10 = 1024 := by decide\n",
    "theorem c9 (x y : Int) : (x + y) ^ 2 = x ^ 2 + 2 * x * y + y ^ 2 := by grind\n",
    "theorem c10 (n : Nat) (h : n < 5) : n < 10 := by omega\n",
    "theorem c11 (a b : Nat) (h : a \xE2\x89\xA4 b) : a < b + 1 := by omega\n",
    "theorem c12 (p : Prop) : \xC2\xAC\xC2\xAC(p \xE2\x88\xA8 \xC2\xACp) := fun h \xE2\x86\xA6 h (Or.inr fun hp \xE2\x86\xA6 h (Or.inl hp))\n",
    "theorem c13 (n : Nat) : n * 0 = 0 := by simp\n",
    "theorem c14 (xs : List Nat) : (xs ++ []).length = xs.length := by simp\n",
    "theorem c15 (a b : Nat) : a + b + 0 = b + a := by\n  have h : a + b = b + a := Nat.add_comm a b\n  show a + b = b + a\n  exact h\n",
    "theorem c16 (a b c : Nat) (h1 : a = b) (h2 : b = c) : a = c := by\n  calc a = b := h1\n    _ = c := h2\n",
    "theorem c17 (n : Nat) : n + 1 \xE2\x89\xA0 0 := by omega\n",
    "theorem c18 : \xE2\x88\x80 n : Nat, n = 0 \xE2\x88\xA8 n \xE2\x89\xA5 1 := by\n  intro n\n  match n with\n  | 0 => exact Or.inl rfl\n  | k + 1 => exact Or.inr (by omega)\n",
    "section Foo\ntheorem bar : True := trivial\nend Foo\n",
    "section\nvariable (n : Nat)\ntheorem c20 : n = n := rfl\nend\n",
    "universe u\ntheorem c21 {\xCE\xB1 : Sort u} (a : \xCE\xB1) : a = a := rfl\n",
    "theorem c22 : Nat.succ 1 = 2 := rfl\n",
    "set_option maxHeartbeats 100000 in\ntheorem c23 : 3 * 7 = 21 := by decide\n",
    "theorem c24 (f : Nat \xE2\x86\x92 Nat) (h : \xE2\x88\x80 x, f x = x) : f (f 2) = 2 := by rw [h, h]\n",
    "theorem c25 (a : Int) : a - a = 0 := by grind\n",
    "theorem c26 (s : List Nat) (h : 3 \xE2\x88\x88 s) : s \xE2\x89\xA0 [] := by\n  intro e\n  simp [e] at h\n",
    "theorem c27 (x : Nat) : x \xE2\x88\xA3 x * 5 := Nat.dvd_mul_right x 5\n",
    "theorem c28 (a b : Nat) : max a b \xE2\x89\xA5 a := Nat.le_max_left a b\n",
    "theorem c29 (h' : 1 = 1) : h' = h' := rfl\n",
    "theorem c30 (\xCE\xB1 \xCE\xB2 : Type) (f : \xCE\xB1 \xE2\x86\x92 \xCE\xB2) (g : \xCE\xB2 \xE2\x86\x92 \xCE\xB1) : (g \xE2\x88\x98 f) = fun x \xE2\x86\xA6 g (f x) := rfl\n",
    "theorem c31 (n : Nat) : (List.range n).length = n := by simp\n",
    "noncomputable section\ntheorem c32 : True := trivial\nend\n",
    "theorem c33 (x : Nat \xC3\x97 Nat) : x = (x.1, x.2) := rfl\n",
    "theorem c34 (a b : Nat) (h : a = b) : b = a := by\n  cases h\n  rfl\n",
    "theorem c35 (p q : Prop) (hq : q) : p \xE2\x86\x92 q := fun _ \xE2\x86\xA6 hq\n",
    NULL
};

/* Inserted with a leading space where a word or a quote could otherwise glue
 * onto a neighbouring identifier (`ab` + `1'` is the prime of `ab1`). */
static const char* const LEAN_REFUSED[] = {
    "\"", "\"x\"", "--", "/-", "-/", "/--", "`", "`(x)", "#", "#eval 1", " #print axioms c1",
    "$", "%", "eval%", "\\", "\\fun", "@[simp]", "\xC2\xAB" "x" "\xC2\xBB", "\xC2\xAB", "\xC2\xBB",
    "\t", "\r", " 1'", " 'a'", " x!'", "\xE2\x86\xA6'", " \xCE\xB1'", "\xE2\x8A\x95",
    "\xE2\x88\x9A", "\xF0\x9D\x94\xB8", "\xE2\x80\x8B", " sorry ", " admit ", " native_decide ",
    " native ", " by_elab ", " run_tac ", " run_cmd ", " eval ", " attribute ", " macro ",
    " elab ", " def ", " instance ", " theorem ", " example ", " open ", " end ", " section ",
    " set_option ", " variable ", " universe ", " namespace ", " syntax ", " initialize ",
    " unsafe ", " IO ", " debug ", " 2theorem ", " 0xftheorem ", " import ", " axiom ",
    " notation ", " lemma ", " deriving ", " mutual ", " private ", " noncomputable ",
    " implemented_by ", " extern ", " ofReduceBool ", " plausible ", " sudo ", " alias ",
    " irreducible_def ", " whatsnew ", " declare_syntax_cat ", " binder_predicate ",
    NULL
};

/* A Mathilda string literal for the bytes s (the caller frees it). */
static char* lean_mstring(const char* s) {
    size_t n = strlen(s);
    char* r = malloc(4 * n + 3);
    size_t w = 0;
    if (!r) abort();
    r[w++] = '"';
    for (size_t i = 0; i < n; i++) {
        char c = s[i];
        if (c == '"' || c == '\\') { r[w++] = '\\'; r[w++] = c; }
        else if (c == '\n') { r[w++] = '\\'; r[w++] = 'n'; }
        else if (c == '\t') { r[w++] = '\\'; r[w++] = 't'; }
        else if (c == '\r') { r[w++] = '\\'; r[w++] = 'r'; }
        else r[w++] = c;
    }
    r[w++] = '"';
    r[w] = '\0';
    return r;
}

/* /tmp/rv4/x1.lean */
static const char* const REVIEW4_X1 =
    "theorem foo : True := trivial\n"
    "  1'\"'\n"
    "#eval IO.FS.writeFile (String.mk [(Char.ofNat 47),(Char.ofNat 100),(Char.ofNat 101),(Char.ofNat 118),(Char.ofNat 47),(Char.ofNat 110),(Char.ofNat 117),(Char.ofNat 108),(Char.ofNat 108)]) (toString 42)\n"
    "-- \"\n";
/* /tmp/rv4/x2.lean */
static const char* const REVIEW4_X2 =
    "theorem warm : True := (fun (_ : Char) \xE2\x86\xA6 trivial) ((fun (_ : Nat) \xE2\x86\xA6'\"') 0)\n"
    "#eval show IO Unit from do let o \xE2\x86\x90 IO.Process.output {cmd := String.ofList ['/','u','s','r','/','b','i','n','/','t','r','u','e'], args := #[]}; IO.print o.stdout; (\xE2\x86\x90 IO.getStdout).flush; IO.Process.exit 0\n"
    "-- \"\n"
    "theorem bad : False := by\n"
    "  decide\n";
/* /tmp/rv4/x3.lean */
static const char* const REVIEW4_X3 =
    "theorem t : True := trivial\n"
    "  x!'\"' \"\n"
    "#eval IO.FS.writeFile (System.FilePath.mk (String.ofList ['/','d','e','v','/','n','u','l','l'])) (toString 3)\n"
    "-- \"\n";
/* /tmp/rv4/x4.lean */
static const char* const REVIEW4_X4 =
    "theorem t : True := trivial\n"
    "/--/ \"\n"
    "-/\n"
    "theorem u : True := trivial\n"
    "#eval IO.FS.writeFile (System.FilePath.mk (String.ofList ['/','d','e','v','/','n','u','l','l'])) (toString 4)\n"
    "-- \"\n";

#endif /* MATHILDA_LEAN_CHECK_CORPUS_H */
