/-
  The oracle's answer as data.

  Mathilda's answer is never elaborated as Lean syntax.  It is parsed here, by
  a small hand-written parser (a pure function on `String`), into a data AST
  that can only describe

  * rational-number arithmetic: natural-number literals, the variables
    `v1, v2, …` (the goal's free variables, in the order the serializer sent
    them), `+ - * /`, unary `-`, and `^` with a natural-number literal
    exponent;
  * lists `[a, b, …]` and tuples `(a, b, …)` of those.

  Anything else -- another identifier, `by`, `fun`, `sorry`, a string, a
  backquote, any other character -- is a parse error.  Sizes are bounded
  (answer length, digits per numeral, exponent, nesting depth), so a hostile
  answer cannot make the parser or the later `grind` check blow up.

  `Mathilda.Tactics` turns the AST into `Syntax` built only from fixed
  quotations (`$a + $b`, …), numerals and the goal's own fvars, elaborates it at
  the goal's type, and then re-checks the certificate in the kernel.
-/
import Lean

namespace Mathilda.Answer

/-- Bounds on what the parser accepts. -/
def maxLength : Nat := 1000000
def maxDigits : Nat := 2000
def maxExponent : Nat := 10000
def maxDepth : Nat := 100

/-- A polynomial / rational-function term. -/
inductive Poly where
  | num (n : Nat)
  | var (i : Nat)          -- `v<i>`, 1-based
  | add (a b : Poly)
  | sub (a b : Poly)
  | mul (a b : Poly)
  | div (a b : Poly)
  | neg (a : Poly)
  | pow (a : Poly) (k : Nat)
  deriving Inhabited, Repr

/-- An answer: a term, a list or a tuple. -/
inductive Ans where
  | poly (p : Poly)
  | list (xs : Array Ans)
  | tuple (xs : Array Ans)
  deriving Inhabited, Repr

inductive Tok where
  | num (n : Nat)
  | var (i : Nat)
  | sym (c : Char)
  deriving Inhabited, BEq, Repr

def Tok.show : Tok → String
  | .num n => toString n
  | .var i => s!"v{i}"
  | .sym c => c.toString

/-- Split the answer into tokens. -/
def tokenize (s : String) : Except String (Array Tok) := Id.run do
  if s.length > maxLength then
    return .error s!"the answer is longer than {maxLength} characters"
  let cs := s.toList.toArray
  let mut out : Array Tok := #[]
  let mut i := 0
  while h : i < cs.size do
    let c := cs[i]
    if c == ' ' || c == '\n' || c == '\r' || c == '\t' then
      i := i + 1
    else if c.isDigit then
      let mut j := i
      let mut n := 0
      while h2 : j < cs.size do
        if cs[j].isDigit && j - i ≤ maxDigits then
          n := 10 * n + (cs[j].toNat - '0'.toNat)
          j := j + 1
        else break
      if j < cs.size && cs[j]!.isDigit then
        return .error s!"a numeral has more than {maxDigits} digits"
      out := out.push (.num n)
      i := j
    else if c.isAlpha || c == '_' then
      let mut j := i
      while h2 : j < cs.size do
        if cs[j].isAlphanum || cs[j] == '_' || cs[j] == '\'' || cs[j] == '.' then j := j + 1
        else break
      let w := cs.extract i j
      let digits := w.extract 1 w.size
      if w[0]? == some 'v' && 0 < digits.size && digits.size ≤ 6 && digits.all Char.isDigit then
        out := out.push (.var (digits.foldl (fun n c => 10 * n + (c.toNat - '0'.toNat)) 0))
      else
        return .error s!"the identifier `{String.ofList w.toList}` is not a variable v1, v2, …"
      i := j
    else if "+-*/^()[],".contains c then
      out := out.push (.sym c)
      i := i + 1
    else
      return .error s!"unexpected character {repr c}"
  return .ok out

/-- Parser state: the tokens and the position. -/
abbrev P := StateT Nat (ReaderT (Array Tok) (Except String))

def peek : P (Option Tok) := do return (← read)[← get]?
def advance : P Unit := modify (· + 1)
def fail {α} (msg : String) : P α := do
  let t := match (← peek) with | some t => s!"`{t.show}`" | none => "the end"
  throw s!"{msg} at {t}"

def expect (c : Char) : P Unit := do
  if (← peek) == some (.sym c) then advance else fail s!"expected `{c}`"

def asPoly : Ans → P Poly
  | .poly p => pure p
  | _ => fail "a list or tuple cannot be an operand"

/-- Binding powers: `+ -` 65, `* /` 70, unary `-` 75, `^` 75 (right). -/
def binPrec : Char → Option Nat
  | '+' => some 65 | '-' => some 65 | '*' => some 70 | '/' => some 70 | '^' => some 75
  | _ => none

mutual
  /-- A value: a list, a tuple or a term. -/
  partial def value (d : Nat) : P Ans := do
    if d > maxDepth then fail s!"nesting deeper than {maxDepth}"
    if (← peek) == some (.sym '[') then
      advance
      let mut xs := #[]
      if (← peek) == some (.sym ']') then advance; return .list xs
      repeat
        xs := xs.push (← value (d + 1))
        if (← peek) == some (.sym ',') then advance else break
      expect ']'
      return .list xs
    expr 0 d

  /-- A term whose binary operators bind at least `minPrec`. -/
  partial def expr (minPrec d : Nat) : P Ans := do
    let mut lhs ← unary d
    repeat
      let some (.sym c) ← peek | break
      let some p := binPrec c | break
      if p < minPrec then break
      advance
      let a ← asPoly lhs
      if c == '^' then
        let some (.num k) ← peek | fail "the exponent must be a natural-number literal"
        if k > maxExponent then fail s!"exponent {k} exceeds {maxExponent}"
        advance
        lhs := .poly (.pow a k)
      else
        let b ← asPoly (← expr (p + 1) d)
        lhs := .poly <| match c with
          | '+' => .add a b | '-' => .sub a b | '*' => .mul a b | _ => .div a b
    return lhs

  /-- Unary minus, else an atom. -/
  partial def unary (d : Nat) : P Ans := do
    if d > maxDepth then fail s!"nesting deeper than {maxDepth}"
    if (← peek) == some (.sym '-') then
      advance
      return .poly (.neg (← asPoly (← expr 75 (d + 1))))
    atom d

  /-- A numeral, a variable, `( value )` or a tuple `(a, b, …)`. -/
  partial def atom (d : Nat) : P Ans := do
    match (← peek) with
    | some (.num n) => advance; return .poly (.num n)
    | some (.var i) => advance; return .poly (.var i)
    | some (.sym '(') =>
      advance
      let mut xs := #[← value (d + 1)]
      while (← peek) == some (.sym ',') do
        advance
        xs := xs.push (← value (d + 1))
      expect ')'
      if h : xs.size = 1 then return xs[0] else return .tuple xs
    | some (.sym '[') => value d
    | _ => fail "expected a number, a variable v<i>, `(` or `[`"
end

/-- Parse a whole answer. -/
def parse (s : String) : Except String Ans := do
  let toks ← tokenize s
  let (a, pos) ← (value 0).run 0 |>.run toks
  if pos < toks.size then
    throw s!"unexpected `{toks[pos]!.show}` after the answer"
  return a

/-- The largest variable index in `p` (0 if none). -/
def Poly.maxVar : Poly → Nat
  | .num _ => 0
  | .var i => i
  | .add a b | .sub a b | .mul a b | .div a b => max a.maxVar b.maxVar
  | .neg a | .pow a _ => a.maxVar

end Mathilda.Answer
