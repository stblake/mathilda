/* check-selection.mjs — unit checks for the notebook tokenizer (mathildaLex) and
 * the bottom-up structural selection (structSelect). Pure logic, no DOM, no
 * CodeMirror; node strips the TS types natively. Run: npm run check:selection */

import { tokenize } from '../src/lib/mathildaLex.ts';
import { buildSpans, ladderAt, chooseExpand, pickAtom, clickGesture, emptyGesture, stepGesture, ladderFrom } from '../src/lib/structSelect.ts';

let failures = 0;
function check(cond, msg) {
  if (cond) { console.log('  ok  ' + msg); }
  else { console.error('  FAIL ' + msg); failures++; }
}
function kinds(src) { return tokenize(src).map((t) => t.kind); }
function texts(src) { return tokenize(src).map((t) => src.slice(t.start, t.end) + ':' + t.kind); }
function hasSpan(spans, from, to) { return spans.some((s) => s.from === from && s.to === to); }
function at(src, sub) { return src.indexOf(sub); }

console.log('tokenizer:');
{
  const src = 'Sin[x] + 2.5 + myVar';
  check(JSON.stringify(texts(src)) ===
        JSON.stringify(['Sin:symbolBuiltin','[:bracket','x:symbolVar',']:bracket',
                        '+:operator','2.5:number','+:operator','myVar:symbolVar']),
        'builtin vs var, bracket, number, operator: ' + JSON.stringify(texts(src)));
}
{
  const src = 'a (* c1 (* nested *) c2 *) b';
  const toks = tokenize(src);
  const comment = toks.find((t) => t.kind === 'comment');
  check(comment && src.slice(comment.start, comment.end) === '(* c1 (* nested *) c2 *)',
        'nested comment is one token');
}
{
  const src = 'x = "line one\nline two" ;';
  const toks = tokenize(src);
  const str = toks.find((t) => t.kind === 'string');
  check(str && src.slice(str.start, str.end) === '"line one\nline two"', 'multiline string one token');
  check(toks.some((t) => t.kind === 'semicolon'), 'lone ; is a semicolon');
}
{
  check(JSON.stringify(kinds('x_ + _Integer + x__h + f[y_.]')) ===
        JSON.stringify(['pattern','operator','pattern','operator','pattern','operator',
                        'symbolVar','bracket','pattern','bracket']),
        'patterns: x_  _Integer  x__h  y_.');
}
{
  check(JSON.stringify(kinds('#1 + ## + #name + %%  + %3')) ===
        JSON.stringify(['slot','operator','slot','operator','slot','operator','out','operator','out']),
        'slots and Out');
}
{
  // 1.. is Repeated of 1, not 1.0 then .  → number '1', operator '..'
  const src = '1..';
  check(JSON.stringify(texts(src)) === JSON.stringify(['1:number','..:operator']),
        '1.. is number then Repeated, not a real');
  check(tokenize('1.5`20')[0].end === 6, 'precision backtick suffix consumed');
  check(tokenize('6.02*^23')[0].end === 8, 'scaled scientific *^ consumed');
}

console.log('selection ladder:');
{
  // a + b*c : clicking b -> b -> whole (NO b*c step, by design)
  const src = 'a + b*c';
  const spans = buildSpans(tokenize(src), src);
  const chain = ladderAt(tokenize(src), spans, at(src, 'b'));
  const asStr = chain.map((s) => src.slice(s.from, s.to));
  check(JSON.stringify(asStr) === JSON.stringify(['b', 'a + b*c']),
        'a + b*c  ladder from b = [b, a + b*c]: ' + JSON.stringify(asStr));
}
{
  // f[a, b*c] : clicking b -> b -> contents "a, b*c" -> group "f[a, b*c]"
  const src = 'f[a, b*c]';
  const spans = buildSpans(tokenize(src), src);
  const b = at(src, 'b');
  const chain = ladderAt(tokenize(src), spans, b).map((s) => src.slice(s.from, s.to));
  check(chain[0] === 'b', 'f[...] ladder bottom is b');
  check(hasSpan(spans, at(src, 'a'), at(src, ']')), 'f[...] contents span "a, b*c" exists');
  check(hasSpan(spans, 0, src.length), 'f[...] whole group span includes head f and brackets');
  check(chain[chain.length - 1] === 'f[a, b*c]', 'f[...] ladder top is f[a, b*c]: ' + JSON.stringify(chain));
}
{
  // nested brackets: {g[x]} clicking x climbs x -> g[x] -> {g[x]}
  const src = '{g[x]}';
  const spans = buildSpans(tokenize(src), src);
  const chain = ladderAt(tokenize(src), spans, at(src, 'x')).map((s) => src.slice(s.from, s.to));
  check(chain.includes('g[x]') && chain[chain.length - 1] === '{g[x]}',
        'nested {g[x]} ladder reaches g[x] then whole: ' + JSON.stringify(chain));
}
{
  // expand/shrink are inverse on the ladder
  const src = 'f[a, b*c]';
  const spans = buildSpans(tokenize(src), src);
  let cur = ladderAt(tokenize(src), spans, at(src, 'b'))[0];
  const up1 = chooseExpand(spans, cur);
  const up2 = up1 && chooseExpand(spans, up1);
  check(up1 && up2 && up2.from <= up1.from && up2.to >= up1.to && (up2.to - up2.from) > (up1.to - up1.from),
        'chooseExpand climbs strictly outward');
  check(chooseExpand(spans, { from: 0, to: src.length }) === null, 'chooseExpand returns null at top');
}
{
  // two statements: spans stay within each half, with absolute offsets
  const src = 'x = 1; y = 2';
  const spans = buildSpans(tokenize(src), src);
  check(hasSpan(spans, 0, at(src, ';')), 'first statement span "x = 1"');
  check(hasSpan(spans, at(src, 'y'), src.length), 'second statement span "y = 2"');
}
{
  // tolerant of incomplete input (no throw, some spans)
  const src = 'f[a, ';
  const spans = buildSpans(tokenize(src), src);
  check(spans.length >= 1, 'incomplete input still yields spans');
}

{
  // clicking the RIGHT edge of a symbol selects the symbol, not the adjacent [
  const src = 'Sin[x]';
  const spans = buildSpans(tokenize(src), src);
  const toks = tokenize(src);
  const atom = pickAtom(toks, 3); // boundary between 'Sin' and '['
  check(atom && src.slice(atom.from, atom.to) === 'Sin', 'boundary click prefers word Sin over [');
}

console.log('click gesture (first-click-selects, precedence-aware):');
{
  // Kernel-style precedence spans for  2 x + 1  (from the real kernel: 2, x, 2 x, 1, 2 x+1)
  const src = '2 x + 1';
  const spans = [[0,1],[2,3],[0,3],[6,7],[0,7]].map(([from,to]) => ({from,to}));
  const toks = tokenize(src);
  const atomX = pickAtom(toks, src.indexOf('x'));
  const txt = (sp) => sp ? src.slice(sp.from, sp.to) : '(none)';

  let g = emptyGesture(), r;
  r = clickGesture(g, spans, atomX); g = r.state;
  check(txt(r.selection) === 'x', 'click 1 on x selects x (first click selects): ' + txt(r.selection));
  r = clickGesture(g, spans, atomX); g = r.state;
  check(txt(r.selection) === '2 x', 'click 2 climbs to 2 x: ' + txt(r.selection));
  r = clickGesture(g, spans, atomX); g = r.state;
  check(txt(r.selection) === '2 x + 1', 'click 3 climbs to 2 x + 1: ' + txt(r.selection));
  r = clickGesture(g, spans, atomX); g = r.state;
  check(txt(r.selection) === 'x', 'click 4 wraps back to x (restartable): ' + txt(r.selection));
}
{
  // THE REPORTED BUG: a second attempt on a different token must select that
  // token, never jump to the whole expression.
  const src = 'f[a, b]';
  const spans = buildSpans(tokenize(src), src);
  const toks = tokenize(src);
  const txt = (sp) => sp ? src.slice(sp.from, sp.to) : '(none)';
  let g = emptyGesture(), r;
  // climb fully on a
  for (let i = 0; i < 6; i++) { r = clickGesture(g, spans, pickAtom(toks, src.indexOf('a'))); g = r.state; }
  // now click a different token b
  r = clickGesture(g, spans, pickAtom(toks, src.indexOf('b'))); g = r.state;
  check(txt(r.selection) === 'b', 'clicking a different token selects THAT token, not the whole expr: ' + txt(r.selection));
}
{
  // local → kernel ladder upgrade mid-gesture: click 1 with the local
  // bracket-level spans selects x; click 2 after the precedence spans arrive
  // climbs to 2 x (not the bracket-level step).
  const src = '2 x + 1';
  const local = buildSpans(tokenize(src), src);         // no "2 x" step
  const kernel = [[0,1],[2,3],[0,3],[6,7],[0,7]].map(([from,to]) => ({from,to}));
  const toks = tokenize(src);
  const atomX = pickAtom(toks, src.indexOf('x'));
  const txt = (sp) => sp ? src.slice(sp.from, sp.to) : '(none)';
  let g = emptyGesture(), r;
  r = clickGesture(g, local, atomX); g = r.state;      // local spans
  check(txt(r.selection) === 'x', 'upgrade: click 1 (local) selects x');
  r = clickGesture(g, kernel, atomX); g = r.state;     // kernel spans now
  check(txt(r.selection) === '2 x', 'upgrade: click 2 (kernel) climbs to 2 x: ' + txt(r.selection));
}
{
  // keyboard stepGesture grows then shrinks along the same ladder
  const src = '2 x + 1';
  const spans = [[0,1],[2,3],[0,3],[6,7],[0,7]].map(([from,to]) => ({from,to}));
  const toks = tokenize(src);
  const atomX = pickAtom(toks, src.indexOf('x'));
  const txt = (sp) => sp ? src.slice(sp.from, sp.to) : '(none)';
  let g = { anchor: atomX, level: 0 }, r;
  r = stepGesture(g, spans, atomX, { from: atomX.from, to: atomX.to }, true); g = r.state;
  check(txt(r.selection) === '2 x', 'Alt-Up grows x -> 2 x');
  r = stepGesture(g, spans, atomX, r.selection, false); g = r.state;
  check(txt(r.selection) === 'x', 'Alt-Down shrinks 2 x -> x');
}

console.log(failures === 0 ? '\nALL PASS' : `\n${failures} FAILURE(S)`);
process.exit(failures === 0 ? 0 : 1);
