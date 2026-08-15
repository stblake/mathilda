// check-prose.mjs — behavioural checks for Markdown text cells.
//
// Compiles src/lib/prose.ts with the project's own tsc and imports the RESULT, so
// these run the shipped renderProse, extractMath and marker constants rather than a
// paraphrase of them. An earlier version re-implemented the marked options here and
// read the constants out with a regex; importing the real thing is strictly better,
// because a paraphrase can agree with the test and disagree with the app.
//
// Covers the PURE half. The DOM half (wrapSelection's execCommand, the
// click-to-edit swap) needs a real pointer in the WKWebView.
//
//     npm run check:prose

import { execFileSync } from 'child_process';
import { mkdtempSync, rmSync } from 'fs';
import { join } from 'path';
import { pathToFileURL } from 'url';

/* Emitted INSIDE the project: prose.ts imports marked and katex, and node resolves
   bare specifiers by walking up from the importing file. */
const out = mkdtempSync(join('node_modules', '.mathilda-prose-'));
try {
  execFileSync('node_modules/.bin/tsc', [
    'src/lib/prose.ts',
    '--target', 'es2020', '--module', 'esnext', '--moduleResolution', 'bundler',
    '--skipLibCheck', '--ignoreConfig', '--outDir', out,
  ], { stdio: 'inherit' });
} catch {
  console.error('tsc failed to compile src/lib/prose.ts');
  process.exit(1);
}
const P = await import(pathToFileURL(join(out, 'prose.js')).href);
const { renderProse, extractMath, PROSE_BOLD, PROSE_ITALIC, PROSE_CODE, PROSE_LINK,
        PROSE_MATH } = P;
/* Imported directly for ONE check: that `breaks: true` is load-bearing. Proving the
   option does something requires calling marked without it. */
const { marked } = await import('marked');

let fail = 0;
const ok = (name, cond, got) => {
  console.log(`${cond ? 'PASS' : 'FAIL'}  ${name}${cond ? '' : `   got: ${JSON.stringify(got)}`}`);
  if (!cond) fail++;
};

// ---- Markdown ----------------------------------------------------------------
// The claim `breaks: true` makes: a single newline is a line break, because a cell
// is typed like prose. The negative control is the point -- without the option
// there is no <br> at all, so the option is load-bearing rather than decorative.
ok('single newline becomes <br>', renderProse('one\ntwo').includes('<br>'),
   renderProse('one\ntwo'));
ok('and would NOT without breaks:true',
   !marked.parse('one\ntwo', { async: false, gfm: true }).includes('<br>'),
   marked.parse('one\ntwo', { async: false, gfm: true }));
ok('bold',        renderProse('**a**').includes('<strong>a</strong>'), renderProse('**a**'));
ok('italic',      renderProse('*a*').includes('<em>a</em>'), renderProse('*a*'));
ok('inline code', renderProse('`a`').includes('<code>a</code>'), renderProse('`a`'));
ok('link',        renderProse('[t](u)').includes('href="u"'), renderProse('[t](u)'));
ok('table (gfm)', renderProse('|a|\n|-|\n|b|').includes('<table>'), 'no table');
ok('empty is exactly the empty string', renderProse('   \n  ') === '', renderProse('   \n  '));

// ---- Math extraction ---------------------------------------------------------
const hasKatex   = (s) => s.includes('class="katex');
const hasDisplay = (s) => s.includes('katex-display');

ok('$x$ renders as math',        hasKatex(renderProse('$x^2$')), renderProse('$x^2$').slice(0, 80));
ok('$$x$$ renders as DISPLAY math', hasDisplay(renderProse('$$x^2$$')), 'not display');
ok('$x$ is inline, not display', !hasDisplay(renderProse('$x^2$')), 'was display');

// THE ordering property: Markdown must not see inside the math. Handed to marked
// first, `$a_1 * b_2$` loses `_1 * b_2` to emphasis and nothing downstream can
// recover it.
const emph = renderProse('$a_1 * b_2$');
ok('markdown does not eat _ and * inside math', hasKatex(emph) && !emph.includes('<em>'), emph.slice(0, 120));

// Context: a dollar inside code is a dollar.
ok('$ inside an inline code span is literal',
   !hasKatex(renderProse('`$x$`')) && renderProse('`$x$`').includes('$x$'), renderProse('`$x$`'));
ok('$ inside a fenced block is literal',
   !hasKatex(renderProse('```\n$x$\n```')), renderProse('```\n$x$\n```'));
ok('an escaped \\$ is literal', !hasKatex(renderProse('\\$5')), renderProse('\\$5'));

// The money heuristic, which is why the delimiters are not a bare regex.
ok('"it cost $5 and $6" stays prose',
   !hasKatex(renderProse('it cost $5 and $6')), renderProse('it cost $5 and $6'));
ok('an unclosed $ is left as typed',
   !hasKatex(renderProse('a $ b')), renderProse('a $ b'));
ok('an unclosed $$ is left as typed',
   !hasKatex(renderProse('a $$ b')), renderProse('a $$ b'));
ok('math does not cross a blank line',
   !hasKatex(renderProse('$a\n\nb$')), renderProse('$a\n\nb$'));

// Math nests inside Markdown structure rather than replacing it.
const boldMath = renderProse('**$x$**');
ok('math inside bold keeps both', boldMath.includes('<strong>') && hasKatex(boldMath), boldMath.slice(0, 100));

// Bad TeX must not take the cell down.
const bad = renderProse('$\\notacommand{x}$');
ok('unparseable TeX still returns a string', typeof bad === 'string' && bad.length > 0, bad);

// The sentinel cannot collide with prose that discusses the sentinel.
const collide = extractMath('@@MATHILDA-MATH-0@@ and $x$');
ok('the math token is chosen so the source cannot already contain it',
   collide.spans.length === 1 && collide.text.includes('@@MATHILDA-MATH-0@@') &&
   collide.text.includes(collide.token(0)) && collide.token(0) !== '@@MATHILDA-MATH-0@@',
   collide.token(0));
ok('and such a cell still renders its math',
   hasKatex(renderProse('@@MATHILDA-MATH-0@@ and $x$')), 'no katex');

// Offsets/counts of extraction itself.
const two = extractMath('$a$ text $$b$$');
ok('extracts both spans with the right display flags',
   two.spans.length === 2 && two.spans[0].display === false && two.spans[0].tex === 'a' &&
   two.spans[1].display === true && two.spans[1].tex === 'b',
   two.spans);

// ---- Marker constants --------------------------------------------------------
for (const [name, spec, sel, want, tag] of [
  ['PROSE_BOLD',   PROSE_BOLD,   'word', '**word**', '<strong>'],
  ['PROSE_ITALIC', PROSE_ITALIC, 'word', '*word*',   '<em>'],
  ['PROSE_CODE',   PROSE_CODE,   'word', '`word`',   '<code>'],
  ['PROSE_LINK',   PROSE_LINK,   'word', '[word]()', null],
  ['PROSE_MATH',   PROSE_MATH,   'x',    '$x$',      null],
]) {
  const wrapped = spec.before + sel + spec.after;
  ok(`${name} wraps a selection as ${want}`, wrapped === want, wrapped);
  const inserted = spec.before + spec.after;
  ok(`${name} collapsed caret sits after '${spec.before}'`,
     inserted.length - spec.caret.collapsed === spec.before.length,
     `${inserted.length} - ${spec.caret.collapsed}`);
  if (tag) ok(`${name} round-trips through the renderer as ${tag}`,
              renderProse(wrapped).includes(tag), renderProse(wrapped));
}
ok('PROSE_MATH round-trips through the renderer as KaTeX',
   renderProse(PROSE_MATH.before + 'x^2' + PROSE_MATH.after).includes('class="katex'),
   renderProse(PROSE_MATH.before + 'x^2' + PROSE_MATH.after));
ok('PROSE_LINK selected caret sits between the parens',
   PROSE_LINK.after.length - PROSE_LINK.caret.selected === PROSE_LINK.after.indexOf('(') + 1,
   `${PROSE_LINK.after} back ${PROSE_LINK.caret.selected}`);

rmSync(out, { recursive: true, force: true });
console.log(fail === 0 ? '\nall prose checks passed' : `\n${fail} FAILED`);
process.exit(fail === 0 ? 0 : 1);
