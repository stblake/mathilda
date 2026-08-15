// check-prose.mjs — the frontend's behavioural checks for Markdown text cells.
//
// Scope, stated plainly: this covers the PURE half of prose.ts -- what the
// renderer does to a string, and what the four marker constants wrap a selection
// with. It does not and cannot cover the DOM half (wrapSelection's execCommand,
// the click-to-edit swap), because synthetic events do not reach the WKWebView
// this app actually runs in; that last mile needs a real pointer.
//
// The marker constants are READ OUT OF prose.ts rather than duplicated here, so
// this file cannot quietly drift from the thing it is checking.
//
//     npm run check:prose

import { marked } from 'marked';
import { readFileSync } from 'fs';

const OPTS = { async: false, breaks: true, gfm: true };
const render = (s) => (!s.trim() ? '' : marked.parse(s, OPTS));

let fail = 0;
const ok = (name, cond, got) => {
  console.log(`${cond ? 'PASS' : 'FAIL'}  ${name}${cond ? '' : `   got: ${JSON.stringify(got)}`}`);
  if (!cond) fail++;
};

// The claim the docstring makes about `breaks: true`: a single newline is a line
// break, because a prose cell is typed like prose, not like a Markdown file.
const br = render('one\ntwo');
ok('single newline becomes <br>', br.includes('<br>'), br);
const noBreaks = marked.parse('one\ntwo', { async: false, gfm: true });
ok('and would NOT without breaks:true', !noBreaks.includes('<br>'), noBreaks);

ok('bold',        render('**a**').includes('<strong>a</strong>'), render('**a**'));
ok('italic',      render('*a*').includes('<em>a</em>'), render('*a*'));
ok('inline code', render('`a`').includes('<code>a</code>'), render('`a`'));
ok('link',        render('[t](u)').includes('href="u"'), render('[t](u)'));
ok('table (gfm)', render('|a|\n|-|\n|b|').includes('<table>'), 'no table');
ok('empty is exactly the empty string', render('   \n  ') === '', render('   \n  '));

// The marker constants, read out of the source so the test cannot drift from it:
// applying each to a selection must produce the source text it claims, and the
// caret offsets must land where the comments say.
const src = readFileSync('src/lib/prose.ts', 'utf8');
const spec = (name) => {
  const m = src.match(new RegExp(`${name}\\s*=\\s*\\{[^}]*before:\\s*'([^']*)'\\s*,\\s*after:\\s*'([^']*)'\\s*,\\s*caret:\\s*\\{\\s*collapsed:\\s*(\\d+),\\s*selected:\\s*(\\d+)`));
  if (!m) throw new Error(`could not read ${name} from prose.ts`);
  return { before: m[1], after: m[2], collapsed: +m[3], selected: +m[4] };
};
for (const [name, sel, want] of [
  ['PROSE_BOLD',   'word', '**word**'],
  ['PROSE_ITALIC', 'word', '*word*'],
  ['PROSE_CODE',   'word', '`word`'],
  ['PROSE_LINK',   'word', '[word]()'],
]) {
  const s = spec(name);
  const out = s.before + sel + s.after;
  ok(`${name} wraps a selection as ${want}`, out === want, out);
  // Collapsed: the caret must end up immediately after `before`, i.e. between the
  // markers -- so collapsed must equal the length of everything after that point.
  const inserted = s.before + s.after;
  ok(`${name} collapsed caret sits after '${s.before}'`,
     inserted.length - s.collapsed === s.before.length,
     `${inserted.length} - ${s.collapsed}`);
  // And that the rendered result of the wrap is the markup intended.
  if (name !== 'PROSE_LINK') {
    const tag = { PROSE_BOLD: '<strong>', PROSE_ITALIC: '<em>', PROSE_CODE: '<code>' }[name];
    ok(`${name} round-trips through the renderer as ${tag}`, render(out).includes(tag), render(out));
  }
}
// The link's caret goes between the parens either way: ']()' minus 1 lands before ')'.
const link = spec('PROSE_LINK');
ok('PROSE_LINK selected caret sits between the parens',
   link.after.length - link.selected === link.after.indexOf('(') + 1,
   `${link.after} back ${link.selected}`);

console.log(fail === 0 ? '\nall prose checks passed' : `\n${fail} FAILED`);
process.exit(fail === 0 ? 0 : 1);
