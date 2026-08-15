// check-snippets.mjs — the Insert palette's templates, judged by Mathilda's parser.
//
// The machinery is CodeMirror's; the templates are ours, and they are the risk. A
// template with an unbalanced bracket produces a broken cell every time the button
// is pressed, and reading `{{${a}, ${b}}, {${c}, ${d}}}` is a poor way to notice.
//
// So each template is expanded (every ${field} left at its placeholder name) and fed
// to the real Mathilda binary inside Hold[...], which parses without evaluating --
// Plot[] must not open a window, and a definition must not define anything. A
// template whose expansion does not parse fails here. That is stronger than counting
// brackets: it is the language's own parser agreeing the insertion is valid.
//
//     npm run check:snippets

import { execFileSync } from 'child_process';
import { mkdtempSync, rmSync, writeFileSync, existsSync } from 'fs';
import { join } from 'path';
import { pathToFileURL } from 'url';

const out = mkdtempSync(join('node_modules', '.mathilda-snip-'));
try {
  execFileSync('node_modules/.bin/tsc', [
    'src/lib/snippets.ts',
    '--target', 'es2020', '--module', 'esnext', '--moduleResolution', 'bundler',
    '--skipLibCheck', '--ignoreConfig', '--outDir', out,
  ], { stdio: 'inherit' });
} catch {
  console.error('tsc failed to compile src/lib/snippets.ts');
  process.exit(1);
}
const { SNIPPETS, expandedText } = await import(pathToFileURL(join(out, 'snippets.js')).href);

let fail = 0;
const ok = (name, cond, got) => {
  console.log(`${cond ? 'PASS' : 'FAIL'}  ${name}${cond ? '' : `   got: ${JSON.stringify(got)}`}`);
  if (!cond) fail++;
};

ok('the palette is not empty', SNIPPETS.length > 0, SNIPPETS.length);
ok('every entry has a label and a template',
   SNIPPETS.every(s => s.label && s.template), SNIPPETS);
ok('labels are unique', new Set(SNIPPETS.map(s => s.label)).size === SNIPPETS.length,
   SNIPPETS.map(s => s.label));
ok('every template has at least one field',
   SNIPPETS.every(s => /\$\{[^}]*\}/.test(s.template)),
   SNIPPETS.filter(s => !/\$\{[^}]*\}/.test(s.template)).map(s => s.label));

for (const s of SNIPPETS) {
  const e = expandedText(s.template);
  ok(`${s.label}: expansion leaves no field syntax`, !e.includes('${'), e);
  // Brackets, counted here as well as parsed below: the count says WHICH kind is
  // unbalanced, where the parser only says the line is wrong.
  for (const [o, c] of [['[', ']'], ['{', '}'], ['(', ')']]) {
    const a = e.split(o).length - 1, b = e.split(c).length - 1;
    ok(`${s.label}: ${o}${c} balance`, a === b, `${a} vs ${b} in ${e}`);
  }
}

// ---- the parser -------------------------------------------------------------
const bin = '../Mathilda';
if (!existsSync(bin)) {
  console.log('\nSKIP  Mathilda binary not built, so the parser half did not run (make -j)');
} else {
  /* One file, one Hold[...] per template. A syntax error stops the whole script
     before anything evaluates and is reported as file:LINE: syntax error, so the
     line number names the template -- which is why they go one per line in order. */
  const lines = SNIPPETS.map(s => `Hold[${expandedText(s.template)}];`);
  const file = join(out, 'snippets.m');
  writeFileSync(file, lines.join('\n') + '\n');
  let res = '', code = 0;
  try {
    res = execFileSync(bin, ['-file', file], { encoding: 'utf8', stdio: ['ignore', 'pipe', 'pipe'] });
  } catch (err) {
    code = err.status ?? 1;
    res = String(err.stdout ?? '') + String(err.stderr ?? '');
  }
  const m = res.match(/:(\d+):\s*syntax error/);
  if (m) {
    const which = SNIPPETS[+m[1] - 1];
    ok(`every expansion parses (failed on ${which ? which.label : 'line ' + m[1]})`, false, res.trim());
  } else {
    ok(`every expansion parses in Mathilda (${SNIPPETS.length} templates)`, code === 0, res.trim());
  }
}

rmSync(out, { recursive: true, force: true });
console.log(fail === 0 ? '\nall snippet checks passed' : `\n${fail} FAILED`);
process.exit(fail === 0 ? 0 : 1);
