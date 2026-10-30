// Checks the in-browser engine (static/engine.js) under node.
//
//   node validate_js.mjs [N] [server_url] [run]
//
// 1. per-row ones / new chains / active pyramids / centre bit vs rule30.c rows.i64
// 2. chains of length >= 18 vs xfiles.c trace_1m/chains.i64
// 3. tiles byte-for-byte vs the C server's /api/tile for the same run (body and histogram), when a server URL
//    is given; b >= the run's base level is compared on the channels both have (no pattern there)
import { readFileSync } from 'node:fs';
import { LiveEngine, T } from './static/engine.js';

const N = +(process.argv[2] || 10000);
const url = process.argv[3];
const run = process.argv[4] || 'test20k';
let fails = 0;
const check = (name, ok, detail = '') => { fails += !ok; console.log(`${ok ? 'PASS' : 'FAIL'}  ${name}  ${detail}`); };

const t0 = Date.now();
const e = new LiveEngine(N);
console.log(`engine: ${N} rows in ${Date.now() - t0} ms, ${e.cs.length} chains`);

const rows = new BigInt64Array(readFileSync('/Volumes/BLIPS/rule30/2026-10-22_stats_1m/rows.i64').buffer.slice(0, N * 20 * 8));
for (const [name, arr, f] of [['centre bit', e.cbit, 1], ['ones', e.ones, 9], ['new chains', e.newp, 6], ['active', e.active, 7]]) {
  let bad = 0;
  for (let n = 0; n < N; n++) if (BigInt(arr[n]) !== rows[n * 20 + f]) bad++;
  check(`JS engine vs rule30.c rows.i64: ${name}`, bad === 0, `${N} rows, ${bad} differ`);
}
const tr = new BigInt64Array(readFileSync('/Volumes/BLIPS/rule30/2026-10-23_xfiles/trace_1m/chains.i64').buffer);
const ref = [];
for (let i = 0; i < tr.length; i += 3) if (tr[i] < BigInt(N)) ref.push(`${tr[i]},${tr[i + 1]},${tr[i + 2]}`);
const mine = e.chainsIn(0, N, -N, N + 1, 18, 255, 1e9).chains.map((c) => c.join(','));
check('JS engine vs xfiles.c chains.i64 (length >= 18)', mine.join(';') === ref.join(';'), `${mine.length} vs ${ref.length} chains`);

if (url) {
  const meta = await (await fetch(`${url}/api/meta`)).json();
  const base = meta.runs[run].base;
  const cases = [];
  for (const pat of ['', '11011', '0001000']) {
    for (const [b, ty, tx] of [[1, 0, 0], [1, 0, -1], [1, 3, -4], [1, 7, 6], [1, 20, -21], [1, 38, 0], [2, 1, -2], [4, 2, 1],
                               [8, 1, -1], [16, 0, 0], [16, 1, -2], [64, 0, -1], [128, 0, 0], [256, 0, -1]])
      if (ty * T * b < N) cases.push([b, ty, tx, pat]);
  }
  let bad = 0;
  for (const [b, ty, tx, pat] of cases) {
    const theirs = new Uint8Array(await (await fetch(`${url}/api/tile?run=${run}&b=${b}&ty=${ty}&tx=${tx}&pat=${pat}`)).arrayBuffer());
    const ours = new Uint8Array(e.tile(b, ty, tx, pat || null));
    // a JS engine of N rows clips tiles at N; the server run is longer, so compare only rows < N
    // (only whole block rows: a block row straddling row N is partial in the engine but full on the server)
    const rowsIn = Math.min(T, Math.floor((N - ty * T * b) / b));
    let ok = true;
    if (b === 1) {
      for (let i = 0; i < rowsIn * T && ok; i++) ok = ours[528 + i] === theirs[528 + i] && ours[528 + T * T + i] === theirs[528 + T * T + i];
    } else {
      const o = new Uint32Array(ours.buffer, 528), t = new Uint32Array(theirs.buffer, 528);
      const nch = b >= base ? 5 : 6;
      for (let i = 0; i < rowsIn * T && ok; i++) for (let c = 0; c < nch; c++) if (o[i * 6 + c] !== t[i * 6 + c]) ok = false;
    }
    const full = (ty + 1) * T * b <= N;
    if (ok && full && b < base) {  // histograms only match when the whole tile is inside both runs
      const ho = new Float64Array(ours.buffer, 16, 64), ht = new Float64Array(theirs.buffer, 16, 64);
      ok = ho.every((x, i) => x === ht[i]);
    }
    if (!ok) { bad++; console.log('   tile mismatch', b, ty, tx, pat); }
  }
  check('JS engine tiles vs C server tiles (cells, labels, chain cells, pattern, blocks, histograms)', bad === 0, `${cases.length} tiles`);
}
console.log(fails ? `${fails} FAILED` : 'ALL PASS');
process.exit(fails ? 1 : 0);
