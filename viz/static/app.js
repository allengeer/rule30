// Rule 30 Observatory: tiled triangle view (live JS engine or precomputed server runs) and a 10^9-bit
// centre-column strip view. All geometry is in cell units: cell (n, s) occupies x in [s, s+1), y in [n, n+1).
import { LiveEngine, T, wedgeDiagonals } from './engine.js';

const $ = (id) => document.getElementById(id);
const fmt = (x) => (Math.abs(x) >= 1e4 ? Math.round(x).toLocaleString('en-US') : (+x.toFixed(3)).toString());
const BG = [13, 17, 23], OUT = [7, 9, 12], LIGHT = [201, 209, 217], MAGENTA = [255, 79, 216];
const STRIP = 120;

// ---- colour helpers -------------------------------------------------------------------------------------
function turbo(t) {  // Google Turbo, polynomial approximation
  t = Math.min(1, Math.max(0, t));
  const r = 0.13572138 + t * (4.6153926 + t * (-42.66032258 + t * (132.13108234 + t * (-152.94239396 + t * 59.28637943))));
  const g = 0.09140261 + t * (2.19418839 + t * (4.84296658 + t * (-14.18503333 + t * (4.27729857 + t * 2.82956604))));
  const b = 0.1066733 + t * (12.64194608 + t * (-60.58204836 + t * (110.36276771 + t * (-89.90310912 + t * 27.34824973))));
  return [r, g, b].map((v) => Math.round(255 * Math.min(1, Math.max(0, v))));
}
const PAL = Array.from({ length: 64 }, (_, L) => turbo(Math.log2(Math.max(2, L) / 2) / Math.log2(32)));  // 2..64 on a log scale
const mix = (a, b, t) => [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t];
function diverge(z) {
  const t = Math.max(-1, Math.min(1, z / 4));
  return t < 0 ? mix([22, 27, 34], [59, 130, 246], -t) : mix([22, 27, 34], [239, 68, 68], t);
}
const css = (c, a = 1) => `rgba(${c[0] | 0},${c[1] | 0},${c[2] | 0},${a})`;

// ---- state ----------------------------------------------------------------------------------------------
const st = {
  meta: null, src: null, live: null, N: 0, maxB: 1, base: 1 << 30, LMIN: 2,
  view: { x0: -100, y0: -10, cpp: 1 },
  tiles: new Map(), imgs: new Map(), pending: 0, diag: [], chainDots: null, dotsKey: '',
  strip: null, stripKey: '', dirty: true,
};
const L = () => ({
  cells: $('L_cells').checked, tomb: $('L_tomb').checked, dots: $('L_dots').checked, pat: $('L_pat').checked,
  wedge: $('L_wedge').checked, pow2: $('L_pow2').checked, centre: $('L_centre').checked, coin: $('L_coin').checked,
  Lmin: +$('Lmin').value, Lmax: +$('Lmax').value, k: +$('pow2k').value, patv: cleanPat($('pat').value),
});
function cleanPat(s) { s = (s || '').trim(); return /^[01]{1,63}$/.test(s) ? s : ''; }

const tri = $('tri'), ctx = tri.getContext('2d');
function size(cv) {
  const r = cv.getBoundingClientRect(), d = devicePixelRatio || 1;
  if (cv.width !== Math.round(r.width * d) || cv.height !== Math.round(r.height * d)) {
    cv.width = Math.round(r.width * d); cv.height = Math.round(r.height * d);
  }
  return { w: r.width, h: r.height, d };
}
const plotW = () => tri.getBoundingClientRect().width - (L().centre ? STRIP : 0);

// ---- sources --------------------------------------------------------------------------------------------
async function init() {
  st.meta = await (await fetch('/api/meta')).json();
  const sel = $('source');
  sel.innerHTML = '';
  for (const [name, r] of Object.entries(st.meta.runs).sort((a, b) => b[1].N - a[1].N))
    sel.add(new Option(`Precomputed ${name} (${r.N.toLocaleString('en-US')} rows, C tiles)`, name));
  sel.add(new Option('Live JS engine (in browser)', 'live'));
  const h = readHash();
  sel.value = h.src && [...sel.options].some((o) => o.value === h.src) ? h.src : sel.options[0].value;
  await setSource(sel.value, h);
  initCentre();
}

async function setSource(src, h = null) {
  st.src = src;
  st.tiles.clear(); st.imgs.clear(); st.queue = []; st.chainDots = null; st.dotsKey = '';
  $('liveNWrap').style.display = src === 'live' ? '' : 'none';
  if (src === 'live') {
    const N = Math.max(64, Math.min(16384, +$('liveN').value | 0));
    const t0 = performance.now();
    st.live = new LiveEngine(N);
    st.N = N; st.base = 1 << 30; st.LMIN = 2;
    st.maxB = 1 << Math.ceil(Math.log2(Math.max(1, (2 * N) / T)));
    st.diag = wedgeDiagonals(st.live.chainsIn(0, N, -N, N + 1, 18, 255, 1e9).chains);
    $('srcInfo').textContent = `${N.toLocaleString('en-US')} rows evolved in ${(performance.now() - t0).toFixed(0)} ms · ${st.live.cs.length.toLocaleString('en-US')} chains`;
  } else {
    st.live = null;
    const r = st.meta.runs[src];
    st.N = r.N; st.base = r.base; st.LMIN = r.LMIN; st.maxB = Math.max(...r.levels);
    st.diag = await (await fetch(`/api/diagonals?run=${src}`)).json();
    $('srcInfo').textContent = `on-demand below ${r.base} cells/px, precomputed levels ${r.base}…${st.maxB} · ${r.chains_logged.toLocaleString('en-US')} chains ≥ ${r.LMIN} indexed`;
  }
  scaleInfo();
  if (h && h.cpp) Object.assign(st.view, { x0: h.x0, y0: h.y0, cpp: h.cpp });
  else preset('fit');
  dirty();
}

function scaleInfo() {
  const c = st.meta.centre;
  $('scaleInfo').innerHTML = [
    `<b>Live JS engine</b>: up to 16,384 rows (all rows held bit-packed in memory; tiles computed in the page).`,
    ...Object.values(st.meta.runs).map((r) => `<b>${r.name}</b>: ${r.N.toLocaleString('en-US')} rows. Full resolution anywhere via on-demand C evolution from checkpoints every ${r.K} rows; precomputed block levels ${r.base}–${Math.max(...r.levels)} cells/px; every chain ≥ ${r.LMIN} indexed (${r.chains_logged.toLocaleString('en-US')}).`),
    c ? `<b>Centre column</b>: ${c.nbits.toLocaleString('en-US')} bits (Wolfram Data Repository), ${c.ones.toLocaleString('en-US')} ones, ${c.ties} ties, last at n = ${c.last_tie.toLocaleString('en-US')}.` : '',
    `A full 10⁹-row triangle is 10¹⁸ cells and cannot be computed; beyond the triangle runs only the centre column exists.`,
  ].join('<br><br>');
}

// ---- tiles ----------------------------------------------------------------------------------------------
function parseTile(buf) {
  const h = new Uint32Array(buf, 0, 4);
  const t = { kind: h[0], b: h[2], flags: h[3], hist: new Float64Array(buf, 16, 64) };
  if (t.kind === 1) { t.cell = new Uint8Array(buf, 528, T * T); t.aux = new Uint8Array(buf, 528 + T * T, T * T); }
  else t.blk = new Uint32Array(buf, 528, T * T * 6);
  return t;
}
const tileKey = (b, ty, tx, pat) => `${st.src}|${b}|${ty}|${tx}|${pat}`;
function needsPat(b) { return st.live || b < st.base || b <= st.meta.runs[st.src].pat_max; }

// Server tiles go through a small queue: at most MAXFETCH requests in flight, only for tiles wanted by the
// latest frame, and a failed tile waits before it is asked for again (no retry storms).
const MAXFETCH = 6;
st.queue = [];
st.wanted = new Set();
function getTile(b, ty, tx, pat) {
  if (pat && !needsPat(b)) pat = '';
  const key = tileKey(b, ty, tx, pat);
  st.wanted.add(key);
  let e = st.tiles.get(key);
  if (e) {
    st.tiles.delete(key); st.tiles.set(key, e);
    if (!e.tile && !e.inflight && !e.queued && performance.now() >= e.retryAt) enqueue(key, e, b, ty, tx, pat);
    return e.tile;
  }
  e = { tile: null, inflight: false, queued: false, retryAt: 0 };
  st.tiles.set(key, e);
  if (st.tiles.size > 500)
    for (const [k, v] of st.tiles) {
      if (v.inflight || v.queued || st.wanted.has(k)) continue;
      st.tiles.delete(k); st.imgs.delete(k);
      if (st.tiles.size <= 400) break;
    }
  if (st.live) e.tile = parseTile(st.live.tile(b, ty, tx, pat || null));
  else enqueue(key, e, b, ty, tx, pat);
  return e.tile;
}
function enqueue(key, e, b, ty, tx, pat) {
  e.queued = true;
  st.queue.push({ key, e, url: `/api/tile?run=${st.src}&b=${b}&ty=${ty}&tx=${tx}&pat=${pat}`, src: st.src });
  pump();
}
function pump() {
  while (st.pending < MAXFETCH && st.queue.length) {
    const j = st.queue.shift();
    j.e.queued = false;
    if (j.src !== st.src || !st.wanted.has(j.key)) continue;  // scrolled away: ask again if it comes back
    j.e.inflight = true;
    st.pending++;
    fetch(j.url)
      .then((r) => { if (!r.ok) throw new Error(r.status); return r.arrayBuffer(); })
      .then((buf) => { j.e.tile = parseTile(buf); })
      .catch(() => { j.e.retryAt = performance.now() + 3000; })
      .finally(() => { j.e.inflight = false; st.pending--; busy(); pump(); dirty(); });
  }
  busy();
}
function busy() {
  const n = st.pending + st.queue.length;
  $('busy').textContent = n ? `computing ${n} tile${n > 1 ? 's' : ''}…` : '';
}

function layerSig(l) { return [l.cells, l.tomb, l.Lmin, l.Lmax, l.pat && l.patv, l.coin].join(','); }

function colorize(t, ty, tx, l) {
  const cv = document.createElement('canvas');
  cv.width = cv.height = T;
  const g = cv.getContext('2d'), im = g.createImageData(T, T), px = im.data;
  const b = t.b, n0 = ty * T * b, s0 = tx * T * b;
  const put = (i, c) => { px[i * 4] = c[0]; px[i * 4 + 1] = c[1]; px[i * 4 + 2] = c[2]; px[i * 4 + 3] = 255; };
  const inRange = (x) => x >= l.Lmin && x <= l.Lmax;
  if (t.kind === 1) {
    const c = t.cell;
    let zs = null;
    if (l.coin) {  // z over 16x16 groups
      zs = new Float32Array(16 * 16);
      for (let gy = 0; gy < 16; gy++) for (let gx = 0; gx < 16; gx++) {
        let ones = 0, cells = 0;
        for (let y = gy * 16; y < gy * 16 + 16; y++) {
          const n = n0 + y;
          for (let x = gx * 16; x < gx * 16 + 16; x++) {
            const s = s0 + x;
            if (n < st.N && s >= -n && s <= n) { cells++; ones += c[y * T + x] & 1; }
          }
        }
        zs[gy * 16 + gx] = cells ? (ones - cells / 2) / Math.sqrt(cells / 4) : 0;
      }
    }
    for (let y = 0; y < T; y++) {
      const n = n0 + y;
      let runStart = -1;
      for (let x = 0; x < T; x++) {
        const i = y * T + x, s = s0 + x, v = c[i];
        if (n >= st.N || s < -n || s > n) { put(i, OUT); continue; }
        let col = l.coin ? diverge(zs[(y >> 4) * 16 + (x >> 4)]) : BG;
        if (l.cells && (v & 1)) col = l.coin ? mix(col, LIGHT, 0.55) : LIGHT;
        if (l.tomb) {
          const lab = v >> 2;
          if (lab && inRange(lab)) col = mix(col, PAL[lab], 0.8);
          if (v & 2) {
            if (runStart < 0) runStart = x;
            // the cap's length is the label stamped one row below (falls back to the run inside this tile)
            let len = y + 1 < T ? c[i + T] >> 2 : 0;
            if (!len) { let e = x; while (e + 1 < T && (c[y * T + e + 1] & 2)) e++; len = e - runStart + 1; }
            if (inRange(len)) col = mix(PAL[Math.min(63, len)], [255, 255, 255], 0.35);
          } else runStart = -1;
        }
        if (l.pat && l.patv && (t.aux[i] & 2)) col = mix(col, MAGENTA, 0.85);
        put(i, col);
      }
    }
  } else {
    const k = l.patv.length, pexp = Math.pow(2, -k);
    for (let i = 0; i < T * T; i++) {
      const q = t.blk.subarray(i * 6, i * 6 + 6), cells = q[0];
      if (!cells) { put(i, OUT); continue; }
      const d = q[1] / cells;
      let col = BG;
      if (l.coin) col = diverge((q[1] - cells / 2) / Math.sqrt(cells / 4));
      if (l.cells) col = l.coin ? mix(col, LIGHT, d * 0.5) : mix(BG, LIGHT, d * 0.9);
      if (l.tomb && q[3] && inRange(q[4])) col = mix(col, PAL[Math.min(63, q[4])], 0.85);
      if (l.pat && k && (t.flags & 1)) {
        const r = q[5] / cells / pexp;
        if (q[5]) col = mix(col, MAGENTA, Math.min(1, r / 2) * 0.9);
      }
      put(i, col);
    }
  }
  g.putImageData(im, 0, 0);
  return cv;
}

function tileImage(b, ty, tx, l) {
  const pat = l.pat ? l.patv : '';
  const t = getTile(b, ty, tx, pat);
  if (!t) return null;
  const key = tileKey(b, ty, tx, pat && needsPat(b) ? pat : '') + '#' + layerSig(l);
  let im = st.imgs.get(key);
  if (!im) { im = colorize(t, ty, tx, l); st.imgs.set(key, im); if (st.imgs.size > 300) st.imgs.delete(st.imgs.keys().next().value); }
  return im;
}

function levelFor(cpp) {
  let b = 1;
  while (b * 2 <= cpp && b < st.maxB) b *= 2;
  return b;
}

// ---- drawing --------------------------------------------------------------------------------------------
function dirty() { st.dirty = true; }
function frame() {
  if (st.dirty && st.N) { st.dirty = false; draw(); writeHash(); }
  requestAnimationFrame(frame);
}

function visibleTiles(b) {
  const { x0, y0, cpp } = st.view, w = plotW(), h = tri.getBoundingClientRect().height;
  const span = T * b;
  const tx0 = Math.floor(x0 / span), tx1 = Math.floor((x0 + w * cpp) / span);
  const ty0 = Math.max(0, Math.floor(y0 / span)), ty1 = Math.min(Math.floor((y0 + h * cpp) / span), Math.floor((st.N - 1) / span));
  const out = [];
  for (let ty = ty0; ty <= ty1; ty++) for (let tx = tx0; tx <= tx1; tx++) {
    const n0 = ty * span, n1 = Math.min(n0 + span, st.N);
    if (tx * span > n1 || (tx + 1) * span <= -n1) continue;
    out.push([ty, tx]);
  }
  return out;
}

function draw() {
  const { w, h, d } = size(tri), l = L(), { x0, y0, cpp } = st.view, pw = plotW();
  st.wanted = new Set();
  ctx.setTransform(d, 0, 0, d, 0, 0);
  ctx.fillStyle = css(OUT); ctx.fillRect(0, 0, w, h);
  ctx.imageSmoothingEnabled = false;
  const b = levelFor(cpp), span = T * b, scr = span / cpp;
  ctx.save();
  ctx.beginPath(); ctx.rect(0, 0, pw, h); ctx.clip();
  for (const [ty, tx] of visibleTiles(b)) {
    // snap to whole device pixels so neighbouring tiles share an edge without seams
    const Xa = Math.round(((tx * span - x0) / cpp) * d) / d, Ya = Math.round(((ty * span - y0) / cpp) * d) / d;
    const Xb = Math.round((((tx + 1) * span - x0) / cpp) * d) / d, Yb = Math.round((((ty + 1) * span - y0) / cpp) * d) / d;
    const X = Xa, Y = Ya;
    const im = tileImage(b, ty, tx, l);
    if (im) { ctx.drawImage(im, Xa, Ya, Xb - Xa, Yb - Ya); continue; }
    for (let pb = b * 2; pb <= Math.min(st.maxB, b * 64); pb *= 2) {  // a coarser cached tile as a stand-in
      const ps = T * pb, pty = Math.floor((ty * span) / ps), ptx = Math.floor((tx * span) / ps);
      const pk = tileKey(pb, pty, ptx, '');
      if (!st.tiles.get(pk)?.tile) continue;
      const pim = tileImage(pb, pty, ptx, { ...l, pat: false });
      const sx = ((tx * span - ptx * ps) / ps) * T, sy = ((ty * span - pty * ps) / ps) * T, sw = (span / ps) * T;
      ctx.drawImage(pim, sx, sy, sw, sw, X, Y, Xb - Xa, Yb - Ya);
      break;
    }
  }
  drawOverlays(l, b);
  ctx.restore();
  if (l.centre) drawStrip(l, pw, w, h);
  updateStats(l, b);
}

const W2S = (x, y) => [(x - st.view.x0) / st.view.cpp, (y - st.view.y0) / st.view.cpp];

function drawOverlays(l, b) {
  const { x0, y0, cpp } = st.view, h = tri.getBoundingClientRect().height, pw = plotW();
  const yTop = Math.max(0, y0), yBot = Math.min(st.N, y0 + h * cpp);
  ctx.lineWidth = 1;
  // triangle edges
  ctx.strokeStyle = 'rgba(139,148,158,0.35)';
  ctx.beginPath(); ctx.moveTo(...W2S(0, 0)); ctx.lineTo(...W2S(-st.N, st.N)); ctx.moveTo(...W2S(1, 0)); ctx.lineTo(...W2S(st.N + 1, st.N)); ctx.stroke();
  if (l.wedge) {
    ctx.strokeStyle = 'rgba(126,231,135,0.9)'; ctx.setLineDash([6, 4]);
    ctx.beginPath(); ctx.moveTo(...W2S(0, 0)); ctx.lineTo(...W2S(-st.N / 4, st.N)); ctx.stroke();
    ctx.setLineDash([]);
    ctx.strokeStyle = 'rgba(126,231,135,0.55)';
    for (const g of st.diag) {
      const a = W2S(g.dl - g.n0 + 0.5, g.n0), z = W2S(g.dl - g.n1 + 0.5, g.n1 + 1);
      if (Math.max(a[1], z[1]) < 0 || Math.min(a[1], z[1]) > h) continue;
      ctx.beginPath(); ctx.moveTo(...a); ctx.lineTo(...z); ctx.stroke();
    }
    const [lx, ly] = W2S(-yBot / 4, yBot);
    ctx.fillStyle = 'rgba(126,231,135,0.9)'; ctx.font = '11px ui-sans-serif';
    if (lx > 0 && lx < pw) ctx.fillText('s = −n/4', lx + 4, Math.min(ly, h) - 6);
  }
  if (l.pow2) {
    const p = 2 ** l.k, m0 = Math.ceil((yTop + 1) / p), m1 = Math.floor((yBot) / p);
    if (m1 - m0 < h / 3) {
      ctx.strokeStyle = 'rgba(242,204,96,0.6)'; ctx.fillStyle = 'rgba(242,204,96,0.95)';
      for (let m = m0; m <= m1; m++) {
        const n = m * p - 1;
        const [ax, ay] = W2S(-n, n + 0.5), [bx] = W2S(n + 1, n + 0.5);
        if (1 / cpp >= 3) {  // zoomed in: a translucent band over the whole row, labelled with v = v2(n+1)
          ctx.fillStyle = 'rgba(242,204,96,0.35)';
          ctx.fillRect(ax, ay - 0.5 / cpp, bx - ax, 1 / cpp);
          ctx.strokeRect(ax, ay - 0.5 / cpp, bx - ax, 1 / cpp);
          ctx.fillStyle = 'rgba(242,204,96,0.95)';
          ctx.font = '12px ui-sans-serif';
          ctx.fillText(`n+1 = 2^${v2(n + 1)}·odd`, Math.min(bx + 6, pw - 110), ay - 0.5 / cpp - 3);
        } else {
          ctx.beginPath(); ctx.moveTo(ax, ay); ctx.lineTo(bx, ay); ctx.stroke();
        }
        ctx.beginPath(); ctx.arc(bx, ay, 3, 0, 7); ctx.fill();
      }
    } else {
      ctx.fillStyle = 'rgba(242,204,96,0.95)'; ctx.font = '12px ui-sans-serif';
      ctx.fillText(`2^${l.k} rows too dense to draw here: raise k or zoom in`, 10, 18);
    }
  }
  if (l.centre) {
    const [cx] = W2S(0, 0), [cx1] = W2S(1, 0);
    ctx.fillStyle = 'rgba(88,166,255,0.25)';
    const ya = Math.max(0, -y0 / cpp), yb = Math.min(h, (st.N - y0) / cpp);
    if (yb > ya) ctx.fillRect(cx, ya, Math.max(1, cx1 - cx), yb - ya);
  }
  if (l.tomb && l.dots && b >= 4) drawDots(l, b);
}

let dotsTimer = 0;
function drawDots(l, b) {
  const { x0, y0, cpp } = st.view, h = tri.getBoundingClientRect().height, pw = plotW();
  const lmin = Math.max(l.Lmin, st.LMIN, b >= 64 ? 14 : 10), key = [st.src, lmin, l.Lmax, ...[x0, y0, cpp].map((v) => Math.round(v))].join(',');
  if (key !== st.dotsKey) {
    st.dotsKey = key;
    clearTimeout(dotsTimer);
    dotsTimer = setTimeout(async () => {
      const q = [Math.floor(y0), Math.ceil(y0 + h * cpp), Math.floor(x0), Math.ceil(x0 + pw * cpp)];
      if (st.live) st.chainDots = { lmin, ...st.live.chainsIn(q[0], q[1], q[2], q[3], lmin, l.Lmax, 6000) };
      else {
        const buf = await (await fetch(`/api/chains?run=${st.src}&n0=${q[0]}&n1=${q[1]}&s0=${q[2]}&s1=${q[3]}&lmin=${lmin}&lmax=${l.Lmax}&limit=6000`)).arrayBuffer();
        const a = new Int32Array(buf), ch = [];
        for (let i = 2; i < a.length; i += 3) ch.push([a[i], a[i + 1], a[i + 2]]);
        st.chainDots = { lmin, total: a[0], chains: ch };
      }
      dirty();
    }, 120);
  }
  if (!st.chainDots) return;
  for (const [n, s, len] of st.chainDots.chains) {
    const [X, Y] = W2S(s + len / 2, n + 0.5);
    if (X < -4 || X > pw + 4 || Y < -4 || Y > h + 4) continue;
    ctx.fillStyle = css(PAL[Math.min(63, len)], 0.95);
    ctx.beginPath(); ctx.arc(X, Y, len >= 22 ? 3.2 : 2.2, 0, 7); ctx.fill();
  }
}

let stripTimer = 0;
function drawStrip(l, pw, w, h) {
  const { y0, cpp } = st.view;
  ctx.fillStyle = '#0b0f14'; ctx.fillRect(pw, 0, STRIP, h);
  ctx.strokeStyle = '#30363d'; ctx.beginPath(); ctx.moveTo(pw + 0.5, 0); ctx.lineTo(pw + 0.5, h); ctx.stroke();
  const lo = Math.max(0, Math.floor(y0)), hi = Math.min(st.meta.centre?.nbits || 0, Math.ceil(y0 + h * cpp));
  ctx.fillStyle = '#8b949e'; ctx.font = '11px ui-sans-serif';
  ctx.fillText('centre excess E(n)', pw + 6, 13);
  if (!st.meta.centre || hi <= lo + 1) return;
  const key = `${lo},${hi},${Math.round(h)}`;
  if (key !== st.stripKey) {
    st.stripKey = key;
    clearTimeout(stripTimer);
    stripTimer = setTimeout(async () => {
      st.strip = await (await fetch(`/api/centre/series?lo=${lo}&hi=${hi}&px=${Math.round(h)}`)).json();
      dirty();
    }, 80);
  }
  const s = st.strip;
  if (!s) return;
  let m = 1;
  for (let i = 0; i < s.emin.length; i++) m = Math.max(m, Math.abs(s.emin[i]), Math.abs(s.emax[i]));
  const X = (e) => pw + STRIP / 2 + (e / m) * (STRIP / 2 - 8), Y = (p) => (p - y0) / cpp;
  ctx.strokeStyle = 'rgba(139,148,158,0.5)'; ctx.beginPath(); ctx.moveTo(X(0), 18); ctx.lineTo(X(0), h); ctx.stroke();
  ctx.fillStyle = 'rgba(88,166,255,0.35)';
  for (let i = 0; i < s.emin.length; i++) {
    const ya = Y(s.p0[i]), yb = Math.max(ya + 1, Y(s.p1[i]));
    ctx.fillRect(X(s.emin[i]), ya, Math.max(1, X(s.emax[i]) - X(s.emin[i])), yb - ya);
  }
  ctx.strokeStyle = '#58a6ff'; ctx.beginPath();
  for (let i = 0; i < s.eend.length; i++) ctx.lineTo(X(s.eend[i]), Y(s.p1[i]));
  ctx.stroke();
  ctx.fillStyle = '#8b949e';
  ctx.fillText(`±${fmt(m)}`, pw + 6, h - 6);
}

// ---- stats ----------------------------------------------------------------------------------------------
function updateStats(l, b) {
  const { x0, y0, cpp } = st.view, h = tri.getBoundingClientRect().height, pw = plotW();
  const vx0 = x0, vx1 = x0 + pw * cpp, vy0 = Math.max(0, y0), vy1 = Math.min(st.N, y0 + h * cpp);
  const span = T * b;
  let cells = 0, ones = 0, pcells = 0, heads = 0, pat = 0, loaded = 0, missing = 0, exactHist = b === 1;
  const hist = new Float64Array(64);
  for (const [ty, tx] of visibleTiles(b)) {
    const t = getTile(b, ty, tx, l.pat ? l.patv : '');
    if (!t) { missing++; continue; }
    loaded++;
    const n0 = ty * span, s0 = tx * span;
    const r0 = Math.max(0, Math.floor((vy0 - n0) / b)), r1 = Math.min(T, Math.ceil((vy1 - n0) / b));
    const c0 = Math.max(0, Math.floor((vx0 - s0) / b)), c1 = Math.min(T, Math.ceil((vx1 - s0) / b));
    if (t.kind === 1) {
      for (let y = r0; y < r1; y++) {
        const n = n0 + y;
        for (let x = c0; x < c1; x++) {
          const s = s0 + x;
          if (s < -n || s > n) continue;
          const v = t.cell[y * T + x];
          cells++; ones += v & 1; pcells += v >> 2 ? 1 : 0; pat += t.aux[y * T + x] & 1;
          if ((v & 2) && !(x > 0 && (t.cell[y * T + x - 1] & 2))) {
            let e = x; while (e + 1 < T && (t.cell[y * T + e + 1] & 2)) e++;
            const len = e - x + 1;
            heads++; hist[Math.min(63, len)]++;
          }
        }
      }
    } else {
      for (let y = r0; y < r1; y++) for (let x = c0; x < c1; x++) {
        const q = (y * T + x) * 6;
        cells += t.blk[q]; ones += t.blk[q + 1]; pcells += t.blk[q + 2]; heads += t.blk[q + 3]; pat += t.blk[q + 5];
      }
      for (let i = 0; i < 64; i++) hist[i] += t.hist[i];
    }
  }
  const z = cells ? (ones - cells / 2) / Math.sqrt(cells / 4) : 0;
  const k = l.patv.length, pe = cells * Math.pow(2, -k);
  const rows = [
    ['rows', `${fmt(Math.floor(vy0))} – ${fmt(Math.ceil(vy1) - 1)}`],
    ['columns', `${fmt(Math.floor(vx0))} – ${fmt(Math.ceil(vx1) - 1)}`],
    ['resolution', b === 1 ? '1 cell per pixel (exact)' : `${b}×${b}-cell blocks${b >= st.base ? ' (precomputed)' : ' (on demand)'}`],
    ['cells in triangle', fmt(cells)],
    ['black density', cells ? `${(ones / cells).toFixed(6)}  (z = ${z.toFixed(2)})` : '–'],
    ['in a pyramid', cells ? `${(pcells / cells).toFixed(5)}` : '–'],
    ['chains (pyramids born)', `${fmt(heads)}${cells ? `  · ${(heads / cells).toFixed(5)}/cell (coin: 1/32 = 0.03125)` : ''}`],
  ];
  if (l.pat && k) rows.push([`pattern ${l.patv}`, (b === 1 || (b < st.base || (!st.live && b <= st.meta.runs[st.src]?.pat_max)) || st.live)
    ? `${fmt(pat)} starts · coin ≈ ${fmt(pe)} · z = ${pe ? ((pat - pe) / Math.sqrt(pe)).toFixed(2) : '–'}` : `zoom in to ≤ ${st.meta.runs[st.src]?.pat_max} cells/px`]);
  if (missing) rows.push(['', `${missing} tile(s) still loading`]);
  if (!exactHist) rows.push(['', 'histogram summed over whole tiles in view']);
  $('stats').innerHTML = '<table>' + rows.map(([a, v]) => `<tr><td class="muted">${a}</td><td>${v}</td></tr>`).join('') + '</table>';
  $('patNote').textContent = l.pat && !l.patv ? 'pattern must be 1–63 characters of 0/1' : '';
  drawHist($('hist'), [{ h: hist, color: '#c9d1d9' }], 2, l);
}

function drawHist(cv, series, Lfrom, l = null) {
  const g = cv.getContext('2d'), W = cv.width, H = cv.height;
  g.fillStyle = '#0d1117'; g.fillRect(0, 0, W, H);
  let Lmax = 16;
  for (const s of series) for (let i = Lfrom; i < 64; i++) if (s.h[i] > 0) Lmax = Math.max(Lmax, i);
  Lmax = Math.min(63, Lmax + 1);
  let top = 1;
  for (const s of series) for (let i = Lfrom; i <= Lmax; i++) top = Math.max(top, s.h[i]);
  const lt = Math.log10(top) + 0.2, padL = 26, padB = 16;
  const X = (i) => padL + ((i - Lfrom) / (Lmax - Lfrom + 1)) * (W - padL - 4), bw = (W - padL - 4) / (Lmax - Lfrom + 1);
  const Y = (v) => (v <= 0 ? H - padB : H - padB - ((Math.log10(v) + 0.3) / (lt + 0.3)) * (H - padB - 6));
  g.strokeStyle = '#30363d'; g.fillStyle = '#8b949e'; g.font = '10px ui-sans-serif';
  for (let p = 0; p <= lt; p++) { const y = Y(10 ** p); g.beginPath(); g.moveTo(padL, y); g.lineTo(W, y); g.stroke(); g.fillText(p ? `1e${p}` : '1', 2, y + 3); }
  for (let i = Lfrom; i <= Lmax; i += Lmax > 40 ? 8 : 4) g.fillText(i, X(i) + 1, H - 4);
  series.forEach((s, si) => {
    for (let i = Lfrom; i <= Lmax; i++) {
      if (!(s.h[i] > 0)) continue;
      const inR = !l || (i >= l.Lmin && i <= l.Lmax);
      g.fillStyle = l ? css(PAL[i], inR ? 0.9 : 0.25) : s.color;
      const x = X(i) + (series.length > 1 ? si * bw / 2 : 0);
      g.fillRect(x, Y(s.h[i]), Math.max(1, (series.length > 1 ? bw / 2 : bw) - 1), H - padB - Y(s.h[i]));
    }
    const tot = s.h.reduce((a, v, i) => a + (i >= Lfrom ? v : 0), 0);
    if (tot > 0) {  // geometric null: P(L) = 2^-(L - Lfrom + 1)
      g.strokeStyle = s.null || '#f0883e'; g.beginPath();
      for (let i = Lfrom; i <= Lmax; i++) {
        const v = tot * 2 ** -(i - Lfrom + 1);
        if (v < 0.3) break;
        g.lineTo(X(i) + bw / 2, Y(v));
      }
      g.stroke();
    }
  });
}

// ---- hover ----------------------------------------------------------------------------------------------
function v2(x) { let v = 0; while (x > 0 && x % 2 === 0) { x /= 2; v++; } return v; }
tri.addEventListener('mousemove', (ev) => {
  const r = tri.getBoundingClientRect(), mx = ev.clientX - r.left, my = ev.clientY - r.top;
  if (drag) {
    st.view.x0 = drag.x0 - (mx - drag.mx) * st.view.cpp; st.view.y0 = drag.y0 - (my - drag.my) * st.view.cpp; dirty();
  }
  const x = st.view.x0 + mx * st.view.cpp, y = st.view.y0 + my * st.view.cpp;
  const n = Math.floor(y), s = Math.floor(x), l = L();
  if (mx > plotW()) { $('hover').textContent = ''; return; }
  const lines = [`n ${fmt(n)}   s ${fmt(s)}`];
  if (n >= 0 && n < st.N && s >= -n && s <= n) {
    lines.push(`dl ${fmt(s + n)} (left edge)   dr ${fmt(n - s)} (right edge)`);
    lines.push(`s/n ${(s / Math.max(1, n)).toFixed(4)}${4 * s < -n ? '   in left wedge' : ''}   v2(n+1) = ${v2(n + 1)}`);
    const b = levelFor(st.view.cpp), span = T * b, ty = Math.floor(n / span), tx = Math.floor(s / span);
    const t = st.tiles.get(tileKey(b, ty, tx, l.pat && needsPat(b) ? l.patv : ''))?.tile;
    if (t && t.kind === 1) {
      const i = (n - ty * span) * T + (s - tx * span), v = t.cell[i];
      let chain = '';
      if (v & 2) {
        let a = s - tx * span, e = a; const y = n - ty * span;
        while (a > 0 && (t.cell[y * T + a - 1] & 2)) a--;
        while (e + 1 < T && (t.cell[y * T + e + 1] & 2)) e++;
        chain = `   chain cell, length ${e - a + 1}`;
      }
      lines.push(`cell ${v & 1 ? 'black' : 'white'}${v >> 2 ? `   in pyramid of chain length ${v >> 2}` : ''}${chain}${t.aux[i] & 2 ? '   pattern' : ''}`);
    } else if (t) {
      const i = (Math.floor((n - ty * span) / b) * T + Math.floor((s - tx * span) / b)) * 6, q = t.blk;
      if (q[i]) lines.push(`block ${b}×${b}: density ${(q[i + 1] / q[i]).toFixed(4)}  z ${((q[i + 1] - q[i] / 2) / Math.sqrt(q[i] / 4)).toFixed(2)}`,
        `chains ${q[i + 3]}  longest ${q[i + 4]}  pyramid cells ${(q[i + 2] / q[i]).toFixed(3)}`);
    }
  }
  if (st.strip && n >= 0) {
    const sp = st.strip; let j = 0;
    while (j < sp.p1.length - 1 && sp.p1[j] < n + 1) j++;
    if (sp.p1.length) lines.push(`centre E(${fmt(sp.p1[j])}) = ${fmt(sp.eend[j])}${sp.mode === 'bits' && sp.p1[j] - sp.p0[j] === 1 ? '' : '  (pixel end)'}`);
  }
  $('hover').textContent = lines.join('\n');
});
tri.addEventListener('mouseleave', () => { $('hover').textContent = ''; drag = null; });
let drag = null;
tri.addEventListener('mousedown', (ev) => { const r = tri.getBoundingClientRect(); drag = { mx: ev.clientX - r.left, my: ev.clientY - r.top, x0: st.view.x0, y0: st.view.y0 }; });
addEventListener('mouseup', () => { drag = null; });
tri.addEventListener('wheel', (ev) => {
  ev.preventDefault();
  const r = tri.getBoundingClientRect(), mx = ev.clientX - r.left, my = ev.clientY - r.top;
  const f = Math.exp(ev.deltaY * 0.0015), v = st.view;
  const ncpp = Math.min(Math.max(v.cpp * f, 1 / 64), st.maxB * 4);
  const wx = v.x0 + mx * v.cpp, wy = v.y0 + my * v.cpp;
  v.cpp = ncpp; v.x0 = wx - mx * ncpp; v.y0 = wy - my * ncpp; dirty();
}, { passive: false });

// ---- navigation -----------------------------------------------------------------------------------------
function centreOn(s, n, cpp) {
  const h = tri.getBoundingClientRect().height;
  st.view = { cpp, x0: s - (plotW() / 2) * cpp, y0: n - (h / 2) * cpp };
  dirty();
}
function preset(p) {
  const N = st.N, h = tri.getBoundingClientRect().height || 800, w = plotW() || 1000;
  if (p === 'fit') { const cpp = Math.max((2 * N + 2) / w, (N + 2) / h) * 1.04; st.view = { cpp, x0: -cpp * w / 2, y0: -N * 0.02 }; dirty(); }
  if (p === 'apex') st.view = { cpp: 1 / 8, x0: -w / 16, y0: -2 };
  if (p === 'edge4095') centreOn(4095 - 40, 4095, 1 / 6);
  if (p === 'edge65535') centreOn(65535 - 40, 65535, 1 / 6);
  if (p === 'wedge') centreOn(-0.8 * N / 4, 0.8 * N, Math.max(1, N / 4000));
  if (p === 'diag' && st.diag.length) { const g = st.diag[0]; centreOn(g.dl - g.n0 - 20, g.n0 + 60, 1 / 3); }
  if (p === 'tie') centreOn(0, 195112, 1 / 3);
  dirty();
}
document.querySelectorAll('[data-preset]').forEach((bt) => bt.addEventListener('click', () => preset(bt.dataset.preset)));
$('go').addEventListener('click', () => centreOn(+$('goS').value, +$('goN').value, Math.max(1 / 64, +$('goZ').value || 1)));
$('source').addEventListener('change', (e) => setSource(e.target.value));
$('liveN').addEventListener('change', () => st.src === 'live' && setSource('live'));
for (const id of ['L_cells', 'L_tomb', 'L_dots', 'L_pat', 'L_wedge', 'L_pow2', 'L_centre', 'L_coin', 'Lmin', 'Lmax', 'pow2k', 'pat'])
  $(id).addEventListener('input', () => {
    if (+$('Lmin').value > +$('Lmax').value) (id === 'Lmin' ? $('Lmax') : $('Lmin')).value = $(id).value;
    $('LminV').textContent = $('Lmin').value; $('LmaxV').textContent = $('Lmax').value; $('pow2kV').textContent = $('pow2k').value;
    st.dotsKey = ''; legend(); dirty();
  });
function legend() {
  const a = +$('Lmin').value, z = +$('Lmax').value;
  $('tombLegend').innerHTML = Array.from({ length: 62 }, (_, i) => i + 2)
    .map((x) => `<span style="background:${css(PAL[x], x >= a && x <= z ? 1 : 0.15)}" title="${x}"></span>`).join('');
  $('tombLegend').insertAdjacentHTML('afterend', '');
  let lab = document.getElementById('tombLabels');
  if (!lab) { lab = document.createElement('div'); lab.id = 'tombLabels'; lab.className = 'legendLabels'; $('tombLegend').after(lab); }
  lab.innerHTML = '<span>2</span><span>8</span><span>16</span><span>32</span><span>63</span>';
}
addEventListener('resize', () => { dirty(); cdirty(); });

function writeHash() {
  if (!st.src) return;
  const v = st.view;
  history.replaceState(null, '', `#src=${st.src}&x0=${v.x0.toFixed(2)}&y0=${v.y0.toFixed(2)}&cpp=${v.cpp.toPrecision(5)}`);
}
function readHash() {
  const p = new URLSearchParams(location.hash.slice(1));
  return { src: p.get('src'), x0: +p.get('x0'), y0: +p.get('y0'), cpp: +p.get('cpp') };
}

// ---- tabs & export --------------------------------------------------------------------------------------
let tab = 'tri';
document.querySelectorAll('.tab').forEach((bt) => bt.addEventListener('click', () => {
  tab = bt.dataset.tab;
  document.querySelectorAll('.tab').forEach((x) => x.classList.toggle('active', x === bt));
  $('viewTri').classList.toggle('hidden', tab !== 'tri'); $('panelTri').classList.toggle('hidden', tab !== 'tri');
  $('viewCentre').classList.toggle('hidden', tab !== 'centre'); $('panelCentre').classList.toggle('hidden', tab !== 'centre');
  dirty(); cdirty();
}));
$('png').addEventListener('click', () => {
  const cv = tab === 'tri' ? tri : cen;
  const a = document.createElement('a');
  a.download = tab === 'tri'
    ? `rule30_${st.src}_n${Math.round(st.view.y0)}_s${Math.round(st.view.x0)}_cpp${st.view.cpp.toPrecision(3)}.png`
    : `rule30_centre_${Math.round(cv_.lo)}_${Math.round(cv_.hi)}.png`;
  a.href = cv.toDataURL('image/png');
  a.click();
});

// ---- centre column view ---------------------------------------------------------------------------------
const cen = $('cen'), cctx = cen.getContext('2d');
const cv_ = { lo: 0, hi: 1e9, data: null, key: '', bits: null, bitsKey: '', ties: [], find: null, runs: null, runsKey: '', dirty: true };
function cdirty() { cv_.dirty = true; }
async function initCentre() {
  if (!st.meta.centre) return;
  cv_.ties = await (await fetch('/api/centre/ties')).json();
  cframe();
}
function cframe() {
  if (cv_.dirty && tab === 'centre') { cv_.dirty = false; cdraw(); }
  requestAnimationFrame(cframe);
}
let cTimer = 0;
function cdraw() {
  const { w, h, d } = size(cen), NB = st.meta.centre.nbits;
  cv_.lo = Math.max(0, cv_.lo); cv_.hi = Math.min(NB, Math.max(cv_.hi, cv_.lo + 16));
  const lo = Math.floor(cv_.lo), hi = Math.ceil(cv_.hi), px = Math.round(w - 70);
  const key = `${lo},${hi},${px}`;
  if (key !== cv_.key) {
    cv_.key = key;
    clearTimeout(cTimer);
    $('cbusy').textContent = 'loading…';
    cTimer = setTimeout(async () => {
      const [s, r] = await Promise.all([
        fetch(`/api/centre/series?lo=${lo}&hi=${hi}&px=${px}`).then((x) => x.json()),
        fetch(`/api/centre/runs?lo=${lo}&hi=${hi}`).then((x) => x.json()),
      ]);
      cv_.data = s; cv_.runs = r;
      if (hi - lo <= 4 * px) cv_.bits = { lo, bits: new Uint8Array(await (await fetch(`/api/centre/bits?lo=${lo}&hi=${hi}`)).arrayBuffer()) };
      else cv_.bits = null;
      $('cbusy').textContent = '';
      cdirty();
    }, 60);
  }
  cctx.setTransform(d, 0, 0, d, 0, 0);
  cctx.fillStyle = '#0d1117'; cctx.fillRect(0, 0, w, h);
  const s = cv_.data;
  if (!s) return;
  const padL = 60, padR = 10, bitsH = cv_.bits ? 40 : 0, findH = cv_.find ? 50 : 0;
  const top = 20, bot = h - 28 - bitsH - findH;
  const X = (p) => padL + ((p - cv_.lo) / (cv_.hi - cv_.lo)) * (w - padL - padR);
  let mn = Math.min(0, ...s.emin), mx = Math.max(0, ...s.emax);
  if (mx - mn < 2) { mx += 1; mn -= 1; }
  const Y = (e) => bot - ((e - mn) / (mx - mn)) * (bot - top);
  // grid
  cctx.strokeStyle = '#30363d'; cctx.fillStyle = '#8b949e'; cctx.font = '11px ui-sans-serif';
  const step = niceStep((mx - mn) / 6);
  for (let e = Math.ceil(mn / step) * step; e <= mx; e += step) {
    cctx.beginPath(); cctx.moveTo(padL, Y(e)); cctx.lineTo(w - padR, Y(e)); cctx.stroke(); cctx.fillText(fmt(e), 4, Y(e) + 4);
  }
  const xs = niceStep((cv_.hi - cv_.lo) / 8);
  for (let p = Math.ceil(cv_.lo / xs) * xs; p <= cv_.hi; p += xs) cctx.fillText(fmtN(p), X(p) - 10, h - 10 - bitsH - findH);
  cctx.strokeStyle = '#8b949e'; cctx.beginPath(); cctx.moveTo(padL, Y(0)); cctx.lineTo(w - padR, Y(0)); cctx.stroke();
  // envelope + line
  cctx.fillStyle = 'rgba(88,166,255,0.35)';
  for (let i = 0; i < s.emin.length; i++) {
    const xa = X(s.p0[i]), xb = Math.max(xa + 1, X(s.p1[i]));
    cctx.fillRect(xa, Y(s.emax[i]), xb - xa, Math.max(1, Y(s.emin[i]) - Y(s.emax[i])));
  }
  cctx.strokeStyle = '#58a6ff'; cctx.lineWidth = 1.2; cctx.beginPath();
  for (let i = 0; i < s.eend.length; i++) cctx.lineTo(X(s.p1[i]), Y(s.eend[i]));
  cctx.stroke(); cctx.lineWidth = 1;
  // ties
  cctx.fillStyle = '#f0883e';
  for (const t of cv_.ties) if (t >= cv_.lo && t <= cv_.hi) cctx.fillRect(X(t) - 0.5, Y(0) - 5, 1.5, 10);
  const last = st.meta.centre.last_tie;
  if (last >= cv_.lo && last <= cv_.hi) { cctx.fillText(`last tie ${fmtN(last)}`, X(last) + 4, Y(0) - 8); }
  cctx.fillStyle = '#e6edf3'; cctx.fillText(`E(n) = black − n/2 over the first n centre bits  ·  ${s.mode}`, padL, 13);
  // bits
  if (cv_.bits) {
    const y = h - 22 - bitsH - findH + 10, bb = cv_.bits;
    for (let i = 0; i < bb.bits.length; i++) {
      const p = bb.lo + i, xa = X(p), xb = X(p + 1);
      cctx.fillStyle = bb.bits[i] ? '#c9d1d9' : '#21262d';
      cctx.fillRect(xa, y, Math.max(1, xb - xa - (xb - xa > 4 ? 1 : 0)), bitsH - 10);
    }
    if (cv_.find) {
      cctx.fillStyle = 'rgba(255,79,216,0.85)';
      const k = cv_.find.pattern.length;
      for (const p of cv_.find.positions) if (p + k >= cv_.lo && p <= cv_.hi) cctx.fillRect(X(p), y - 4, Math.max(2, X(p + k) - X(p)), 3);
    }
  }
  if (cv_.find && cv_.find.lo === lo && cv_.find.hi === hi) {
    const f = cv_.find, y0 = h - 10 - findH, bmax = Math.max(1, ...f.bins);
    cctx.fillStyle = 'rgba(255,79,216,0.8)';
    const bw = (w - padL - padR) / f.bins.length;
    f.bins.forEach((c, i) => { const hh = (c / bmax) * (findH - 14); cctx.fillRect(padL + i * bw, y0 + findH - 10 - hh, Math.max(1, bw), hh); });
    cctx.fillText(`${f.pattern}: matches per pixel (max ${bmax})`, padL, y0 + 2);
  }
  // window stats
  const e0 = s.eend.length ? s.emin[0] : 0, eEnd = s.eend[s.eend.length - 1];
  const n = hi - lo;
  const ones = (n + 2 * (eEnd - (lo ? s.eend[0] - 0 : 0))) / 2;
  const Elo = lo === 0 ? 0 : null;
  $('cStats').innerHTML = '<table>' + [
    ['bits', `${fmtN(lo)} – ${fmtN(hi)}  (${fmt(n)})`],
    ['E at end', fmt(eEnd)],
    ['min / max E in view', `${fmt(Math.min(...s.emin))} / ${fmt(Math.max(...s.emax))}`],
    ['ties in view', fmt(cv_.ties.filter((t) => t > lo && t <= hi).length)],
    ['resolution', s.mode],
  ].map(([a, v]) => `<tr><td class="muted">${a}</td><td>${v}</td></tr>`).join('') + '</table>';
  void e0; void ones; void Elo;
  if (cv_.runs) drawHist($('cHist'), [{ h: cv_.runs.black, color: '#c9d1d9', null: '#f0883e' }, { h: cv_.runs.white, color: '#58a6ff', null: '#58a6ff' }], 1);
}
function niceStep(x) { const p = 10 ** Math.floor(Math.log10(Math.max(x, 1e-9))), m = x / p; return (m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10) * p; }
function fmtN(p) { return p >= 1e6 && p % 1e5 === 0 ? `${(p / 1e6).toLocaleString('en-US')}M` : Math.round(p).toLocaleString('en-US'); }

let cdrag = null;
cen.addEventListener('mousedown', (ev) => { cdrag = { x: ev.clientX, lo: cv_.lo, hi: cv_.hi }; });
addEventListener('mouseup', () => { cdrag = null; });
cen.addEventListener('mousemove', (ev) => {
  const r = cen.getBoundingClientRect(), mx = ev.clientX - r.left, W = r.width - 70;
  const per = (cv_.hi - cv_.lo) / W;
  if (cdrag) {
    const dp = (ev.clientX - cdrag.x) * per;
    cv_.lo = cdrag.lo - dp; cv_.hi = cdrag.hi - dp;
    if (cv_.lo < 0) { cv_.hi -= cv_.lo; cv_.lo = 0; }
    cdirty();
  }
  const p = cv_.lo + (mx - 60) * per, s = cv_.data;
  if (!s || mx < 60) return;
  let j = 0;
  while (j < s.p1.length - 1 && s.p1[j] < p) j++;
  const lines = [`n ≈ ${fmtN(p)}`, `E(${fmtN(s.p1[j])}) = ${fmt(s.eend[j])}`, `pixel min/max ${fmt(s.emin[j])} / ${fmt(s.emax[j])}`];
  if (cv_.bits) { const i = Math.floor(p) - cv_.bits.lo; if (i >= 0 && i < cv_.bits.bits.length) lines.push(`bit ${fmtN(Math.floor(p))} = ${cv_.bits.bits[i]}`); }
  $('chover').textContent = lines.join('\n');
});
cen.addEventListener('wheel', (ev) => {
  ev.preventDefault();
  const r = cen.getBoundingClientRect(), mx = ev.clientX - r.left, W = r.width - 70;
  const p = cv_.lo + ((mx - 60) / W) * (cv_.hi - cv_.lo), f = Math.exp(ev.deltaY * 0.0015);
  let lo = p - (p - cv_.lo) * f, hi = p + (cv_.hi - p) * f;
  if (hi - lo < 32) { const m = (lo + hi) / 2; lo = m - 16; hi = m + 16; }
  cv_.lo = Math.max(0, lo); cv_.hi = Math.min(st.meta.centre.nbits, hi); cdirty();
}, { passive: false });
const parseN = (s) => Number(String(s).replace(/[, _]/g, ''));
$('cGo').addEventListener('click', () => { cv_.lo = parseN($('cLo').value); cv_.hi = parseN($('cHi').value); cdirty(); });
document.querySelectorAll('[data-cpreset]').forEach((bt) => bt.addEventListener('click', () => {
  const p = bt.dataset.cpreset;
  const r = { all: [0, 1e9], first1m: [0, 1e6], lasttie: [150000, 240000], min3m: [3.4e6, 3.46e6], bits: [0, 400] }[p];
  [cv_.lo, cv_.hi] = r; cdirty();
}));
$('cFind').addEventListener('click', async () => {
  const pat = cleanPat($('cPat').value);
  if (!pat) { $('cFindOut').textContent = 'pattern must be 1–63 characters of 0/1'; return; }
  const lo = Math.floor(cv_.lo), hi = Math.ceil(cv_.hi), px = Math.round(cen.getBoundingClientRect().width - 70);
  $('cFindOut').textContent = 'scanning…';
  const f = await (await fetch(`/api/centre/find?lo=${lo}&hi=${hi}&pat=${pat}&px=${px}`)).json();
  cv_.find = f;
  $('cFindOut').innerHTML = `<table><tr><td class="muted">${pat} in ${fmtN(lo)}–${fmtN(hi)}</td><td>${fmt(f.count)}</td></tr>`
    + `<tr><td class="muted">fair coin expects</td><td>${fmt(f.expected)}  (z = ${f.z.toFixed(2)})</td></tr>`
    + `<tr><td class="muted">first match</td><td>${f.positions.length ? fmtN(f.positions[0]) : '–'}</td></tr></table>`;
  cdirty();
});
$('cClear').addEventListener('click', () => { cv_.find = null; $('cFindOut').textContent = ''; cdirty(); });

legend();
init().then(() => requestAnimationFrame(frame));
