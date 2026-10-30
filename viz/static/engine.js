// Live in-browser Rule 30 engine, bit-packed in 32-bit words. Definitions match rule30.c / xfiles.c:
// rows start at n = 0 with one black cell, column s = -n..n lives at bit MID + s, bit i's left neighbour is
// bit i-1. A chain is a maximal run of >= 2 adjacent 111 neighbourhoods; its seed cells are the chain's cells
// plus the 110 cell to its right, and pyramid membership P persists through 000 neighbourhoods.
//
// tile(b, ty, tx, pat) returns the same bytes as the server's /api/tile (see server.py / r30lib.c).

export const T = 256;

const pc = (x) => {
  x -= (x >>> 1) & 0x55555555;
  x = (x & 0x33333333) + ((x >>> 2) & 0x33333333);
  return (Math.imul((x + (x >>> 4)) & 0x0f0f0f0f, 0x01010101) >>> 24);
};
const span = (w, lo, hi) => {  // absolute cells lo..hi inside word w
  const a = w * 32, b = a + 31;
  if (hi < a || lo > b) return 0;
  if (lo < a) lo = a;
  if (hi > b) hi = b;
  const len = hi - lo + 1;
  return (len === 32 ? 0xffffffff : ((1 << len) - 1) >>> 0) << (lo - a) >>> 0;
};
const lsh = (x, o, w) => ((x[o + w] << 1) | (x[o + w - 1] >>> 31)) >>> 0;  // left neighbour's value
const rsh = (x, o, w) => ((x[o + w] >>> 1) | (x[o + w + 1] << 31)) >>> 0;  // right neighbour's value
const ahead = (x, o, w, j) => (j ? ((x[o + w] >>> j) | (x[o + w + 1] << (32 - j))) >>> 0 : x[o + w]);  // cell i+j

export class LiveEngine {
  constructor(N) {
    this.N = N;
    const MID = (this.MID = Math.ceil((N + 64) / 32) * 32);
    const WN = (this.WN = ((MID + N + 64) >> 5) + 2);
    this.rows = new Uint32Array(N * WN);
    this.P = new Uint32Array(N * WN);
    this.ones = new Int32Array(N);
    this.newp = new Int32Array(N);
    this.active = new Int32Array(N);
    this.cbit = new Uint8Array(N);
    this.hist = new Float64Array(64);
    let cap = Math.max(1024, (N * N) / 24) | 0;
    let cs = new Int32Array(cap), cl = new Uint8Array(cap), nc = 0;
    this.rowOff = new Int32Array(N + 1);

    let row = new Uint32Array(WN), nxt = new Uint32Array(WN), P = new Uint32Array(WN), Pn = new Uint32Array(WN);
    const is111 = new Uint32Array(WN), is000 = new Uint32Array(WN), in2 = new Uint32Array(WN);
    row[MID >> 5] |= 1 << (MID & 31);
    for (let n = 0; n < N; n++) {
      const wlo = (MID - n - 1) >> 5, whi = (MID + n + 1) >> 5;
      this.rows.set(row, n * WN);
      this.P.set(P, n * WN);
      this.rowOff[n] = nc;
      let ones = 0, newp = 0, active = 0;
      for (let w = wlo; w <= whi; w++) {
        const L = lsh(row, 0, w), C = row[w], R = rsh(row, 0, w);
        is111[w] = L & C & R;
        is000[w] = ~(L | C | R) >>> 0;
        nxt[w] = (L ^ (C | R)) >>> 0;
        ones += pc(C);
      }
      for (let w = wlo; w <= whi; w++) in2[w] = (is111[w] & (lsh(is111, 0, w) | rsh(is111, 0, w))) >>> 0;
      let cur = 0;
      for (let w = wlo; w <= whi; w++) {
        const l = lsh(in2, 0, w);
        Pn[w] = (in2[w] | l | (is000[w] & P[w])) >>> 0;
        const starts = (in2[w] & ~l) >>> 0, ends = (in2[w] & ~rsh(in2, 0, w)) >>> 0;
        let ev = (starts | ends) >>> 0;
        newp += pc(starts);
        while (ev) {
          const j = 31 - Math.clz32(ev & -ev);
          ev = (ev & (ev - 1)) >>> 0;
          const pos = w * 32 + j;
          if ((starts >>> j) & 1) cur = pos;
          if ((ends >>> j) & 1) {
            const len = pos - cur + 1;
            if (nc === cap) {
              cap *= 2;
              const a = new Int32Array(cap), b = new Uint8Array(cap);
              a.set(cs); b.set(cl); cs = a; cl = b;
            }
            cs[nc] = cur - MID; cl[nc] = len > 255 ? 255 : len; nc++;
            this.hist[len < 63 ? len : 63]++;
          }
        }
      }
      for (let w = wlo; w <= whi; w++) active += pc((Pn[w] & ~lsh(Pn, 0, w)) >>> 0);
      this.ones[n] = ones; this.newp[n] = newp; this.active[n] = active;
      this.cbit[n] = (row[MID >> 5] >>> (MID & 31)) & 1;
      let t = row; row = nxt; nxt = t;
      t = P; P = Pn; Pn = t;
    }
    this.rowOff[N] = nc;
    this.cs = cs.subarray(0, nc);
    this.cl = cl.subarray(0, nc);
  }

  bit(n, s) {
    const i = this.MID + s;
    return (this.rows[n * this.WN + (i >> 5)] >>> (i & 31)) & 1;
  }

  // chains (n, s, len) inside a window, longest first when over the limit
  chainsIn(n0, n1, s0, s1, lmin, lmax, limit = 20000) {
    const out = [];
    for (let n = Math.max(0, n0); n < Math.min(n1, this.N); n++)
      for (let k = this.rowOff[n]; k < this.rowOff[n + 1]; k++) {
        const s = this.cs[k], L = this.cl[k];
        if (s >= s0 && s < s1 && L >= lmin && L <= lmax) out.push([n, s, L]);
      }
    const total = out.length;
    if (total > limit) out.sort((a, b) => b[2] - a[2]).length = limit;
    return { total, chains: out };
  }

  tile(b, ty, tx, pat) {
    const N = this.N, MID = this.MID, WN = this.WN;
    const n0 = ty * T * b, s0 = tx * T * b, n1 = Math.min(n0 + T * b, N), s1 = s0 + T * b;
    const k = pat ? pat.length : 0;
    const pbits = pat ? [...pat].map((c) => c === '1') : null;
    const bodyBytes = b === 1 ? 2 * T * T : T * T * 6 * 4;
    const buf = new ArrayBuffer(16 + 512 + bodyBytes);
    new Uint32Array(buf, 0, 4).set([b === 1 ? 1 : 2, T, b, pat ? 1 : 0]);
    const hist = new Float64Array(buf, 16, 64);
    if (ty < 0 || n0 >= N || s0 > n1 || s1 <= -n1) return buf;
    const a0 = MID + s0, a1 = MID + s1 - 1, W = s1 - s0;
    const rows = this.rows, Pr = this.P;

    if (b === 1) {
      const cell = new Uint8Array(buf, 528, T * T), aux = new Uint8Array(buf, 528 + T * T, T * T);
      let lab = new Uint8Array(W), labn = new Uint8Array(W);
      const bitAt = (o, i) => (rows[o + (i >> 5)] >>> (i & 31)) & 1;
      for (let r = Math.max(0, n0 - 128); r < n1; r++) {
        const o = r * WN;
        for (let i = 0; i < W; i++) {
          const a = a0 + i;
          labn[i] = bitAt(o, a - 1) | bitAt(o, a) | bitAt(o, a + 1) ? 0 : lab[i];
        }
        const y = r - n0;
        for (let q = this.rowOff[r]; q < this.rowOff[r + 1]; q++) {
          const st = MID + this.cs[q], len = this.cl[q], L = len < 63 ? len : 63;
          for (let i = Math.max(st, a0); i <= Math.min(st + len, a1); i++) labn[i - a0] = L;
          if (r >= n0) {
            for (let i = Math.max(st, a0); i <= Math.min(st + len - 1, a1); i++) cell[y * T + i - a0] |= 2;
            if (st >= a0 && st <= a1) hist[L]++;
          }
        }
        if (r >= n0) {
          for (let i = 0; i < W; i++) {
            const a = a0 + i;
            let v = bitAt(o, a);
            if ((Pr[o + (a >> 5)] >>> (a & 31)) & 1) v |= lab[i] << 2;
            cell[y * T + i] |= v;
          }
          if (k) {
            for (let s = Math.max(s0 - k + 1, -r); s < s1 && s <= r - k + 1; s++) {
              let ok = true;
              for (let j = 0; j < k && ok; j++) ok = bitAt(o, MID + s + j) === +pbits[j];
              if (!ok) continue;
              if (s >= s0) aux[y * T + s - s0] |= 1;
              for (let j = Math.max(0, s0 - s); j < k && s + j < s1; j++) aux[y * T + s + j - s0] |= 2;
            }
          }
        }
        const t = lab; lab = labn; labn = t;
      }
      return buf;
    }

    const blk = new Uint32Array(buf, 528, T * T * 6);
    const wa = a0 >> 5, wb = a1 >> 5;
    for (let r = n0; r < n1; r++) {
      const o = r * WN, y = r - n0, brow = ((y / b) | 0) * T * 6;
      for (let q = this.rowOff[r]; q < this.rowOff[r + 1]; q++) {
        const st = MID + this.cs[q], len = this.cl[q];
        if (st < a0 || st > a1) continue;
        hist[len < 63 ? len : 63]++;
        const p = brow + (((st - a0) / b) | 0) * 6;
        blk[p + 3]++;
        if (len > blk[p + 4]) blk[p + 4] = len;
      }
      for (let w = wa; w <= wb; w++) {
        const inside = (span(w, MID - r, MID + r) & span(w, a0, a1)) >>> 0;
        if (!inside) continue;
        const C = rows[o + w], Pw = Pr[o + w];
        let M = 0;
        if (k) {
          M = span(w, MID - r, MID + r - k + 1);
          for (let j = 0; j < k && M; j++) M = (M & (pbits[j] ? ahead(rows, o, w, j) : ~ahead(rows, o, w, j))) >>> 0;
        }
        const step = b >= 32 ? 32 : b;
        for (let off = 0; off < 32; off += step) {
          const sub = step === 32 ? inside : (inside & (((1 << step) - 1) << off)) >>> 0;
          if (!sub) continue;
          const p = brow + (((w * 32 + off - a0) / b) | 0) * 6;
          blk[p] += pc(sub);
          blk[p + 1] += pc((C & sub) >>> 0);
          blk[p + 2] += pc((Pw & sub) >>> 0);
          blk[p + 5] += pc((M & sub) >>> 0);
        }
      }
    }
    return buf;
  }
}

// Persistent left-wedge structures: chains of length >= 18 with s < -n/4 recurring every 32 rows at a fixed
// distance dl = s + n from the left edge (same rule as server.py).
export function wedgeDiagonals(chains) {
  const by = new Map();
  for (const [n, s, L] of chains) {
    if (L < 18 || 4 * s >= -n) continue;
    const dl = s + n;
    if (!by.has(dl)) by.set(dl, []);
    by.get(dl).push([n, L]);
  }
  const out = [];
  for (const [dl, v] of by) {
    if (v.length < 3) continue;
    v.sort((a, b) => a[0] - b[0]);
    if (v.every((x, i) => i === 0 || (x[0] - v[i - 1][0]) % 32 === 0))
      out.push({ dl, n0: v[0][0], n1: v[v.length - 1][0], count: v.length, phase: v[0][0] % 32,
                 lengths: [...new Set(v.map((x) => x[1]))].sort((a, b) => a - b) });
  }
  return out;
}
