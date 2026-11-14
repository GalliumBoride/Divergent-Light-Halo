// 曲面网格贴边测试. 原则: 精确边界线 / 脊线 (特征线) 必须是网格的边.
//   路标 = 层线落在特征线上的点 (独立于缝合器另算一遍). 沿每条特征线把路标按弧长排:
//     相邻两个路标在相邻两层上 → 网格里必须有这条边 (缺顶点 missV / 缺边 missE), 两者之间特征线上的点到网格的距离应只有弦高 (否则 teeth);
//     在同一层上 → 特征线在这一层外鼓出去一段 (极值帽 cap), 网格够不到, 单列;
//     隔了层 → skip (中间的层在这里没有线), 单列.
// node test/mesh2/edges.mjs code mode el D surf fix   (例: 4-1-8 lowitz 3.5 15 caustic 0; mode = column | lowitz)
//   env: STITCH (换一个 stitch2.js 来比), FEAT=0 (不给缝合器特征线), DUMP=n [GRP=cusp] (打印明细), KINK=n, NOCACHE=1
//   层与特征线算一次后存在 test/mesh2/cache/ (已 gitignore)
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
const HERE = path.dirname(fileURLToPath(import.meta.url));
const B = process.env.ENG || new URL('../../src/engine/', import.meta.url).href;
const CR = await import(B + 'crystal.js'), LPm = await import(B + 'lp.js'), EP = await import(B + 'epdiv.js'), ST = await import(process.env.STITCH || (B + 'stitch2.js'));
const PG = await import('./pagehs.mjs');
const [code = '4-1-8', mode = 'lowitz', el = '3.5', Ds = '15', surfName = 'caustic', fixAng = '0'] = process.argv.slice(2);
const D = +Ds;
CR.setApex(56, 56); LPm.LP.partial = true; LPm.LP.upperApexDeg = 56; LPm.LP.lowerApexDeg = 56; LPm.setIndex(1.31);
const lp = CR.parsePath(code).map(CR.lpFaceFromCompact), surf = EP.resolveSurface(lp, surfName);
const low = mode === 'lowitz', z = low ? 0 : +fixAng, th = low ? +fixAng : 0;
const L = [D * Math.cos(el * Math.PI / 180), 0, D * Math.sin(el * Math.PI / 180)];
// 场景的层与特征线算一次存盘 (cache/), 换缝合器重跑时直接读
fs.mkdirSync(path.join(HERE, 'cache'), { recursive: true });
const cname = `${HERE}/cache/${[code, mode, el, Ds, surfName, fixAng].join('_')}${process.env.HF ? '_' + process.env.HF.replace(/\W/g, '') : ''}.json`;
const enc = a => Buffer.from(a.buffer, a.byteOffset, a.byteLength).toString('base64'), dec = s => { const b = Buffer.from(s, 'base64'); return new Float32Array(b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength)); };
let feats = [], layers = [];
if (fs.existsSync(cname) && !process.env.NOCACHE) {
  const j = JSON.parse(fs.readFileSync(cname, 'utf8'));
  feats = j.feats.map(f => ({ g: f.g, p: dec(f.p) }));
  layers = j.layers.map(l => ({ h: l.h, curves: l.curves.map(c => ({ p: dec(c.p), d: dec(c.d), x: dec(c.x), b: c.b })) }));
} else {
  const e = EP.edges2Dof(lp, +el, D, 0, 80, low, z, th, surf, null, true);
  const rimsC = surf === 'caustic' ? (EP.causticRims2Dof(lp, +el, D, 0, 80, low, z, th, 101) || []) : [];
  const groups = { entry: e.entry, tir: e.tir, seam: e.seam, turn: e.turn, cusp: e.cusp, rim: rimsC };
  for (const g in groups) for (const s of groups[g] || []) if ((s.p || s).length >= 6) feats.push({ g, p: s.p || s });
  if (process.env.HF) { for (const h of fs.readFileSync(process.env.HF, 'utf8').trim().split(',').map(Number)) { const l = EP.lightCone2DofLayers(lp, +el, D, 0, 1, 101, 80, low, z, th, surf, 0, 1, h, true)[0]; if (l.curves.length) layers.push({ h, curves: l.curves }); } }
  else layers = PG.pageLayers(EP, lp, +el, D, low, z, th, surf);
  fs.writeFileSync(cname, JSON.stringify({ feats: feats.map(f => ({ g: f.g, p: enc(f.p) })), layers: layers.map(l => ({ h: l.h, curves: l.curves.map(c => ({ p: enc(c.p), d: enc(c.d), x: enc(c.x), b: c.b || 0 })) })) }));
}
const useFeat = process.env.FEAT !== '0';
const t0 = Date.now();
const res = PG.pageMesh(EP, ST, lp, +el, D, low, z, th, surf, L, layers, useFeat ? feats.map(f => ({ p: f.p })) : null, !!process.env.HF);
const m = res.m, tS = Date.now() - t0;
layers = res.layers;

// ---- 网格: 顶点按位置焊接, 边表, 三角形分桶
const P = m.pos, nT = P.length / 9;
const vkey = (x, y, z) => Math.round(x * 2e4) + ',' + Math.round(y * 2e4) + ',' + Math.round(z * 2e4);
const vid = new Map(), edges = new Set();
const vpos = [], adj = [];
const idOf = o => { const k = vkey(P[o], P[o + 1], P[o + 2]); let i = vid.get(k); if (i === undefined) { vid.set(k, i = vid.size); vpos.push([P[o], P[o + 1], P[o + 2]]); adj.push([]); } return i; };
for (let t = 0; t < nT; t++) { const a = idOf(9 * t), b = idOf(9 * t + 3), c = idOf(9 * t + 6); for (const [u, v] of [[a, b], [b, c], [c, a]]) { const k = u < v ? u * 4194304 + v : v * 4194304 + u; if (!edges.has(k)) { edges.add(k); adj[u].push(v); adj[v].push(u); } } }
const hasEdge = (pa, pb) => { const a = vid.get(vkey(...pa)), b = vid.get(vkey(...pb)); if (a === undefined || b === undefined) return -1; if (a === b) return 1; return edges.has(a < b ? a * 4194304 + b : b * 4194304 + a) ? 1 : 0; };
const CELL = Math.max(0.25, D / 15), grid = new Map(), key = (i, j, k) => (i + 512) * 1048576 + (j + 512) * 1024 + (k + 512);
for (let t = 0; t < nT; t++) {
  const o = 9 * t, lo = [0, 1, 2].map(a => Math.floor(Math.min(P[o + a], P[o + 3 + a], P[o + 6 + a]) / CELL)), hi = [0, 1, 2].map(a => Math.floor(Math.max(P[o + a], P[o + 3 + a], P[o + 6 + a]) / CELL));
  for (let i = lo[0]; i <= hi[0]; i++) for (let j = lo[1]; j <= hi[1]; j++) for (let k = lo[2]; k <= hi[2]; k++) { const q = key(i, j, k); let a = grid.get(q); if (!a) grid.set(q, a = []); a.push(t); }
}
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]], dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
function ptTri(p, a, b, c) {      // Ericson 最近点
  const ab = sub(b, a), ac = sub(c, a), ap = sub(p, a), d1 = dot(ab, ap), d2 = dot(ac, ap);
  if (d1 <= 0 && d2 <= 0) return Math.hypot(...ap);
  const bp = sub(p, b), d3 = dot(ab, bp), d4 = dot(ac, bp);
  if (d3 >= 0 && d4 <= d3) return Math.hypot(...bp);
  const vc = d1 * d4 - d3 * d2;
  if (vc <= 0 && d1 >= 0 && d3 <= 0) { const v = d1 / (d1 - d3); return Math.hypot(ap[0] - v * ab[0], ap[1] - v * ab[1], ap[2] - v * ab[2]); }
  const cp = sub(p, c), d5 = dot(ab, cp), d6 = dot(ac, cp);
  if (d6 >= 0 && d5 <= d6) return Math.hypot(...cp);
  const vb = d5 * d2 - d1 * d6;
  if (vb <= 0 && d2 >= 0 && d6 <= 0) { const w = d2 / (d2 - d6); return Math.hypot(ap[0] - w * ac[0], ap[1] - w * ac[1], ap[2] - w * ac[2]); }
  const va = d3 * d6 - d5 * d4;
  if (va <= 0 && d4 - d3 >= 0 && d5 - d6 >= 0) { const w = (d4 - d3) / ((d4 - d3) + (d5 - d6)); return Math.hypot(bp[0] - w * (c[0] - b[0]), bp[1] - w * (c[1] - b[1]), bp[2] - w * (c[2] - b[2])); }
  const den = 1 / (va + vb + vc), v = vb * den, w = vc * den;
  return Math.hypot(ap[0] - v * ab[0] - w * ac[0], ap[1] - v * ab[1] - w * ac[1], ap[2] - v * ab[2] - w * ac[2]);
}
function query(p) {               // → [d_mesh, d_vert]
  const c = p.map(x => Math.floor(x / CELL));
  let dm = Infinity, dv = Infinity;
  for (let R = 1; R <= 6; R++) {
    for (let i = c[0] - R; i <= c[0] + R; i++) for (let j = c[1] - R; j <= c[1] + R; j++) for (let k = c[2] - R; k <= c[2] + R; k++) {
      if (R > 1 && Math.max(Math.abs(i - c[0]), Math.abs(j - c[1]), Math.abs(k - c[2])) < R) continue;
      const a = grid.get(key(i, j, k)); if (!a) continue;
      for (const t of a) {
        const o = 9 * t, A = [P[o], P[o + 1], P[o + 2]], Bq = [P[o + 3], P[o + 4], P[o + 5]], C = [P[o + 6], P[o + 7], P[o + 8]];
        const d = ptTri(p, A, Bq, C); if (d < dm) dm = d;
        for (const V of [A, Bq, C]) { const dd = Math.hypot(p[0] - V[0], p[1] - V[1], p[2] - V[2]); if (dd < dv) dv = dd; }
      }
    }
    if (dv < (R - 0.5) * CELL) break;
  }
  return [dm, dv];
}

// ---- 路标 (测试自己算): 每层每条曲线的每个顶点到各特征线的距离, 容差内连着的一串取最近的
const rl = q => Math.hypot(q[0] - L[0], q[1] - L[1], q[2] - L[2]);
const tolAt = r => 1.5e-4 * D + 3e-4 * r;
const FC = Math.max(0.05 * D, 0.03), fgrid = new Map(), fkey = (i, j, k) => ((i + 4096) * 8192 + (j + 4096)) * 8192 + (k + 4096);
const FA = feats.map(f => { const p = f.p, n = p.length / 3, t = new Float64Array(n); for (let i = 1; i < n; i++) t[i] = t[i - 1] + Math.hypot(p[3 * i] - p[3 * i - 3], p[3 * i + 1] - p[3 * i - 2], p[3 * i + 2] - p[3 * i - 1]); return t; });
feats.forEach((f, fi) => {
  const p = f.p;
  for (let i = 0; i + 1 < p.length / 3; i++) {
    const o = 3 * i, pad = 0.04 + 4e-4 * 100;
    const lo = [0, 1, 2].map(a => Math.floor((Math.min(p[o + a], p[o + 3 + a]) - pad) / FC)), hi = [0, 1, 2].map(a => Math.floor((Math.max(p[o + a], p[o + 3 + a]) + pad) / FC));
    for (let a = lo[0]; a <= hi[0]; a++) for (let b = lo[1]; b <= hi[1]; b++) for (let c = lo[2]; c <= hi[2]; c++) { const k = fkey(a, b, c); let q = fgrid.get(k); if (!q) fgrid.set(k, q = []); q.push(fi, i); }
  }
});
const nearF = (q, tl) => {
  const a = fgrid.get(fkey(Math.floor(q[0] / FC), Math.floor(q[1] / FC), Math.floor(q[2] / FC)));
  if (!a) return null;
  const out = new Map();
  for (let s = 0; s < a.length; s += 2) {
    const fi = a[s], i = a[s + 1], p = feats[fi].p, o = 3 * i;
    const b = [p[o + 3] - p[o], p[o + 4] - p[o + 1], p[o + 5] - p[o + 2]], w = [q[0] - p[o], q[1] - p[o + 1], q[2] - p[o + 2]], l2 = dot(b, b);
    let u = l2 > 0 ? dot(w, b) / l2 : 0; u = Math.max(0, Math.min(1, u));
    const d = Math.hypot(w[0] - u * b[0], w[1] - u * b[1], w[2] - u * b[2]);
    if (d > tl) continue;
    const t = FA[fi][i] + u * Math.sqrt(l2), old = out.get(fi);
    if (!old || d < old.d) out.set(fi, { d, t });
  }
  return out.size ? out : null;
};
const marks = feats.map(() => []);       // 每条特征线: [{ t, row, pos }]
layers.forEach((l, row) => {
  for (const c of l.curves) {
    const p = c.p, n = p.length / 3, hit = new Array(n);
    for (let j = 0; j < n; j++) { const q = [p[3 * j], p[3 * j + 1], p[3 * j + 2]]; hit[j] = nearF(q, tolAt(rl(q))); }
    const lm = [];
    for (let j = 0; j < n; j++) if (hit[j]) for (const [fi, q] of hit[j]) {
      if (j > 0 && hit[j - 1] && hit[j - 1].has(fi)) continue;
      let bj = j, bd = q.d, bt = q.t;
      for (let k = j + 1; k < n && hit[k] && hit[k].has(fi); k++) { const w = hit[k].get(fi); if (w.d < bd) { bd = w.d; bj = k; bt = w.t; } }
      lm.push({ j: bj, fi, t: bt });
    }
    // 并在一起走的两条特征线: 各自最近的顶点差一两个, 相距不到两倍容差的算同一个路标 (取序号小的那个的位置)
    lm.sort((a, b) => a.j - b.j);
    let ref = -1;
    for (const w of lm) {
      const q = [p[3 * w.j], p[3 * w.j + 1], p[3 * w.j + 2]];
      if (ref >= 0 && Math.hypot(q[0] - p[3 * ref], q[1] - p[3 * ref + 1], q[2] - p[3 * ref + 2]) < 2 * tolAt(rl(q))) w.j = ref; else ref = w.j;
      marks[w.fi].push({ t: w.t, row, pos: [p[3 * w.j], p[3 * w.j + 1], p[3 * w.j + 2]] });
    }
  }
});
const featPt = (fi, t) => {            // 特征线上弧长 t 处的点
  const A = FA[fi], p = feats[fi].p; let lo = 0, hi = A.length - 1;
  while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (A[mid] <= t) lo = mid; else hi = mid; }
  const u = A[hi] > A[lo] ? (t - A[lo]) / (A[hi] - A[lo]) : 0;
  return [0, 1, 2].map(a => p[3 * lo + a] + u * (p[3 * hi + a] - p[3 * lo + a]));
};
const stat = {}, detail = [];
const S = g => stat[g] || (stat[g] = { marks: 0, links: 0, missV: 0, missE: 0, samp: 0, teeth: 0, caps: 0, capLen: 0, skips: 0 });
feats.forEach((f, fi) => {
  const a = marks[fi].sort((x, y) => x.t - y.t), g = [], st = S(f.g);
  // 同一层同一位置的几个 (几条曲线在这里相交 / 收头; 位置差在毫米以内但不逐位相等) 并成一个路标, 各自的位置都记着
  for (const mk of a) { const q = g[g.length - 1]; if (q && q.row === mk.row && mk.t - q.t1 < 8 * tolAt(rl(mk.pos))) { q.all.push(mk.pos); q.t1 = mk.t; continue; } mk.all = [mk.pos]; mk.t1 = mk.t; g.push(mk); }
  st.marks += g.length;
  const prs = [];
  for (let i = 0; i + 1 < g.length; i++) {
    prs.push([g[i], g[i + 1]]);
    // 两层共用的捏合点: 两个身份都与前后的路标成对
    if (g[i].row !== g[i + 1].row && Math.hypot(g[i].pos[0] - g[i + 1].pos[0], g[i].pos[1] - g[i + 1].pos[1], g[i].pos[2] - g[i + 1].pos[2]) < 2e-4 * D) {
      if (i + 2 < g.length && Math.abs(g[i].row - g[i + 2].row) === 1) prs.push([g[i], g[i + 2]]);
      if (i > 0 && Math.abs(g[i - 1].row - g[i + 1].row) === 1) prs.push([g[i - 1], g[i + 1]]);
    }
  }
  for (const [A, Bm] of prs) {
    const dt = Bm.t - A.t1, r = rl(A.pos);
    if (r < 0.02 * D || rl(Bm.pos) < 0.02 * D || Math.hypot(...A.pos) < 0.02 * D || Math.hypot(...Bm.pos) < 0.02 * D) continue;      // 灯 / 眼旁的捏合点不算
    if (A.row === Bm.row) { if (dt > 2 * tolAt(r)) { st.caps++; st.capLen += dt; } continue; }
    if (Math.abs(A.row - Bm.row) > 1) { st.skips++; continue; }
    if (dt > 0.5 * Math.max(r, rl(Bm.pos)) + 0.05 * D) { st.skips++; continue; }
    st.links++;
    let he = -1;
    for (const pa of A.all) for (const pb of Bm.all) he = Math.max(he, hasEdge(pa, pb));
    // 没有直接的边: 沿特征线插了点 (边界线在两个路标之间是弯的) —— 从 A 沿网格的边走, 只走落在这条特征线这一段上的顶点, 走得到 B 也算
    if (he === 0) {
      const tgt = new Set(Bm.all.map(q => vid.get(vkey(...q))).filter(x => x !== undefined)), seen = new Set(A.all.map(q => vid.get(vkey(...q))).filter(x => x !== undefined));
      let fr = [...seen];
      for (let it = 0; it < 400 && fr.length && he === 0; it++) {
        const nx = [];
        for (const u of fr) for (const v of adj[u]) {
          if (seen.has(v)) continue;
          if (tgt.has(v)) { he = 1; break; }
          const hq = nearF(vpos[v], 1.5 * tolAt(rl(vpos[v]))), w = hq && hq.get(fi);
          if (w && w.t > A.t - 1e-6 && w.t < Bm.t1 + 1e-6) { seen.add(v); nx.push(v); }
        }
        fr = nx;
      }
    }
    if (he < 0) { st.missV++; detail.push([f.g, 'missV', A, Bm]); } else if (he === 0) { st.missE++; detail.push([f.g, 'missE', A, Bm]); }
    for (const u of [0.25, 0.5, 0.75]) {
      const q = featPt(fi, A.t1 + u * dt), [dm, dv] = query(q);
      st.samp++;
      if (dm > tolAt(rl(q)) + 0.02 * dt && dm > 0.15 * dv) { st.teeth++; detail.push([f.g, 'tooth', A, Bm, dm, dv, q]); }
    }
  }
});
// 折叠 / 串面的三角形: 自身法向与顶点法向 (相邻三角形的面积加权平均) 夹角 > 60° —— 画出来是面上一块块暗斑
let kinks = 0, kArea = 0; const kPos = [];
for (let t = 0; t < nT; t++) {
  const o = 9 * t, u = [P[o + 3] - P[o], P[o + 4] - P[o + 1], P[o + 5] - P[o + 2]], v = [P[o + 6] - P[o], P[o + 7] - P[o + 1], P[o + 8] - P[o + 2]];
  const n = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]], l = Math.hypot(...n);
  if (!(l > 1e-12)) continue;
  let bad = false;
  for (let k = 0; k < 3; k++) { const q = o + 3 * k; if (Math.abs(n[0] * m.nrm[q] + n[1] * m.nrm[q + 1] + n[2] * m.nrm[q + 2]) / l < 0.5) bad = true; }
  if (bad && l / 2 > 1e-5 * D * D) { kinks++; kArea += l / 2; kPos.push([l / 2, P[o], P[o + 1], P[o + 2]]); }
}
const out = [`${code} ${mode} Σ${el} D${D} ${surf} fix${fixAng}: layers ${layers.length} tris ${m.nTri} bad ${m.bad.length} stitch ${tS} ms kinks ${kinks} (${kArea.toFixed(2)} m2)`];
if (process.env.KINK) { kPos.sort((a, b) => b[0] - a[0]); for (const k of kPos.slice(0, +process.env.KINK)) console.log('kink area', k[0].toFixed(4), 'at', k.slice(1).map(x => x.toFixed(2)).join(',')); }
for (const g in stat) { const s = stat[g]; if (s.marks) out.push(`  ${g.padEnd(5)} marks ${String(s.marks).padStart(5)} links ${String(s.links).padStart(5)} missV ${String(s.missV).padStart(4)} missE ${String(s.missE).padStart(4)} teeth ${String(s.teeth).padStart(4)}/${String(s.samp).padEnd(5)} caps ${String(s.caps).padStart(3)} (${s.capLen.toFixed(2)} m) skips ${s.skips}`); }
console.log(out.join('\n'));
const fmtP = q => q.map(x => x.toFixed(3)).join(',');
if (process.env.DUMP) {
  const w = detail.filter(x => !process.env.GRP || x[0] === process.env.GRP), stp = Math.max(1, Math.floor(w.length / +process.env.DUMP));
  for (let i = 0; i < w.length; i += stp) { const x = w[i]; console.log('   ', x[0], x[1], 'rows', x[2].row, x[3].row, 'h', layers[x[2].row].h.toFixed(3), layers[x[3].row].h.toFixed(3), 'A', fmtP(x[2].pos), 'B', fmtP(x[3].pos), x[4] !== undefined ? `d_mesh ${x[4].toFixed(4)} d_vert ${x[5].toFixed(4)} at ${fmtP(x[6])}` : ''); }
}
