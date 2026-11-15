// 曲面网格封闭性测试. 原则: 网格的开口边 (只属于一个三角形的边) 必须落在特征线上 (入射掠射 / 出射临界 / 外缘 …),
//   或在最大距离截断处、灯 / 眼的捏合点旁. 其余的开口边都是缺陷:
//     layer —— 沿一条层线的开口边: 这条层线的一侧没有三角形 (极值帽: 边界线在最末一层之外鼓出去的一段; 或面中间的洞);
//     cross —— 跨两层的开口边, 不在特征线上: 楔形缺口 / 洞的侧边.
//   同一位置的顶点先焊接. 输出两类缺陷的条数与总长, 及最长的几个连通分量的位置 (三维坐标 + 从眼看的方位 / 仰角).
// node test/mesh2/open.mjs code mode el D surf fix      env: STITCH, FEAT=0, N=明细条数 (默认 6)
//   读 test/mesh2/cache/ 里 edges.mjs 存下的场景 (没有就先跑 edges.mjs)
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
const HERE = path.dirname(fileURLToPath(import.meta.url));
const B = process.env.ENG || new URL('../../src/engine/', import.meta.url).href;
const ST = await import(process.env.STITCH || (B + 'stitch2.js'));
const [code = '4-1-8', mode = 'lowitz', el = '3.5', Ds = '15', surfName = 'caustic', fixAng = '0'] = process.argv.slice(2);
const D = +Ds, MAXD = 80;
const cname = `${process.env.CACHE || HERE + '/cache'}/${[code, mode, el, Ds, surfName, fixAng].join('_')}.json`;
if (!fs.existsSync(cname)) { console.log('no cache: run edges.mjs first', cname); process.exit(2); }
const dec = s => { const b = Buffer.from(s, 'base64'); return new Float32Array(b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength)); };
const j = JSON.parse(fs.readFileSync(cname, 'utf8'));
const feats = j.feats.map(f => ({ g: f.g, p: dec(f.p) }));
let layers = j.layers.map(l => ({ h: l.h, curves: l.curves.map(c => ({ p: dec(c.p), d: dec(c.d), x: dec(c.x), b: c.b })) }));
const L = [D * Math.cos(el * Math.PI / 180), 0, D * Math.sin(el * Math.PI / 180)];
const useFeat = process.env.FEAT !== '0';
let m = ST.stitchLayers2(layers, L, D, useFeat ? feats.map(f => ({ p: f.p, g: f.g })) : null, { maxD: MAXD });
if (m.bad.length && process.env.GROW) {       // 页面那一轮补层 (要重算层, 慢; 默认只缝缓存里的层)
  const CR = await import(B + 'crystal.js'), LPm = await import(B + 'lp.js'), EP = await import(B + 'epdiv.js');
  CR.setApex(56, 56); LPm.LP.partial = true; LPm.LP.upperApexDeg = 56; LPm.LP.lowerApexDeg = 56; LPm.setIndex(1.31);
  const lp = CR.parsePath(code).map(CR.lpFaceFromCompact), surf = EP.resolveSurface(lp, surfName), low = mode === 'lowitz', z = low ? 0 : +fixAng, th = low ? +fixAng : 0;
  for (const [a, b] of m.bad) { const h = 0.5 * (a + b), c = EP.lightCone2DofLayers(lp, +el, D, 0, 1, 101, MAXD, low, z, th, surf, 0, 1, h, true)[0].curves; if (c.length) layers.push({ h, curves: c }); }
  layers.sort((a, b) => a.h - b.h);
  m = ST.stitchLayers2(layers, L, D, useFeat ? feats.map(f => ({ p: f.p, g: f.g })) : null, { maxD: MAXD });
}
const P = m.pos, H = m.h, nT = P.length / 9;

// ---- 焊接顶点 + 边的三角形计数
const vkey = (x, y, z) => Math.round(x * 2e4) + ',' + Math.round(y * 2e4) + ',' + Math.round(z * 2e4);
const vid = new Map(), vpos = [], vh = [], ecount = new Map();
const idOf = (o, t, k) => { const q = vkey(P[o], P[o + 1], P[o + 2]); let i = vid.get(q); if (i === undefined) { vid.set(q, i = vid.size); vpos.push([P[o], P[o + 1], P[o + 2]]); vh.push(H[3 * t + k]); } return i; };
for (let t = 0; t < nT; t++) {
  const a = idOf(9 * t, t, 0), b = idOf(9 * t + 3, t, 1), c = idOf(9 * t + 6, t, 2);
  if (a === b || b === c || c === a) continue;
  for (const [u, v] of [[a, b], [b, c], [c, a]]) { const k = u < v ? u * 4194304 + v : v * 4194304 + u; ecount.set(k, (ecount.get(k) || 0) + 1); }
}

// ---- 特征线分桶 (同 edges.mjs)
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
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

// ---- 开口边分类
const open = [];
const cnt = { feat: 0, trunc: 0, pinch: 0, layer: 0, cross: 0 }, len = { feat: 0, trunc: 0, pinch: 0, layer: 0, cross: 0 };
const onTrunc = q => Math.max(rl(q), Math.hypot(...q)) > 0.97 * MAXD;
const atPinch = q => rl(q) < 0.02 * D || Math.hypot(...q) < 0.02 * D;
for (const [k, c] of ecount) {
  if (c !== 1) continue;
  const u = Math.floor(k / 4194304), v = k % 4194304, a = vpos[u], b = vpos[v], l = Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
  const mid = [0, 1, 2].map(i => 0.5 * (a[i] + b[i]));
  let kind;
  const fa = nearF(a, 2 * tolAt(rl(a))), fb = nearF(b, 2 * tolAt(rl(b)));
  let onF = false;
  // 两端在同一条特征线上, 沿线相距与边长相当 (不是恰好跨在特征线两处之间的弦)
  if (fa && fb) for (const [fi, w] of fa) { const w2 = fb.get(fi); if (w2 && Math.abs(w2.t - w.t) < 1.5 * l + 4 * tolAt(rl(a))) { onF = true; break; } }
  if (onF) kind = 'feat';
  else if (atPinch(a) || atPinch(b)) kind = 'pinch';
  else if (onTrunc(a) && onTrunc(b)) kind = 'trunc';
  else if (vh[u] === vh[v]) kind = 'layer';
  else kind = 'cross';
  cnt[kind]++; len[kind] += l;
  if (kind === 'layer' || kind === 'cross') open.push({ u, v, l, kind, mid });
}
// 缺陷开口边的连通分量
const adj = new Map();
for (const e of open) { (adj.get(e.u) || adj.set(e.u, []).get(e.u)).push(e); (adj.get(e.v) || adj.set(e.v, []).get(e.v)).push(e); }
const seen = new Set(), comps = [];
for (const [s0] of adj) {
  if (seen.has(s0)) continue;
  const q = [s0]; seen.add(s0); const es = new Set();
  for (let i = 0; i < q.length; i++) for (const e of adj.get(q[i])) { es.add(e); const w = e.u === q[i] ? e.v : e.u; if (!seen.has(w)) { seen.add(w); q.push(w); } }
  let tl = 0, big = null, nl = 0, nc = 0; const hs = new Set();
  for (const e of es) { tl += e.l; if (!big || e.l > big.l) big = e; if (e.kind === 'layer') nl++; else nc++; hs.add(vh[e.u]); hs.add(vh[e.v]); }
  comps.push({ tl, n: es.size, nl, nc, at: big.mid, hmin: Math.min(...hs), hmax: Math.max(...hs) });
}
comps.sort((a, b) => b.tl - a.tl);
const sky = q => { const r = Math.hypot(...q); return `az ${(Math.atan2(q[1], q[0]) * 180 / Math.PI).toFixed(1)} el ${(Math.asin(q[2] / r) * 180 / Math.PI).toFixed(1)}`; };
console.log(`${code} ${mode} Σ${el} D${D} ${surfName} fix${fixAng}: layers ${layers.length} tris ${m.nTri} open ${cnt.feat + cnt.trunc + cnt.pinch + cnt.layer + cnt.cross}`);
console.log(`  ok   feat ${cnt.feat} (${len.feat.toFixed(1)} m)  trunc ${cnt.trunc} (${len.trunc.toFixed(1)} m)  pinch ${cnt.pinch}`);
console.log(`  BAD  layer ${cnt.layer} (${len.layer.toFixed(2)} m)  cross ${cnt.cross} (${len.cross.toFixed(2)} m)  comps ${comps.length}`);
for (const c of comps.slice(0, +(process.env.N || 6))) console.log(`    comp ${c.tl.toFixed(3)} m, ${c.n} edges (layer ${c.nl} cross ${c.nc}) h ${c.hmin.toFixed(2)}..${c.hmax.toFixed(2)} at ${c.at.map(x => x.toFixed(3)).join(',')}  ${sky(c.at)}  r_lamp ${rl(c.at).toFixed(2)}`);
