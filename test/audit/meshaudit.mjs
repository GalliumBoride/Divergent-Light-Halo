// 贴边审计 —— 球壳法网格的整体验收, 走页面的整条管线 (worker.js 原样: 边界线 → 按半径分段 → 各段球壳法网格 → 拼起来).
//   node test/audit/meshaudit.mjs code mode el D [surf] [fix] [NW=12]
//   mode: plate | parry | custom (一自由度, fix = "ζ,θ") | column | lowitz (二自由度, fix = ζ 或 θ);  surf: caustic / entry / exit / deg1 / blue1 …
//   env: SRC = src 目录的 file URL (默认本仓库的 src; 可指向另一份源码做对比), Q = 精度 0/1/2, MAXD = 截断距离,
//        TOP = 列出最大的几处缺口, BLTOP = 列出最长的几段没贴上的边界线, DUMP = 存缺口边 (json), VERBOSE = 逐段打印
//
// 解决什么: 回归测试 (test/mesh1, test/mesh2) 只测单段网格的覆盖, 测不到页面实际用的分段管线, 也不看网格的边是不是贴在边界上.
//   这里把各段的三角形按顶点位置 (逐位相同) 焊在一起, 找开口边 (只属于一个三角形的边), 逐条判:
//     在引擎直接算的精确边界线上 / 截断球面上 / 最里最外一层 / 灯旁 → 合法;
//     其余 = 缺口: 面里的洞、两段之间没对上的缝 (两端都在段的交界半径上: 段间缝)、没贴到边界线的边.
//   反向: 精确边界线上的点离网格多远 (边界线没贴上的长度).
//   降自由度的光路 (没有边界线可比) 看网格器自己的自检 (stat.openDeg: 不在末端点 / 截断点 / 边界线点之间的开口边).
// 输出: 一行摘要 + 最大的几处缺口 + 一行 "AUDIT {json}" (test/audit/run.mjs 批量跑时读它)
import fs from 'node:fs';
const SRC = process.env.SRC || new URL('../../src/', import.meta.url).href;
const [code, mode, els, Ds, surfArg = 'caustic', fixArg = '0', NWs = '12'] = process.argv.slice(2);
globalThis.onmessage = null;
await import(SRC + 'worker.js');
let response = null;
globalThis.postMessage = m => { response = m; };
const run = (kind, args) => { globalThis.onmessage({ data: { kind, args } }); const r = response; response = null; if (!r.ok) throw new Error(kind + ': ' + r.error); return r.data; };
const CR = await import(SRC + 'engine/crystal.js'), MS = await import(SRC + 'engine/mesh2s.js');
const { shellChunks } = await import(SRC + 'app/shellparts.js');

const DEG = Math.PI / 180, el = +els, dist = +Ds, maxD = +(process.env.MAXD || 80), q = +(process.env.Q || 0), NW = +NWs;
const er = el * DEG, L = [dist * Math.cos(er), 0, dist * Math.sin(er)], E = [0, 0, 0], D = dist;
const dof = mode === 'column' || mode === 'lowitz' ? 2 : 1;
let zeta = 0, theta = 0;
if (mode === 'parry') zeta = -90;
else if (mode === 'custom') [zeta, theta] = fixArg.split(',').map(Number);
else if (mode === 'column') zeta = +fixArg;
else if (mode === 'lowitz') theta = +fixArg;
const fp = CR.parsePath(code);
const base = { n: 1.31, apexUp: 56, apexLo: 56, lang: 'zh', fp, el, dist, eye: 0, maxD, zeta, theta, lowitz: mode === 'lowitz' };
const bud = [60000, 90000, 140000][q];
const t00 = Date.now();
let bl = [], chunks, kindM, margs, note = '';
if (dof === 2) {
  const nG = [101, 141, 181][q];
  const e = run('edges2', { ...base, surf: surfArg, nG });
  const rims = run('rims2', { ...base, surf: surfArg, nG }) || [];
  for (const k of ['entry', 'tir', 'seam', 'turn', 'cusp']) for (const c of e[k] || []) { const p = c.p || c; if (p.length >= 6) bl.push({ p, kind: k }); }
  for (const c of rims) { const p = c.p || c; if (p.length >= 6) bl.push({ p, kind: 'rim' }); }
  if (e.rank1) note = ' [降自由度]';
  const crit = MS.criticalRadii(bl.map(l => l.p), L, E);
  let rHi = 0;
  for (const l of bl) for (let i = 0; i + 2 < l.p.length; i += 3) rHi = Math.max(rHi, Math.hypot(l.p[i], l.p[i + 1], l.p[i + 2]));
  chunks = shellChunks(L, E, crit, NW, rHi, maxD);
  kindM = 'mesh2s';
  const sf = surfArg === 'caustic' ? null : surfArg.replace('deg', 'seam');
  margs = { ...base, L, E, surf: surfArg, surfFence: sf, budgetAll: bud, maxMs: 60000, lines: bl };
} else {
  const b = run('bdry1', base);
  for (const [k, arr] of [['entry', b.entry], ['mirror', b.mirror], ['tir', b.tir]]) for (const c of arr || []) { const p = c.p || c; if (p.length >= 6) bl.push({ p, kind: k }); }
  chunks = shellChunks(L, E, [], NW, 0, maxD);
  kindM = 'mesh1s';
  margs = { ...base, budgetAll: bud, maxMs: 30000 };
}
const tEdges = Date.now() - t00;

// ---- 各段网格
const parts = [];
let nTri = 0, tMax = 0, tSum = 0, skipped = 0, loose = 0, failed = 0, selfOpen = 0, selfFilled = 0;
for (const [ci, c] of chunks.entries()) {
  const t = Date.now();
  let m;
  try { m = run(kindM, { ...margs, rMin: c[0], rMax: c[1], rows: c[2], budget: Math.max(8000, Math.round(2 * bud * c[3])) }); }
  catch (err) { m = { ok: false, why: String(err.message).split('\n')[0] }; }
  const ms = Date.now() - t;
  tMax = Math.max(tMax, ms); tSum += ms;
  if (m.ok) { nTri += m.nTri; skipped += m.skipped || 0; loose += m.loose || 0; selfOpen += m.open || 0; selfFilled += m.filled || 0; parts.push({ ci, c, m }); }
  else if (m.why !== 'empty') failed++;
  if (process.env.VERBOSE) console.log(`  段 ${ci} [${c[0].toFixed(4)}, ${c[1].toFixed(4)}]: ${m.ok ? m.nTri + ' 三角形, ' + (m.rows || 0) + ' 层' : m.why}, ${ms} ms` + (m.ok ? `, 没缝上 ${(m.loose || 0).toFixed(1)}°, 丢超长边 ${m.skipped === undefined ? '?' : m.skipped}` : ''));
}

// ---- 焊接 (逐位相同的顶点) → 开口边
const vid = new Map(), VP = [];
const idOf = (x, y, z) => { const k = x + ',' + y + ',' + z; let i = vid.get(k); if (i === undefined) { i = VP.length / 3; vid.set(k, i); VP.push(x, y, z); } return i; };
const edge = new Map();                          // "a,b" (a<b) → 次数
const T = [];
for (const { m } of parts) {
  const p = m.pos;
  for (let t = 0; t < p.length; t += 9) {
    const a = idOf(p[t], p[t + 1], p[t + 2]), b = idOf(p[t + 3], p[t + 4], p[t + 5]), c = idOf(p[t + 6], p[t + 7], p[t + 8]);
    if (a === b || b === c || c === a) continue;
    T.push(a, b, c);
    for (const [u, v] of [[a, b], [b, c], [c, a]]) { const k = u < v ? u * 4194304 + v : v * 4194304 + u; edge.set(k, (edge.get(k) || 0) + 1); }
  }
}
const P3 = i => [VP[3 * i], VP[3 * i + 1], VP[3 * i + 2]];
const rOf = p => Math.hypot(p[0] - E[0], p[1] - E[1], p[2] - E[2]);
// 边界线: 线段的哈希格 (每段登记到它包围盒碰到的格子), 点到线段的精确距离
const CELL = 0.25, LH = new Map(), SEG = [];
for (const ln of bl) {
  const p = ln.p;
  for (let i = 0; i + 5 < p.length; i += 3) {
    const si = SEG.length / 6;
    SEG.push(p[i], p[i + 1], p[i + 2], p[i + 3], p[i + 4], p[i + 5]);
    const lo = [0, 1, 2].map(k => Math.floor(Math.min(p[i + k], p[i + 3 + k]) / CELL)), hi = [0, 1, 2].map(k => Math.floor(Math.max(p[i + k], p[i + 3 + k]) / CELL));
    if ((hi[0] - lo[0] + 1) * (hi[1] - lo[1] + 1) * (hi[2] - lo[2] + 1) > 4096) continue;
    for (let x = lo[0]; x <= hi[0]; x++) for (let y = lo[1]; y <= hi[1]; y++) for (let z = lo[2]; z <= hi[2]; z++) { const k = x + ',' + y + ',' + z; let l = LH.get(k); if (!l) LH.set(k, l = []); l.push(si); }
  }
}
const lineDist = p => {                          // 到边界线的最近距离, 只看周围 27 格 (≥ CELL 的算远)
  const cx = Math.floor(p[0] / CELL), cy = Math.floor(p[1] / CELL), cz = Math.floor(p[2] / CELL);
  let best = Infinity;
  for (let x = -1; x <= 1; x++) for (let y = -1; y <= 1; y++) for (let z = -1; z <= 1; z++) {
    const l = LH.get((cx + x) + ',' + (cy + y) + ',' + (cz + z));
    if (!l) continue;
    for (const si of l) {
      const ax = SEG[6 * si], ay = SEG[6 * si + 1], az = SEG[6 * si + 2], ex = SEG[6 * si + 3] - ax, ey = SEG[6 * si + 4] - ay, ez = SEG[6 * si + 5] - az;
      const qx = p[0] - ax, qy = p[1] - ay, qz = p[2] - az, l2 = ex * ex + ey * ey + ez * ez;
      let u = l2 > 0 ? (qx * ex + qy * ey + qz * ez) / l2 : 0; u = u < 0 ? 0 : u > 1 ? 1 : u;
      const d = Math.hypot(qx - u * ex, qy - u * ey, qz - u * ez);
      if (d < best) best = d;
    }
  }
  return best;
};
const tolB = r => Math.min(0.01 + 2e-3 * r, 0.02 * r) + 1e-6;      // 贴边容差 (三维), 近眼按角度 (~1°)
const rMinAll = chunks.length ? chunks[0][0] : 0, rMaxAll = chunks.length ? chunks[chunks.length - 1][1] : 0;
const cuts = new Set(); for (const c of chunks) { cuts.add(c[0]); cuts.add(c[1]); }
const nearCut = r => { for (const c of cuts) if (Math.abs(r / c - 1) < 2e-6) return c; return 0; };
const stats = { open: 0, line: 0, cut: 0, rim: 0, lamp: 0, defect: 0, seam: 0 };
const defects = [];
let defLen = 0, defSky = 0, seamLen = 0, seamSky = 0;
const ang = (a, b) => { const la = Math.hypot(...a), lb = Math.hypot(...b); return Math.acos(Math.max(-1, Math.min(1, (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) / (la * lb || 1)))); };
for (const [k, cnt] of edge) {
  if (cnt !== 1) continue;
  stats.open++;
  const a = Math.floor(k / 4194304), b = k % 4194304, pa = P3(a), pb = P3(b), pm = [0.5 * (pa[0] + pb[0]), 0.5 * (pa[1] + pb[1]), 0.5 * (pa[2] + pb[2])];
  const ra = rOf(pa), rb = rOf(pb), rm = rOf(pm);
  const dla = Math.hypot(pa[0] - L[0], pa[1] - L[1], pa[2] - L[2]), dlb = Math.hypot(pb[0] - L[0], pb[1] - L[1], pb[2] - L[2]);
  if (dla >= maxD * (1 - 2e-3) && dlb >= maxD * (1 - 2e-3)) { stats.cut++; continue; }
  if ((ra <= rMinAll * (1 + 1e-5) && rb <= rMinAll * (1 + 1e-5)) || (ra >= rMaxAll * (1 - 1e-5) && rb >= rMaxAll * (1 - 1e-5))) { stats.rim++; continue; }
  if (Math.max(dla, dlb) < 3e-3 * D) { stats.lamp++; continue; }
  if (lineDist(pa) <= tolB(ra) && lineDist(pb) <= tolB(rb) && lineDist(pm) <= tolB(rm)) { stats.line++; continue; }
  const len = Math.hypot(pa[0] - pb[0], pa[1] - pb[1], pa[2] - pb[2]), sky = ang(pa, pb) / DEG;
  const sc = nearCut(ra) && nearCut(ra) === nearCut(rb);
  if (sc) { stats.seam++; seamLen += len; seamSky += sky; } else { stats.defect++; defLen += len; defSky += sky; }
  defects.push({ a, b, len, sky, seam: !!sc, r: rm, ld: Math.min(lineDist(pa), lineDist(pb)) });
}
// 缺口聚成团 (共顶点的连在一起), 列最大的几团
const comp = new Map();
const find = x => { let r = x; while (comp.get(r) !== r) r = comp.get(r); while (comp.get(x) !== r) { const n = comp.get(x); comp.set(x, r); x = n; } return r; };
for (const d of defects) { for (const v of [d.a, d.b]) if (!comp.has(v)) comp.set(v, v); comp.set(find(d.a), find(d.b)); }
const groups = new Map();
for (const d of defects) { const g = find(d.a); let o = groups.get(g); if (!o) groups.set(g, o = { len: 0, sky: 0, n: 0, seam: 0, c: [0, 0, 0], rlo: Infinity, rhi: 0, ld: 0 }); o.len += d.len; o.sky += d.sky; o.n++; if (d.seam) o.seam++; const pa = P3(d.a); o.c[0] += pa[0]; o.c[1] += pa[1]; o.c[2] += pa[2]; o.rlo = Math.min(o.rlo, d.r); o.rhi = Math.max(o.rhi, d.r); o.ld = Math.max(o.ld, d.ld); }
const top = [...groups.values()].sort((x, y) => y.sky - x.sky);

// ---- 反向: 边界线 → 网格 (三角形的哈希格)
const C3 = 0.5, H3 = new Map();
for (let t = 0; t < T.length; t += 3) {
  const A = P3(T[t]), B = P3(T[t + 1]), C = P3(T[t + 2]);
  const lo = [0, 1, 2].map(k => Math.floor(Math.min(A[k], B[k], C[k]) / C3)), hi = [0, 1, 2].map(k => Math.floor(Math.max(A[k], B[k], C[k]) / C3));
  if ((hi[0] - lo[0] + 1) * (hi[1] - lo[1] + 1) * (hi[2] - lo[2] + 1) > 3000) continue;
  for (let x = lo[0]; x <= hi[0]; x++) for (let y = lo[1]; y <= hi[1]; y++) for (let z = lo[2]; z <= hi[2]; z++) { const k = x + ',' + y + ',' + z; let l = H3.get(k); if (!l) H3.set(k, l = []); l.push(t); }
}
const ptTri = (p, a, b, c) => {
  const sub = (u, v) => [u[0] - v[0], u[1] - v[1], u[2] - v[2]], dot = (u, v) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
  const ab = sub(b, a), ac = sub(c, a), ap = sub(p, a);
  const d1 = dot(ab, ap), d2 = dot(ac, ap); if (d1 <= 0 && d2 <= 0) return Math.hypot(...ap);
  const bp = sub(p, b), d3 = dot(ab, bp), d4 = dot(ac, bp); if (d3 >= 0 && d4 <= d3) return Math.hypot(...bp);
  const vc = d1 * d4 - d3 * d2; if (vc <= 0 && d1 >= 0 && d3 <= 0) { const v = d1 / (d1 - d3); return Math.hypot(ap[0] - v * ab[0], ap[1] - v * ab[1], ap[2] - v * ab[2]); }
  const cp = sub(p, c), d5 = dot(ab, cp), d6 = dot(ac, cp); if (d6 >= 0 && d5 <= d6) return Math.hypot(...cp);
  const vb = d5 * d2 - d1 * d6; if (vb <= 0 && d2 >= 0 && d6 <= 0) { const w = d2 / (d2 - d6); return Math.hypot(ap[0] - w * ac[0], ap[1] - w * ac[1], ap[2] - w * ac[2]); }
  const va = d3 * d6 - d5 * d4; if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) { const w = (d4 - d3) / ((d4 - d3) + (d5 - d6)); return Math.hypot(p[0] - (b[0] + w * (c[0] - b[0])), p[1] - (b[1] + w * (c[1] - b[1])), p[2] - (b[2] + w * (c[2] - b[2]))); }
  const den = 1 / (va + vb + vc), v = vb * den, w = vc * den;
  return Math.hypot(ap[0] - ab[0] * v - ac[0] * w, ap[1] - ab[1] * v - ac[1] * w, ap[2] - ab[2] * v - ac[2] * w);
};
const meshDist = (p, lim) => {
  const c = p.map(x => Math.floor(x / C3)); let best = Infinity;
  for (let x = -1; x <= 1; x++) for (let y = -1; y <= 1; y++) for (let z = -1; z <= 1; z++) for (const t of H3.get((c[0] + x) + ',' + (c[1] + y) + ',' + (c[2] + z)) || []) {
    const d = ptTri(p, P3(T[t]), P3(T[t + 1]), P3(T[t + 2])); if (d < best) { best = d; if (best <= lim) return best; }
  }
  return best;
};
let blLen = 0, blMiss = 0, blSky = 0, blMissSky = 0, worstL = null;
const byKind = {};
const MISS = [];           // BLTOP: 没贴上的线段连成的段
// 网格给出的特殊线 (焦线 / 折叠线): 边界引擎的线与焦线重合时 (一整条取向曲线都把光射向同一点, 那里没有三角形) 也算贴上
const SPEC = [];
for (const { m } of parts) for (const s of m.special || []) SPEC.push(s.p);
const specDist = (q, lim) => {
  let bd = Infinity;
  for (const p of SPEC) for (let i = 0; i + 5 < p.length; i += 3) {
    const ex = p[i + 3] - p[i], ey = p[i + 4] - p[i + 1], ez = p[i + 5] - p[i + 2], l2 = ex * ex + ey * ey + ez * ez;
    let t = l2 > 0 ? ((q[0] - p[i]) * ex + (q[1] - p[i + 1]) * ey + (q[2] - p[i + 2]) * ez) / l2 : 0;
    t = t < 0 ? 0 : t > 1 ? 1 : t;
    const d = Math.hypot(q[0] - p[i] - t * ex, q[1] - p[i + 1] - t * ey, q[2] - p[i + 2] - t * ez);
    if (d < bd) { bd = d; if (bd <= lim) return bd; }
  }
  return bd;
};
for (const ln of bl) {
  const p = ln.p;
  for (let i = 0; i + 5 < p.length; i += 3) {
    const a = [p[i], p[i + 1], p[i + 2]], b = [p[i + 3], p[i + 4], p[i + 5]], m = [0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2])];
    const r = rOf(m), dl = Math.hypot(m[0] - L[0], m[1] - L[1], m[2] - L[2]);
    if (r < rMinAll * 1.5 || r > rMaxAll || dl > maxD * (1 - 2e-3) || dl < 3e-3 * D) continue;
    const len = Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]), sky = ang(a, b) / DEG;
    blLen += len; blSky += sky;
    const o = byKind[ln.kind] || (byKind[ln.kind] = { len: 0, miss: 0 }); o.len += len;
    let d = meshDist(m, tolB(r));
    if (d > tolB(r) && SPEC.length) d = Math.min(d, specDist(m, tolB(r)));
    if (d > tolB(r) && process.env.BLTOP) { const last = MISS[MISS.length - 1]; if (last && last.ln === ln && last.i1 === i - 3) { last.i1 = i; last.len += len; last.rlo = Math.min(last.rlo, r); last.rhi = Math.max(last.rhi, r); last.d = Math.max(last.d, d); last.b = m; } else MISS.push({ ln, i0: i, i1: i, len, rlo: r, rhi: r, d, a: m, b: m }); }
    if (d > tolB(r)) { blMiss += len; blMissSky += sky; o.miss += len; if (!worstL || sky > worstL.sky) worstL = { sky, d, r, kind: ln.kind, az: Math.atan2(m[1], m[0]) / DEG, el: Math.asin(m[2] / r) / DEG }; }
  }
}

const sceneName = [code, mode, els, Ds, dof === 2 ? surfArg : '', fixArg].filter(s => s !== '').join(' ');
console.log(`${sceneName}${note}: ${nTri} 三角形 (预算 ${bud}), 边界 ${tEdges} ms, 网格合计 ${tSum} ms (最慢一段 ${tMax} ms)${failed ? ', 失败 ' + failed + ' 段' : ''}, 丢超长边 ${skipped}, 没缝上 ${loose.toFixed(1)}°, 自检: 开口 ${selfOpen.toFixed(1)}° 补洞 ${selfFilled}`);
console.log(`  开口边 ${stats.open}: 边界线上 ${stats.line}, 截断面 ${stats.cut}, 最里/最外层 ${stats.rim}, 灯旁 ${stats.lamp}; 缺口 ${stats.defect} 条 ${defLen.toFixed(2)} m / 天空 ${defSky.toFixed(1)}°, 段间缝 ${stats.seam} 条 ${seamLen.toFixed(2)} m / 天空 ${seamSky.toFixed(1)}°`);
console.log(`  边界线 ${blLen.toFixed(1)} m / 天空 ${blSky.toFixed(0)}°: 没贴上网格 ${blMiss.toFixed(2)} m / ${blMissSky.toFixed(1)}°` + (worstL ? ` (最长一段 ${worstL.kind} 离网格 ${worstL.d === Infinity ? '>0.5' : worstL.d.toFixed(3)} m, r ${worstL.r.toFixed(2)}, 天空 ${worstL.az.toFixed(1)}/${worstL.el.toFixed(1)})` : '') + ' | ' + Object.entries(byKind).map(([k, o]) => `${k} ${o.miss.toFixed(1)}/${o.len.toFixed(1)}`).join(' '));
const fsky = q => { const r = Math.hypot(...q); return `${(Math.atan2(q[1], q[0]) / DEG).toFixed(1)}/${(Math.asin(q[2] / r) / DEG).toFixed(1)}`; };
for (const u of MISS.sort((x, y) => y.len - x.len).slice(0, +(process.env.BLTOP || 0))) console.log(`    没贴上的边界线: ${u.ln.kind} ${u.len.toFixed(2)} m, r ${u.rlo.toFixed(3)}–${u.rhi.toFixed(3)}, 天空 ${fsky(u.a)} → ${fsky(u.b)}, 离网格最远 ${u.d === Infinity ? '>0.5' : u.d.toFixed(3)} m`);
for (const g of top.slice(0, +(process.env.TOP || 0))) {
  const c = g.c.map(x => x / g.n), r = Math.hypot(...c);
  console.log(`    缺口团: ${g.n} 条边 ${g.len.toFixed(3)} m / 天空 ${g.sky.toFixed(2)}°${g.seam ? ' (段间缝 ' + g.seam + ' 条)' : ''}, r ${g.rlo.toFixed(4)}–${g.rhi.toFixed(4)}, 天空 ${(Math.atan2(c[1], c[0]) / DEG).toFixed(2)}/${(Math.asin(c[2] / r) / DEG).toFixed(2)}, 离边界线 ≥ ${g.ld === Infinity ? '>0.25' : g.ld.toFixed(3)} m`);
}
if (process.env.DUMP) fs.writeFileSync(process.env.DUMP, JSON.stringify(defects.map(d => ({ a: P3(d.a), b: P3(d.b), seam: d.seam }))));
if (process.env.SAVE) { const fd = fs.openSync(process.env.SAVE, 'w'); for (const { m } of parts) fs.writeSync(fd, Buffer.from(m.dir.buffer, m.dir.byteOffset, m.dir.byteLength)); fs.closeSync(fd); }
console.log(`AUDIT ${JSON.stringify({ scene: sceneName, rank1: note !== '', nLines: bl.length, nTri, bud, tEdges, tSum, tMax, failed, skipped, selfOpen: +selfOpen.toFixed(2), selfFilled, defN: stats.defect, defLen: +defLen.toFixed(3), defSky: +defSky.toFixed(2), seamN: stats.seam, seamSky: +seamSky.toFixed(2), blMiss: +blMiss.toFixed(3), blMissSky: +blMissSky.toFixed(2), blLen: +blLen.toFixed(1) })}`);
