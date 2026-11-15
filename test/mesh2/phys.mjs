// 二自由度网格的几何校验: 网格与真实晕面互相对照, 三维里量一遍, 观测者天空里再量一遍.
//   node test/mesh2/phys.mjs code mode el D surf fix      env: NT = 各查多少个三角形 (默认 500), NL = 另算多少层 (默认 40), TOP = 明细条数, STITCH
// 晕面的定义与网格怎么缝无关: 晕面 = 各 h 的层线 (灯处仰角 h 的等值线) 的并集, 边上是精确边界线. 所以
//   (1) 网格 → 晕面 (有没有连错 / 补错的三角形 —— "面里必须有层线"):
//       三维: 取面积最大的 NT/2 个三角形 + 等距抽 NT/2 个, 重心 X 的 h 由位置直接算出 (h = asin((L_z − X_z)/|X − L|)),
//             用引擎另算这一层的层线, 量 X 到层线的三维距离. 面与光锥相切处层线随 h 挪得快, 超差的再算 h ± 0.02°、± 0.06° 取最小.
//       天空: 取天空里面积 (立体角) 最大的 NT/2 个 + 等距抽 NT/2 个, 三个顶点视方向的平均方向 d, h 在三个顶点的 h 之间取几个值,
//             量 d 到这些层线 (各点的视方向连成的线) 的夹角, 取最小. 贴着眼睛的三角形三维里只有几厘米、天空里却可以横跨几十度,
//             三维那一遍看不见, 要靠这一遍.
//       真在面上的三角形只差弦高; 横跨空洞、把两张面连在一起的, 三维里差分米到米, 天空里差几度到几十度.
//   (2) 晕面 → 网格 (有没有漏的面 —— 层线和边界线都得被网格盖到):
//       在网格的层之间另算 NL 层 (相邻两层 h 的中点, 沿层序等距取), 层线上每个点 + 每条精确边界线 / 脊线上每个点:
//       三维里量到网格 (点到三角形) 的距离; 天空里量该点的视方向到网格 (按顶点视方向画的三角形) 的夹角.
//   容差: 三维 3 cm + 0.3% r 记一档, 10 cm + 1% r 记一档 (r = 到灯、眼较近者的距离), 灯 / 眼 0.02 D 之内不算 (捏合点, 交给天空那一遍);
//         天空 0.3° 记一档, 1° 记一档.
//   读 test/mesh2/cache/ 里 edges.mjs 存下的场景 (没有就先跑 edges.mjs); 网格按页面的流程缝 (含补层那一轮).
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
const HERE = path.dirname(fileURLToPath(import.meta.url));
const B = process.env.ENG || new URL('../../src/engine/', import.meta.url).href;
const CR = await import(B + 'crystal.js'), LPm = await import(B + 'lp.js'), EP = await import(B + 'epdiv.js'), ST = await import(process.env.STITCH || (B + 'stitch2.js'));
const PG = await import('./pagehs.mjs');
const MAXD = 80, DEG = Math.PI / 180;

// mesh: 不给 → 按页面流程用光锥层缝 (stitch2); 给 { pos, dir } → 校验这张网格 (levels: 另算层线取哪些 h, 不给就在网格的 h 范围里等距取 NL 个)
export function physCheck2(code, mode, el, Ds, surfName, fixAng, NT = 500, NL = 40, mesh = null, levels = null) {
  const D = +Ds, cname = `${process.env.CACHE || HERE + '/cache'}/${[code, mode, el, Ds, surfName, fixAng].join('_')}.json`;
  if (!fs.existsSync(cname)) return null;
  CR.setApex(56, 56); LPm.LP.partial = true; LPm.LP.upperApexDeg = 56; LPm.LP.lowerApexDeg = 56; LPm.setIndex(1.31);
  const lp = CR.parsePath(code).map(CR.lpFaceFromCompact), surf = EP.resolveSurface(lp, surfName), low = mode === 'lowitz', z = low ? 0 : +fixAng, th = low ? +fixAng : 0;
  const dec = s => { const b = Buffer.from(s, 'base64'); return new Float32Array(b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength)); };
  const j = JSON.parse(fs.readFileSync(cname, 'utf8'));
  const feats = j.feats.map(f => ({ g: f.g, p: dec(f.p) }));
  const layers0 = j.layers.map(l => ({ h: l.h, curves: l.curves.map(c => ({ p: dec(c.p), d: dec(c.d), x: dec(c.x), b: c.b })) }));
  const er = +el * DEG, L = [D * Math.cos(er), 0, D * Math.sin(er)];
  let m = mesh, layers = null;
  if (!m) ({ m, layers } = PG.pageMesh(EP, ST, lp, +el, D, low, z, th, surf, L, layers0, feats.map(f => ({ p: f.p, g: f.g }))));
  const P = m.pos, Dr = m.dir, nT = P.length / 9;
  const hAt = o => Math.asin((L[2] - P[o + 2]) / (Math.hypot(P[o] - L[0], P[o + 1] - L[1], P[o + 2] - L[2]) || 1)) / DEG;
  const mh = new Float32Array(3 * nT);
  for (let i = 0; i < 3 * nT; i++) mh[i] = hAt(3 * i);
  const lcache = new Map();
  const layerAt = h => { const k = Math.round(h * 1e5); let c = lcache.get(k); if (!c) { c = EP.lightCone2DofLayers(lp, +el, D, 0, 1, 101, MAXD, low, z, th, surf, 0, 1, h, true)[0].curves; lcache.set(k, c); if (lcache.size > 64) lcache.delete(lcache.keys().next().value); } return c; };
  const rOf = X => Math.min(Math.hypot(X[0] - L[0], X[1] - L[1], X[2] - L[2]), Math.hypot(X[0], X[1], X[2]));
  // 点到折线组的最近距离 (q = 各曲线的坐标数组: 三维位置, 或单位视方向 —— 后者的弦长当夹角用)
  const dPoly = (X, cs, key) => {
    let bd = Infinity;
    for (const c of cs) { const q = c[key]; for (let k = 0; k + 1 < q.length / 3; k++) {
      const ax = q[3 * k], ay = q[3 * k + 1], az = q[3 * k + 2], bx = q[3 * k + 3] - ax, by = q[3 * k + 4] - ay, bz = q[3 * k + 5] - az, px = X[0] - ax, py = X[1] - ay, pz = X[2] - az, l2 = bx * bx + by * by + bz * bz;
      let u = l2 > 0 ? (px * bx + py * by + pz * bz) / l2 : 0; u = u < 0 ? 0 : u > 1 ? 1 : u;
      const d = Math.hypot(px - u * bx, py - u * by, pz - u * bz); if (d < bd) bd = d; } }
    return bd;
  };
  const pick = (score, n) => {                           // 分数最大的 n/2 个 + 其余等距抽 n/2 个
    const a = score.map((s, t) => [s, t]).sort((p, q) => q[0] - p[0]), nBig = Math.min(a.length, n >> 1), step = Math.max(1, Math.floor((a.length - nBig) / Math.max(1, n - nBig))), out = [];
    for (let r = 0; r < a.length; r++) if (r < nBig || (r - nBig) % step === 0) out.push(a[r]);
    return out;
  };
  // ---- (1a) 网格 → 晕面, 三维
  const area3 = new Float64Array(nT), areaS = new Float64Array(nT);
  let area = 0;
  for (let t = 0; t < nT; t++) {
    const o = 9 * t, cr = Q => { const u = [Q[o + 3] - Q[o], Q[o + 4] - Q[o + 1], Q[o + 5] - Q[o + 2]], w = [Q[o + 6] - Q[o], Q[o + 7] - Q[o + 1], Q[o + 8] - Q[o + 2]]; return 0.5 * Math.hypot(u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]); };
    area3[t] = cr(P); areaS[t] = cr(Dr) / (DEG * DEG); area += area3[t];
  }
  const A = { checked: 0, areaChecked: 0, bad3: 0, bad10: 0, noRef: 0, worst: [] };
  for (const [a, t] of pick(Array.from(area3), NT)) {
    const o = 9 * t, X = [(P[o] + P[o + 3] + P[o + 6]) / 3, (P[o + 1] + P[o + 4] + P[o + 7]) / 3, (P[o + 2] + P[o + 5] + P[o + 8]) / 3], rr = rOf(X);
    if (rr < 0.02 * D) continue;
    const hX = Math.asin((L[2] - X[2]) / Math.hypot(X[0] - L[0], X[1] - L[1], X[2] - L[2])) / DEG;
    let miss = dPoly(X, layerAt(hX), 'p');
    if (miss > 0.03 + 0.003 * rr) for (const dh of [0.02, -0.02, 0.06, -0.06]) { miss = Math.min(miss, dPoly(X, layerAt(hX + dh), 'p')); if (miss <= 0.03 + 0.003 * rr) break; }
    // 参照引擎 (光锥法) 在这个 h 附近一条层线都没有 (例: Lowitz 族 ζ → 0 的边上, 光锥法不出层线): 无从比较, 单列 "无参照", 不算假面
    if (miss === Infinity) { A.noRef += a; continue; }
    A.checked++; A.areaChecked += a;
    if (miss > 0.03 + 0.003 * rr) A.bad3 += a;
    if (miss > 0.10 + 0.01 * rr) { A.bad10 += a; A.worst.push({ miss, a, X, hs: [mh[3 * t], mh[3 * t + 1], mh[3 * t + 2]], v: [0, 1, 2].map(k => [P[o + 3 * k], P[o + 3 * k + 1], P[o + 3 * k + 2]]) }); }
  }
  A.worst.sort((p, q) => q.a - p.a);
  // ---- (1b) 网格 → 晕面, 天空
  const azel = d => [Math.atan2(d[1], d[0]) / DEG, Math.asin(Math.max(-1, Math.min(1, d[2] / (Math.hypot(d[0], d[1], d[2]) || 1)))) / DEG];
  const S = { checked: 0, areaChecked: 0, bad03: 0, bad1: 0, worst: [], total: 0 };
  for (let t = 0; t < nT; t++) S.total += areaS[t];
  for (const [a, t] of pick(Array.from(areaS), NT)) {
    const o = 9 * t, d = [Dr[o] + Dr[o + 3] + Dr[o + 6], Dr[o + 1] + Dr[o + 4] + Dr[o + 7], Dr[o + 2] + Dr[o + 5] + Dr[o + 8]], dl = Math.hypot(d[0], d[1], d[2]);
    if (!(dl > 1e-6) || !(a > 0)) continue;
    d[0] /= dl; d[1] /= dl; d[2] /= dl;
    const hs = [mh[3 * t], mh[3 * t + 1], mh[3 * t + 2]], h0 = Math.min(...hs), h1 = Math.max(...hs), hm = (hs[0] + hs[1] + hs[2]) / 3;
    let miss = Infinity;
    for (const h of [hm, h0, h1, 0.5 * (hm + h0), 0.5 * (hm + h1)]) { miss = Math.min(miss, dPoly(d, layerAt(h), 'd') / DEG); if (miss <= 0.3) break; }
    S.checked++; S.areaChecked += a;
    if (miss > 0.3) S.bad03 += a;
    if (miss > 1) { S.bad1 += a; S.worst.push({ miss, a, dir: azel(d), hs, v: [0, 1, 2].map(k => [P[o + 3 * k], P[o + 3 * k + 1], P[o + 3 * k + 2]]), vd: [0, 1, 2].map(k => azel([Dr[o + 3 * k], Dr[o + 3 * k + 1], Dr[o + 3 * k + 2]])) }); }
  }
  S.worst.sort((p, q) => q.a - p.a);
  // ---- (2) 晕面 → 网格: 三维分桶 + 天空分桶
  const CELL = 0.5, grid = new Map(), key = (i, jj, k) => (i + 512) * 1048576 + (jj + 512) * 1024 + (k + 512);
  const bigT = [], sgrid = new Map(), skey = (ia, ie) => ((ia % 360) + 360) % 360 * 256 + (ie + 90);
  let giant = 0;
  for (let t = 0; t < nT; t++) {
    const o = 9 * t, lo = [0, 1, 2].map(c => Math.floor(Math.min(P[o + c], P[o + 3 + c], P[o + 6 + c]) / CELL)), hi = [0, 1, 2].map(c => Math.floor(Math.max(P[o + c], P[o + 3 + c], P[o + 6 + c]) / CELL));
    if ((hi[0] - lo[0] + 1) * (hi[1] - lo[1] + 1) * (hi[2] - lo[2] + 1) > 4096) bigT.push(t);
    else for (let i = lo[0]; i <= hi[0]; i++) for (let jj = lo[1]; jj <= hi[1]; jj++) for (let k = lo[2]; k <= hi[2]; k++) { const q = key(i, jj, k); let a = grid.get(q); if (!a) grid.set(q, a = []); a.push(t); }
    // 天空: (方位, 仰角) 1° 的格; 跨度超过 30° 的三角形不算覆盖 (本身就是错的)
    const ae = [0, 1, 2].map(k => azel([Dr[o + 3 * k], Dr[o + 3 * k + 1], Dr[o + 3 * k + 2]]));
    let a0 = Math.min(ae[0][0], ae[1][0], ae[2][0]), a1 = Math.max(ae[0][0], ae[1][0], ae[2][0]);
    if (a1 - a0 > 180) { const w = ae.map(q => q[0] < 0 ? q[0] + 360 : q[0]); a0 = Math.min(...w); a1 = Math.max(...w); }
    const e0 = Math.min(ae[0][1], ae[1][1], ae[2][1]), e1 = Math.max(ae[0][1], ae[1][1], ae[2][1]);
    const edge = (i, k) => Math.acos(Math.max(-1, Math.min(1, Dr[o + 3 * i] * Dr[o + 3 * k] + Dr[o + 3 * i + 1] * Dr[o + 3 * k + 1] + Dr[o + 3 * i + 2] * Dr[o + 3 * k + 2]))) / DEG;
    if (Math.max(edge(0, 1), edge(1, 2), edge(2, 0)) > 30) { giant++; continue; }
    if (e1 > 80 || e0 < -80) { a0 = -180; a1 = 179.999; }          // 天顶 / 天底附近方位角失真: 整圈都登记
    for (let ia = Math.floor(a0); ia <= Math.floor(a1); ia++) for (let ie = Math.floor(e0); ie <= Math.floor(e1); ie++) { const q = skey(ia, ie); let a = sgrid.get(q); if (!a) sgrid.set(q, a = []); a.push(t); }
  }
  const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]], dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const ptTri = (p, a, b, c) => {
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
    if (va <= 0 && d4 - d3 >= 0 && d5 - d6 >= 0) { const w = (d4 - d3) / ((d4 - d3) + (d5 - d6)), bc = sub(c, b); return Math.hypot(bp[0] - w * bc[0], bp[1] - w * bc[1], bp[2] - w * bc[2]); }
    const den = 1 / (va + vb + vc), v = vb * den, w = vc * den;
    return Math.hypot(ap[0] - v * ab[0] - w * ac[0], ap[1] - v * ab[1] - w * ac[1], ap[2] - v * ab[2] - w * ac[2]);
  };
  const dMesh = p => {
    const c = p.map(x => Math.floor(x / CELL)); let dm = Infinity;
    for (let i = c[0] - 1; i <= c[0] + 1; i++) for (let jj = c[1] - 1; jj <= c[1] + 1; jj++) for (let k = c[2] - 1; k <= c[2] + 1; k++) {
      const a = grid.get(key(i, jj, k)); if (!a) continue;
      for (const t of a) { const o = 9 * t, d = ptTri(p, [P[o], P[o + 1], P[o + 2]], [P[o + 3], P[o + 4], P[o + 5]], [P[o + 6], P[o + 7], P[o + 8]]); if (d < dm) dm = d; }
    }
    for (const t of bigT) { const o = 9 * t, d = ptTri(p, [P[o], P[o + 1], P[o + 2]], [P[o + 3], P[o + 4], P[o + 5]], [P[o + 6], P[o + 7], P[o + 8]]); if (d < dm) dm = d; }     // 大三角形单列
    return dm;
  };
  // 视方向 q 到网格天空像的夹角 (度): 三角形的三个视方向按心射投影到 q 的切平面上, 量原点到平面三角形的距离
  const dSky = q => {
    const [qa, qe] = azel(q); let dm = Infinity;
    for (let ia = Math.floor(qa) - 1; ia <= Math.floor(qa) + 1; ia++) for (let ie = Math.floor(qe) - 1; ie <= Math.floor(qe) + 1; ie++) {
      const a = sgrid.get(skey(ia, ie)); if (!a) continue;
      for (const t of a) {
        const o = 9 * t, pr = k => { const x = Dr[o + 3 * k], y = Dr[o + 3 * k + 1], zz = Dr[o + 3 * k + 2], w = x * q[0] + y * q[1] + zz * q[2]; return w > 0.2 ? [x / w - q[0], y / w - q[1], zz / w - q[2]] : null; };
        const a3 = pr(0), b3 = pr(1), c3 = pr(2);
        if (!a3 || !b3 || !c3) continue;
        const d = ptTri([0, 0, 0], a3, b3, c3) / DEG; if (d < dm) dm = d;
      }
    }
    return dm;
  };
  const mk = () => ({ len: 0, bad3: 0, bad10: 0, worst: [], ang: 0, bad03: 0, bad1: 0, worstS: [] });
  const scan = (R, q, dq, tag) => {              // 一条折线 (q: 位置; dq: 视方向, 没有就按 眼 → 点)
    const n = q.length / 3;
    let run = null, runS = null;
    for (let k = 0; k < n; k++) {
      const X = [q[3 * k], q[3 * k + 1], q[3 * k + 2]], rr = rOf(X), rl = Math.hypot(X[0] - L[0], X[1] - L[1], X[2] - L[2]), re = Math.hypot(X[0], X[1], X[2]);
      const seg = k + 1 < n ? Math.hypot(q[3 * k + 3] - X[0], q[3 * k + 4] - X[1], q[3 * k + 5] - X[2]) : 0;
      if (rl > MAXD - 0.5 || rl < 1e-3 * D || re < 1e-3 * D) { run = runS = null; continue; }
      const dv = dq ? [dq[3 * k], dq[3 * k + 1], dq[3 * k + 2]] : [X[0] / re, X[1] / re, X[2] / re];
      const dn = k + 1 < n ? (dq ? [dq[3 * k + 3], dq[3 * k + 4], dq[3 * k + 5]] : (() => { const r2 = Math.hypot(q[3 * k + 3], q[3 * k + 4], q[3 * k + 5]) || 1; return [q[3 * k + 3] / r2, q[3 * k + 4] / r2, q[3 * k + 5] / r2]; })()) : dv;
      const segS = Math.acos(Math.max(-1, Math.min(1, dot(dv, dn)))) / DEG;
      if (rr >= 0.02 * D) {
        const d = dMesh(X);
        R.len += seg;
        if (d > 0.03 + 0.003 * rr) R.bad3 += seg;
        if (d > 0.10 + 0.01 * rr) { R.bad10 += seg; if (!run) R.worst.push(run = { len: 0, d: 0, X, tag }); run.len += seg; if (d > run.d) { run.d = Math.min(d, 9.99); run.X = X; } } else run = null;
      } else run = null;
      const ds = dSky(dv);
      R.ang += segS;
      if (ds > 0.3) R.bad03 += segS;
      if (ds > 1) { R.bad1 += segS; if (!runS) R.worstS.push(runS = { len: 0, d: 0, dir: azel(dv), X, tag }); runS.len += segS; if (ds > runS.d) { runS.d = Math.min(ds, 99.9); runS.dir = azel(dv); runS.X = X; } } else runS = null;
    }
  };
  const Bq = mk(), Fq = mk();
  let nl = 0;
  let hsX = levels;
  if (!hsX && layers) {
    hsX = [];
    const nPair = layers.length - 1, stepL = Math.max(1, Math.floor(nPair / Math.max(1, NL)));
    for (let i = stepL >> 1; i < nPair; i += stepL) if (layers[i + 1].h - layers[i].h >= 1e-4) hsX.push(0.5 * (layers[i].h + layers[i + 1].h));
  }
  if (!hsX) {                    // 外来的网格: 场景缓存里的光锥层 h 范围内等距取 (错开整层)
    const a = layers0[0].h, b = layers0[layers0.length - 1].h;
    hsX = Array.from({ length: NL }, (_, i) => a + (b - a) * (i + 0.37) / NL);
  }
  for (const h of hsX) { nl++; for (const c of layerAt(h)) scan(Bq, c.p, c.d, 'h ' + h.toFixed(3)); }
  feats.forEach((f, fi) => scan(Fq, f.p, null, f.g + fi));
  for (const R of [Bq, Fq]) { R.worst.sort((p, q) => q.len - p.len); R.worstS.sort((p, q) => q.len - p.len); }
  return { nT, area, rows: layers ? layers.length : 0, A, S, B: Bq, F: Fq, nl, giant };
}

export function printCheck(title, r, top = 5) {
  const f = X => X.map(x => x.toFixed(2)).join(','), f1 = X => X.map(x => x.toFixed(1)).join(',');
  console.log(`${title}: ${r.nT} 三角形 ${r.area.toFixed(1)} m², ${r.rows} 层`);
  console.log(`  网格→晕面 三维: 查 ${r.A.checked} 个 (${r.A.areaChecked.toFixed(1)} m²); 假面 >3cm ${r.A.bad3.toFixed(2)} m², >10cm ${r.A.bad10.toFixed(2)} m²${r.A.noRef ? `; 参照引擎在该 h 没有层线 (无从比较) ${r.A.noRef.toFixed(2)} m²` : ''}`);
  for (const w of r.A.worst.slice(0, top)) console.log(`     偏差 ${w.miss.toFixed(2)} m 面积 ${w.a.toFixed(3)} m² 重心 ${f(w.X)}; 顶点 ${w.v.map((q, k) => `(${f(q)} h${w.hs[k].toFixed(3)})`).join(' ')}`);
  console.log(`  网格→晕面 天空: 查 ${r.S.checked} 个 (${r.S.areaChecked.toFixed(0)} 平方度, 全网格 ${r.S.total.toFixed(0)}); 假面 >0.3° ${r.S.bad03.toFixed(1)} 平方度, >1° ${r.S.bad1.toFixed(1)} 平方度; 跨度 >30° 的三角形 ${r.giant} 个`);
  for (const w of r.S.worst.slice(0, top)) console.log(`     偏差 ${w.miss.toFixed(1)}° 面积 ${w.a.toFixed(1)} 平方度 方向 (方位 ${w.dir[0].toFixed(1)}, 仰角 ${w.dir[1].toFixed(1)}); 顶点 ${w.v.map((q, k) => `(${f(q)} h${w.hs[k].toFixed(3)} 天空 ${f1(w.vd[k])})`).join(' ')}`);
  for (const [name, R] of [[`层线 (另算 ${r.nl} 层)`, r.B], ['边界线 / 脊线', r.F]]) {
    console.log(`  晕面→网格 ${name}: 三维 ${R.len.toFixed(1)} m 里漏 >3cm ${R.bad3.toFixed(2)} m, >10cm ${R.bad10.toFixed(2)} m; 天空 ${R.ang.toFixed(0)}° 里漏 >0.3° ${R.bad03.toFixed(1)}°, >1° ${R.bad1.toFixed(1)}°`);
    for (const w of R.worst.slice(0, top)) console.log(`     三维: ${w.tag} 有 ${w.len.toFixed(2)} m 离网格最远 ${w.d.toFixed(2)} m, 在 ${f(w.X)}`);
    for (const w of R.worstS.slice(0, top)) console.log(`     天空: ${w.tag} 有 ${w.len.toFixed(1)}° 离网格最远 ${w.d.toFixed(1)}°, 方向 (方位 ${w.dir[0].toFixed(1)}, 仰角 ${w.dir[1].toFixed(1)}), 位置 ${f(w.X)}`);
  }
}

if (process.argv[1] && path.resolve(process.argv[1]) === fileURLToPath(import.meta.url)) {
  const [code = '4-1-8', mode = 'lowitz', el = '3.5', Ds = '15', surfName = 'caustic', fixAng = '0'] = process.argv.slice(2);
  const t0 = Date.now(), r = physCheck2(code, mode, el, Ds, surfName, fixAng, +(process.env.NT || 500), +(process.env.NL || 40));
  if (!r) { console.log('no cache: run edges.mjs first'); process.exit(2); }
  printCheck(`${[code, mode, el, Ds, surfName, fixAng].join(' ')} (${Date.now() - t0} ms)`, r, +(process.env.TOP || 5));
}
