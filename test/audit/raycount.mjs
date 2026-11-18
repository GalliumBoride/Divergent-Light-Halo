// 独立的光线追踪核对 (不用任何引擎的链式公式): 空间一点 P 处, 有多少个取向能把灯光经给定光路折进眼睛.
//   node test/audit/raycount.mjs code mode el D fix  x,y,z [x,y,z ...]      (二自由度: mode = lowitz | column)
//   取向 R = Rz(α)·Ry(ζ)·Rz(θ); lowitz: θ 定 (fix), α、ζ 自由; column: ζ 定 (fix), α、θ 自由.
//   对每个 P: 入射 u = (P − L)/|…|, 应有的出射 o = (E − P)/|…|. 逐面按 Snell 定律追 (入射面折入, 各反射面镜面反射, 出射面折出),
//   在两个自由角上扫 720×720 网格, 夹角 < 1° 的点逐个加密到 < 1e-4°, 数互不相同的解.
//   怎么读: 二自由度取向族照亮的是空间里的一块体积, 一般的点上有 1、2、3 … 个孤立的解.
//     穿过一张面时解的个数差 2 → 焦散 (折叠);  差 1 → 一张边界像 (一个解在这里走到了有效区的边上);
//     一个点上有一整条取向曲线都是解 → 焦线上的点.
//   用途: 判断引擎给的一条 "边界线" 是不是真在焦散面上 (见 README: 3-5-6-8 Lowitz Σ10.5 焦散的 8 条伪线)
const SRC = new URL('../../src/', import.meta.url).href;
const CR = await import(SRC + 'engine/crystal.js');
CR.setApex(56, 56);
const [code, mode, els, Ds, fixs, ...pts] = process.argv.slice(2);
const DEG = Math.PI / 180, N = 1.31, el = +els * DEG, D = +Ds, L = [D * Math.cos(el), 0, D * Math.sin(el)], E = [0, 0, 0], fix = +fixs * DEG;
const nr = CR.parsePath(code).map(CR.baseNormalFromCompact).map(n => { const l = Math.hypot(...n); return n.map(x => x / l); });
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const Rz = a => [[Math.cos(a), -Math.sin(a), 0], [Math.sin(a), Math.cos(a), 0], [0, 0, 1]];
const Ry = a => [[Math.cos(a), 0, Math.sin(a)], [0, 1, 0], [-Math.sin(a), 0, Math.cos(a)]];
const mul = (A, B) => A.map((r, i) => [0, 1, 2].map(j => r[0] * B[0][j] + r[1] * B[1][j] + r[2] * B[2][j]));
const app = (A, v) => A.map(r => r[0] * v[0] + r[1] * v[1] + r[2] * v[2]);
const appT = (A, v) => [0, 1, 2].map(j => A[0][j] * v[0] + A[1][j] * v[1] + A[2][j] * v[2]);
// 折射: d 传播方向 (单位), n 面的外法向; eta = n1/n2. 返回 null = 全反射 / 方向不对
const refr = (d, n, eta, entering) => {
  const c = dot(d, n);
  if (entering ? c >= 0 : c <= 0) return null;
  const nn = entering ? n : n.map(x => -x), ci = -dot(d, nn), k = 1 - eta * eta * (1 - ci * ci);
  if (k < 0) return null;
  return d.map((x, i) => eta * x + (eta * ci - Math.sqrt(k)) * nn[i]);
};
const trace = dIn => {                     // 晶体系入射传播方向 → 出射传播方向 (null = 走不通)
  let d = refr(dIn, nr[0], 1 / N, true);
  if (!d) return null;
  for (let k = 1; k < nr.length - 1; k++) { const c = dot(d, nr[k]); if (c <= 0) return null; d = d.map((x, i) => x - 2 * c * nr[k][i]); }
  return refr(d, nr[nr.length - 1], N, false);
};
const low = mode === 'lowitz';
const Rof = (a, b) => low ? mul(mul(Rz(a), Ry(b)), Rz(fix)) : mul(mul(Rz(a), Ry(fix)), Rz(b));
for (const ps of pts) {
  const P = ps.split(',').map(Number), u0 = P.map((x, i) => x - L[i]), ul = Math.hypot(...u0), u = u0.map(x => x / ul), o0 = E.map((x, i) => x - P[i]), ol = Math.hypot(...o0), o = o0.map(x => x / ol);
  const f = (a, b) => { const R = Rof(a, b), d = trace(appT(R, u)); if (!d) return Infinity; const w = app(R, d); return Math.acos(Math.max(-1, Math.min(1, dot(w, o)))); };
  const NA = 720, NB = 720;
  let best = Infinity, ba = 0, bb = 0, nHit = 0;
  const hits = [];
  for (let i = 0; i < NA; i++) for (let j = 0; j < NB; j++) { const a = 2 * Math.PI * i / NA, b = 2 * Math.PI * j / NB, v = f(a, b); if (v < best) { best = v; ba = a; bb = b; } if (v < 1.0 * DEG) hits.push([a, b, v]); }
  // 局部加密
  let st = 2 * Math.PI / NA;
  for (let it = 0; it < 40; it++) { st *= 0.5; for (const [da, db] of [[1, 0], [-1, 0], [0, 1], [0, -1], [1, 1], [-1, -1], [1, -1], [-1, 1]]) { const v = f(ba + da * st, bb + db * st); if (v < best) { best = v; ba += da * st; bb += db * st; } } }
  // 网格上夹角 < 1° 的点各自加密, 数能收敛到 < 1e-4° 的有多少个 (互相隔开 > 0.5° 的算不同的解)
  const sols = [];
  for (const [a0, b0] of hits) {
    let a = a0, b = b0, v = f(a, b), s = 2 * Math.PI / NA;
    for (let it = 0; it < 40; it++) { s *= 0.5; for (const [da, db] of [[1, 0], [-1, 0], [0, 1], [0, -1], [1, 1], [-1, -1], [1, -1], [-1, 1]]) { const w = f(a + da * s, b + db * s); if (w < v) { v = w; a += da * s; b += db * s; } } }
    if (v < 1e-4 * DEG && !sols.some(q => Math.hypot(Math.sin(q[0] - a), Math.sin(q[1] - b), Math.cos(q[0] - a) - 1, Math.cos(q[1] - b) - 1) < 0.5 * DEG)) sols.push([a, b]);
  }
  console.log(`P (${P.map(x => x.toFixed(3)).join(', ')}) r ${Math.hypot(...P).toFixed(3)}: 最小夹角 ${(best / DEG).toExponential(2)}° 在 (${(ba / DEG).toFixed(2)}, ${(bb / DEG).toFixed(2)}); 网格上 < 1° 的 ${hits.length} 点, 精确解 (夹角 < 1e-4°, 互隔 > 0.5°) ${sols.length} 个${sols.length > 12 ? ' —— 一整条取向曲线 (焦线)' : ''}`);
}
