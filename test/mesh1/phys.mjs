// 一自由度网格的物理校验: 网格上的点是不是真的在晕面上.
//   node test/mesh1/phys.mjs <光路> <Σ> <D> <ζ> <θ> [预算]
// 晕面的定义与网格怎么构造无关: 点 X 在晕面上 ⟺ 存在晶体方位 α, 使灯 → X 的光线经光路偏折后正好射向眼睛.
//   对每个三角形的重心 X: û = (X − L)/|X − L|, ê = (E − X)/|E − X|, 扫 α 求出射方向 d(α) 与 ê 的最小夹角 δ (黄金分割细化).
//   换成长度: 偏差 ≈ δ·t·s/(t + s) (t、s = X 到灯、眼的距离; X 横移 ε 时 û、ê 各转 ε/t、ε/s).
//   这个换算在出射临界 (全反射) 边界旁不成立: 那里出射方向 ∝ √余量, X 偏出边界 ε 时夹角是 √(ε/t) 量级, 几毫米的弦高
//   会算成十几厘米. 所以偏差超过 2 cm 的再查一遍: 在 X 周围 1 cm、4 cm 的 26 个邻点上各算一次, 取 (挪动距离 + 该点偏差)
//   的最小值 —— X 离晕面不到 1 cm 时总有一个邻点落在有效一侧, 偏差回到线性.
//   真在面上的三角形偏差只有弦高那么大 (毫米级); 连错的三角形 (跨过空洞、两片之间乱连) 偏差是分米到米.
// 输出: 偏差超过 3 cm + 0.3% r / 10 cm + 1% r 的三角形面积 (r = 到灯、眼较近者的距离: 几十米外的弦高本来就大), 以及最差的几个. 环境变量 M1 = 换一份 mesh1 实现 (对比用).
import { fileURLToPath, pathToFileURL } from 'node:url';
import path from 'node:path';

const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../..');
const eng = f => pathToFileURL(path.join(ROOT, 'src/engine', f)).href;
const CR = await import(eng('crystal.js')), LPm = await import(eng('lp.js'));
const M1 = await import(process.env.M1 ? pathToFileURL(path.resolve(process.env.M1)).href : eng('mesh1.js'));

// mesh: 不给 → 校验 mesh1.js 的网格; 给 { pos } → 校验这张网格 (球壳法的网格用它)
export function physCheck(code, el, D, zeta, theta, budget = 160000, maxD = 80, mesh = null) {
  CR.setApex(56, 56); LPm.LP.partial = true; LPm.setIndex(1.31);
  const fp = CR.parsePath(code), N = 1.31, eta = 1 / N;
  const m = mesh || M1.surfaceMesh1(fp, el, D, 0, N, maxD, zeta, theta, { budget });
  const canon = CR.buildCanonNormals(zeta, theta), nf = fp.length, cn = fp.map(f => canon[f]);
  const er = el * Math.PI / 180, L = [D * Math.cos(er), 0, D * Math.sin(er)];
  const AS = 720, TWO = 2 * Math.PI, wn = new Float64Array(3 * nf);
  // 出射方向与 ê 的夹角余弦 (无效 → −2)
  const cosAt = (a, U, e) => {
    const ca = Math.cos(a), sa = Math.sin(a);
    for (let k = 0; k < nf; k++) { const q = cn[k]; wn[3 * k] = ca * q[0] - sa * q[1]; wn[3 * k + 1] = sa * q[0] + ca * q[1]; wn[3 * k + 2] = q[2]; }
    const cosI = -(U[0] * wn[0] + U[1] * wn[1] + U[2] * wn[2]);
    if (!(cosI > 0)) return -2;
    const c1 = eta * cosI - Math.sqrt(1 - eta * eta * (1 - cosI * cosI));
    let ix = eta * U[0] + c1 * wn[0], iy = eta * U[1] + c1 * wn[1], iz = eta * U[2] + c1 * wn[2];
    for (let k = 1; k < nf - 1; k++) {
      const dn = ix * wn[3 * k] + iy * wn[3 * k + 1] + iz * wn[3 * k + 2];
      if (!(dn > 0)) return -2;
      ix -= 2 * dn * wn[3 * k]; iy -= 2 * dn * wn[3 * k + 1]; iz -= 2 * dn * wn[3 * k + 2];
    }
    const j = 3 * (nf - 1), c2i = ix * wn[j] + iy * wn[j + 1] + iz * wn[j + 2];
    if (!(c2i > 0)) return -2;
    const tir = 1 - N * N * (1 - c2i * c2i);
    if (!(tir >= 0)) return -2;
    const c2 = Math.sqrt(tir) - N * c2i;
    const dx = N * ix + c2 * wn[j], dy = N * iy + c2 * wn[j + 1], dz = N * iz + c2 * wn[j + 2], dl = Math.hypot(dx, dy, dz);
    return (dx * e[0] + dy * e[1] + dz * e[2]) / dl;
  };
  const missAt = X => {
    let r = miss0(X);
    if (!(r > 0.02)) return r;
    for (const eps of [0.01, 0.04]) {
      for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) for (let dz = -1; dz <= 1; dz++) {
        if (!dx && !dy && !dz) continue;
        const q = miss0([X[0] + eps * dx, X[1] + eps * dy, X[2] + eps * dz]) + eps * Math.hypot(dx, dy, dz);
        if (q < r) r = q;
      }
      if (r <= 3 * eps) break;
    }
    return r;
  };
  const miss0 = X => {
    const t = Math.hypot(X[0] - L[0], X[1] - L[1], X[2] - L[2]), s = Math.hypot(X[0], X[1], X[2]);
    if (!(t > 1e-9 && s > 1e-9)) return 0;
    const U = [(X[0] - L[0]) / t, (X[1] - L[1]) / t, (X[2] - L[2]) / t], e = [-X[0] / s, -X[1] / s, -X[2] / s];
    let best = -2, bi = 0;
    for (let i = 0; i < AS; i++) { const c = cosAt(TWO * i / AS, U, e); if (c > best) { best = c; bi = i; } }
    if (best < -1) return Infinity;
    let lo = TWO * (bi - 1) / AS, hi = TWO * (bi + 1) / AS;
    const g = 0.381966011250105;
    let x1 = lo + g * (hi - lo), x2 = hi - g * (hi - lo), f1 = cosAt(x1, U, e), f2 = cosAt(x2, U, e);
    for (let it = 0; it < 40; it++) {
      if (f1 < f2) { lo = x1; x1 = x2; f1 = f2; x2 = hi - g * (hi - lo); f2 = cosAt(x2, U, e); }
      else { hi = x2; x2 = x1; f2 = f1; x1 = lo + g * (hi - lo); f1 = cosAt(x1, U, e); }
    }
    best = Math.max(best, f1, f2);
    return Math.acos(Math.min(1, best)) * t * s / (t + s);
  };
  const P = m.pos, nT = P.length / 9, tris = [];
  let area = 0;
  for (let t = 0; t < nT; t++) {
    const o = 9 * t, u = [P[o + 3] - P[o], P[o + 4] - P[o + 1], P[o + 5] - P[o + 2]], w = [P[o + 6] - P[o], P[o + 7] - P[o + 1], P[o + 8] - P[o + 2]];
    const a = 0.5 * Math.hypot(u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]);
    area += a; tris.push([a, t]);
  }
  // 全查太慢时: 面积最大的 4000 个全查 (连错的三角形都很大), 其余等距抽样
  tris.sort((p, q) => q[0] - p[0]);
  const cap = +(process.env.CAP || 30000), step = Math.max(1, Math.ceil((nT - 4000) / Math.max(1, cap - 4000)));
  let bad2 = 0, bad10 = 0, checked = 0, worst = [];
  for (let r = 0; r < tris.length; r++) {
    if (r >= 4000 && (r - 4000) % step) continue;
    const [a, t] = tris[r], o = 9 * t, wgt = r < 4000 ? 1 : step;
    const X = [(P[o] + P[o + 3] + P[o + 6]) / 3, (P[o + 1] + P[o + 4] + P[o + 7]) / 3, (P[o + 2] + P[o + 5] + P[o + 8]) / 3];
    const miss = missAt(X);
    checked++;
    const rr = Math.min(Math.hypot(X[0] - L[0], X[1] - L[1], X[2] - L[2]), Math.hypot(X[0], X[1], X[2]));
    if (miss > 0.03 + 0.003 * rr) bad2 += a * wgt;
    if (miss > 0.10 + 0.01 * rr) { bad10 += a * wgt; worst.push([miss, a, X, m.h ? [m.h[3 * t], m.h[3 * t + 1], m.h[3 * t + 2]] : [0, 0, 0]]); }
  }
  worst.sort((p, q) => q[1] - p[1]);
  return { nT, area, checked, bad2, bad10, worst, ms: m.ms, nh: m.nh };
}

if (process.argv[1] && path.resolve(process.argv[1]) === fileURLToPath(import.meta.url)) {
  const [code = '3-5', el = '10.5', D = '15', zeta = '-138', theta = '60', budget = '160000'] = process.argv.slice(2);
  const r = physCheck(code, +el, +D, +zeta, +theta, +budget);
  console.log(`${code} Σ${el} D${D} ζ${zeta} θ${theta}: ${r.nT} 三角形 ${r.area.toFixed(2)} m², 查了 ${r.checked} 个; 偏差 > 3 cm + 0.3% r 的面积 ${r.bad2.toFixed(3)} m², > 10 cm + 1% r 的 ${r.bad10.toFixed(3)} m²`);
  for (const [miss, a, X, hs] of r.worst.slice(0, +(process.env.TOP || 6))) console.log(`   偏差 ${miss.toFixed(2)} m 面积 ${a.toFixed(3)} m² 重心 ${X.map(x => x.toFixed(2)).join(',')}  三顶点 h ${hs.map(x => x.toFixed(5)).join(' ')}`);
}
