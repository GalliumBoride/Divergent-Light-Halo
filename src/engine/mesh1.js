// 一自由度晕面的三角网格: 按 h 行求精确区间 + 区间内余弦 (Chebyshev) 取点 + 整行自适应加行.
//   h = 灯处入射仰角, α = 晶体方位. 每个节点: 形式链 (无效区也连续延拓) → 共面闭式 ψ 双支 → 晶体位置 P
//   与观测者视方向 (−d_div, 精确; 不用 P−E 相减, 靠近眼时不失真).
//   每支一个"有效余量" M = min(入射 cosθ1, 镜面 i·N_k, 出射 cosθ1', 全反射余量, 根余量 1−|rhs|, t, s, 截断余量):
//   M > 0 ⟺ 该支有点.
// 为什么不用 (h, α) 均匀网格: 映射在有效区的边上是奇异的 — 出射临界 (全反射) 旁出射方向 ∝ √余量, ψ 两支在根余量 = 0
//   处以 √ 汇合; 有效区收口的地方 (两条 √ 边夹成的细楔) 参数上只有一条细缝, 三维里却铺开几米到十几米. 均匀网格
//   在那里只落进去寥寥几个格点, 交点之间用直线连, 晕面就碎成锯齿翼片. 加密 (暴力) 治标不治本.
// 这里换坐标让映射变光滑:
//   1) 每个 h 行先把该支的有效 α 区间精确求出来 (余量变号处二分到 1e-13; 比取样还窄的细区间用局部极大值找出来);
//   2) 区间内取点 α = α_a + (α_b − α_a)·(1 − cos πτ)/2. 边上位置 ∝ √(α − α_a), 而 α − α_a ∝ τ², 所以位置对 τ 光滑 —
//      √ 奇点被坐标吸收, 边界点就是精确根, 网格边缘严丝合缝落在边界上; 两支在根余量 = 0 处用的是同一个根, 自然接上;
//   3) 行: 相邻两行在三维里离得远 (或天空张角大), 或两行的区间对不上 (区间出生 / 消失 / 合并 / 分裂), 就在中间插一整行,
//      细楔收口处行自动变密. 整行插入, 不会有 T 形接头;
//   4) 每条"区间链" (相邻行一一对应的区间串) 用同一个点数, 相邻两行的点按 τ 对齐连成四边形带.
// 节点总数不超预算. 坐标与光锥法完全一致: û = Rz(ψ)·(−cos h, 0, −sin h), 世界法向 = Rz(α)·Ry(ζ)·Rz(θ)·n_f.

import { buildCanonNormals } from './crystal.js';

const DEG = Math.PI / 180;

function makeCtx(facePath, lampElev, lampDist, eyeHeight, refN, maxDist, zetaDeg, thetaDeg) {
  const e = lampElev * DEG;
  const Lx = lampDist * Math.cos(e), Lz = lampDist * Math.sin(e);
  const canon = buildCanonNormals(zetaDeg, thetaDeg);
  return {
    fp: facePath, nf: facePath.length, canon: facePath.map(f => canon[f]),
    N: refN, eta: 1 / refN, L: [Lx, 0, Lz], Dx: Lx, Dy: 0, Dz: Lz - eyeHeight, maxDist,
    wn: new Float64Array(3 * facePath.length),
  };
}

// 光路各面在方位 α 下的世界法向 → cx.wn
function setAlpha(cx, a) {
  const ca = Math.cos(a), sa = Math.sin(a), wn = cx.wn;
  for (let k = 0; k < cx.nf; k++) {
    const m = cx.canon[k];
    wn[3 * k] = ca * m[0] - sa * m[1]; wn[3 * k + 1] = sa * m[0] + ca * m[1]; wn[3 * k + 2] = m[2];
  }
}

// 节点记录布局: [M_A, M_B, PA(3), PB(3), dirA(3), dirB(3)] — 无点的支 P/dir 为 NaN
export const REC = 14;

/** 节点求值 (法向已由 setAlpha 设好), 写入 out (长度 REC). */
function evalNode(cx, h, out) {
  const ch = Math.cos(h), sh = Math.sin(h);
  const Ux = -ch, Uy = 0, Uz = -sh;
  const wn = cx.wn, nf = cx.nf, eta = cx.eta, N = cx.N;
  for (let k = 2; k < REC; k++) out[k] = NaN;
  // 形式链
  const ax = wn[0], ay = wn[1], az = wn[2];
  const cosI = -(Ux * ax + Uy * ay + Uz * az);
  let mask = cosI;
  const cc = cosI > 0 ? cosI : 0;
  const c1 = eta * cc - Math.sqrt(Math.max(0, 1 - eta * eta * (1 - cc * cc)));
  let ix = eta * Ux + c1 * ax, iy = eta * Uy + c1 * ay, iz = eta * Uz + c1 * az;
  for (let k = 1; k < nf - 1; k++) {
    const mx = wn[3 * k], my = wn[3 * k + 1], mz = wn[3 * k + 2];
    const dn = ix * mx + iy * my + iz * mz;
    if (dn < mask) mask = dn;
    ix -= 2 * dn * mx; iy -= 2 * dn * my; iz -= 2 * dn * mz;
  }
  const j = 3 * (nf - 1);
  const bx = wn[j], by = wn[j + 1], bz = wn[j + 2];
  const cosI2 = ix * bx + iy * by + iz * bz;
  if (cosI2 < mask) mask = cosI2;
  const tir = 1 - N * N * (1 - cosI2 * cosI2);
  if (tir < mask) mask = tir;
  const c2 = Math.sqrt(Math.max(0, tir)) - N * (cosI2 > 0 ? cosI2 : 0);
  let dx = N * ix + c2 * bx, dy = N * iy + c2 * by, dz = N * iz + c2 * bz;
  const dl = Math.sqrt(dx * dx + dy * dy + dz * dz);
  if (!(dl > 1e-12)) { out[0] = out[1] = -1; return; }
  dx /= dl; dy /= dl; dz /= dl;
  // 共面闭式: D·(Rz(ψ)U × Rz(ψ)d) = 0
  const cxv = Uy * dz - Uz * dy, cyv = Uz * dx - Ux * dz, czv = Ux * dy - Uy * dx;
  const A = cx.Dy * cxv - cx.Dx * cyv, B = cx.Dx * cxv + cx.Dy * cyv, C = -cx.Dz * czv;
  const mag = Math.sqrt(A * A + B * B);
  if (mag < 1e-12) { out[0] = out[1] = -1; return; }
  const rhs = C / mag, rootM = 1 - Math.abs(rhs);
  if (rootM < 0) { out[0] = out[1] = Math.min(rootM, mask); return; }
  const phi = Math.atan2(B, A), asr = Math.asin(Math.max(-1, Math.min(1, rhs)));
  for (let b = 0; b < 2; b++) {
    const psi = b === 0 ? asr - phi : Math.PI - asr - phi;
    const cp = Math.cos(psi), sp = Math.sin(psi);
    const ux = Ux * cp - Uy * sp, uy = Ux * sp + Uy * cp, uz = Uz;
    const ex = dx * cp - dy * sp, ey = dx * sp + dy * cp, ez = dz;
    // 两条线 (灯 L + t·û, 眼 E − s·ê) 的交点. 分母用 |û×ê|² 而不是 1 − (û·ê)²: 两者相等, 但 û、ê 接近平行
    //   (晶体趋于无穷远) 时后者有相消误差, 叉积保持完整的相对精度; 分子同样写成叉积 (D = L − E):
    //   t = −(D×ê)·(û×ê) / |û×ê|²,  s = (D×û)·(û×ê) / |û×ê|²
    const wx = uy * ez - uz * ey, wy = uz * ex - ux * ez, wz = ux * ey - uy * ex, det = wx * wx + wy * wy + wz * wz;
    let pm = -1, t = NaN;
    if (det >= 1e-12) {
      const Dx = cx.Dx, Dy = cx.Dy, Dz = cx.Dz;
      t = -((Dy * ez - Dz * ey) * wx + (Dz * ex - Dx * ez) * wy + (Dx * ey - Dy * ex) * wz) / det;
      const s = ((Dy * uz - Dz * uy) * wx + (Dz * ux - Dx * uz) * wy + (Dx * uy - Dy * ux) * wz) / det;
      pm = Math.min(t - 1e-6, s - 1e-6, cx.maxDist - Math.abs(t));
    }
    const M = Math.min(mask, rootM, pm);
    out[b] = M;
    if (M > 0) {
      const o = 2 + 3 * b, q = 8 + 3 * b;
      out[o] = cx.L[0] + t * ux; out[o + 1] = cx.L[1] + t * uy; out[o + 2] = cx.L[2] + t * uz;
      out[q] = -ex; out[q + 1] = -ey; out[q + 2] = -ez;          // 观测者看到它的方向 = 出射传播方向反向
    }
  }
}

/**
 * 网格: 返回 { pos: Float32Array (三角形顶点, 不共享), dir: Float32Array (每顶点视方向), h: Float32Array (每顶点 h, 度),
 *   nrm (每顶点法向), nTri, ms, nh (行数), na (每行平均点数), nDrop (缺角的四边形), nSub (插入的行) }.
 * opts: budget (节点总数上限).
 */
export function surfaceMesh1(facePath, lampElev, lampDist, eyeHeight, refN, maxDist, zetaDeg, thetaDeg, opts = {}) {
  const t0 = Date.now();
  const budget = opts.budget || 250000;
  const L3 = Math.max(0.05, Math.min(lampDist, maxDist) / 120), A3 = 0.5 * DEG;
  const cx = makeCtx(facePath, lampElev, lampDist, eyeHeight, refN, maxDist, zetaDeg, thetaDeg);
  const own = cx.wn, rec = new Float64Array(REC), nf3 = 3 * cx.nf, TWO = 2 * Math.PI;
  const empty = { pos: new Float32Array(0), dir: new Float32Array(0), h: new Float32Array(0), nrm: new Float32Array(0), nTri: 0, ms: 0, nh: 0, na: 0, nDrop: 0, nSub: 0 };
  // α 取样 (每行找区间用) 的法向预先算好
  const AS = 1440, WN = new Float64Array(AS * nf3);
  for (let i = 0; i < AS; i++) { cx.wn = WN.subarray(nf3 * i, nf3 * (i + 1)); setAlpha(cx, TWO * i / AS); }
  cx.wn = own;
  const at = (h, a) => { setAlpha(cx, a); evalNode(cx, h, rec); return rec; };
  const Mof = (h, a, b) => at(h, a)[b];

  // ---- 1) h 范围: 粗扫 (h 0.25°, α 0.5°)
  const HC = 720;
  let i0 = -1, i1 = -1;
  for (let ih = 0; ih < HC; ih++) {
    const h = (-90 + 180 * (ih + 0.5) / HC) * DEG;
    for (let i = 0; i < AS; i += 2) {
      cx.wn = WN.subarray(nf3 * i, nf3 * (i + 1)); evalNode(cx, h, rec);
      if (rec[0] > 0 || rec[1] > 0) { if (i0 < 0) i0 = ih; i1 = ih; break; }
    }
  }
  cx.wn = own;
  if (i0 < 0) return { ...empty, ms: Date.now() - t0 };
  const hLo = Math.max(-89.99, -90 + 180 * (i0 - 2) / HC) * DEG, hHi = Math.min(89.99, -90 + 180 * (i1 + 3) / HC) * DEG;

  // ---- 2) 一行: 两支各自的有效区间 [{ a, b (> a, 可越过 2π), full }]
  const root = (h, b, aIn, aOut) => {               // aIn 有效, aOut 无效 → 边上最后一个有效的 α
    let lo = aIn, hi = aOut;
    for (let it = 0; it < 52; it++) { const md = 0.5 * (lo + hi); if (Mof(h, md, b) > 0) lo = md; else hi = md; }
    return lo;
  };
  const gmax = (h, b, lo, hi) => {                   // 黄金分割求 M 的最大值
    const g = 0.381966011250105;
    let x1 = lo + g * (hi - lo), x2 = hi - g * (hi - lo), f1 = Mof(h, x1, b), f2 = Mof(h, x2, b);
    for (let it = 0; it < 40; it++) {
      if (f1 < f2) { lo = x1; x1 = x2; f1 = f2; x2 = hi - g * (hi - lo); f2 = Mof(h, x2, b); }
      else { hi = x2; x2 = x1; f2 = f1; x1 = lo + g * (hi - lo); f1 = Mof(h, x1, b); }
    }
    return f1 > f2 ? [x1, f1] : [x2, f2];
  };
  const intervalsOf = (h, M, b) => {
    const da = TWO / AS;
    let all = true;
    for (let i = 0; i < AS; i++) if (!(M[i] > 0)) { all = false; break; }
    if (all) return [{ a: 0, b: TWO, full: true }];
    const ev = [];                                    // [α, +1 起 / −1 止]
    for (let i = 0; i < AS; i++) {
      const j = (i + 1) % AS, ai = da * i, aj = ai + da, pi = M[i] > 0, pj = M[j] > 0;
      if (pi && !pj) ev.push([root(h, b, ai, aj), -1]);
      else if (!pi && pj) ev.push([root(h, b, aj, ai), 1]);
      else if (!pi && !pj) {
        // 比取样还窄的区间: 两个无效样本之间 M 的局部极大若够近 0, 求出最大值, 大于 0 就是一个细区间
        const k = (i + AS - 1) % AS;
        if (M[i] > -0.05 && M[i] >= M[k] && M[i] >= M[j]) {
          const [am, fm] = gmax(h, b, ai - da, aj);
          if (fm > 0) { ev.push([root(h, b, am, ai - da), 1]); ev.push([root(h, b, am, aj), -1]); }
        }
      }
    }
    if (!ev.length) return [];
    ev.sort((p, q) => p[0] - q[0]);
    let s0 = ev.findIndex(e => e[1] === 1);
    if (s0 < 0) return [];
    const out = [];
    for (let n = 0; n < ev.length; n++) {
      const e = ev[(s0 + n) % ev.length];
      if (e[1] !== 1) continue;
      const f = ev[(s0 + n + 1) % ev.length];
      if (!f || f[1] !== -1) continue;
      let bb = f[0];
      while (bb <= e[0]) bb += TWO;
      out.push({ a: e[0], b: bb, full: false });
    }
    return out;
  };
  // 一行. hints = 相邻已有的行: 它们有而这一行取样没抓到的区间 (宽度小于取样间隔, 有效区正在收口), 在预测位置
  //   用黄金分割找 M 的最大值, 大于 0 就补上 —— 细区间不会这一行有、下一行没有 (那样每对行都像拓扑变化, 会无止境地插行)
  const makeRow = (h, hints = []) => {
    const MA = new Float64Array(AS), MB = new Float64Array(AS);
    for (let i = 0; i < AS; i++) { cx.wn = WN.subarray(nf3 * i, nf3 * (i + 1)); evalNode(cx, h, rec); MA[i] = rec[0]; MB[i] = rec[1]; }
    cx.wn = own;
    const ints = [intervalsOf(h, MA, 0), intervalsOf(h, MB, 1)];
    for (let b = 0; b < 2; b++) for (const H of hints) for (const I of H.ints[b]) {
      if (I.full || ints[b].some(K => overlapShift(I, K) !== null)) continue;
      const w = Math.max(I.b - I.a, 0.5 * DEG), lo = I.a - w, hi = I.b + w;
      const [am, fm] = gmax(h, b, lo, hi);
      if (!(fm > 0) || Mof(h, lo, b) > 0 || Mof(h, hi, b) > 0) continue;
      const K = { a: root(h, b, am, lo), b: root(h, b, am, hi), full: false };
      if (!ints[b].some(J => overlapShift(J, K) !== null)) ints[b].push(K);
    }
    for (const L of ints) L.sort((p, q) => p.a - q.a);
    return { h, ints, pilot: [new Map(), new Map()] };
  };

  // 区间内 τ 处的 α (余弦分布; 整圈的用均匀分布)
  const alphaAt = (I, tau) => I.full ? I.a + TWO * tau : tau <= 0 ? I.a : tau >= 1 ? I.b : I.a + (I.b - I.a) * 0.5 * (1 - Math.cos(Math.PI * tau));
  // 两行区间的对应: 放宽 d 后重叠 (考虑 2π 平移). 返回平移量或 null.
  //   端点是 h 的光滑函数, 相邻两行之间最多挪 "速度 × Δh"; 很细的区间 (有效区收口处) 挪得比自己还宽, 只看重叠会误判成
  //   "一个消失、一个出生". 放宽 20·Δh; 放宽后若一对多, 照样当事件插行夹紧, 不会错连.
  function overlapShift(I, K, d = 0) {
    if (I.full || K.full) return 0;
    for (const sh of [0, -TWO, TWO]) if (Math.max(I.a - d, K.a + sh) < Math.min(I.b + d, K.b + sh)) return sh;
    return null;
  }
  // 两行某支的配对: pairs = 互相只重叠对方的区间对 [[i, k, 平移]] (连成链); loose = 没有一一对应的区间 [[行 0/1, 号]]
  //   (出生 / 消失 / 合并 / 分裂 — 要在这两行之间插行把事件夹紧)
  const pairUp = (R0, R1, b) => {
    const A = R0.ints[b], B = R1.ints[b], all = [];
    const degA = new Array(A.length).fill(0), degB = new Array(B.length).fill(0);
    const d = 20 * Math.abs(R1.h - R0.h);
    A.forEach((I, i) => B.forEach((K, k) => { const sh = overlapShift(I, K, d); if (sh !== null) { all.push([i, k, sh]); degA[i]++; degB[k]++; } }));
    const pairs = all.filter(([i, k]) => degA[i] === 1 && degB[k] === 1 && A[i].full === B[k].full);
    const loose = [];
    A.forEach((_, i) => { if (!pairs.some(p => p[0] === i)) loose.push([0, i]); });
    B.forEach((_, k) => { if (!pairs.some(p => p[1] === k)) loose.push([1, k]); });
    return { pairs, loose };
  };
  // 试探点 (9 个 τ) 的位置与视方向, 按行缓存
  const NP = 9;
  const pilot = (R, b, i) => {
    let v = R.pilot[b].get(i);
    if (v) return v;
    const I = R.ints[b][i];
    v = new Float64Array(NP * 7);
    for (let j = 0; j < NP; j++) {
      const r = at(R.h, alphaAt(I, j / (NP - (I.full ? 0 : 1)))), o = 2 + 3 * b, q = 8 + 3 * b;
      v.set([r[b], r[o], r[o + 1], r[o + 2], r[q], r[q + 1], r[q + 2]], 7 * j);
    }
    R.pilot[b].set(i, v);
    return v;
  };
  // 区间在三维里的长度与天空张角 (按试探点)
  const extentOf = (R, b, i) => {
    const P = pilot(R, b, i);
    let l = 0, a = 0;
    for (let j = 1; j < NP; j++) {
      const o = 7 * j, q = o - 7;
      if (!(P[o] > 0 && P[q] > 0)) continue;
      l += Math.hypot(P[o + 1] - P[q + 1], P[o + 2] - P[q + 2], P[o + 3] - P[q + 3]);
      a += Math.acos(Math.min(1, P[o + 4] * P[q + 4] + P[o + 5] * P[q + 5] + P[o + 6] * P[q + 6]));
    }
    return [l, a];
  };
  // 相邻两行之间该不该再插一行: 有对不上的区间 (且它在三维里不是一个小点) → 是, 直到行距 HEV;
  //   对得上的区间看对应 τ 点的三维距离 / 天空张角
  const HEV = 2e-4 * DEG, HGAP = 1e-3 * DEG;
  let Lrow = 2 * L3, Arow = 2 * A3;
  const gapNeed = (R0, R1) => {
    const dh = R1.h - R0.h;
    let worst = 0;
    for (let b = 0; b < 2; b++) {
      const { pairs, loose } = pairUp(R0, R1, b);
      if (dh > HEV && loose.some(([w, i]) => { const [l, a] = extentOf(w ? R1 : R0, b, i); return l > 0.25 * L3 || a > 0.25 * A3; })) return Infinity;
      for (const [i, k] of pairs) {
        const P0 = pilot(R0, b, i), P1 = pilot(R1, b, k);
        for (let j = 0; j < NP; j++) {
          const o = 7 * j;
          if (!(P0[o] > 0 && P1[o] > 0)) continue;
          const d = Math.hypot(P0[o + 1] - P1[o + 1], P0[o + 2] - P1[o + 2], P0[o + 3] - P1[o + 3]);
          const c = P0[o + 4] * P1[o + 4] + P0[o + 5] * P1[o + 5] + P0[o + 6] * P1[o + 6];
          worst = Math.max(worst, d / Lrow, Math.acos(Math.min(1, c)) / Arow);
        }
      }
    }
    return dh > HGAP ? worst : 0;
  };

  // ---- 3) 行: 先均匀 48 行, 再按需插行 (区间对不上 / 太远). 行数受预算约束
  let rows = [];
  for (let k = 0; k <= 48; k++) rows.push(makeRow(hLo + (hHi - hLo) * k / 48, k ? [rows[k - 1]] : []));
  // 粗扫看不见比取样还细的楔尖 (有效区在 h 的两头常常收成一条细缝): 最外一行还有区间, 就借它带路继续往外加行, 直到真空
  const stepH = (hHi - hLo) / 48, HLIM = 89.99 * DEG;
  const hasInts = R => R.ints[0].length + R.ints[1].length > 0;
  for (let n = 0; n < 400 && hasInts(rows[0]) && rows[0].h > -HLIM; n++) rows.unshift(makeRow(Math.max(-HLIM, rows[0].h - stepH), [rows[0]]));
  for (let n = 0; n < 400 && hasInts(rows[rows.length - 1]) && rows[rows.length - 1].h < HLIM; n++) rows.push(makeRow(Math.min(HLIM, rows[rows.length - 1].h + stepH), [rows[rows.length - 1]]));
  const nInts = () => rows.reduce((s, R) => s + R.ints[0].length + R.ints[1].length, 0);
  let nSub = 0;
  for (let pass = 0; pass < 40; pass++) {
    const mark = new Uint8Array(rows.length);
    for (let k = 0; k + 1 < rows.length; k++) if (gapNeed(rows[k], rows[k + 1]) > 1) mark[k] = 1;
    // 边界的弯曲: 一行的区间端点离上下两行端点连线 (按 h 插值) 的偏差 = 折线边缘的弦高; 超过 0.04·L3 (或从眼看 0.05°)
    //   就把上下两个间隔都再分 (只用已算好的试探点, 不另求值)
    for (let k = 1; k + 1 < rows.length; k++) {
      if (mark[k - 1] && mark[k]) continue;
      const h1 = rows[k].h - rows[k - 1].h, h2 = rows[k + 1].h - rows[k].h;
      if (h1 < HGAP || h2 < HGAP) continue;
      for (let b = 0; b < 2 && !(mark[k - 1] && mark[k]); b++) {
        const lo = pairUp(rows[k - 1], rows[k], b).pairs, hi = pairUp(rows[k], rows[k + 1], b).pairs;
        for (const [i, k1] of lo) {
          const p2 = hi.find(p => p[0] === k1);
          if (!p2) continue;
          const A = pilot(rows[k - 1], b, i), Bp = pilot(rows[k], b, k1), C = pilot(rows[k + 1], b, p2[1]);
          for (const j of [0, NP - 1]) {
            const o = 7 * j;
            if (!(A[o] > 0 && Bp[o] > 0 && C[o] > 0)) continue;
            let dev = 0, dd = 0;
            for (let c = 1; c <= 3; c++) dev = Math.hypot(dev, Bp[o + c] - (h2 * A[o + c] + h1 * C[o + c]) / (h1 + h2));
            for (let c = 4; c <= 6; c++) dd = Math.hypot(dd, Bp[o + c] - (h2 * A[o + c] + h1 * C[o + c]) / (h1 + h2));
            if (dev > 0.04 * L3 * Lrow / (2 * L3) || dd > 0.05 * DEG * Arow / (2 * A3)) { mark[k - 1] = mark[k] = 1; break; }
          }
        }
      }
    }
    const add = [];
    for (let k = 0; k + 1 < rows.length; k++) if (mark[k]) add.push(k);
    if (!add.length) break;
    // 预算: 每个区间至少 ~12 个点; 超了就放宽距离标准
    if ((nInts() + add.length * 2) * 12 > budget) {
      Lrow *= 1.5; Arow *= 1.5;
      if (Lrow > 200 * L3) break;
      continue;
    }
    const next = [];
    let ai = 0;
    for (let k = 0; k < rows.length; k++) {
      next.push(rows[k]);
      if (ai < add.length && add[ai] === k) { next.push(makeRow(0.5 * (rows[k].h + rows[k + 1].h), [rows[k], rows[k + 1]])); ai++; nSub++; }
    }
    rows = next;
  }

  // ---- 4) 区间链: 相邻行一一对应的区间串成链; 每条链一个点数 (按最长那一行的三维长度 / 天空张角)
  const chains = [];
  for (let b = 0; b < 2; b++) {
    let open = new Map();                                  // 上一行区间号 → 链
    for (let k = 0; k < rows.length; k++) {
      const R = rows[k], next = new Map();
      const pairs = k > 0 ? pairUp(rows[k - 1], R, b).pairs : null;
      R.ints[b].forEach((I, i) => {
        let ch = null, sh = 0;
        if (pairs) { const p = pairs.find(q => q[1] === i); if (p && open.has(p[0])) { ch = open.get(p[0]); sh = p[2]; } }
        if (!ch) { ch = { b, rows: [], ints: [], full: I.full }; chains.push(ch); sh = 0; }
        const prev = ch.ints.length ? ch.ints[ch.ints.length - 1] : null;
        // 平移到与上一行连续 (α 展开)
        let J = { a: I.a, b: I.b, full: I.full };
        if (prev && !I.full) { const s2 = sh + TWO * Math.round((prev.a - (I.a + sh)) / TWO); J = { a: I.a + s2, b: I.b + s2, full: false }; }
        ch.rows.push(k); ch.ints.push(J); ch.src = ch.src || []; ch.src.push(i);
        next.set(i, ch);
      });
      open = next;
    }
  }
  const lenOf = ch => {                                   // 链上最长一行的三维长度与天空张角
    let L = 0, A = 0;
    ch.rows.forEach((k, m) => { const [l, a] = extentOf(rows[k], ch.b, ch.src[m]); L = Math.max(L, l); A = Math.max(A, a); });
    return [L, A];
  };
  const want = chains.map(ch => { const [L, A] = lenOf(ch); return Math.max(L / L3, A / A3); });
  const nodesFor = sc => chains.reduce((s, ch, c) => s + ch.rows.length * Math.min(600, Math.max(6, Math.ceil(want[c] / sc) + 1)), 0);
  let sc = 1;
  if (nodesFor(1) > budget) { let lo = 1, hi = 2; while (nodesFor(hi) > budget && hi < 1e6) hi *= 2; for (let it = 0; it < 30; it++) { const md = Math.sqrt(lo * hi); if (nodesFor(md) > budget) lo = md; else hi = md; } sc = hi; }

  // ---- 5) 节点 + 四边形带
  const pos = [], dirv = [], hv = [], vids = [];
  let nid = 0, nDrop = 0, nodes = 0;
  for (let c = 0; c < chains.length; c++) {
    const ch = chains[c], b = ch.b, o = 2 + 3 * b, q = 8 + 3 * b;
    const N = Math.min(600, Math.max(6, Math.ceil(want[c] / sc) + 1));
    const ncol = ch.full ? N : N;                          // 整圈: N 个点首尾相接; 区间: N 个点含两端
    const grid = ch.rows.map((k, m) => {
      const I = ch.ints[m], h = rows[k].h, row = [];
      for (let j = 0; j < ncol; j++) {
        let r = at(h, alphaAt(I, ch.full ? j / N : j / (N - 1)));
        // 端点就是根, 舍入 (α 平移 2π 等) 可能让它落到无效那一侧: 往区间里挪一点点, 别把边界节点丢了
        if (!(r[b] > 0) && !ch.full && (j === 0 || j === ncol - 1))
          for (let e = 1e-13, sg = j === 0 ? 1 : -1; e < 1e-6 && !(r[b] > 0); e *= 2) r = at(h, (j === 0 ? I.a : I.b) + sg * e);
        row.push(r[b] > 0 ? { id: nid++, P: [r[o], r[o + 1], r[o + 2]], d: [r[q], r[q + 1], r[q + 2]], h: h / DEG } : null);
      }
      nodes += ncol;
      return row;
    });
    const put = (A, B, C) => {
      pos.push(...A.P, ...B.P, ...C.P); dirv.push(...A.d, ...B.d, ...C.d); hv.push(A.h, B.h, C.h);
      vids.push(A.id, B.id, C.id);
    };
    const d2 = (A, B) => (A.P[0] - B.P[0]) ** 2 + (A.P[1] - B.P[1]) ** 2 + (A.P[2] - B.P[2]) ** 2;
    for (let m = 0; m + 1 < grid.length; m++) {
      const r0 = grid[m], r1 = grid[m + 1];
      for (let j = 0; j < (ch.full ? ncol : ncol - 1); j++) {
        const j1 = (j + 1) % ncol;
        const A = r0[j], B = r0[j1], C = r1[j1], D = r1[j];   // (τ, h) 平面里逆时针
        const n = (A ? 1 : 0) + (B ? 1 : 0) + (C ? 1 : 0) + (D ? 1 : 0);
        if (n === 4) { if (d2(A, C) <= d2(B, D)) { put(A, B, C); put(A, C, D); } else { put(A, B, D); put(B, C, D); } }
        else if (n === 3) { nDrop++; if (!A) put(B, C, D); else if (!B) put(A, C, D); else if (!C) put(A, B, D); else put(A, B, C); }
        else if (n) nDrop++;
      }
    }
  }
  return {
    pos: Float32Array.from(pos), dir: Float32Array.from(dirv), h: Float32Array.from(hv), nrm: vertexNormals(pos, vids, nid),
    nTri: hv.length / 3, ms: Date.now() - t0, nh: rows.length, na: Math.round(nodes / Math.max(1, rows.length)), nDrop, nSub,
  };
}

// 每顶点法向 (3D 视图平滑明暗用): 同一个顶点 (同一支) 在各三角形里是同一个身份, 累加面积加权的面法向.
//   同一支的三角形在参数平面里都是逆时针 (格子两三角形、裁剪扇形、对分子三角形都保持方向), 映到三维后朝向一致,
//   所以法向可以带符号累加; 两支分开 (在转折处它们朝向相反, 混在一起会抵消).
function vertexNormals(pos, vids, nid) {
  const nV = vids.length, out = new Float32Array(3 * nV), acc = new Float64Array(6 * nid);
  for (let t = 0; 3 * t < nV; t++) {
    const o = 9 * t;
    const ux = pos[o + 3] - pos[o], uy = pos[o + 4] - pos[o + 1], uz = pos[o + 5] - pos[o + 2];
    const vx = pos[o + 6] - pos[o], vy = pos[o + 7] - pos[o + 1], vz = pos[o + 8] - pos[o + 2];
    const nx = uy * vz - uz * vy, ny = uz * vx - ux * vz, nz = ux * vy - uy * vx;
    for (let k = 0; k < 3; k++) { const i = 3 * vids[3 * t + k]; acc[i] += nx; acc[i + 1] += ny; acc[i + 2] += nz; }
  }
  for (let v = 0; v < nV; v++) {
    const i = 3 * vids[v], x = acc[i], y = acc[i + 1], z = acc[i + 2], l = Math.hypot(x, y, z) || 1;
    out[3 * v] = x / l; out[3 * v + 1] = y / l; out[3 * v + 2] = z / l;
  }
  return out;
}
