// 一自由度晕面的三角网格: 参数 (h, α) 网格三角化 + 边界二分拟合.
//   h = 灯处入射仰角, α = 晶体方位. 每个节点: 形式链 (无效区也连续延拓) → 共面闭式 ψ 双支 → 晶体位置 P.
//   每支一个"有效余量" M = min(入射 cosθ1, 镜面 i·N_k, 出射 cosθ1', 全反射余量, 根余量 1−|rhs|, t, s, 截断余量):
//   M > 0 ⟺ 该支有点. 三角形按 M 的符号裁剪, 跨零的边沿精确 M 二分 30 次 → 网格边缘精确落在
//   入射掠射 / 镜面掠射 / 出射临界 / ψ 转折 / 捏合 / 截断 上 (与边界引擎 v4 的曲线重合).
// 网格范围: 先用粗网格 (h 0.25°, α 0.5°) 找出有点的 h 区间与 α 区间, 再在其中按预算铺细网格 (格子在度数上近似正方).
// 坐标与光锥法完全一致: û = Rz(ψ)·(−cos h, 0, −sin h), 世界法向 = Rz(α)·Ry(ζ)·Rz(θ)·n_f.

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

/**
 * 节点求值 (法向已由 setAlpha 设好): 写入 out = [M_A, M_B, PA(3), PB(3)]. 无点的支 P 为 NaN.
 */
function evalNode(cx, h, out) {
  const ch = Math.cos(h), sh = Math.sin(h);
  const Ux = -ch, Uy = 0, Uz = -sh;
  const wn = cx.wn, nf = cx.nf, eta = cx.eta, N = cx.N;
  // 形式链
  let ax = wn[0], ay = wn[1], az = wn[2];
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
  out[2] = out[3] = out[4] = out[5] = out[6] = out[7] = NaN;
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
    const c = ux * ex + uy * ey + uz * ez, det = 1 - c * c;
    let pm = -1, t = NaN;
    if (Math.abs(det) >= 1e-12) {
      const uD = ux * cx.Dx + uy * cx.Dy + uz * cx.Dz, dD = ex * cx.Dx + ey * cx.Dy + ez * cx.Dz;
      t = (c * dD - uD) / det;
      const s = (c * uD - dD) / det;
      pm = Math.min(t - 1e-6, s - 1e-6, cx.maxDist - Math.abs(t));
    }
    const M = Math.min(mask, rootM, pm);
    out[b] = M;
    if (M > 0) {
      const o = 2 + 3 * b;
      out[o] = cx.L[0] + t * ux; out[o + 1] = cx.L[1] + t * uy; out[o + 2] = cx.L[2] + t * uz;
    }
  }
}

// 有点的区间: 布尔占用 (可周期) → 膨胀 pad 格 → 区间列表 [i0, i1] (下标)
function intervals(occ, pad, periodic) {
  const n = occ.length, dil = new Uint8Array(n);
  for (let i = 0; i < n; i++) if (occ[i])
    for (let k = -pad; k <= pad; k++) {
      let j = i + k;
      if (periodic) j = ((j % n) + n) % n; else if (j < 0 || j >= n) continue;
      dil[j] = 1;
    }
  const out = [];
  let i = 0;
  while (i < n) {
    if (!dil[i]) { i++; continue; }
    let j = i;
    while (j + 1 < n && dil[j + 1]) j++;
    out.push([i, j]);
    i = j + 1;
  }
  // 周期: 首尾相接的两段合并 (末段延伸过 n)
  if (periodic && out.length > 1 && out[0][0] === 0 && out[out.length - 1][1] === n - 1) {
    const last = out.pop();
    out[0] = [last[0], out[0][1] + n];
  }
  if (periodic && out.length === 1 && out[0][0] === 0 && out[0][1] === n - 1) out[0] = [0, n];
  return out;
}

/**
 * 网格: 返回 { pos: Float32Array (三角形顶点, 不共享), h: Float32Array (每顶点 h, 度), nTri, ms }.
 * budget = 细网格节点总数 (约).
 */
export function surfaceMesh1(facePath, lampElev, lampDist, eyeHeight, refN, maxDist, zetaDeg, thetaDeg, budget = 250000) {
  const t0 = Date.now();
  const cx = makeCtx(facePath, lampElev, lampDist, eyeHeight, refN, maxDist, zetaDeg, thetaDeg);
  const rec = new Float64Array(8);
  // ---- 1) 粗网格找范围
  const HC = 720, AC = 720;                 // h: −89.875…89.875 每 0.25°, α: 每 0.5°
  const hOcc = new Uint8Array(HC), aOcc = new Uint8Array(AC);
  const hAt = i => (-90 + 180 * (i + 0.5) / HC) * DEG;
  for (let ia = 0; ia < AC; ia++) {
    setAlpha(cx, 2 * Math.PI * ia / AC);
    for (let ih = 0; ih < HC; ih++) {
      evalNode(cx, hAt(ih), rec);
      if (rec[0] > 0 || rec[1] > 0) { hOcc[ih] = 1; aOcc[ia] = 1; }
    }
  }
  const hInt = intervals(hOcc, 3, false), aInt = intervals(aOcc, 3, true);
  if (!hInt.length || !aInt.length) return { pos: new Float32Array(0), h: new Float32Array(0), nTri: 0, ms: Date.now() - t0 };
  const h0 = Math.max(-89.99, -90 + 180 * hInt[0][0] / HC), h1 = Math.min(89.99, -90 + 180 * (hInt[hInt.length - 1][1] + 1) / HC);
  const aSpans = aInt.map(([i0, i1]) => [360 * i0 / AC, 360 * Math.min(i1 + 1, i0 + AC) / AC]);
  const aTot = aSpans.reduce((s, [a, b]) => s + (b - a), 0);
  // ---- 2) 细网格 (格子近似正方)
  const cell = Math.max(0.004, Math.sqrt((h1 - h0) * aTot / budget));
  const NH = Math.min(4000, Math.max(8, Math.ceil((h1 - h0) / cell) + 1));
  const hs = new Float64Array(NH);
  for (let i = 0; i < NH; i++) hs[i] = (h0 + (h1 - h0) * i / (NH - 1)) * DEG;
  const pos = [], hv = [];
  // 二分: 从有点端 (M>0) 到无点端, 返回边上最后一个有点的参数处的 P 与 h (度)
  const tmp = new Float64Array(8);
  const crossing = (b, ha, aa, hb, ab) => {
    let lo = 0, hi = 1;
    for (let it = 0; it < 30; it++) {
      const mid = 0.5 * (lo + hi);
      setAlpha(cx, aa + (ab - aa) * mid);
      evalNode(cx, ha + (hb - ha) * mid, tmp);
      if (tmp[b] > 0) lo = mid; else hi = mid;
    }
    setAlpha(cx, aa + (ab - aa) * lo);
    const hh = ha + (hb - ha) * lo;
    evalNode(cx, hh, tmp);
    const o = 2 + 3 * b;
    return [tmp[o], tmp[o + 1], tmp[o + 2], hh / DEG];
  };
  for (const [a0, a1] of aSpans) {
    const NA = Math.min(20000, Math.max(4, Math.ceil((a1 - a0) / cell) + 1));
    const as = new Float64Array(NA);
    for (let j = 0; j < NA; j++) as[j] = (a0 + (a1 - a0) * j / (NA - 1)) * DEG;
    // 一列一列 (α 外层): 每列 NH 个节点, 存 M[2] 与 P[2]
    let prev = null;
    for (let j = 0; j < NA; j++) {
      setAlpha(cx, as[j]);
      const col = { M: new Float64Array(2 * NH), P: new Float64Array(6 * NH) };
      for (let i = 0; i < NH; i++) {
        evalNode(cx, hs[i], rec);
        col.M[2 * i] = rec[0]; col.M[2 * i + 1] = rec[1];
        for (let k = 0; k < 6; k++) col.P[6 * i + k] = rec[2 + k];
      }
      if (prev) {
        const ap = as[j - 1], ac = as[j];
        for (let i = 0; i + 1 < NH; i++) {
          // 格子四角: (i, j−1) (i+1, j−1) (i+1, j) (i, j); 分两个三角形
          const V = [[prev, i, ap], [prev, i + 1, ap], [col, i + 1, ac], [col, i, ac]];
          for (let b = 0; b < 2; b++) {
            const m = V.map(([c, ii]) => c.M[2 * ii + b]);
            if (m[0] <= 0 && m[1] <= 0 && m[2] <= 0 && m[3] <= 0) continue;
            for (const tri of [[0, 1, 2], [0, 2, 3]]) clipTri(tri.map(k => V[k]), b, hs, crossing, pos, hv, cx.maxDist);
          }
        }
      }
      prev = col;
    }
  }
  return { pos: Float32Array.from(pos), h: Float32Array.from(hv), nTri: hv.length / 3, ms: Date.now() - t0 };
}

// 三角形按 M_b > 0 裁剪 (Sutherland–Hodgman), 跨零边二分; 输出扇形三角形
function clipTri(T, b, hs, crossing, pos, hv, maxDist) {
  const pts = [];
  const inside = v => v[0].M[2 * v[1] + b] > 0;
  const vert = v => { const o = 6 * v[1] + 3 * b; return [v[0].P[o], v[0].P[o + 1], v[0].P[o + 2], hs[v[1]] / DEG]; };
  for (let k = 0; k < 3; k++) {
    const p = T[k], q = T[(k + 1) % 3], ip = inside(p), iq = inside(q);
    if (ip) pts.push(vert(p));
    if (ip !== iq) {
      // 规范方向: 总是从有点端二分, 两个相邻三角形算出同一点 (网格不开缝)
      const [a, c] = ip ? [p, q] : [q, p];
      pts.push(crossing(b, hs[a[1]], a[2], hs[c[1]], c[2]));
    }
  }
  if (pts.length < 3) return;
  const lim = 0.5 * maxDist, lim2 = lim * lim;
  for (let k = 1; k + 1 < pts.length; k++) {
    const A = pts[0], B = pts[k], C = pts[k + 1];
    if (!(Number.isFinite(A[0]) && Number.isFinite(B[0]) && Number.isFinite(C[0]))) continue;
    if (d2(A, B) > lim2 || d2(B, C) > lim2 || d2(A, C) > lim2) continue;     // 奇异拉伸保护
    pos.push(A[0], A[1], A[2], B[0], B[1], B[2], C[0], C[1], C[2]);
    hv.push(A[3], B[3], C[3]);
  }
}
const d2 = (a, b) => (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2;
