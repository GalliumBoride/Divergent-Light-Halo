// 一自由度晕面的三角网格: 参数 (h, α) 网格三角化 + 边界二分拟合.
//   h = 灯处入射仰角, α = 晶体方位. 每个节点: 形式链 (无效区也连续延拓) → 共面闭式 ψ 双支 → 晶体位置 P
//   与观测者视方向 (−d_div, 精确; 不用 P−E 相减, 靠近眼时不失真).
//   每支一个"有效余量" M = min(入射 cosθ1, 镜面 i·N_k, 出射 cosθ1', 全反射余量, 根余量 1−|rhs|, t, s, 截断余量):
//   M > 0 ⟺ 该支有点. 三角形按 M 的符号裁剪, 跨零的边沿精确 M 二分 → 网格边缘精确落在
//   入射掠射 / 镜面掠射 / 出射临界 / ψ 转折 / 捏合 / 截断 上 (与边界引擎 v4 的曲线重合).
// 网格: 先用粗网格 (h 0.25°, α 0.5°) 找出有点的 h 区间与 α 区间; 再在其中铺张量积网格并【各向异性加密】:
//   每个 h 行间隔 / α 列间隔按它在三维里的长度 (≤ L3) 与天空张角 (≤ A) 细分 — 板晶这类整张面只占 0.7° h 的
//   极端情形, 行会自动加密而列不会, 不会出现细长条三角形. 节点总数严格不超预算.
// 三角化时太大的三角形 (灯处捏合 / 近水平入射这类参数化奇异的细条) 在参数空间递归四分, 不再丢成破洞.
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
      const o = 2 + 3 * b, q = 8 + 3 * b;
      out[o] = cx.L[0] + t * ux; out[o + 1] = cx.L[1] + t * uy; out[o + 2] = cx.L[2] + t * uz;
      out[q] = -ex; out[q + 1] = -ey; out[q + 2] = -ez;          // 观测者看到它的方向 = 出射传播方向反向
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
  if (periodic && out.length > 1 && out[0][0] === 0 && out[out.length - 1][1] === n - 1) {
    const last = out.pop();
    out[0] = [last[0], out[0][1] + n];
  }
  if (periodic && out.length === 1 && out[0][0] === 0 && out[0][1] === n - 1) out[0] = [0, n];
  return out;
}

const dist3 = (A, o, B, q) => Math.hypot(A[o] - B[q], A[o + 1] - B[q + 1], A[o + 2] - B[q + 2]);
const angle3 = (A, o, B, q) => {
  const c = A[o] * B[q] + A[o + 1] * B[q + 1] + A[o + 2] * B[q + 2];
  return Math.acos(c > 1 ? 1 : c < -1 ? -1 : c);
};

/**
 * 网格: 返回 { pos: Float32Array (三角形顶点, 不共享), dir: Float32Array (每顶点视方向), h: Float32Array (每顶点 h, 度), nTri, ms, nh, na }.
 * opts: budget (节点总数上限), L3 (三维边长目标, m), A (天空张角目标, 度).
 */
export function surfaceMesh1(facePath, lampElev, lampDist, eyeHeight, refN, maxDist, zetaDeg, thetaDeg, opts = {}) {
  const t0 = Date.now();
  const budget = opts.budget || 250000;
  let L3 = opts.L3 || Math.max(0.05, Math.min(lampDist, maxDist) / 120);
  let Aang = (opts.A || 0.5) * DEG;
  const cx = makeCtx(facePath, lampElev, lampDist, eyeHeight, refN, maxDist, zetaDeg, thetaDeg);
  const rec = new Float64Array(REC);
  // ---- 1) 粗网格找范围
  const HC = 720, AC = 720;
  const hOcc = new Uint8Array(HC), aOcc = new Uint8Array(AC);
  const hAt = i => (-90 + 180 * (i + 0.5) / HC) * DEG;
  for (let ia = 0; ia < AC; ia++) {
    setAlpha(cx, 2 * Math.PI * ia / AC);
    for (let ih = 0; ih < HC; ih++) {
      evalNode(cx, hAt(ih), rec);
      if (rec[0] > 0 || rec[1] > 0) { hOcc[ih] = 1; aOcc[ia] = 1; }
    }
  }
  const empty = { pos: new Float32Array(0), dir: new Float32Array(0), h: new Float32Array(0), nTri: 0, ms: Date.now() - t0, nh: 0, na: 0 };
  const hInt = intervals(hOcc, 3, false), aInt = intervals(aOcc, 3, true);
  if (!hInt.length || !aInt.length) return empty;
  const h0 = Math.max(-89.99, -90 + 180 * hInt[0][0] / HC), h1 = Math.min(89.99, -90 + 180 * (hInt[hInt.length - 1][1] + 1) / HC);
  const aSpans = aInt.map(([i0, i1]) => [360 * i0 / AC, 360 * Math.min(i1 + 1, i0 + AC) / AC]);
  const aTot = aSpans.reduce((s, [a, b]) => s + (b - a), 0);

  // ---- 2) 张量积网格 + 各向异性加密
  // 初始: h 64 份; α 按跨度比例共 ~128 份
  let hs = linspace(h0 * DEG, h1 * DEG, 65);
  let as = aSpans.map(([a, b]) => linspace(a * DEG, b * DEG, Math.max(9, Math.round(128 * (b - a) / aTot) + 1)));
  const evalGrid = () => {
    const NH = hs.length, cols = [];
    for (const span of as) {
      const g = [];
      for (let j = 0; j < span.length; j++) {
        setAlpha(cx, span[j]);
        const col = new Float64Array(REC * NH);
        for (let i = 0; i < NH; i++) { evalNode(cx, hs[i], rec); col.set(rec, REC * i); }
        g.push(col);
      }
      cols.push(g);
    }
    return cols;
  };
  // 一对节点 (同支都有效) 的"太大"倍数
  const need = (A, ia, B, ib) => {
    let m = 0;
    for (let b = 0; b < 2; b++) {
      if (!(A[ia + b] > 0 && B[ib + b] > 0)) continue;
      const o = 2 + 3 * b, q = 8 + 3 * b;
      m = Math.max(m, dist3(A, ia + o, B, ib + o) / L3, angle3(A, ia + q, B, ib + q) / Aang);
    }
    return m;
  };
  let grid = evalGrid();
  let relax = 1;                       // 预算不够时目标放宽的倍数 (后面的细分判据随之放宽)
  for (let pass = 0; pass < 3; pass++) {
    const NH = hs.length;
    const rowNeed = new Float64Array(NH - 1);
    const colNeed = as.map(s => new Float64Array(s.length - 1));
    grid.forEach((g, si) => {
      for (let j = 0; j < g.length; j++) {
        const col = g[j];
        for (let i = 0; i + 1 < NH; i++) rowNeed[i] = Math.max(rowNeed[i], need(col, REC * i, col, REC * (i + 1)));
        if (j + 1 < g.length) {
          const nx = g[j + 1];
          for (let i = 0; i < NH; i++) colNeed[si][j] = Math.max(colNeed[si][j], need(col, REC * i, nx, REC * i));
        }
      }
    });
    const cap = (x) => Math.min(24, Math.ceil(x));
    const parts = (v, sc) => (v / sc > 1 ? cap(v / sc) : 1);
    // 按放宽倍数 sc 细分后的节点总数
    const count = sc => {
      let nh = 1, na = 0;
      for (const v of rowNeed) nh += parts(v, sc);
      for (const cn of colNeed) { na += 1; for (const v of cn) na += parts(v, sc); }
      return nh * na;
    };
    if (!rowNeed.some(v => v > 1) && !colNeed.some(cn => cn.some(v => v > 1))) break;
    // 超预算: 二分出最小的放宽倍数, 使节点数严格不超预算 (不细分时的节点数 = 现有的, 本就在预算内)
    let scale = 1;
    if (count(1) > budget) {
      let lo = 1, hi = 2;
      while (count(hi) > budget && hi < 1e6) hi *= 2;
      for (let it = 0; it < 24; it++) { const mid = Math.sqrt(lo * hi); if (count(mid) > budget) lo = mid; else hi = mid; }
      scale = hi;
    }
    relax = Math.max(relax, scale);
    const sub = (arr, needs) => {
      const out = [arr[0]];
      for (let i = 0; i + 1 < arr.length; i++) {
        const k = parts(needs[i], scale);
        for (let q = 1; q <= k; q++) out.push(arr[i] + (arr[i + 1] - arr[i]) * q / k);
      }
      return out;
    };
    const hs2 = sub(hs, rowNeed), as2 = as.map((s, si) => sub(s, colNeed[si]));
    if (hs2.length === hs.length && as2.every((s, si) => s.length === as[si].length)) break;
    hs = hs2; as = as2;
    grid = evalGrid();
    if (scale > 1) break;
  }

  // ---- 3) 三角化 + 裁剪; 太大的三角形在参数空间里递归四分
  //   参数化在两处奇异: 灯处捏合 (t→0) 附近与近水平入射 (晕面伸向远处) — 那里 (h, α) 走一小步, 晶体位置能跳几十米.
  //   张量积网格照顾不到这种细条, 于是对每个裁剪后的多边形检查: 边长 > Lbig 或天空张角 > Abig 就四分 (最多 6 层),
  //   直到都够小. 递归到底仍奇大的才当作真正的不连续丢掉 (旧版是一律丢, 留下破洞).
  const pos = [], dirv = [], hv = [], tb = [];      // tb: 每个三角形属于哪一支
  const tmp = new Float64Array(REC);
  let nDrop = 0, nSub = 0, evalLeft = 1500000;
  const Lbig = 8 * L3 * relax, cosAbig = Math.cos(Math.min(Math.PI / 3, 8 * Aang * relax));
  const lim2 = (0.5 * cx.maxDist) ** 2;
  // 通用顶点 { h, a, r: 节点记录 }
  const node = (h, a) => { setAlpha(cx, a); const r = new Float64Array(REC); evalNode(cx, h, r); evalLeft--; return { h, a, r }; };
  const vert = (v, b) => { const r = v.r, o = 2 + 3 * b, q = 8 + 3 * b; return [r[o], r[o + 1], r[o + 2], r[q], r[q + 1], r[q + 2], v.h / DEG]; };
  // 二分: 从有点端 A (M>0) 到无点端 B, 返回边上最后一个有点的参数处的 [P(3), dir(3), h°]
  //   规范方向 (总是从有点端出发): 相邻三角形算出同一点, 网格不开缝
  const crossing = (b, A, B) => {
    let lo = 0, hi = 1;
    for (let it = 0; it < 34; it++) {
      const mid = 0.5 * (lo + hi);
      setAlpha(cx, A.a + (B.a - A.a) * mid);
      evalNode(cx, A.h + (B.h - A.h) * mid, tmp);
      if (tmp[b] > 0) lo = mid; else hi = mid;
    }
    evalLeft -= 35;
    setAlpha(cx, A.a + (B.a - A.a) * lo);
    const hh = A.h + (B.h - A.h) * lo;
    evalNode(cx, hh, tmp);
    const o = 2 + 3 * b, q = 8 + 3 * b;
    return [tmp[o], tmp[o + 1], tmp[o + 2], tmp[q], tmp[q + 1], tmp[q + 2], hh / DEG];
  };
  // 三角形按 M_b > 0 裁剪 (Sutherland–Hodgman), 跨零边二分 → 多边形顶点
  const clipPoly = (T, b) => {
    const pts = [];
    for (let k = 0; k < 3; k++) {
      const p = T[k], q = T[(k + 1) % 3], ip = p.r[b] > 0, iq = q.r[b] > 0;
      if (ip) pts.push(vert(p, b));
      if (ip !== iq) pts.push(ip ? crossing(b, p, q) : crossing(b, q, p));
    }
    return pts;
  };
  const tooBig = pts => {
    for (let k = 0; k < pts.length; k++) {
      const A = pts[k], B = pts[(k + 1) % pts.length];
      if (d2(A, B) > Lbig * Lbig || A[3] * B[3] + A[4] * B[4] + A[5] * B[5] < cosAbig) return true;
    }
    return false;
  };
  const emit = (T, b, dep) => {
    const pts = clipPoly(T, b);
    if (pts.length < 3) return;
    if (dep < 6 && evalLeft > 0 && tooBig(pts)) {
      nSub++;
      const m = (u, v) => node(0.5 * (u.h + v.h), 0.5 * (u.a + v.a));
      const m01 = m(T[0], T[1]), m12 = m(T[1], T[2]), m20 = m(T[2], T[0]);
      emit([T[0], m01, m20], b, dep + 1); emit([m01, T[1], m12], b, dep + 1);
      emit([m20, m12, T[2]], b, dep + 1); emit([m01, m12, m20], b, dep + 1);
      return;
    }
    for (let k = 1; k + 1 < pts.length; k++) {           // 扇形三角化
      const A = pts[0], B = pts[k], C = pts[k + 1];
      if (d2(A, B) > lim2 || d2(B, C) > lim2 || d2(A, C) > lim2) { nDrop++; continue; }
      pos.push(A[0], A[1], A[2], B[0], B[1], B[2], C[0], C[1], C[2]);
      tb.push(b);
      dirv.push(A[3], A[4], A[5], B[3], B[4], B[5], C[3], C[4], C[5]);
      hv.push(A[6], B[6], C[6]);
    }
  };
  const NH = hs.length;
  grid.forEach((g, si) => {
    const A = as[si];
    const V = g.map((col, j) => {
      const arr = new Array(NH);
      for (let i = 0; i < NH; i++) arr[i] = { h: hs[i], a: A[j], r: col.subarray(REC * i, REC * i + REC) };
      return arr;
    });
    for (let j = 0; j + 1 < g.length; j++) {
      for (let i = 0; i + 1 < NH; i++) {
        // 格子四角: 0 (i, j) 1 (i+1, j) 2 (i+1, j+1) 3 (i, j+1)
        const Q = [V[j][i], V[j][i + 1], V[j + 1][i + 1], V[j + 1][i]];
        for (let b = 0; b < 2; b++) {
          const in0 = Q[0].r[b] > 0, in1 = Q[1].r[b] > 0, in2 = Q[2].r[b] > 0, in3 = Q[3].r[b] > 0;
          if (!in0 && !in1 && !in2 && !in3) continue;
          // 对角线取三维更短的一条 (四角都有效时), 三角形更匀称
          let diag02 = true;
          if (in0 && in1 && in2 && in3) { const o = 2 + 3 * b; diag02 = dist3(Q[0].r, o, Q[2].r, o) <= dist3(Q[1].r, o, Q[3].r, o); }
          if (diag02) { emit([Q[0], Q[1], Q[2]], b, 0); emit([Q[0], Q[2], Q[3]], b, 0); }
          else { emit([Q[0], Q[1], Q[3]], b, 0); emit([Q[1], Q[2], Q[3]], b, 0); }
        }
      }
    }
  });
  return {
    pos: Float32Array.from(pos), dir: Float32Array.from(dirv), h: Float32Array.from(hv), nrm: vertexNormals(pos, tb),
    nTri: hv.length / 3, ms: Date.now() - t0, nh: NH, na: as.reduce((s, a) => s + a.length, 0), nDrop, nSub,
  };
}

// 每顶点法向 (3D 视图平滑明暗用): 同一支里位置相同 (1 µm 网格) 的顶点焊在一起, 累加面积加权的面法向.
//   同一支的三角形在参数平面里都是逆时针 (格子两三角形、裁剪扇形、四分子三角形都保持方向), 映到三维后朝向一致,
//   所以法向可以带符号累加; 两支分开焊 (在转折处它们朝向相反, 混在一起会抵消).
function vertexNormals(pos, tb) {
  const nT = tb.length, out = new Float32Array(9 * nT);
  const ids = new Map(), acc = [], vid = new Int32Array(3 * nT);
  const q = v => Math.round(v * 1e6);
  for (let t = 0; t < nT; t++) {
    const o = 9 * t;
    const ux = pos[o + 3] - pos[o], uy = pos[o + 4] - pos[o + 1], uz = pos[o + 5] - pos[o + 2];
    const vx = pos[o + 6] - pos[o], vy = pos[o + 7] - pos[o + 1], vz = pos[o + 8] - pos[o + 2];
    const nx = uy * vz - uz * vy, ny = uz * vx - ux * vz, nz = ux * vy - uy * vx;
    for (let k = 0; k < 3; k++) {
      const key = tb[t] + ',' + q(pos[o + 3 * k]) + ',' + q(pos[o + 3 * k + 1]) + ',' + q(pos[o + 3 * k + 2]);
      let id = ids.get(key);
      if (id === undefined) { id = acc.length / 3; ids.set(key, id); acc.push(0, 0, 0); }
      acc[3 * id] += nx; acc[3 * id + 1] += ny; acc[3 * id + 2] += nz;
      vid[3 * t + k] = id;
    }
  }
  for (let v = 0; v < vid.length; v++) {
    const i = 3 * vid[v], x = acc[i], y = acc[i + 1], z = acc[i + 2], l = Math.hypot(x, y, z) || 1;
    out[3 * v] = x / l; out[3 * v + 1] = y / l; out[3 * v + 2] = z / l;
  }
  return out;
}

function linspace(a, b, n) {
  const out = new Array(n);
  for (let i = 0; i < n; i++) out[i] = a + (b - a) * i / (n - 1);
  return out;
}
const d2 = (a, b) => (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2;
