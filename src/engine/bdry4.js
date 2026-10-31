// 一自由度晕面的物理边界 — 边界引擎 v4.1 (Divergentlightform.cs BdryV4* 的 JS 移植, 与 C# 逐句对应).
//
// 边界族 = 判据取等号 ⟹ 链上某点被钉在一个固定圆上 ⟹ 每族一维闭式构造 (参数 τ):
//   入射掠射 (m1=0):  入射方向 U 躺在入射面切平面的单位圆上;
//   镜面掠射 (λ_k=0): 面 k 前的内部光线 g (|g|=n) 躺在该面平面内半径 n 的圆上, 逆链反射回入射面后由切向连续性重构 U;
//   出射临界 (m3=0):  出射前内部光线 g = √(n²−1)·N_out + (cosτ e1 + sinτ e2), 同样逆链重构 U.
// 每个 τ: 形式链 → 出射 d 与全部 margins (自身那条不进掩码) → 共面闭式 ψ 双支 → P.
// 事件 (沿 τ 二分 48 次): 角点 = 掩码翻转; 转折 = ψ 双根合并; 捏合/截断 = t,s 过零或出截断球.
// 随后按 3D 距离递归密化, 两支就近配对装配为折线.

import { buildCanonNormals } from './crystal.js';

const f32 = Math.fround;
const fin = x => Number.isFinite(x);
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
function dist2(a, b) { const x = a[0] - b[0], y = a[1] - b[1], z = a[2] - b[2]; return x * x + y * y + z * z; }
function reflect(v, m) { const d = dot(v, m); return [v[0] - 2 * d * m[0], v[1] - 2 * d * m[1], v[2] - 2 * d * m[2]]; }
function orthoBasis(v) {
  const t = Math.abs(v[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0];
  const e1 = cross(v, t);
  let l = Math.sqrt(dot(e1, e1));
  if (l < 1e-12) l = 1;
  e1[0] /= l; e1[1] /= l; e1[2] /= l;
  return [e1, cross(v, e1)];
}

// 形式链 (规范系, 不提前返回): m = [cosθ1, i·N_k (各镜面), cosθ1(出), TIR margin]; 返回归一化出射方向 d (退化 null)
function formalChain(cx, U, m) {
  const fp = cx.faces, last = fp.length - 1;
  const eta = 1.0 / cx.N;
  const a = cx.canon[fp[0]];
  const cosI = -dot(U, a);
  m[0] = cosI;
  const cc = Math.max(cosI, 0);
  const c1 = eta * cc - Math.sqrt(Math.max(0, 1 - eta * eta * (1 - cc * cc)));
  let ix = eta * U[0] + c1 * a[0], iy = eta * U[1] + c1 * a[1], iz = eta * U[2] + c1 * a[2];
  for (let k = 1; k < last; k++) {
    const mm = cx.canon[fp[k]];
    const dn = ix * mm[0] + iy * mm[1] + iz * mm[2];
    m[k] = dn;
    ix -= 2 * dn * mm[0]; iy -= 2 * dn * mm[1]; iz -= 2 * dn * mm[2];
  }
  const b = cx.canon[fp[last]];
  const cosI2 = ix * b[0] + iy * b[1] + iz * b[2];
  m[last] = cosI2;
  const tir = 1 - cx.N * cx.N * (1 - cosI2 * cosI2);
  m[last + 1] = tir;
  const c2 = Math.sqrt(Math.max(0, tir)) - cx.N * Math.max(cosI2, 0);
  const dx = cx.N * ix + c2 * b[0], dy = cx.N * iy + c2 * b[1], dz = cx.N * iz + c2 * b[2];
  const l = Math.sqrt(dx * dx + dy * dy + dz * dz);
  if (l < 1e-12) return null;
  return [dx / l, dy / l, dz / l];
}

// 通用锥点定位: u_world=Rz(ψ)·u, d_div=Rz(ψ)·d
function conePoint(cx, u, d, psi) {
  const cp = Math.cos(psi), sp = Math.sin(psi);
  const ux = u[0] * cp - u[1] * sp, uy = u[0] * sp + u[1] * cp, uz = u[2];
  const dx = d[0] * cp - d[1] * sp, dy = d[0] * sp + d[1] * cp, dz = d[2];
  const c = ux * dx + uy * dy + uz * dz, det = 1 - c * c;
  if (Math.abs(det) < 1e-12) return null;
  const uD = ux * cx.Dx + uy * cx.Dy + uz * cx.Dz, dD = dx * cx.Dx + dy * cx.Dy + dz * cx.Dz;
  const t = (c * dD - uD) / det, s = (c * uD - dD) / det;
  if (t < 1e-6 || s < 1e-6) return null;
  const Px = cx.L[0] + t * ux, Py = cx.L[1] + t * uy, Pz = cx.L[2] + t * uz;
  const rx = Px - cx.L[0], ry = Py - cx.L[1], rz = Pz - cx.L[2];
  if (rx * rx + ry * ry + rz * rz > cx.maxDist * cx.maxDist) return null;
  return [Px, Py, Pz];
}

// 点存在余量 min(t−1e-6, s−1e-6, maxD−|t|): >0 ⇔ conePoint 有点
function coneMargin(cx, u, d, psi) {
  const cp = Math.cos(psi), sp = Math.sin(psi);
  const ux = u[0] * cp - u[1] * sp, uy = u[0] * sp + u[1] * cp, uz = u[2];
  const dx = d[0] * cp - d[1] * sp, dy = d[0] * sp + d[1] * cp, dz = d[2];
  const c = ux * dx + uy * dy + uz * dz, det = 1 - c * c;
  if (Math.abs(det) < 1e-12) return -1;
  const uD = ux * cx.Dx + uy * cx.Dy + uz * cx.Dz, dD = dx * cx.Dx + dy * cx.Dy + dz * cx.Dz;
  const t = (c * dD - uD) / det, s = (c * uD - dD) / det;
  return Math.min(Math.min(t - 1e-6, s - 1e-6), cx.maxDist - Math.abs(t));
}

// 共面闭式 (通用 u): rhs = C/mag, φ; 根 ψ = asin(rhs)−φ, π−asin(rhs)−φ
function psiRoots(cx, u, d) {
  const c = cross(u, d);
  const A = cx.Dy * c[0] - cx.Dx * c[1], B = cx.Dx * c[0] + cx.Dy * c[1], C = -cx.Dz * c[2];
  const mag = Math.sqrt(A * A + B * B);
  if (mag < 1e-12) { const collapse = Math.abs(C) < 1e-12; return { rhs: collapse ? 0 : 1e9, phi: 0, collapse }; }
  return { rhs: C / mag, phi: Math.atan2(B, A), collapse: false };
}

// 一维闭式构造: 族 kind (0 入射掠射 / 1 镜面掠射 mi / 2 出射临界) 在参数 τ 处的边界态
function evalState(cx, kind, mi, tau) {
  const st = { tau, mask: -1, rootM: NaN, pA: null, pB: null, pm: NaN, k: 0 };
  const fp = cx.faces, last = fp.length - 1;
  const a = cx.canon[fp[0]];
  const ct = Math.cos(tau), sn = Math.sin(tau);
  let U;
  let consMin = 1;
  if (kind === 0) {
    const e1 = cx.basE[0], e2 = cx.basE[1];
    U = [ct * e1[0] + sn * e2[0], ct * e1[1] + sn * e2[1], ct * e1[2] + sn * e2[2]];
  } else {
    let g;
    if (kind === 2) {
      const e1 = cx.basT[0], e2 = cx.basT[1], bN = cx.canon[fp[last]];
      g = [cx.K * bN[0] + ct * e1[0] + sn * e2[0], cx.K * bN[1] + ct * e1[1] + sn * e2[1], cx.K * bN[2] + ct * e1[2] + sn * e2[2]];
      for (let j = last - 1; j >= 1; j--) g = reflect(g, cx.canon[fp[j]]);
    } else {
      const e1 = cx.basM[mi][0], e2 = cx.basM[mi][1];
      g = [cx.N * (ct * e1[0] + sn * e2[0]), cx.N * (ct * e1[1] + sn * e2[1]), cx.N * (ct * e1[2] + sn * e2[2])];
      for (let j = mi - 1; j >= 1; j--) g = reflect(g, cx.canon[fp[j]]);
    }
    const ga = dot(g, a);
    const gt = [g[0] - ga * a[0], g[1] - ga * a[1], g[2] - ga * a[2]];
    const gt2 = dot(gt, gt);
    consMin = Math.min(-ga - 1e-12, 1 - gt2);
    if (consMin <= 0) { st.mask = consMin; return st; }
    const ua = -Math.sqrt(Math.max(0, 1 - gt2));
    U = [gt[0] + ua * a[0], gt[1] + ua * a[1], gt[2] + ua * a[2]];
  }
  const m = new Float64Array(fp.length + 1);          // [入射, 镜1..镜k, 出面cosI2, TIR]
  const d = formalChain(cx, U, m);
  if (d === null) { st.mask = -1; return st; }
  const own = kind === 0 ? 0 : (kind === 2 ? m.length - 1 : mi);
  let mask = consMin;
  for (let i = 0; i < m.length; i++) { if (i === own) continue; if (m[i] < mask) mask = m[i]; }
  st.mask = mask;
  const pr = psiRoots(cx, U, d);
  st.rootM = pr.collapse ? 1 : (1 - Math.abs(pr.rhs));
  if (!pr.collapse && st.rootM >= -1e-12) {
    const rh = Math.max(-1, Math.min(1, pr.rhs));
    const asr = Math.asin(rh);
    let P = conePoint(cx, U, d, asr - pr.phi); if (P) st.pA = P;
    P = conePoint(cx, U, d, Math.PI - asr - pr.phi); if (P) st.pB = P;
    st.pm = Math.max(coneMargin(cx, U, d, asr - pr.phi), coneMargin(cx, U, d, Math.PI - asr - pr.phi));
  }
  return st;
}

// 事件二分: f 在 [t0,t1] 变号 → 精确事件 τ*, 返回严格 f>0 一侧的 state
function bisectEvent(evalS, f, t0, t1, f0) {
  let a, b;
  if (f0 > 0) { a = t0; b = t1; } else { a = t1; b = t0; }
  for (let i = 0; i < 48; i++) {
    const mid = 0.5 * (a + b);
    const fm = f(evalS(mid));
    if (fm > 0) a = mid; else b = mid;
  }
  return evalS(a);
}

// 稳定排序 (按 τ, 同 τ 保持原次序) — 与 C# 的 idx.Sort((x,y) => τ 比较, 再比下标) 一致
function stableSortByTau(arr) {
  const idx = arr.map((_, i) => i);
  idx.sort((x, y) => { const c = arr[x].tau - arr[y].tau; return c !== 0 ? (c < 0 ? -1 : 1) : x - y; });
  return idx.map(i => arr[i]);
}

// 构造一族边界: 基网格 → 窄窗口搜索 → 事件精化 → 掩码过滤 → 拉伸密化 → 两片就近配对装配
function buildFamily(cx, kind, mi, NT) {
  const evalS = t => evalState(cx, kind, mi, t);

  // 1) 基网格
  let samples = [];
  for (let i = 0; i <= NT; i++) samples.push(evalS(2 * Math.PI * i / NT));

  // 1b) 窄窗口搜索 (v4.1)
  {
    const score = s => Math.min(s.mask, Math.min(fin(s.rootM) ? s.rootM : -1, fin(s.pm) ? s.pm : -1));
    const extra = [];
    const NS = 8;
    const gr = (Math.sqrt(5) - 1) / 2;
    for (let i = 0; i < NT; i++) {
      const s0 = samples[i], s1 = samples[i + 1];
      if (score(s0) >= 0 || score(s1) >= 0) continue;
      const subs = [];
      let bestV = -Infinity;
      for (let k = 1; k < NS; k++) {
        const st = evalS(s0.tau + (s1.tau - s0.tau) * k / NS);
        subs.push(st);
        const v = score(st); if (v > bestV) bestV = v;
      }
      if (bestV >= 0) { extra.push(...subs); continue; }
      if (bestV <= Math.max(score(s0), score(s1))) continue;
      let a = s0.tau, b = s1.tau;
      let c = b - gr * (b - a), d = a + gr * (b - a);
      let fc = score(evalS(c)), fd = score(evalS(d));
      for (let it = 0; it < 40; it++) {
        if (fc > fd) { b = d; d = c; fd = fc; c = b - gr * (b - a); fc = score(evalS(c)); }
        else { a = c; c = d; fc = fd; d = a + gr * (b - a); fd = score(evalS(d)); }
      }
      const tp = 0.5 * (a + b);
      const sp = evalS(tp);
      if (score(sp) < 0) continue;
      const edge = (tOk, tBad) => {
        let lo = tOk, hi = tBad;
        for (let q = 0; q < 48; q++) { const mm = 0.5 * (lo + hi); if (score(evalS(mm)) >= 0) lo = mm; else hi = mm; }
        return lo;
      };
      const tL = edge(tp, s0.tau), tR = edge(tp, s1.tau);
      extra.push(evalS(tL)); extra.push(evalS(tL + (tp - tL) * 0.5)); extra.push(sp);
      extra.push(evalS(tp + (tR - tp) * 0.5)); extra.push(evalS(tR));
    }
    if (extra.length > 0) samples = stableSortByTau(samples.concat(extra));
  }

  // 2) 事件精化 (递归 ≤4 层): 角点 / 转折 / 捏合·截断 / 通用"有点↔无点"翻转
  const okS = s => s.mask >= -1e-9 && (s.pA !== null || s.pB !== null);
  const pinchable = s => s.mask >= -1e-9 && fin(s.rootM) && s.rootM >= -1e-12;
  const rk = (k, s) => { s.k = k; return s; };       // 0 转折 > 1 捏合 > 2 角点 > 3 通用
  const findEvents = (s0, s1) => {
    const evs = [];
    if (fin(s0.mask) && fin(s1.mask) && (s0.mask > 0) !== (s1.mask > 0)
      && Math.max(Math.abs(s0.mask), Math.abs(s1.mask)) > 1e-12)
      evs.push(rk(2, bisectEvent(evalS, s => s.mask, s0.tau, s1.tau, s0.mask)));            // 角点
    if (fin(s0.rootM) && fin(s1.rootM) && (s0.rootM > 0) !== (s1.rootM > 0)
      && Math.max(Math.abs(s0.rootM), Math.abs(s1.rootM)) > 1e-12) {
      const ev = rk(0, bisectEvent(evalS, s => s.rootM, s0.tau, s1.tau, s0.rootM));         // ψ 转折
      if (ev.pA !== null && ev.pB !== null) {
        const mpt = [0.5 * (ev.pA[0] + ev.pB[0]), 0.5 * (ev.pA[1] + ev.pB[1]), 0.5 * (ev.pA[2] + ev.pB[2])];
        ev.pA = mpt; ev.pB = mpt;
      } else if (ev.pA !== null && ev.pB === null) ev.pB = ev.pA;
      else if (ev.pB !== null && ev.pA === null) ev.pA = ev.pB;
      evs.push(ev);
    }
    if (pinchable(s0) && pinchable(s1)) {                                                   // 捏合/截断
      if ((s0.pA !== null) !== (s1.pA !== null))
        evs.push(rk(1, bisectEvent(evalS, s => s.pA !== null ? 1 : -1, s0.tau, s1.tau, s0.pA !== null ? 1 : -1)));
      if ((s0.pB !== null) !== (s1.pB !== null))
        evs.push(rk(1, bisectEvent(evalS, s => s.pB !== null ? 1 : -1, s0.tau, s1.tau, s0.pB !== null ? 1 : -1)));
    }
    if (okS(s0) !== okS(s1))
      evs.push(rk(3, bisectEvent(evalS, s => okS(s) ? 1 : -1, s0.tau, s1.tau, okS(s0) ? 1 : -1)));
    // 去重 (TOL=1e-7 rad): 优先留有点的, 再按类型等级
    const TOL = 1e-7;
    const order = stableSortByTau(evs);
    const kept = [];
    const better = (e, l) => { const oe = okS(e), ol = okS(l); if (oe !== ol) return oe; return e.k < l.k; };
    for (const e of order) {
      if (Math.abs(e.tau - s0.tau) < TOL || Math.abs(e.tau - s1.tau) < TOL) continue;
      if (kept.length > 0 && Math.abs(e.tau - kept[kept.length - 1].tau) < TOL) {
        if (better(e, kept[kept.length - 1])) kept[kept.length - 1] = e;
        continue;
      }
      kept.push(e);
    }
    return kept;
  };
  const refine = (s0, s1, dep, acc) => {
    const evs = findEvents(s0, s1);
    if (evs.length === 0) return;
    if (dep >= 4) { acc.push(...evs); return; }
    let prev = s0;
    for (const e of evs) { refine(prev, e, dep + 1, acc); acc.push(e); prev = e; }
    refine(prev, s1, dep + 1, acc);
  };
  const outL = [];
  for (let i = 0; i < samples.length; i++) {
    outL.push(samples[i]);
    if (i === samples.length - 1) break;
    refine(samples[i], samples[i + 1], 0, outL);
  }

  // 掩码过滤
  const vA = s => (s.mask >= -1e-9 && s.pA !== null) ? s.pA : null;
  const vB = s => (s.mask >= -1e-9 && s.pB !== null) ? s.pB : null;

  // 3) 拉伸密化 (真无效才断)
  const Lden = Math.max(0.15, cx.maxDist / 160), Lden2 = Lden * Lden;
  const dense = [];
  let budget = 5000;
  for (let i = 0; i < outL.length; i++) {
    dense.push(outL[i]);
    if (i === outL.length - 1) break;
    const stack = [[outL[i], outL[i + 1], 0]];
    const seg = [];
    while (stack.length > 0 && budget > 0) {
      const [sa, sb, dep] = stack.pop();
      const a0 = vA(sa), a1 = vA(sb), b0 = vB(sa), b1 = vB(sb);
      const needA = a0 !== null && a1 !== null && dist2(a0, a1) > Lden2;
      const needB = b0 !== null && b1 !== null && dist2(b0, b1) > Lden2;
      if ((!needA && !needB) || dep >= 11) continue;
      const sm = evalS(0.5 * (sa.tau + sb.tau));
      budget--;
      seg.push(sm);
      stack.push([sa, sm, dep + 1]);
      stack.push([sm, sb, dep + 1]);
    }
    // C# List.Sort 是不稳定排序, 但 τ 互不相同 (二分中点), 次序唯一
    seg.sort((x, y) => x.tau - y.tau);
    dense.push(...seg);
  }

  // 4) 两片就近配对 + 装配
  const strips = [];
  let sA = [], sB = [];
  let pA0 = null, pB0 = null;
  const brk = Math.max(4, 0.0625 * cx.maxDist * cx.maxDist);
  const flushA = () => { if (sA.length >= 6) strips.push(Float32Array.from(sA)); sA = []; };
  const flushB = () => { if (sB.length >= 6) strips.push(Float32Array.from(sB)); sB = []; };
  for (const st of dense) {
    let pa = vA(st), pb = vB(st);
    if (pA0 !== null && pB0 !== null && pa !== null && pb !== null) {
      if (dist2(pa, pB0) + dist2(pb, pA0) < dist2(pa, pA0) + dist2(pb, pB0)) { const tmp = pa; pa = pb; pb = tmp; }
    }
    if (pa === null) flushA();
    else {
      if (pA0 !== null && dist2(pa, pA0) > brk) flushA();
      sA.push(f32(pa[0]), f32(pa[1]), f32(pa[2]));
    }
    if (pb === null) flushB();
    else {
      if (pB0 !== null && dist2(pb, pB0) > brk) flushB();
      sB.push(f32(pb[0]), f32(pb[1]), f32(pb[2]));
    }
    pA0 = pa; pB0 = pb;
  }
  flushA(); flushB();
  return strips;
}

/**
 * 光锥法物理边界主入口 (ComputeLightConeBoundaries v4): 入射掠射 → entry (红),
 * 镜面掠射 (每个中间反射面一族) → mirror (橙), 出射临界 → tir (蓝). 每条为 Float32Array (x,y,z 连排).
 */
export function lightConeBoundaries(facePath, lampElev, lampDist, eyeHeight, refN, maxDist, zetaDeg, thetaDeg, NT = 900) {
  const elevRad = lampElev * Math.PI / 180.0;
  const Lx = lampDist * Math.cos(elevRad), Lz = lampDist * Math.sin(elevRad);
  const cx = {
    faces: facePath, N: refN, K: Math.sqrt(Math.max(0, refN * refN - 1)),
    L: [Lx, 0, Lz], Dx: Lx, Dy: 0, Dz: Lz - eyeHeight, maxDist,
    canon: buildCanonNormals(zetaDeg, thetaDeg),
  };
  const last = facePath.length - 1;
  cx.basE = orthoBasis(cx.canon[facePath[0]]);
  cx.basT = orthoBasis(cx.canon[facePath[last]]);
  cx.basM = new Array(facePath.length);
  for (let mi = 1; mi < last; mi++) cx.basM[mi] = orthoBasis(cx.canon[facePath[mi]]);
  const entry = buildFamily(cx, 0, 0, NT);
  const mirror = [];
  for (let mi = 1; mi < last; mi++) mirror.push(...buildFamily(cx, 1, mi, NT));
  const tir = buildFamily(cx, 2, 0, NT);
  return { entry, mirror, tir };
}
