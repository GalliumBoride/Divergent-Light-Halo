// 二自由度 (柱晶 Column: α+θ 自由 / Lowitz: α+ζ 自由) 发散光晕面 —
//   (1) Divergentlightform.cs 的 ComputeHaloSurfaceLightCone2Dof / SliceRuns2Dof / RunsFromGrid / ConeMapRuns:
//       每个锥仰角 h 把灯当"太阳", 用光点法在取向平面上追 J=0 焦散 (或选定的一条平行光曲线),
//       焦散点的平行光出射方向 d = −D 代入与一自由度相同的共面闭式 ψ 定位晶体 P = L + t·û;
//   (2) HaloEndpointEngine.cs 的发散光部分 (ComputeDivergentCausticEdges / ComputeDivergentEdges / FenceRuns …):
//       曲面边界线 = 平行光末端点身份的 h-延拓 × 锥闭式.
// 与 C# 逐句对应 (对账见 tools/).

import { LP, OMode, setSun, PathCtx, Chain3, jana3, validAt, axesMeta, traceContours, allFams, famEval3, SqK } from './lp.js';
import { PathC, FamC, Scene, polishPair, fenceSamples, cuspState, minMarginUV, computeEndpoints } from './ep.js';
import { tryConePoint } from './cone.js';

const f32 = Math.fround;
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const nrm = a => Math.sqrt(dot(a, a));
const scl = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];

// ##################################################################  HaloEP 发散光部分

// 转折专用焦散泛函: 可达圆上 G=0 ⇒ Λ ∝ (ê×w), 光滑多项式
function jturnW(sc, w) {
  const { Na, Nb, qE, TE } = sc.P.naNb(w);
  const lam = cross(sc.F.Ec, w);
  return TE * dot(Na, lam) - qE * dot(Nb, lam);
}

function lightScene(reference, fam) {
  const sc = Object.create(Scene.prototype);
  sc.P = reference.P; sc.F = fam; sc.PC = reference.PC;
  sc.Geo = reference.Geo; sc.Vacuous = reference.Vacuous;
  sc.Reach = fam.reachFences();
  sc.DegenCaustic = [];
  return sc;
}

function skyOfW(sc, w, eps) {
  const E = sc.F.lift(w, eps);
  if (E === null) return null;
  const [e1, e2, e3] = E;
  const Dc = sc.P.chainW(w).Dc;
  const D = [
    e1[0] * Dc[0] + e2[0] * Dc[1] + e3[0] * Dc[2],
    e1[1] * Dc[0] + e2[1] * Dc[1] + e3[1] * Dc[2],
    e1[2] * Dc[0] + e2[2] * Dc[1] + e3[2] * Dc[2]];
  const nn = nrm(D);
  if (Number.isNaN(nn) || nn < 1e-12) return null;
  return scl(D, 1.0 / nn);
}

function solveEndW(sc, fen, eps, seed, turn) {
  const F2 = turn ? (x => jturnW(sc, x)) : (x => sc.jhatW(x, eps));
  const w = polishPair(x => sc.P.fval(fen, x), x => sc.P.fgrad(fen, x), F2, null, seed);
  if (w === null) return null;
  if (Math.abs(sc.P.fval(fen, w)) > 1e-8 || Math.abs(F2(w)) > 1e-8) return null;
  if (sc.F.gram(w) < -1e-9) return null;
  if (sc.minMargin(w) < -1e-8) return null;
  return w;
}

function solveCuspUV(sc, seed) {
  let u = seed[0], v = seed[1];
  for (let it = 0; it < 30; it++) {
    const st = cuspState(sc, u, v);
    const r = [st.J0, st.V[0], st.V[1], st.V[2]];
    if (r.some(x => !Number.isFinite(x))) return null;
    if (Math.max(...r.map(Math.abs)) < 1e-11) break;
    const e2 = 1e-6;
    const Jc = [[0, 0], [0, 0], [0, 0], [0, 0]];
    for (let k = 0; k < 2; k++) {
      const du = k === 0 ? e2 : 0, dv = k === 0 ? 0 : e2;
      const p = cuspState(sc, u + du, v + dv);
      const q = cuspState(sc, u - du, v - dv);
      Jc[0][k] = (p.J0 - q.J0) / (2 * e2);
      for (let i = 0; i < 3; i++) Jc[i + 1][k] = (p.V[i] - q.V[i]) / (2 * e2);
    }
    let a11 = 0, a12 = 0, a22 = 0, b1 = 0, b2 = 0;
    for (let i = 0; i < 4; i++) {
      a11 += Jc[i][0] * Jc[i][0]; a12 += Jc[i][0] * Jc[i][1]; a22 += Jc[i][1] * Jc[i][1];
      b1 -= Jc[i][0] * r[i]; b2 -= Jc[i][1] * r[i];
    }
    const det = a11 * a22 - a12 * a12;
    if (Math.abs(det) < 1e-30) return null;
    let dxu = (b1 * a22 - b2 * a12) / det, dxv = (a11 * b2 - a12 * b1) / det;
    const nn = Math.sqrt(dxu * dxu + dxv * dxv);
    if (Number.isNaN(nn)) return null;
    if (nn > 0.2) { dxu *= 0.2 / nn; dxv *= 0.2 / nn; }
    u += dxu; v += dxv;
  }
  const st = cuspState(sc, u, v);
  if (Math.abs(st.J0) > 1e-9 || nrm(st.V) > 1e-6) return null;
  if (minMarginUV(sc, u, v) < -1e-7) return null;
  return [u, v];
}

// 锥闭式: 天球方向 (ψ=0 帧) → 两个 ψ 支的晶体位置; shadow = 锥影 (E4)
function coneBranches(Dsky, hRad, L, E, maxDist) {
  const DDv = sub(L, E);
  const d = scl(Dsky, -1.0);
  const i0 = [-Math.cos(hRad), 0.0, -Math.sin(hRad)];
  const vv = cross(i0, d);
  const A_ = DDv[1] * vv[0] - DDv[0] * vv[1];
  const B_ = DDv[0] * vv[0] + DDv[1] * vv[1];
  const C_ = -DDv[2] * vv[2];
  const Rm = Math.sqrt(A_ * A_ + B_ * B_);
  const outB = [null, null];
  if (Rm < 1e-15 || Rm < Math.abs(C_)) return { br: outB, shadow: true };
  const ph = Math.atan2(B_, A_);
  const bs = Math.asin(Math.max(-1, Math.min(1, C_ / Rm)));
  for (let b = 0; b < 2; b++) {
    const psi = b === 0 ? bs - ph : Math.PI - bs - ph;
    const cp = Math.cos(psi), sp = Math.sin(psi);
    const u = [cp * i0[0] - sp * i0[1], sp * i0[0] + cp * i0[1], i0[2]];
    const dd = [cp * d[0] - sp * d[1], sp * d[0] + cp * d[1], d[2]];
    const c = dot(u, dd), den = 1 - c * c;
    if (den < 1e-12) continue;
    const tR = (c * dot(dd, DDv) - dot(u, DDv)) / den;
    const sR = (c * dot(u, DDv) - dot(dd, DDv)) / den;
    if (tR < 1e-9 || sR < 1e-9 || tR > maxDist) continue;
    outB[b] = [L[0] + tR * u[0], L[1] + tR * u[1], L[2] + tR * u[2]];
  }
  return { br: outB, shadow: false };
}

function newTrack(kind, eps, label) {
  return { kind, label: label || '', eps, hs: [], sols: [], ds: [], pts3: [], nShadow: 0 };
}

// 天球方向 + 锥映射 → 3D 折线段 (两 ψ 支, 有效性断开; end/turn 用 w+lift, cusp 用 (u,v))
function projectTracks(tracks, m, fixedAng, path, L, E, maxDist) {
  const pctx = new PathCtx(path);
  const pc = new PathC(path, LP.partial);
  setSun(10);
  const refScene = new Scene(pc, FamC.make(m, fixedAng), pctx);
  for (const tk of tracks) {
    const segs = [[], []];
    const flush = b => { if (segs[b].length >= 6) tk.pts3.push(Float32Array.from(segs[b])); segs[b] = []; };
    for (let i = 0; i < tk.hs.length; i++) {
      const hd = tk.hs[i];
      setSun(hd);
      const sc = lightScene(refScene, FamC.make(m, fixedAng));
      let D = null;
      if (tk.kind === 'cusp') {
        const c3 = Chain3(tk.sols[i][0], tk.sols[i][1], m, fixedAng, pctx);
        const Dv = [c3.Dx, c3.Dy, c3.Dz];
        const nn = nrm(Dv);
        if (nn > 1e-12) D = scl(Dv, 1.0 / nn);
      } else D = skyOfW(sc, tk.sols[i], tk.eps);
      tk.ds.push(D);
      if (D === null) { flush(0); flush(1); continue; }
      const { br, shadow } = coneBranches(D, hd * Math.PI / 180, L, E, maxDist);
      if (shadow) tk.nShadow++;
      for (let b = 0; b < 2; b++) {
        if (br[b] === null) { flush(b); continue; }
        segs[b].push(f32(br[b][0]), f32(br[b][1]), f32(br[b][2]));
      }
    }
    flush(0); flush(1);
  }
}

// ------------------------------------------------------------------ 修正版 (显示用): 加密 + 端点二分 + 精确视方向
// 锥闭式, 每支给 { P, s: 观测者视方向 (= Rz(ψ)·D̂, 精确), u: 从灯出发的方向 }; 锥影时 null
// 锥闭式的折叠量 |C/R| − 1: > 0 时两支都不存在, = 0 处两支重合 (锥面上的折叠)
function coneFold(Dsky, hRad, L, E) {
  const DDv = sub(L, E), d = scl(Dsky, -1.0), i0 = [-Math.cos(hRad), 0.0, -Math.sin(hRad)], vv = cross(i0, d);
  const A_ = DDv[1] * vv[0] - DDv[0] * vv[1], B_ = DDv[0] * vv[0] + DDv[1] * vv[1], C_ = -DDv[2] * vv[2];
  const Rm = Math.sqrt(A_ * A_ + B_ * B_);
  return Rm < 1e-300 ? (Math.abs(C_) > 0 ? 1 : -1) : Math.abs(C_) / Rm - 1;
}
// 样本序列里, 相邻两点一个在折叠外 (两支都没有) 一个在内: 二分到折叠处插入一个点 (两支在此重合).
//   折叠附近一支的位置随参数呈 √ 变化, 可以在一个样本间隔里从折叠点一路走回灯 —— 那一支在两个样本上都无效,
//   不插这个点就整段漏掉 (层曲线在折叠处断开, 晕面边缘看起来接不上).
function insertFolds(samples, foldAt, mid, iters = 40) {
  const out = [];
  let pf = null;
  for (let i = 0; i < samples.length; i++) {
    const x = samples[i], f = foldAt(x);
    if (i > 0 && pf !== null && f !== null && (pf > 0) !== (f > 0)) {
      let lo = f > 0 ? samples[i - 1] : x, hi = f > 0 ? x : samples[i - 1];      // lo: 折叠内 (有解), hi: 折叠外
      for (let it = 0; it < iters; it++) {
        const m = mid(lo, hi);
        if (m === null) break;
        const fm = foldAt(m);
        if (fm === null) break;
        if (fm > 0) hi = m; else lo = m;
      }
      if (lo !== samples[i - 1] && lo !== x) out.push(lo);
    }
    out.push(x);
    pf = f;
  }
  return out;
}

function coneBranchesFull(Dsky, hRad, L, E, maxDist) {
  const DDv = sub(L, E);
  const d = scl(Dsky, -1.0);
  const i0 = [-Math.cos(hRad), 0.0, -Math.sin(hRad)];
  const vv = cross(i0, d);
  const A_ = DDv[1] * vv[0] - DDv[0] * vv[1], B_ = DDv[0] * vv[0] + DDv[1] * vv[1], C_ = -DDv[2] * vv[2];
  const Rm = Math.sqrt(A_ * A_ + B_ * B_);
  const out = [null, null];
  if (Rm < 1e-15 || Rm < Math.abs(C_)) return out;
  const ph = Math.atan2(B_, A_), bs = Math.asin(Math.max(-1, Math.min(1, C_ / Rm)));
  for (let b = 0; b < 2; b++) {
    const psi = b === 0 ? bs - ph : Math.PI - bs - ph;
    const cp = Math.cos(psi), sp = Math.sin(psi);
    const u = [cp * i0[0] - sp * i0[1], sp * i0[0] + cp * i0[1], i0[2]];
    const dd = [cp * d[0] - sp * d[1], sp * d[0] + cp * d[1], d[2]];
    // 交点分母用 |û×d|² (= 1 − c², 但接近平行时没有相消误差), 分子也写成叉积, 见 mesh1.js evalNode
    const w = cross(u, dd), den = dot(w, w);
    if (den < 1e-12) continue;
    const tR = -dot(cross(DDv, dd), w) / den, sR = dot(cross(DDv, u), w) / den;
    if (tR < 1e-9 || sR < 1e-9 || tR > maxDist) continue;
    out[b] = { P: [L[0] + tR * u[0], L[1] + tR * u[1], L[2] + tR * u[2]], s: [-dd[0], -dd[1], -dd[2]], u };
  }
  return out;
}
const angD = (a, b) => { const c = dot(a, b); return Math.acos(c > 1 ? 1 : c < -1 ? -1 : c) * 180 / Math.PI; };
const dist3 = (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
// 相邻两点是否还要加密: 三维长度 / 从眼看的张角 / 从灯看的张角
function tooFar(a, b, Lden) {
  return dist3(a.P, b.P) > Lden || angD(a.s, b.s) > 0.3 || angD(a.u, b.u) > 1.0;
}

// 点到线段 ab 的距离
function chordDev(p, a, b) {
  const ab = sub(b, a), ap = sub(p, a), l2 = dot(ab, ab);
  const t = l2 > 0 ? Math.max(0, Math.min(1, dot(ab, ap) / l2)) : 0;
  return Math.hypot(ap[0] - t * ab[0], ap[1] - t * ab[1], ap[2] - t * ab[2]);
}

// 一条参数化曲线 (参数 x: 标量或数组, 任意) 的一支 → 加密好的折线们. 通用于 h-层曲线与 h-延拓轨迹.
//   samples: 参数序列; mapAt(x) → 该支状态或 null; mid(xa, xb) → 中间参数 (可含校正) 或 null.
//   有效 ↔ 无效 的转折处二分 (iters 次) 找精确端点. 返回 [{ p, d }].
//   sag > 0 (边界线用): 段长/张角已够细时还要看中点离弦多远, 超过 sag 就再分 —— 弯得急的地方 (√ 区、尖点旁)
//   按长度分段的折线会从真曲线上偏出去, 层线端点明明落在曲线上, 却离画出来的折线一两厘米.
function refineBranch(samples, mapAt, mid, Lden, maxDepth = 9, iters = 40, sag = 0) {
  const out = [];
  let P = [], D = [];
  const push = st => { P.push(st.P[0], st.P[1], st.P[2]); D.push(st.s[0], st.s[1], st.s[2]); };
  // 整条落在一个点上的不要: 退化纤维 (一圈取向的天空方向都一样, J 沿它恒为 0) 映成一个点, 是点状焦散, 画不成线
  const extent = () => { let e = 0; for (let i = 3; i < P.length; i += 3) e = Math.max(e, Math.abs(P[i] - P[0]) + Math.abs(P[i + 1] - P[1]) + Math.abs(P[i + 2] - P[2])); return e; };
  const flush = () => { if (P.length >= 6 && extent() > 1e-6) out.push({ p: Float32Array.from(P), d: Float32Array.from(D) }); P = []; D = []; };
  // 有效 xa 与无效 xb 之间的转折点: 返回最后一个有效的 { x, s }
  const edge = (xa, sa, xb) => {
    let lo = xa, slo = sa, hi = xb;
    for (let k = 0; k < iters; k++) {
      const xm = mid(lo, hi);
      if (xm === null) break;
      const sm = mapAt(xm);
      if (sm) { lo = xm; slo = sm; } else hi = xm;
    }
    return { x: lo, s: slo };
  };
  // 两端都有效的一段: 递归加密 (不推入起点, 推入终点)
  const dens = (xa, sa, xb, sb, dep) => {
    const far = tooFar(sa, sb, Lden);
    if (dep >= maxDepth || (!far && !(sag > 0))) { push(sb); return; }
    const xm = mid(xa, xb);
    if (xm === null) { push(sb); return; }
    const sm = mapAt(xm);
    if (!sm) {                                    // 中间无效: 这一段其实断开了
      const e1 = edge(xa, sa, xm);
      if (e1.x !== xa) dens(xa, sa, e1.x, e1.s, dep + 1);
      flush();
      const e2 = edge(xb, sb, xm);
      push(e2.s);
      if (e2.x !== xb) dens(e2.x, e2.s, xb, sb, dep + 1);
      return;
    }
    if (!far && chordDev(sm.P, sa.P, sb.P) <= sag) { push(sb); return; }
    dens(xa, sa, xm, sm, dep + 1);
    dens(xm, sm, xb, sb, dep + 1);
  };
  let px = null, ps = null;
  for (const x of samples) {
    const s = mapAt(x);
    if (px !== null) {
      if (ps && s) dens(px, ps, x, s, 0);
      else if (ps && !s) { const e = edge(px, ps, x); if (e.x !== px) dens(px, ps, e.x, e.s, 0); flush(); }
      else if (!ps && s) { const e = edge(x, s, px); push(e.s); if (e.x !== x) dens(e.x, e.s, x, s, 0); }
    } else if (s) push(s);
    px = x; ps = s;
  }
  flush();
  return out;
}

// 修正版轨迹投影: 每条轨迹 (h 样本 + 解) 两支各自加密, 中间 h 用 tk.slv 以邻点为种子重解
function projectTracksFine(tracks, m, fixedAng, path, L, E, maxDist) {
  const pctx = new PathCtx(path);
  const pc = new PathC(path, LP.partial);
  setSun(10);
  const refScene = new Scene(pc, FamC.make(m, fixedAng), pctx);
  const Lden = Math.max(1e-3, Math.min(maxDist, dist3(L, E)) / 150), sag = Math.max(2e-4, 1e-4 * dist3(L, E));
  const hSing = Math.atan2(L[2] - E[2], Math.hypot(L[0] - E[0], L[1] - E[1])) * 180 / Math.PI;
  for (const tk of tracks) {
    tk.pts3 = [];
    // 参数 = { h, sol }; 解失败 → null
    const skyAt = (h, sol) => {
      setSun(h);
      if (tk.kind === 'cusp') {
        const c3 = Chain3(sol[0], sol[1], m, fixedAng, pctx);
        const Dv = [c3.Dx, c3.Dy, c3.Dz], nn = nrm(Dv);
        return nn > 1e-12 ? scl(Dv, 1.0 / nn) : null;
      }
      return skyOfW(lightScene(refScene, FamC.make(m, fixedAng)), sol, tk.eps);
    };
    // 灯仰角 (从眼看) 附近锥映射的有效性变得极快 (晶体贴到灯或眼上), 有效窗口可能比 0.5° 的 h 步还窄, 两端都无效就整段漏掉:
    // 距它 1° 以内的 h 间隔再插 8 个样本 (以邻点为种子重解)
    // 弧长延拓的轨迹 (tk.ys): 中间点取在弦上该处垂直于弦的超平面与曲线的交点 (h 不单调、跨折叠也对)
    const RAD = Math.PI / 180, arc = !!tk.ys;
    const fromY = y => ({ h: y[0] / RAD, sol: y.slice(1), y });
    const samples = [];
    for (let i = 0; i < tk.hs.length; i++) {
      samples.push(arc ? fromY(tk.ys[i]) : { h: tk.hs[i], sol: tk.sols[i] });
      if (i + 1 < tk.hs.length && tk.slv) {
        const ha = tk.hs[i], hb = tk.hs[i + 1];
        if (Math.min(Math.abs(ha - hSing), Math.abs(hb - hSing)) < 1.0 || (ha - hSing) * (hb - hSing) < 0) {
          let seed = tk.sols[i];
          for (let k = 1; k < 8; k++) {
            if (arc) { const y = arcMid(tk.arcF, tk.ys[i], tk.ys[i + 1], k / 8, tk.arcOpt); if (y) samples.push(fromY(y)); continue; }
            const hk = ha + (hb - ha) * k / 8, sk = tk.slv(hk, seed);
            if (sk) { samples.push({ h: hk, sol: sk }); seed = sk; }
          }
        }
      }
    }
    // 超平面校正解不出时退回定 h 求解 (两端 h 不同 ⟹ 中间没有折叠)
    const mid = arc ? (a, b) => {
      const y = arcMid(tk.arcF, a.y, b.y, 0.5, tk.arcOpt);
      if (y) return fromY(y);
      if (a.h === b.h || !tk.slv) return null;
      const h = 0.5 * (a.h + b.h), sol = tk.slv(h, a.sol);
      return sol ? fromY([h * RAD, ...sol]) : null;
    } : (a, b) => {
      const h = 0.5 * (a.h + b.h);
      const sol = tk.slv ? tk.slv(h, a.sol) : null;
      return sol ? { h, sol } : null;
    };
    const cache = new Map();
    const branches = x => {
      let v = cache.get(x);
      if (v === undefined) {
        const D = skyAt(x.h, x.sol);
        v = D ? coneBranchesFull(D, x.h * Math.PI / 180, L, E, maxDist) : [null, null];
        cache.set(x, v);
      }
      return v;
    };
    const foldAt = x => { const D = skyAt(x.h, x.sol); return D ? coneFold(D, x.h * Math.PI / 180, L, E) : null; };
    // 轨迹上每个中点都要重解一次末端点 (Newton), 二分 24 次已是 0.5° 步长的 2⁻²⁴
    const smp = insertFolds(samples, foldAt, mid, 24);
    for (let b = 0; b < 2; b++) tk.pts3.push(...refineBranch(smp, x => branches(x)[b], mid, Lden, 9, 24, sag));
  }
}

// ------------------------------------------------------------------ 伪弧长延拓 (修正版轨迹用)
// 边界线 = (h, 取向) 空间里的一条曲线 F(y) = 0 (y ∈ Rⁿ, n−1 个方程, y[0] = h 的弧度). 按 h 一步步走 (自然参数延拓)
//   到 h 方向的折叠点 (两个解汇合、dh/ds = 0) 就走不下去: 折叠另一侧那一支丢了, 停下的端点也只准到 √(h 的二分精度).
//   这里把 h 当未知数之一, 沿弧长走 (Keller): 预测 = 沿切线 (雅可比的零向量) 走一步; 校正 = Newton 解
//   [F(y) = 0; t·(y − y预测) = 0], 即拉回曲线且落在过预测点、垂直切线的超平面上. 折叠处切线垂直于 h 轴, 照走不误, h 自然折返.
//   步长按校正量与转角自适应; 出界 (opt.valid 不过 / F 无定义) 就减半, 减到 sMin 为止 (= 精确到 sMin 的端点).

// n×n 线性方程 (部分主元消去), 奇异返回 null
function solveN(A, b) {
  const n = b.length, M = A.map((r, i) => [...r, b[i]]);
  for (let c = 0; c < n; c++) {
    let pr = c;
    for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[pr][c])) pr = r;
    if (!(Math.abs(M[pr][c]) > 1e-300)) return null;
    [M[c], M[pr]] = [M[pr], M[c]];
    for (let r = c + 1; r < n; r++) { const f = M[r][c] / M[c][c]; for (let k = c; k <= n; k++) M[r][k] -= f * M[c][k]; }
  }
  const x = new Array(n);
  for (let r = n - 1; r >= 0; r--) { let v = M[r][n]; for (let k = r + 1; k < n; k++) v -= M[r][k] * x[k]; x[r] = v / M[r][r]; }
  return x.every(Number.isFinite) ? x : null;
}
const dotN = (a, b) => { let v = 0; for (let i = 0; i < a.length; i++) v += a[i] * b[i]; return v; };
const nrmN = a => Math.sqrt(dotN(a, a));
// 差分雅可比 ((n−1)×n); 一侧无定义 (贴着有效区边) 就用另一侧. fwd: 前向差分 (n 次求值而不是 2n; 光滑的 F 够用)
function jacN(F, y, f0, e, fwd = false) {
  const n = y.length, J = f0.map(() => new Array(n));
  for (let j = 0; j < n; j++) {
    const yp = y.slice(), ym = y.slice(); yp[j] += e; ym[j] -= e;
    const fp = F(yp), fm = fwd && fp ? null : F(ym);
    for (let i = 0; i < f0.length; i++) {
      if (fp && fm) J[i][j] = (fp[i] - fm[i]) / (2 * e);
      else if (fp) J[i][j] = (fp[i] - f0[i]) / e;
      else if (fm) J[i][j] = (f0[i] - fm[i]) / e;
      else return null;
    }
  }
  return J;
}
// (n−1)×n 矩阵的零向量 (广义叉积: 第 j 分量 = (−1)^j · 去掉第 j 列的行列式), 单位化
function nullVec(J) {
  const n = J[0].length, det = M => {
    if (M.length === 1) return M[0][0];
    let d = 0;
    for (let c = 0; c < M.length; c++) d += (c % 2 ? -1 : 1) * M[0][c] * det(M.slice(1).map(r => r.filter((_, k) => k !== c)));
    return d;
  };
  const t = [];
  for (let j = 0; j < n; j++) t.push((j % 2 ? -1 : 1) * det(J.map(r => r.filter((_, k) => k !== j))));
  const l = nrmN(t);
  return l > 0 && Number.isFinite(l) ? t.map(v => v / l) : null;
}
// 校正: 从 y0 出发解 [F(y) = 0; nv·y = c]. 返回 { y, t (切线) } 或 null.
//   y0 取曲线上已有的点 (而不是预测点/弦中点): 曲线贴着有效区的弯边走时 (外缘就在出射临界边内 1e-6), 切线上的预测点
//   或弦中点会落到边外 (F 无定义); 从曲线点出发, Newton 第一步沿切线移向超平面, 出界还能回溯
//   opt.stall: 残差 (以容差为单位) 连续两次迭代没降就放弃 —— 预测点离曲线太远时 Newton 方向一再指向域外, 回溯把步长压到 0
//   原地空转, 不如早点失败让外层减步长 (末端点轨迹上这类空转曾占延拓时间的大半). 外缘贴着出射临界边 (√ 病态), Newton 本来就
//   收敛得慢, 不用这条.
function arcSolve(F, y0, nv, opt, c = dotN(nv, y0)) {
  let y = y0.slice(), r1 = Infinity, r2 = Infinity;
  for (let it = 0; it < 30; it++) {
    const f = F(y);
    if (!f) return null;
    const g = dotN(nv, y) - c, rn = Math.max(...f.map((v, i) => Math.abs(v) / opt.tol[i]), Math.abs(g) / 1e-12);
    if (opt.stall && it >= 2 && rn > 0.9 * r2) return null;
    r2 = r1; r1 = rn;
    const J = jacN(F, y, f, opt.eps, opt.fwd);
    if (!J) return null;
    if (rn < 1) {
      const t = nullVec(J);
      return t && opt.valid(y) ? { y, t } : null;
    }
    let d = solveN([...J, nv], [...f.map(v => -v), -g]);
    if (!d) return null;
    const dl = nrmN(d);
    if (dl > opt.maxCorr) d = d.map(v => v * opt.maxCorr / dl);
    // 回溯: 一步跨出 F 的定义域 (出射临界边外 Ĵ 无定义) 就减半
    let bt = 0;
    for (; bt < 40 && !F(y.map((v, i) => v + d[i])); bt++) d = d.map(v => 0.5 * v);
    if (bt === 40) return null;
    y = y.map((v, i) => v + d[i]);
  }
  return null;
}
// 周期坐标 (opt.wrap 里的下标, 周期 2π) 折回后的差 a − b
function arcDiff(a, b, wrap) {
  const d = a.map((v, i) => v - b[i]);
  if (wrap) for (const i of wrap) d[i] -= 2 * Math.PI * Math.round(d[i] / (2 * Math.PI));
  return d;
}
// 从 y0 (切线 t0) 沿弧长走到头: 返回点列 (不含 y0); closed = 绕回了 y0 (周期坐标按 2π 折回比较 ——
//   取向角转一整圈回到同一取向, 不比折回就一圈圈走到步数上限)
function arcTrace(F, y0, t0, opt) {
  const out = [];
  let y = y0, t = t0, s = opt.s0, len = 0;
  for (let n = 0; n < opt.maxN; n++) {
    if (opt.cancelled && opt.cancelled()) break;
    const yp = y.map((v, i) => v + s * t[i]);
    const r = arcSolve(F, y, t, opt, dotN(t, yp));
    let tn = null;
    if (r) {
      tn = dotN(r.t, t) < 0 ? r.t.map(v => -v) : r.t;
      const dc = nrmN(r.y.map((v, i) => v - yp[i]));
      if (!(dc < 0.5 * s) || dotN(tn, t) < Math.cos(15 * Math.PI / 180)) tn = null;   // 校正跑远 / 转角过大: 可能跳到了别的分支
    }
    if (!tn) { if ((s *= 0.5) < opt.sMin) break; continue; }
    out.push(r.y);
    len += nrmN(r.y.map((v, i) => v - y[i]));
    y = r.y; t = tn;
    if (opt.stop && opt.stop(y)) break;
    if (len > 4 * opt.sMax && nrmN(arcDiff(y, y0, opt.wrap)) < 1.5 * s && dotN(t, t0) > 0) { out.closed = true; break; }
    s = Math.min(opt.sMax, 1.5 * s);
  }
  return out;
}
// 两个样本之间曲线上的点 (frac 处): 过弦上该点、垂直于弦的超平面与曲线的交点 —— 折叠两侧的样本之间也取得对
function arcMid(F, ya, yb, frac, opt) {
  const ch = yb.map((v, i) => v - ya[i]), l = nrmN(ch);
  if (!(l > 0)) return null;
  const nv = ch.map(v => v / l), c = dotN(nv, ya) + frac * l;
  // 解须落在两样本之间 (到两端都比弦短), 否则是 Newton 收敛到了别的分支
  const between = r => r && nrmN(r.y.map((v, i) => v - ya[i])) < l && nrmN(r.y.map((v, i) => v - yb[i])) < l ? r.y : null;
  return between(arcSolve(F, frac <= 0.5 ? ya : yb, nv, opt, c)) || between(arcSolve(F, frac <= 0.5 ? yb : ya, nv, opt, c));
}
// 一条轨迹: 从 y0 两个方向各走到头, 拼成一列 (h 不必单调)
function arcBoth(F, y0, opt) {
  const r0 = arcSolve(F, y0, [1, ...new Array(y0.length - 1).fill(0)], opt);
  if (!r0) return null;
  const fw = arcTrace(F, r0.y, r0.t, opt);
  if (fw.closed) {                                  // 补上闭合点 (起点按周期平移到末点那一圈), 中点插值才连续
    const last = fw[fw.length - 1], d = arcDiff(last, r0.y, opt.wrap);
    return { ys: [r0.y, ...fw, last.map((v, i) => v - d[i])], closed: true };
  }
  const bw = arcTrace(F, r0.y, r0.t.map(v => -v), opt);
  return { ys: [...bw.reverse(), r0.y, ...fw], closed: false };
}

const DIV_SEEDS = [3, 10, 18, 26, 29.4, 36, 50, 62, 68.5, 75, 82];
const DIV_CUSP_SEEDS = [10, 26, 45, 70];
const fmtEps = e => (e >= 0 ? '+' : '-') + Math.abs(Math.round(e));

// h-延拓: 从 h0 出发两个方向各走 dh, 失败处二分 10 次补端点
function extendTrack(tk, h0, sd0, slv, dh, cancelled) {
  tk.slv = slv;                       // 修正版投影要在中间 h 上重解 (加密 / 端点二分)
  const b0 = slv(h0, sd0);
  if (b0 === null) return;
  tk.hs.push(h0); tk.sols.push(b0);
  for (let dir = 0; dir < 2; dir++) {
    const sgn = dir === 0 ? 1 : -1;
    let h = h0, sd = b0;
    while (true) {
      if (cancelled && cancelled()) return;
      const h2 = h + sgn * dh;
      if (Math.abs(h2) > 89.0) break;
      const r = slv(h2, sd);
      if (r === null) {
        let lo = h, hi = h2, rg = null;
        for (let bi = 0; bi < 10; bi++) {
          const mid = 0.5 * (lo + hi);
          const rm = slv(mid, sd);
          if (rm === null) hi = mid; else { lo = mid; rg = rm; sd = rm; }
        }
        if (rg !== null) {
          if (sgn > 0) { tk.hs.push(lo); tk.sols.push(rg); }
          else { tk.hs.unshift(lo); tk.sols.unshift(rg); }
        }
        break;
      }
      if (sgn > 0) { tk.hs.push(h2); tk.sols.push(r); }
      else { tk.hs.unshift(h2); tk.sols.unshift(r); }
      h = h2; sd = r;
    }
  }
}

/** 发散光焦散曲面边界线 (末端点身份 h-延拓 × 锥闭式). 会临时改 LP 太阳, 结束时还原. path = HaloLP 面号. */
export function computeDivergentCausticEdges(path, m, fixedAng, lampElevDeg, lampDist, eyeHeight, maxDist,
  dh = 0.5, withCusps = true, cancelled = null, fine = false) {
  const oldTH = LP.TH * 180 / Math.PI;
  const tracks = [];
  try {
    const se = lampElevDeg * Math.PI / 180;
    const L = [lampDist * Math.cos(se), 0.0, lampDist * Math.sin(se)];
    const E = [0.0, 0.0, eyeHeight];
    const fam0 = FamC.make(m, fixedAng);
    if (fam0 === null) return tracks;
    setSun(10);
    const pctx = new PathCtx(path);
    const pc = new PathC(path, LP.partial);
    const refScene = new Scene(pc, FamC.make(m, fixedAng), pctx);
    const sceneAt = hd => { setSun(hd); return lightScene(refScene, FamC.make(m, fixedAng)); };

    const onTrack = (probe, hd, sol) => {
      for (const tk of tracks) {
        if (probe !== null && tk.kind !== probe.kind) continue;
        if (probe !== null && tk.kind !== 'cusp' && Math.abs(tk.eps - probe.eps) > 0.5) continue;
        let best = -1, bd = 0.8;
        for (let i = 0; i < tk.hs.length; i++)
          if (Math.abs(tk.hs[i] - hd) < bd) { bd = Math.abs(tk.hs[i] - hd); best = i; }
        if (best < 0) continue;
        const a = tk.sols[best];
        let same = a.length === sol.length;
        if (same)
          for (let i = 0; i < a.length; i++)
            if (Math.abs(a[i] - sol[i]) > 3e-2) { same = false; break; }
        if (same) return true;
      }
      return false;
    };

    const extend = (tk, h0, seedW, fen, rname) => {
      const slv = (hd, sd) => {
        const sc = sceneAt(hd);
        if (tk.kind === 'cusp') return solveCuspUV(sc, sd);
        const f = tk.kind === 'end' ? fen : (sc.Reach.find(x => x.name === rname) || null);
        if (f === null) return null;
        return solveEndW(sc, f, tk.eps, sd, tk.kind === 'turn');
      };
      if (fine && tk.kind !== 'cusp') arcExtend(tk, h0, seedW, fen, rname, slv);
      else extendTrack(tk, h0, seedW, slv, dh, cancelled);
    };
    // 修正版: end / turn 轨迹 = (h, w) 里的曲线 { 篱笆 = 0, Ĵ (或转折泛函) = 0, |w|² = 1 }, 伪弧长延拓 ——
    //   两个末端点在某个 h 汇合消失 (h 方向的折叠) 时绕过去接着走另一支, 不再停在折叠点、丢掉另一支
    const DEG = Math.PI / 180;
    // 场景按 h 缓存几份 (差分雅可比在 h ± ε 与 h 之间来回)
    const scCache = new Map();
    const scAt = hd => {
      let sc = scCache.get(hd);
      if (sc) { setSun(hd); return sc; }
      sc = sceneAt(hd);
      if (scCache.size >= 6) scCache.delete(scCache.keys().next().value);
      scCache.set(hd, sc);
      return sc;
    };
    const arcExtend = (tk, h0, seedW, fen, rname, slv) => {
      tk.slv = slv;
      const b0 = slv(h0, seedW);
      if (b0 === null) return;
      const turn = tk.kind === 'turn';
      const F = y => {
        const sc = scAt(y[0] / DEG), w = [y[1], y[2], y[3]];
        const f = turn ? (sc.Reach.find(x => x.name === rname) || null) : fen;
        if (f === null) return null;
        const a = sc.P.fval(f, w), b = turn ? jturnW(sc, w) : sc.jhatW(w, tk.eps);
        return Number.isFinite(a) && Number.isFinite(b) ? [a, b, dot(w, w) - 1] : null;
      };
      const opt = { eps: 1e-7, stall: true, tol: [1e-11, 1e-11, 1e-13], maxCorr: 0.05, s0: 0.005, sMax: 0.02, sMin: 1e-7, maxN: 3000, cancelled,
        valid: y => { if (Math.abs(y[0]) > 89 * DEG) return false; const sc = scAt(y[0] / DEG), w = [y[1], y[2], y[3]]; return sc.F.gram(w) >= -1e-9 && sc.minMargin(w) >= -1e-8; } };
      const r = arcBoth(F, [h0 * DEG, ...b0], opt);
      if (!r) return;
      tk.ys = r.ys; tk.arcF = F; tk.arcOpt = opt;
      tk.hs = r.ys.map(y => y[0] / DEG); tk.sols = r.ys.map(y => [y[1], y[2], y[3]]);
    };

    // (1) end / turn 身份发现 + 延拓
    for (const s0 of DIV_SEEDS)
      for (const hs of [s0, -s0]) {
        if (cancelled && cancelled()) return tracks;
        setSun(hs);
        let recs;
        try { recs = computeEndpoints(path, m, fixedAng, false); } catch (e) { continue; }
        for (const rc of recs) {
          if (!rc.label.includes('焦散')) continue;
          if (rc.kind === 'end') {
            const probe = { kind: 'end', eps: rc.eps };
            if (onTrack(probe, hs, rc.W)) continue;
            const sc0 = sceneAt(hs);
            let fname = null;
            for (const f of sc0.Geo)
              if (rc.label.includes('焦散∩' + f.name)) { fname = f.name; break; }
            if (fname === null) continue;
            const fenC = sc0.Geo.find(f => f.name === fname);
            const tk = newTrack('end', rc.eps, `焦散∩${fname} ε=${fmtEps(rc.eps)}`);
            extend(tk, hs, rc.W, fenC, null);
            if (tk.hs.length >= 3) tracks.push(tk);
          } else if (rc.kind === 'turn') {
            const probe = { kind: 'turn', eps: 1 };
            if (onTrack(probe, hs, rc.W)) continue;
            const rname = rc.label.includes('可达带上界') ? '可达带上界'
                        : rc.label.includes('可达带下界') ? '可达带下界' : null;
            if (rname === null) continue;
            const tk = newTrack('turn', 1, `焦散∩${rname} (缝)`);
            extend(tk, hs, rc.W, null, rname);
            if (tk.hs.length >= 3) tracks.push(tk);
          }
        }
      }
    // (2) cusp 身份
    if (withCusps)
      for (const s0 of DIV_CUSP_SEEDS)
        for (const hs of [s0, -s0]) {
          if (cancelled && cancelled()) return tracks;
          setSun(hs);
          let recs;
          try { recs = computeEndpoints(path, m, fixedAng, true); } catch (e) { continue; }
          for (const rc of recs) {
            if (rc.kind !== 'cusp') continue;
            const uv = [rc.sc * Math.PI / 180, rc.sw * Math.PI / 180];
            const probe = { kind: 'cusp', eps: 0 };
            if (onTrack(probe, hs, uv)) continue;
            const tk = newTrack('cusp', 0, '焦散尖点 (pleat 脊线)');
            extend(tk, hs, uv, null, null);
            if (tk.hs.length >= 3) tracks.push(tk);
          }
        }
    // (2b) 修正版: 闭环焦散上的尖点. HaloEP 只从焦散的 end/turn 端点出发沿线找尖点, 整圈都在有效区里的闭环没有端点,
    //   上面的尖点它看不到 —— 层线在那里折回 (天空曲线的尖点经锥映射仍是尖点), 却没有脊线. 每 2° 一层追修正版焦散线,
    //   沿线看天空方向的弦是否反向 (尖点处天空曲线折回), 候选处用 solveCuspUV 精确解 (Ĵ = 0 且切向在映射的核里, V = 0),
    //   再同样沿 h 延拓. 判重按周期折回比较 (u, v); 退化纤维 (V 沿整条线 ≈ 0) 用两侧 |V| 的对比度剔除.
    if (withCusps && fine) {
      const TWO = 2 * Math.PI, wd = a => a - TWO * Math.round(a / TWO);
      const known = (hd, uv) => tracks.some(tk => tk.kind === 'cusp' && tk.hs.some((h, i) =>
        Math.abs(h - hd) < 0.8 && Math.abs(wd(tk.sols[i][0] - uv[0])) < 3e-2 && Math.abs(wd(tk.sols[i][1] - uv[1])) < 3e-2));
      const skyN = x => { const c = Chain3(x[0], x[1], m, fixedAng, pctx), n = Math.hypot(c.Dx, c.Dy, c.Dz); return [c.Dx / n, c.Dy / n, c.Dz / n]; };
      for (let hs = -88; hs <= 88; hs += 2) {
        if (cancelled && cancelled()) return tracks;
        const sc = sceneAt(hs);
        for (const run of traceCausticFine(pctx, m, fixedAng, 61)) {
          const n = run.length;
          if (n < 7) continue;
          const closed = Math.hypot(wd(run[0][0] - run[n - 1][0]), wd(run[0][1] - run[n - 1][1])) < 0.05;
          const Ds = run.map(skyN), at = i => closed ? (i + n) % n : i;
          for (let i = closed ? 0 : 1; i < (closed ? n : n - 1); i++) {
            const a = Ds[at(i - 1)], b = Ds[at(i)], c = Ds[at(i + 1)], c0 = sub(b, a), c1 = sub(c, b);
            if (!(nrm(c0) > 1e-9 && nrm(c1) > 1e-9 && dot(c0, c1) < 0)) continue;
            const uv = solveCuspUV(sc, run[at(i)]);
            if (uv === null || known(hs, uv)) continue;
            if (!closed && (i < 3 || i + 3 >= n)) continue;
            const vc = nrm(cuspState(sc, uv[0], uv[1]).V);
            const vn = Math.min(nrm(cuspState(sc, ...run[at(i - 3)]).V), nrm(cuspState(sc, ...run[at(i + 3)]).V));
            if (!(vn > Math.max(1e2 * vc, 1e-8))) continue;
            const tk = newTrack('cusp', 0, '焦散尖点 (pleat 脊线)');
            extend(tk, hs, uv, null, null);
            if (tk.hs.length >= 3) tracks.push(tk);
          }
        }
      }
    }
    // (3) 锥映射
    if (fine) projectTracksFine(tracks, m, fixedAng, path, L, E, maxDist);
    else projectTracks(tracks, m, fixedAng, path, L, E, maxDist);
  } finally { setSun(oldTH); }
  return tracks;
}

// ------------------------------------------------------------------ 可选曲面族
export const SURF_IDS = ['caustic', 'entry', 'exit', 'deg1', 'blue1', 'deg2', 'blue2'];

function fenceById(pc, id) {
  if (id === 'entry') return pc.Fences.find(f => f.name === 'entry') || null;
  if (id === 'exit') return pc.Fences.find(f => f.name === 'exit') || null;
  if (id.startsWith('deg') || id.startsWith('blue')) {
    const ks = id.substring(id.startsWith('deg') ? 3 : 4);
    if (!/^\d+$/.test(ks)) return null;
    const k = parseInt(ks, 10);
    const suff = id.startsWith('deg') ? '=0' : '=-√K';
    return pc.Fences.find(f => f.name.startsWith(`q${k}(`) && f.name.endsWith(suff)) || null;
  }
  return null;
}

/** 曲面自己的曲线名 (分类边界线用): caustic → "焦散", 否则篱笆名. */
export function surfaceFenceName(path, surfId) {
  if (surfId == null || surfId === 'caustic') return '焦散';
  const f = fenceById(new PathC(path, LP.partial), surfId);
  return f !== null ? f.name : '焦散';
}

/** 选择的曲面在该光路上不存在时回退焦散. */
export function resolveSurface(path, id) {
  if (id == null || id === 'caustic') return 'caustic';
  return fenceById(new PathC(path, LP.partial), id) !== null ? id : 'caustic';
}

// 每 h 层的篱笆曲面等参线 (闭式取样 → 有效段 → (sc,sw) 运行); 用当前 LP 太阳. 参考场景按 key 缓存.
let _frKey = null, _frScene = null;
// fine = true (显示用): 段端不停在最后一个有效样本上, 而在样本间沿篱笆二分到精确终点 —
//   终点要么是与别的篱笆相交 (落在 h-延拓的边界线上), 要么是可达带边缘 (ε=±1 两段在此汇合, 两段端点重合);
//   该处取向随参数呈 √ 变化, 停在样本上会差出零点几米到两米. 另把周期坐标展开成连续的, 加密时中点不会跨接缝.
export function fenceRuns(path, m, fixedAng, surfId, fine = false) {
  const runs = [];
  const key = path.join('-') + '|' + m + '|' + fixedAng.toFixed(6) + '|' + LP.partial + '|' + LP.upperApexDeg + '|' + LP.lowerApexDeg;
  if (_frKey !== key || _frScene === null) {
    const pc0 = new PathC(path, LP.partial);
    const keep = LP.TH * 180 / Math.PI;
    setSun(10);
    _frScene = new Scene(pc0, FamC.make(m, fixedAng), new PathCtx(path));
    setSun(keep);
    _frKey = key;
  }
  const fam = FamC.make(m, fixedAng);
  if (fam === null) return runs;
  const sc = lightScene(_frScene, fam);
  const fen = fenceById(sc.P, surfId);
  if (fen === null) return runs;
  if (fine) return fenceRunsFine(fam, sc, fen, runs);
  for (const comp of fenceSamples(sc, fen, 700)) {
    for (let e = 0; e < 2; e++) {
      const eps = e === 0 ? 1 : -1;
      let cur = null;
      for (const w of comp) {
        let ok = fam.gram(w) >= -1e-9 && sc.minMargin(w) >= -1e-7;
        let scv = 0, swv = 0;
        if (ok) {
          const E = fam.lift(w, eps);
          ok = E !== null;
          if (ok) [scv, swv] = fam.uvOf(E[0], E[1], E[2]);
        }
        if (!ok) {
          if (cur !== null && cur.length >= 2) runs.push(cur);
          cur = null;
          continue;
        }
        if (cur === null) cur = [];
        cur.push([scv, swv]);
      }
      if (cur !== null && cur.length >= 2) runs.push(cur);
    }
  }
  return runs;
}

function fenceRunsFine(fam, sc, fen, runs) {
  const uvAt = (w, eps) => {
    if (!(fam.gram(w) >= -1e-9 && sc.minMargin(w) >= -1e-7)) return null;
    const E = fam.lift(w, eps);
    return E === null ? null : fam.uvOf(E[0], E[1], E[2]);
  };
  // a 有效, b 无效: 沿弦 (再归一化到单位球, 偏离篱笆是二阶小量) 二分 40 次 → 最后一个有效点
  const edge = (wa, wb, eps) => {
    let lo = wa, hi = wb, uv = uvAt(wa, eps);
    for (let it = 0; it < 40; it++) {
      const mm = [0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]), 0.5 * (lo[2] + hi[2])], n = nrm(mm);
      const wm = n > 1e-15 ? scl(mm, 1 / n) : mm;
      const u = uvAt(wm, eps);
      if (u) { lo = wm; uv = u; } else hi = wm;
    }
    return uv;
  };
  const TWO = 2 * Math.PI;
  const unwrap = run => {
    for (let i = 1; i < run.length; i++)
      for (let k = 0; k < 2; k++) run[i][k] -= TWO * Math.round((run[i][k] - run[i - 1][k]) / TWO);
    return run;
  };
  for (const comp of fenceSamples(sc, fen, 700)) {
    for (let e = 0; e < 2; e++) {
      const eps = e === 0 ? 1 : -1;
      let cur = null, pw = null, puv = null;
      for (const w of comp) {
        const uv = uvAt(w, eps);
        if (uv && !puv && pw) cur = [edge(w, pw, eps)];
        if (!uv && puv) { cur.push(edge(pw, w, eps)); if (cur.length >= 2) runs.push(unwrap(cur)); cur = null; }
        if (uv) { if (!cur) cur = []; cur.push(uv); }
        pw = w; puv = uv;
      }
      if (cur && cur.length >= 2) runs.push(unwrap(cur));
    }
  }
  return runs;
}

function solveFF(sc, fa, fb, eps, seed) {
  const w = polishPair(x => sc.P.fval(fa, x), x => sc.P.fgrad(fa, x),
                       x => sc.P.fval(fb, x), x => sc.P.fgrad(fb, x), seed);
  if (w === null) return null;
  if (Math.abs(sc.P.fval(fa, w)) > 1e-9 || Math.abs(sc.P.fval(fb, w)) > 1e-9) return null;
  if (sc.F.gram(w) < -1e-9) return null;
  if (sc.minMargin(w) < -1e-8) return null;
  return w;
}

/** 任意选定曲面的边界线: surfId="caustic" 走焦散版; 否则该篱笆的末端点 (篱笆∩篱笆 + 篱笆∩可达带) 的 h-扫描. */
export function computeDivergentEdges(path, m, fixedAng, surfId, lampElevDeg, lampDist, eyeHeight, maxDist,
  dh = 0.5, withCusps = true, cancelled = null, fine = false) {
  if (surfId == null || surfId === 'caustic')
    return computeDivergentCausticEdges(path, m, fixedAng, lampElevDeg, lampDist, eyeHeight, maxDist, dh, withCusps, cancelled, fine);
  const oldTH = LP.TH * 180 / Math.PI;
  const tracks = [];
  try {
    const se = lampElevDeg * Math.PI / 180;
    const L = [lampDist * Math.cos(se), 0.0, lampDist * Math.sin(se)];
    const E = [0.0, 0.0, eyeHeight];
    const fam0 = FamC.make(m, fixedAng);
    if (fam0 === null) return tracks;
    setSun(10);
    const pctx = new PathCtx(path);
    const pc = new PathC(path, LP.partial);
    const refScene = new Scene(pc, FamC.make(m, fixedAng), pctx);
    const fenS = fenceById(pc, surfId);
    if (fenS === null) return tracks;
    const fname = fenS.name;
    const sceneAt = hd => { setSun(hd); return lightScene(refScene, FamC.make(m, fixedAng)); };
    const onTrack = (kind, hd, sol, eps) => {
      for (const tk of tracks) {
        if (tk.kind !== kind) continue;
        if (Math.abs(tk.eps - eps) > 0.5) continue;
        let best = -1, bd = 0.8;
        for (let i = 0; i < tk.hs.length; i++)
          if (Math.abs(tk.hs[i] - hd) < bd) { bd = Math.abs(tk.hs[i] - hd); best = i; }
        if (best < 0) continue;
        const a = tk.sols[best];
        let same = true;
        for (let i = 0; i < 3; i++)
          if (Math.abs(a[i] - sol[i]) > 3e-2) { same = false; break; }
        if (same) return true;
      }
      return false;
    };
    const extend = (tk, h0, sd0, otherName) => {
      const slv = (hd, sd) => {
        const sc = sceneAt(hd);
        const fa = sc.Geo.find(f => f.name === fname) || null;
        if (fa === null) return null;
        const fb = tk.kind === 'turn'
          ? (sc.Reach.find(f => f.name === otherName) || null)
          : (sc.Geo.find(f => f.name === otherName) || null);
        if (fb === null) return null;
        return solveFF(sc, fa, fb, tk.eps, sd);
      };
      extendTrack(tk, h0, sd0, slv, dh, cancelled);
    };
    for (const s0 of DIV_SEEDS)
      for (const hs of [s0, -s0]) {
        if (cancelled && cancelled()) return tracks;
        setSun(hs);
        let recs;
        try { recs = computeEndpoints(path, m, fixedAng, false); } catch (e) { continue; }
        for (const rc of recs) {
          if (rc.kind !== 'end' && rc.kind !== 'cross' && rc.kind !== 'turn') continue;
          for (const part of rc.label.split(' | ')) {
            const ab = part.split('∩');
            if (ab.length !== 2) continue;
            const other = ab[0] === fname ? ab[1] : ab[1] === fname ? ab[0] : null;
            if (other === null || other.startsWith('焦散')) continue;
            const isTurn = other.startsWith('可达带');
            const kind = isTurn ? 'turn' : 'end';
            if (!isTurn) {
              const ofen = pc.Fences.find(f => f.name === other) || null;
              if (ofen === null || ofen.side === 0.0) continue;   // 非约束方不终结弧
            }
            const eps = isTurn ? 1.0 : rc.eps;
            if (onTrack(kind, hs, rc.W, eps)) continue;
            const tk = newTrack(kind, eps, `${fname}∩${other}` + (isTurn ? ' (缝)' : ''));
            extend(tk, hs, rc.W, other);
            if (tk.hs.length >= 3) tracks.push(tk);
          }
        }
      }
    if (fine) projectTracksFine(tracks, m, fixedAng, path, L, E, maxDist);
    else projectTracks(tracks, m, fixedAng, path, L, E, maxDist);
  } finally { setSun(oldTH); }
  return tracks;
}

// ##################################################################  窗体二自由度光锥

/** 等值线 → 有效段折线 (宽容差 1e-4, 有效性翻转处 22 次二分出精确端点). */
function runsFromGrid(F, SCs, SWs, perX, perY, Fv, P, m, fx) {
  const refine = (xa, ya, xb, yb, fa, fb) => {
    let lo = 0, hi = 1, flo = fa;
    for (let it = 0; it < 18; it++) {
      const mid = 0.5 * (lo + hi);
      const fv = Fv(xa + (xb - xa) * mid, ya + (yb - ya) * mid);
      if (Number.isNaN(fv)) break;
      if (flo * fv <= 0) hi = mid; else { lo = mid; flo = fv; }
    }
    return 0.5 * (lo + hi);
  };
  const VTOL = 1e-4;
  const runs = [];
  for (const c of traceContours(F, SCs, SWs, perX, perY, refine, Fv)) {
    let cur = [];
    let prev = null, prevOk = false;
    for (const pt of c.pts) {
      const ok = validAt(P, m, fx, pt[0], pt[1], VTOL);
      if (prev !== null && ok !== prevOk) {
        const v = prevOk ? prev : pt, iv = prevOk ? pt : prev;
        let lo = 0, hi = 1;
        for (let it = 0; it < 22; it++) {
          const mid = 0.5 * (lo + hi);
          if (validAt(P, m, fx, v[0] + (iv[0] - v[0]) * mid, v[1] + (iv[1] - v[1]) * mid, VTOL)) lo = mid; else hi = mid;
        }
        cur.push([v[0] + (iv[0] - v[0]) * lo, v[1] + (iv[1] - v[1]) * lo]);
        if (!ok) { if (cur.length >= 2) runs.push(cur); cur = []; }
      }
      if (ok) cur.push(pt);
      prev = pt; prevOk = ok;
    }
    if (cur.length >= 2) runs.push(cur);
  }
  return runs;
}

/** J=0 焦散有效段 (错相网格抗退化 + 宽容差 + 精确端点 + 断头延拓). 用当前 LP 太阳. */
export function sliceRuns2Dof(P, m, fx, nGrid) {
  const ax = axesMeta(m);
  const X0 = ax.x0 * Math.PI / 180, X1 = ax.x1 * Math.PI / 180;
  const Y0 = ax.y0 * Math.PI / 180, Y1 = ax.y1 * Math.PI / 180;
  const perX = ax.perX, perY = ax.perY;
  const dX = (X1 - X0) / (nGrid - 1), dY = (Y1 - Y0) / (nGrid - 1);
  const offX = perX ? dX * 0.1371 : 0;
  const offY = perY ? dY * 0.1371 : 0;
  const SCs = new Float64Array(nGrid), SWs = new Float64Array(nGrid);
  for (let i = 0; i < nGrid; i++)
    SCs[i] = perX ? X0 + offX + dX * i : X0 + 1e-4 + (X1 - X0 - 2e-4) * i / (nGrid - 1);
  for (let i = 0; i < nGrid; i++)
    SWs[i] = perY ? Y0 + offY + dY * i : Y0 + 1e-4 + (Y1 - Y0 - 2e-4) * i / (nGrid - 1);
  const Fj = new Float64Array(nGrid * nGrid);
  for (let iw = 0; iw < nGrid; iw++)
    for (let isc = 0; isc < nGrid; isc++)
      Fj[iw * nGrid + isc] = jana3(SCs[isc], SWs[iw], m, fx, P);
  return runsFromGrid(Fj, SCs, SWs, perX, perY, (xa, ya) => jana3(xa, ya, m, fx, P), P, m, fx);
}

/** 取向平面折线 → Chain3 → d=−D → ψ 闭式解 (两支各自成线) → 锥点. verts 非 null 时同时累积点云. */
function coneMapRuns(runs, cosH, sinH, m, fx, P, Dx, Dy, Dz, L, maxDist, verts, outCurves) {
  for (const run of runs) {
    let cur1 = [], cur2 = [];
    const fl1 = () => { if (cur1.length >= 6) outCurves.push(Float32Array.from(cur1)); cur1 = []; };
    const fl2 = () => { if (cur2.length >= 6) outCurves.push(Float32Array.from(cur2)); cur2 = []; };
    for (const pt of run) {
      const ch = Chain3(pt[0], pt[1], m, fx, P);
      const Dn = Math.sqrt(ch.Dx * ch.Dx + ch.Dy * ch.Dy + ch.Dz * ch.Dz);
      if (Number.isNaN(Dn) || Dn < 1e-12) { fl1(); fl2(); continue; }
      const dxp = -ch.Dx / Dn, dyp = -ch.Dy / Dn, dzp = -ch.Dz / Dn;   // 传播方向 d = −D
      const Pc = sinH * (Dx * dxp + Dy * dyp) - cosH * Dx * dzp;
      const Qc = sinH * (Dx * dyp - Dy * dxp) + cosH * Dy * dzp;
      const Rc = cosH * Dz * dyp;
      const mag = Math.sqrt(Pc * Pc + Qc * Qc);
      if (mag < 1e-12) { fl1(); fl2(); continue; }
      let rhs = Rc / mag;
      if (rhs > 1.0 + 1e-6 || rhs < -1.0 - 1e-6) { fl1(); fl2(); continue; }
      if (rhs > 1.0) rhs = 1.0; else if (rhs < -1.0) rhs = -1.0;
      const phi = Math.atan2(Qc, Pc), asr = Math.asin(rhs);
      let Pv = tryConePoint(asr - phi, cosH, sinH, dxp, dyp, dzp, Dx, Dy, Dz, L, maxDist);
      if (Pv) {
        cur1.push(f32(Pv[0]), f32(Pv[1]), f32(Pv[2]));
        if (verts !== null) verts.push(f32(Pv[0]), f32(Pv[1]), f32(Pv[2]));
      } else fl1();
      Pv = tryConePoint(Math.PI - asr - phi, cosH, sinH, dxp, dyp, dzp, Dx, Dy, Dz, L, maxDist);
      if (Pv) {
        cur2.push(f32(Pv[0]), f32(Pv[1]), f32(Pv[2]));
        if (verts !== null) verts.push(f32(Pv[0]), f32(Pv[1]), f32(Pv[2]));
      } else fl2();
    }
    fl1(); fl2();
  }
}


// 修正版: 取向平面曲线 → 两支各自加密 (中点校正回 J=0), 有效↔无效转折处二分; 每条带精确视方向 { p, d }
// checkValid = false: 段已按篱笆自己的有效性 (可达带 + 其余约束, w 空间) 精确裁好; 再用取向空间的宽容差判一次会在
//   段内别处切断 (出射临界篱笆本身就压在出射约束的边上), 段端偏离真正的边界
// 周期坐标差折回到 (−π, π] (闭环可能绕周期轴整圈, 首尾差 2π)
const wd2 = a => a - 2 * Math.PI * Math.round(a / (2 * Math.PI));
function coneMapRunsFine(runs, hDeg, m, fx, P, L, E, maxDist, corr, outCurves, checkValid = true) {
  const hRad = hDeg * Math.PI / 180;
  const Lden = Math.max(1e-3, Math.min(maxDist, dist3(L, E)) / 150);
  for (const run of runs) {
    const cache = new Map();
    const branches = x => {
      let v = cache.get(x);
      if (v === undefined) {
        v = [null, null];
        if (!checkValid || validAt(P, m, fx, x[0], x[1], 1e-4)) {
          const ch = Chain3(x[0], x[1], m, fx, P);
          const Dn = Math.sqrt(ch.Dx * ch.Dx + ch.Dy * ch.Dy + ch.Dz * ch.Dz);
          if (Dn > 1e-12) v = coneBranchesFull([ch.Dx / Dn, ch.Dy / Dn, ch.Dz / Dn], hRad, L, E, maxDist);
        }
        cache.set(x, v);
      }
      return v;
    };
    const mid = (a, b) => corr(a, b);
    const foldAt = x => {
      if (checkValid && !validAt(P, m, fx, x[0], x[1], 1e-4)) return null;
      const ch = Chain3(x[0], x[1], m, fx, P);
      const Dn = Math.sqrt(ch.Dx * ch.Dx + ch.Dy * ch.Dy + ch.Dz * ch.Dz);
      return Dn > 1e-12 ? coneFold([ch.Dx / Dn, ch.Dy / Dn, ch.Dz / Dn], hRad, L, E) : null;
    };
    const smp = insertFolds(run, foldAt, mid);
    const n = run.length, closed = n > 3 && Math.hypot(wd2(run[0][0] - run[n - 1][0]), wd2(run[0][1] - run[n - 1][1])) < 1e-3;
    for (let b = 0; b < 2; b++) {
      const pcs = refineBranch(smp, x => branches(x)[b], mid, Lden);
      // 闭环: 起点 (追线的种子, 位置任意) 两侧的首段与末段其实是同一段, 接起来 —— 否则层线在种子处断成两截,
      //   种子恰在尖点旁时两截都从这里朝同一方向出发, 看起来像停在半空
      if (closed && pcs.length >= 2) {
        const f = pcs[0], l = pcs[pcs.length - 1], k = l.p.length - 3;
        if (Math.hypot(f.p[0] - l.p[k], f.p[1] - l.p[k + 1], f.p[2] - l.p[k + 2]) < 1e-3) {
          const p = new Float32Array(l.p.length + f.p.length - 3), d = new Float32Array(p.length);
          p.set(l.p); p.set(f.p.subarray(3), l.p.length); d.set(l.d); d.set(f.d.subarray(3), l.d.length);
          pcs[0] = { p, d }; pcs.pop();
        }
      }
      outCurves.push(...pcs);
    }
  }
}

// 修正版 (显示) 用的规范化焦散函数 Ĵ = J·√(出射余量). 出射接近全反射临界时出射方向对取向的导数 ∝ 1/√余量,
//   J 随之发散, 零线在出射临界边界前被"甩开", Newton 走不过去 (段端差出零点几米). Ĵ 一直有限到边界, 零集与 J 相同.
function jHat(m, fx, P) {
  return (u, v) => {
    const J = jana3(u, v, m, fx, P);
    if (!Number.isFinite(J)) return J;
    const mE = -SqK - Chain3(u, v, m, fx, P).qE;
    return mE > 0 ? J * Math.sqrt(mE) : J;
  };
}

// 焦散段端点抛光: 端点 = 焦散 (Ĵ=0) 与当前最紧的有效约束 (入射 g0 / 出射临界 / 中间反射) 的交点.
//   二分只能把点放在一条弦上, 弦偏离零线, √ 区又把偏差放大 → 各层端点忽前忽后 (外缘呈锯齿).
//   这里用两个方程 (Ĵ, 约束余量) 的 Newton (差分雅可比, 步长受限) 直接解交点; 不收敛就退回原点.
// 约束余量: [入射 g0, 出射 (−√K − qE), 各中间反射 (−qs)], 全 ≥ 0 才有效
function marginsAt(u, v, P, m, fx) {
  const c = Chain3(u, v, m, fx, P);
  const out = [c.g0, -SqK - c.qE];
  for (const q of c.qs || []) out.push(-q);
  return out;
}
// 解 Ĵ = 0 且 第 k 个余量 = 0 (焦散与有效区边的交点), 从 x0 出发; 失败返回 null. 用当前 LP 太阳.
function edgeNewton(x0, k, P, m, fx, maxMove = 0.05, target = 1e-6) {
  const Jh = jHat(m, fx, P);
  // 目标余量 1e-6 而不是 0: 焦散在出射临界边旁几乎是贴着边走的, 正好在边上的交点病态; 余量 1e-6 处交点良态,
  //   离真正的边只差约 1/1000 的距离 (20 m 处 2 cm)
  const F = (u, v) => [Jh(u, v), marginsAt(u, v, P, m, fx)[k] - target];
  let u = x0[0], v = x0[1];
  // 起点可能已在有效区外 (h 一变, 边界就挪过来了; 那边 Ĵ 无定义): 先沿该余量的梯度把它推回到余量 2e-6 处 (余量两边都有定义)
  for (let it = 0; it < 20; it++) {
    const mk = marginsAt(u, v, P, m, fx)[k];
    if (!Number.isFinite(mk)) return null;
    if (mk >= target && Number.isFinite(Jh(u, v))) break;
    const e = 1e-7, gu = (marginsAt(u + e, v, P, m, fx)[k] - mk) / e, gv = (marginsAt(u, v + e, P, m, fx)[k] - mk) / e, g2 = gu * gu + gv * gv;
    if (!(g2 > 0) || !Number.isFinite(g2)) return null;
    const t = (2 * target - mk) / g2;
    u += t * gu; v += t * gv;
  }
  // 差分步长: 出射临界边旁 Ĵ 随余量呈 √ 变化, 步长得远小于目标余量, 否则雅可比估不准 (Newton 退化成线性收敛甚至不收敛)
  const eps = k === 1 ? Math.max(1e-11, Math.min(1e-7, 2e-3 * target)) : 1e-7;
  for (let it = 0; it < 40; it++) {
    const f = F(u, v);
    if (!f.every(Number.isFinite)) return null;
    if (Math.abs(f[0]) < 1e-10 && Math.abs(f[1]) < 1e-3 * target) break;
    // 差分: 贴着边界, 向外一侧可能已无定义, 就改用向内一侧
    const diff = (du, dv) => { for (const e of [eps, -eps]) { const g = F(u + du * e, v + dv * e); if (g.every(Number.isFinite)) return [(g[0] - f[0]) / e, (g[1] - f[1]) / e]; } return null; };
    const gu = diff(1, 0), gv = diff(0, 1);
    if (!gu || !gv) return null;
    const a = gu[0], b = gv[0], c = gu[1], d = gv[1];
    const det = a * d - b * c;
    if (!(Math.abs(det) > 1e-300)) return null;
    let du = -(d * f[0] - b * f[1]) / det, dv = -(-c * f[0] + a * f[1]) / det;
    const st = Math.hypot(du, dv);
    if (st > 0.02) { du *= 0.02 / st; dv *= 0.02 / st; }
    // 回溯: 一步跨出有效区边 (那边 Ĵ 无定义) 就减半, 直到落回来
    let ok = false;
    for (let bt = 0; bt < 40 && !ok; bt++) {
      const g = F(u + du, v + dv);
      if (g.every(Number.isFinite) && g[1] > -0.9 * target) ok = true; else { du *= 0.5; dv *= 0.5; }
    }
    if (!ok) return null;
    u += du; v += dv;
  }
  const f = F(u, v);
  if (!(Math.abs(f[0]) < 1e-6 && Math.abs(f[1]) < 1e-2 * target) || Math.hypot(u - x0[0], v - x0[1]) > maxMove) return null;
  return validAt(P, m, fx, u, v, 1e-7) ? [u, v] : null;
}

// 焦散与有效区边的精确交点.
//   出射临界边 (k=1, 全反射临界) 旁, 出射方向随余量呈 √ 变化: 正好在边上 (余量 0) 解是病态的, 停在余量 1e-6 又差 0.15°
//   (天上看得出, 三维里就是贴着精确边界的第二条线). 做法: 在余量 4e-6 / 1e-6 / 2.5e-7 各解一次,
//   以 σ = √余量 为自变量二次外推到 σ = 0 — 与 HaloEP 末端点 (焦散∩出射) 一致到 1e-12.
//   其它约束旁是光滑的, 直接解到余量 1e-10.
function edgeExact(x0, k, P, m, fx, maxMove = 0.05) {
  if (k !== 1) return edgeNewton(x0, k, P, m, fx, maxMove, 1e-10);
  const tg = [4e-6, 1e-6, 2.5e-7], xs = [];
  let seed = x0;
  for (const t of tg) {
    const y = edgeNewton(seed, k, P, m, fx, maxMove, t);
    if (edgeExact.stat) edgeExact.stat[t] = (edgeExact.stat[t] || 0) + (y ? 0 : 1);
    if (!y) return null;
    xs.push(y); seed = y;
  }
  const [s1, s2, s3] = tg.map(Math.sqrt);
  const w1 = s2 * s3 / ((s1 - s2) * (s1 - s3)), w2 = s1 * s3 / ((s2 - s1) * (s2 - s3)), w3 = s1 * s2 / ((s3 - s1) * (s3 - s2));
  return [w1 * xs[0][0] + w2 * xs[1][0] + w3 * xs[2][0], w1 * xs[0][1] + w2 * xs[1][1] + w3 * xs[2][1]];
}

// 焦散段端点抛光: 端点 = 焦散 (Ĵ=0) 与当前最紧的有效约束的交点. 二分只能把点放在一条弦上,
//   弦偏离零线, √ 区又把偏差放大 → 各层端点忽前忽后. 直接解交点; 不收敛就退回原点.
function polishEnd(x, P, m, fx) {
  const m0 = marginsAt(x[0], x[1], P, m, fx);
  if (!m0.every(Number.isFinite)) return x;
  let k = 0;
  for (let i = 1; i < m0.length; i++) if (m0[i] < m0[k]) k = i;
  if (m0[k] > 0.05) return x;                          // 离任何约束都远: 不是有效区边上的端点
  return edgeExact(x, k, P, m, fx) || x;
}

// 修正版 (显示) 的焦散追线: 数值延拓 (预测–校正), 不再拼网格追线的结果.
//   Ĵ = 0 在取向平面里是光滑曲线 (鞍点处是两条光滑曲线交叉). 网格追线 (marching squares) 在这种曲线上会出各种毛病:
//   原路折回、一个圈绕两三遍、一格里两个交点排反成"之"字、进不了有 NaN 角的格子 —— 层曲线因此停在半空, 参差不齐.
//   这里网格只用来找种子 (格边上 J 变号处); 每条线从种子出发: 沿切线 (⊥∇Ĵ) 走一步 → 沿 ∇Ĵ 的 Newton 拉回 Ĵ = 0,
//   步长按转角自适应 (转角 > 12° 或校正量过大就减半, 走顺了再放大). 于是:
//     · 到 X 形交叉处沿切线直接穿过 (两条线各走各的);
//     · 出了有效区: 在最后一个有效点与出界点之间二分 (中点拉回零线) 到边上, 再解精确交点 (polishEnd);
//     · 碰到 Ĵ 的 0/0 奇点 (几条焦散汇于一点, 常映到灯): Newton 不收敛, 步长缩到 2e-6 仍不行就停在那里;
//     · 走回起点 = 闭合.
//   已追过的线上的种子跳过. 坐标沿线连续展开 (周期轴不折回).
function traceCausticFine(P, m, fx, nGrid) {
  const ax = axesMeta(m);
  const X0 = ax.x0 * Math.PI / 180, X1 = ax.x1 * Math.PI / 180, Y0 = ax.y0 * Math.PI / 180, Y1 = ax.y1 * Math.PI / 180;
  const perX = ax.perX, perY = ax.perY;
  const dX = (X1 - X0) / (nGrid - 1), dY = (Y1 - Y0) / (nGrid - 1);
  const SCs = new Float64Array(nGrid), SWs = new Float64Array(nGrid);
  for (let i = 0; i < nGrid; i++) {
    SCs[i] = perX ? X0 + dX * (0.1371 + i) : X0 + 1e-4 + (X1 - X0 - 2e-4) * i / (nGrid - 1);
    SWs[i] = perY ? Y0 + dY * (0.1371 + i) : Y0 + 1e-4 + (Y1 - Y0 - 2e-4) * i / (nGrid - 1);
  }
  const Fg = new Float64Array(nGrid * nGrid);
  for (let iw = 0; iw < nGrid; iw++)
    for (let isc = 0; isc < nGrid; isc++) Fg[iw * nGrid + isc] = jana3(SCs[isc], SWs[iw], m, fx, P);

  const J = jHat(m, fx, P);
  const TWO = 2 * Math.PI, wd = a => a - TWO * Math.round(a / TWO);
  // 有效: 严格; 出射余量 < 1e-6 当作到边 (有的光路 J 在出射临界边上不发散, Ĵ 沿整条边为 0, 不能让线顺着边滑)
  const ok = x => validAt(P, m, fx, x[0], x[1], 0) && -SqK - Chain3(x[0], x[1], m, fx, P).qE > 1e-6;
  const grad = (u, v) => {
    const e = 1e-7;
    let gx = (J(u + e, v) - J(u - e, v)) / (2 * e), gy = (J(u, v + e) - J(u, v - e)) / (2 * e);
    if (!Number.isFinite(gx)) { const j0 = J(u, v); gx = (J(u + e, v) - j0) / e; if (!Number.isFinite(gx)) gx = (j0 - J(u - e, v)) / e; }
    if (!Number.isFinite(gy)) { const j0 = J(u, v); gy = (J(u, v + e) - j0) / e; if (!Number.isFinite(gy)) gy = (j0 - J(u, v - e)) / e; }
    return [gx, gy];
  };
  // 校正: 从 x 沿 ∇Ĵ 的 Newton 回到 Ĵ = 0; 总位移超过 lim 或不收敛 → null
  const correct = (x, lim) => {
    let u = x[0], v = x[1];
    for (let it = 0; it < 12; it++) {
      const j0 = J(u, v);
      if (!Number.isFinite(j0)) return null;
      const [gx, gy] = grad(u, v), g2 = gx * gx + gy * gy;
      if (!(g2 > 0) || !Number.isFinite(g2)) return null;
      const du = -j0 * gx / g2, dv = -j0 * gy / g2;
      u += du; v += dv;
      if (Math.hypot(u - x[0], v - x[1]) > lim) return null;
      if (Math.hypot(du, dv) < 1e-13) break;
    }
    const j1 = J(u, v);
    return Number.isFinite(j1) && Math.abs(j1) < 1e-9 ? [u, v] : null;
  };
  const tangent = (x, ref) => {
    const [gx, gy] = grad(x[0], x[1]), l = Math.hypot(gx, gy);
    if (!(l > 0)) return null;
    let t = [-gy / l, gx / l];
    if (t[0] * ref[0] + t[1] * ref[1] < 0) t = [-t[0], -t[1]];
    return t;
  };
  // 已追的线: 线段按 (折回到一个周期内的) 格子哈希
  const cell = 0.02, nc = Math.ceil(TWO / cell), wrap = a => a - TWO * Math.floor(a / TWO);
  const hash = new Map();
  const hk = (cu, cv) => ((cu % nc) + nc) % nc + ',' + ((cv % nc) + nc) % nc;
  const addRun = run => { for (let i = 0; i + 1 < run.length; i++) {
    const a = run[i], cu = Math.floor(wrap(a[0]) / cell), cv = Math.floor(wrap(a[1]) / cell), k = hk(cu, cv);
    let l = hash.get(k); if (!l) hash.set(k, l = []); l.push([a, run[i + 1]]); } };
  const nearTraced = (x, tol) => {
    const cu = Math.floor(wrap(x[0]) / cell), cv = Math.floor(wrap(x[1]) / cell);
    for (let du = -2; du <= 2; du++) for (let dv = -2; dv <= 2; dv++) for (const [a, b] of hash.get(hk(cu + du, cv + dv)) || []) {
      const px = wd(x[0] - a[0]), py = wd(x[1] - a[1]), sx = b[0] - a[0], sy = b[1] - a[1], l2 = sx * sx + sy * sy;
      const t = l2 > 0 ? Math.max(0, Math.min(1, (px * sx + py * sy) / l2)) : 0;
      if (Math.hypot(px - t * sx, py - t * sy) < tol) return true;
    }
    return false;
  };
  const HMAX = 0.02, HMIN = 2e-6, COS_T = Math.cos(12 * Math.PI / 180), COS_E = Math.cos(4 * Math.PI / 180);
  // 从 x0 沿 d0 走到头: 返回 { pts (不含 x0), closed }
  // 跨过奇点: 梯度为零的点 (两条零线交叉, 或 Ĵ 的 0/0 点) 上 Newton 失效, 步长缩到底也过不去. 沿切线以上一步的
  //   2 / 4 / 8 倍试跳, 校正后方向一致 (< 20°) 就算跨过去了, 沿原来那条线接着走 (文献里的 "try and resume").
  //   跳到有效区外 = 奇点在边界上, 就停在这里.
  const COS_J = Math.cos(20 * Math.PI / 180);
  const hop = (x, dir, hLast) => {
    for (const f of [2, 4, 8]) {
      const hj = f * Math.max(hLast, 1e-3), y = correct([x[0] + hj * dir[0], x[1] + hj * dir[1]], 0.4 * hj);
      if (!y) continue;
      const sx = y[0] - x[0], sy = y[1] - x[1], sl = Math.hypot(sx, sy);
      if (!(sl > 0.5 * hj) || sx * dir[0] + sy * dir[1] < COS_J * sl) continue;
      const t = tangent(y, [sx, sy]);
      if (t && t[0] * dir[0] + t[1] * dir[1] > COS_J) return ok(y) ? { y, t } : null;
    }
    return null;
  };
  // 踩回自己: 本条已走过的线段 (同一格子哈希, 记弧长). 近似 X 形交叉 (近鞍点, 零线是两支几乎相碰的双曲线) 处,
  //   第一次经过可能顺着双曲线拐上一个闭环, 绕一圈回来再直穿, 从此绕着闭环转到步数上限 —— 闭合判据只认起点, 认不出.
  //   这里: 新点落在本条较早一段 (弧长差 > 4·HMAX) 的线段上且同向 (< 8°, 重走同一条线时方向几乎不差; 浅角自交不算), 说明这段已经追过, 停下.
  const COS_R = Math.cos(8 * Math.PI / 180);
  const march = (x0, d0) => {
    const pts = [], own = new Map();
    const ownAdd = (a, b, s0) => { const k = hk(Math.floor(wrap(a[0]) / cell), Math.floor(wrap(a[1]) / cell)); let l = own.get(k); if (!l) own.set(k, l = []); l.push([a, b, s0]); };
    const ownHit = (y, sl, sNow) => {
      const cu = Math.floor(wrap(y[0]) / cell), cv = Math.floor(wrap(y[1]) / cell), tol = Math.max(1e-5, 0.05 * sl);
      for (let du = -1; du <= 1; du++) for (let dv = -1; dv <= 1; dv++) for (const [a, b, s0] of own.get(hk(cu + du, cv + dv)) || []) {
        if (s0 > sNow - 4 * HMAX) continue;
        const px = wd(y[0] - a[0]), py = wd(y[1] - a[1]), sx = b[0] - a[0], sy = b[1] - a[1], l2 = sx * sx + sy * sy;
        if (!(l2 > 0) || (sx * dir[0] + sy * dir[1]) < COS_R * Math.sqrt(l2)) continue;
        const t = Math.max(0, Math.min(1, (px * sx + py * sy) / l2));
        if (Math.hypot(px - t * sx, py - t * sy) < tol) return true;
      }
      return false;
    };
    let x = x0, dir = d0, h = 0.005, run = 0, hLast = 0.005, hops = 0;
    for (let n = 0; n < 6000; n++) {
      const y = correct([x[0] + h * dir[0], x[1] + h * dir[1]], 0.6 * h);
      const sx = y ? y[0] - x[0] : 0, sy = y ? y[1] - x[1] : 0, sl = Math.hypot(sx, sy);
      if (!y || !(sl > 0.3 * h) || sx * dir[0] + sy * dir[1] < COS_T * sl) {
        if ((h *= 0.5) < HMIN) {                                    // 奇点
          const j = hops < 20 ? hop(x, dir, hLast) : null;
          if (!j) {
            // 停下的点若其实贴着有效区边 (零线斜着扎向出射临界边, 边外 Ĵ 无定义, 先失败的是校正而不是有效性判据),
            //   照样抛光到精确交点 —— 出射临界旁位置随 √余量 变化, 余量 1e-5 就差出几厘米
            if (pts.length && Math.min(...marginsAt(x[0], x[1], P, m, fx)) < 1e-3) pts[pts.length - 1] = polishEnd(x, P, m, fx);
            return { pts, closed: false };
          }
          hops++; pts.push(j.y); x = j.y; dir = j.t; h = hLast;
        }
        continue;
      }
      hLast = h;
      if (!ok(y)) {                                                 // 出有效区: 二分到边上, 再解精确交点
        let lo = x, hi = y;
        for (let it = 0; it < 40; it++) {
          const mm = [0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1])], span = Math.hypot(hi[0] - lo[0], hi[1] - lo[1]);
          const mc = correct(mm, span) || mm;
          if (ok(mc)) lo = mc; else hi = mc;
        }
        if (lo !== x) pts.push(polishEnd(lo, P, m, fx));
        return { pts, closed: false };
      }
      run += sl;
      // 闭合: 这一步的线段掠过起点, 且走向与出发时相同 —— 细长的"手指" (两条臂相距 1e-4 级、并排走很远的 U 形零线)
      //   回程那条臂也会从起点旁擦过, 但方向相反, 不是闭合
      if (run > 4 * HMAX && sx * d0[0] + sy * d0[1] > 0.5 * sl) {
        const px = wd(x0[0] - x[0]), py = wd(x0[1] - x[1]), t = Math.max(0, Math.min(1, (px * sx + py * sy) / (sl * sl)));
        if (Math.hypot(px - t * sx, py - t * sy) < 1e-4 * Math.max(1, sl / 1e-3)) { pts.push([x[0] + px, x[1] + py]); return { pts, closed: true }; }
      }
      if (run > 8 * HMAX && ownHit(y, sl, run)) return { pts, closed: false };
      ownAdd(x, y, run - sl);
      pts.push(y);
      const t = tangent(y, [sx, sy]);
      if (!t) return { pts, closed: false };
      if (sx * dir[0] + sy * dir[1] > COS_E * sl && t[0] * dir[0] + t[1] * dir[1] > COS_E) h = Math.min(HMAX, 1.5 * h);
      x = y; dir = t;
    }
    return { pts, closed: false };
  };
  const runs = [];
  const trace = s0 => {
    const t0 = tangent(s0, [1, 0]);
    if (!t0) return;
    const f = march(s0, t0);
    let run;
    if (f.closed) run = [s0, ...f.pts];
    else {
      const b = march(s0, [-t0[0], -t0[1]]);
      run = b.pts.reverse().concat([s0], f.pts);
    }
    if (run.length >= 2) { runs.push(run); addRun(run); }
  };
  // 种子: 格边上 J 变号处 (二分到边上的根, 再拉回 Ĵ = 0)
  const seedOn = (ua, va, fa, ub, vb, fb) => {
    if (!(fa * fb < 0)) return;
    let lo = 0, hi = 1, flo = fa;
    for (let it = 0; it < 40; it++) {
      const md = 0.5 * (lo + hi), fv = jana3(ua + (ub - ua) * md, va + (vb - va) * md, m, fx, P);
      if (!Number.isFinite(fv)) return;
      if (flo * fv <= 0) hi = md; else { lo = md; flo = fv; }
    }
    const x = [ua + (ub - ua) * lo, va + (vb - va) * lo];
    if (!ok(x) || nearTraced(x, 2e-3)) return;
    const y = correct(x, 0.01);
    if (y && ok(y) && !nearTraced(y, 2e-3)) trace(y);
  };
  for (let iw = 0; iw < nGrid; iw++)
    for (let isc = 0; isc < nGrid; isc++) {
      const f0 = Fg[iw * nGrid + isc];
      if (!Number.isFinite(f0)) continue;
      if (isc + 1 < nGrid) seedOn(SCs[isc], SWs[iw], f0, SCs[isc + 1], SWs[iw], Fg[iw * nGrid + isc + 1]);
      if (iw + 1 < nGrid) seedOn(SCs[isc], SWs[iw], f0, SCs[isc], SWs[iw + 1], Fg[(iw + 1) * nGrid + isc]);
    }
  // 比一格还小的闭环 (有效区在 h 的极值附近收缩时, 焦散缩成一个小圈, 格边上不变号): 从 |J| 的局部极小出发,
  //   Newton 能在一格之内走到 Ĵ = 0 就以它为种子
  for (let iw = 1; iw + 1 < nGrid; iw++)
    for (let isc = 1; isc + 1 < nGrid; isc++) {
      const f0 = Fg[iw * nGrid + isc];
      if (!Number.isFinite(f0)) continue;
      let minimal = true, spread = 0;
      for (let a = -1; a <= 1 && minimal; a++) for (let c = -1; c <= 1; c++) {
        if (!a && !c) continue;
        const f = Fg[(iw + a) * nGrid + isc + c];
        if (!Number.isFinite(f) || f * f0 <= 0 || Math.abs(f) <= Math.abs(f0)) { minimal = false; break; }
        spread = Math.max(spread, Math.abs(f - f0));
      }
      if (!minimal || Math.abs(f0) > spread) continue;
      const x = [SCs[isc], SWs[iw]], y = correct(x, 1.5 * Math.max(dX, dY));
      if (y && ok(y) && !nearTraced(y, 2e-3)) trace(y);
    }
  return runs;
}

// 中点校正: 焦散 = 沿 ∇J 的 Newton 回到 J=0 (步长受限, 跑远了就用线性中点); 篱笆曲线 = 线性中点 (取样本身很密)
function midCorrector(P, m, fx, caustic) {
  if (!caustic) return (a, b) => [0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])];
  const J = jHat(m, fx, P);
  return (a, b) => {
    const m0 = [0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])];
    const span = Math.hypot(b[0] - a[0], b[1] - a[1]);
    let u = m0[0], v = m0[1];
    for (let it = 0; it < 5; it++) {
      const J0 = J(u, v);
      if (!Number.isFinite(J0)) return m0;
      const e = 1e-6;
      const gx = (J(u + e, v) - J(u - e, v)) / (2 * e), gy = (J(u, v + e) - J(u, v - e)) / (2 * e);
      const g2 = gx * gx + gy * gy;
      if (!(g2 > 0) || !Number.isFinite(g2)) break;
      let du = -J0 * gx / g2, dv = -J0 * gy / g2;
      const st = Math.hypot(du, dv), lim = 0.25 * span + 1e-9;
      if (st > lim) { du *= lim / st; dv *= lim / st; }
      u += du; v += dv;
      if (st < 1e-11) break;
    }
    return Math.hypot(u - m0[0], v - m0[1]) <= 0.5 * span + 1e-9 ? [u, v] : m0;
  };
}

/** 窗体光路 (紧凑索引) → HaloLP 面号 + 取向模式 + 固定角. */
export function lp2DofSetup(lpFaces, zetaFree, fixedZetaDeg, fixedThetaDeg) {
  const m = zetaFree ? OMode.ZenFree : OMode.ZetaFixed;
  const fx = zetaFree ? (fixedThetaDeg - 90.0) * Math.PI / 180.0 : fixedZetaDeg * Math.PI / 180.0;
  return { m, fx, P: new PathCtx(lpFaces) };
}

/**
 * 2 自由度光锥: h 层 ih ∈ [from, to) (共 nH 层, h = −90..90). 每层 = 选定曲面 (焦散或某条篱笆) 的平行光曲线
 * → 锥映射. 返回 [{ curves: [Float32Array…], verts: Float32Array }] (每层一项); hFixed 给定时只算仰角 hFixed 一层.
 * 调用方先设好 LP.partial / 锥顶角 (窗体: ModePartial = true).
 */
export function lightCone2DofLayers(lpFaces, lampElev, lampDist, eyeHeight, nH, nGrid, maxDist,
  zetaFree, fixedZetaDeg, fixedThetaDeg, surfId, from = 0, to = nH, hFixed = null, fine = false) {
  const { m, fx, P } = lp2DofSetup(lpFaces, zetaFree, fixedZetaDeg, fixedThetaDeg);
  const elevRad = lampElev * Math.PI / 180.0;
  const Lx = lampDist * Math.cos(elevRad), Lz = lampDist * Math.sin(elevRad);
  const L = [Lx, 0, Lz];
  const Dx = Lx, Dy = 0, Dz = Lz - eyeHeight;
  const oldTH = LP.TH * 180 / Math.PI;
  const layers = [];
  try {
    for (let ih = from; ih < to; ih++) {
      const hDeg = hFixed !== null ? hFixed : -90.0 + 180.0 * ih / Math.max(1, nH - 1);   // hFixed: 单层 (锥切片)
      const hRad = hDeg * Math.PI / 180.0;
      const cosH = Math.cos(hRad), sinH = Math.sin(hRad);
      setSun(hDeg);
      const caus = surfId === 'caustic' ? sliceRuns2Dof(P, m, fx, nGrid) : fenceRuns(lpFaces, m, fx, surfId);
      const verts = [], curves = [];
      coneMapRuns(caus, cosH, sinH, m, fx, P, Dx, Dy, Dz, L, maxDist, verts, curves);
      if (fine) {                                   // 修正版曲线 (显示); 采样点仍是原窗体的
        const fc = [], corr = midCorrector(P, m, fx, surfId === 'caustic');
        const fr = surfId === 'caustic' ? traceCausticFine(P, m, fx, nGrid) : fenceRuns(lpFaces, m, fx, surfId, true);
        coneMapRunsFine(fr, hDeg, m, fx, P, L, [0, 0, eyeHeight], maxDist, corr, fc, surfId === 'caustic');
        layers.push({ curves: fc, verts: Float32Array.from(verts) });
      } else layers.push({ curves, verts: Float32Array.from(verts) });
    }
  } finally { setSun(oldTH); }
  return layers;
}

/** 平行光弧 = 灯仰角处的平行光曲线, 投到远天球 (E + 200 m·D̂). */
export function parallelArc2Dof(lpFaces, lampElev, eyeHeight, nGrid, zetaFree, fixedZetaDeg, fixedThetaDeg, surfId, fine = false) {
  const { m, fx, P } = lp2DofSetup(lpFaces, zetaFree, fixedZetaDeg, fixedThetaDeg);
  const oldTH = LP.TH * 180 / Math.PI;
  const pArc = [];
  try {
    setSun(lampElev);
    let runs = surfId === 'caustic' ? sliceRuns2Dof(P, m, fx, nGrid) : fenceRuns(lpFaces, m, fx, surfId, fine);
    if (fine) {
      const corr = midCorrector(P, m, fx, surfId === 'caustic');
      if (surfId === 'caustic') runs = traceCausticFine(P, m, fx, nGrid);
      for (const run of runs) {
        const mapAt = x => {
          if (!validAt(P, m, fx, x[0], x[1], 1e-4)) return null;
          const ch = Chain3(x[0], x[1], m, fx, P);
          const Dn = Math.sqrt(ch.Dx * ch.Dx + ch.Dy * ch.Dy + ch.Dz * ch.Dz);
          if (!(Dn > 1e-12)) return null;
          const D = [ch.Dx / Dn, ch.Dy / Dn, ch.Dz / Dn];
          return { P: [200 * D[0], 200 * D[1], eyeHeight + 200 * D[2]], s: D, u: D };
        };
        pArc.push(...refineBranch(run, mapAt, corr, 3.0));
      }
      return pArc;
    }
    for (const run of runs) {
      let seg = [];
      const fl = () => { if (seg.length >= 6) pArc.push(Float32Array.from(seg)); seg = []; };
      for (const pt of run) {
        const ch = Chain3(pt[0], pt[1], m, fx, P);
        const Dn = Math.sqrt(ch.Dx * ch.Dx + ch.Dy * ch.Dy + ch.Dz * ch.Dz);
        if (Number.isNaN(Dn) || Dn < 1e-12) { fl(); continue; }
        seg.push(f32(200.0 * ch.Dx / Dn), f32(200.0 * ch.Dy / Dn), f32(eyeHeight + 200.0 * ch.Dz / Dn));
      }
      fl();
    }
  } finally { setSun(oldTH); }
  return pArc;
}

/**
 * 二自由度边界线 (精确): 末端点身份 h-延拓 × 锥闭式, 按"对方曲线"分五桶
 * (entry 掠入 / tir 出射临界·蓝边 / seam 接缝 / turn 转折缝 / cusp 尖点脊线), 同窗体.
 */
export function edges2Dof(lpFaces, lampElev, lampDist, eyeHeight, maxDist, zetaFree, fixedZetaDeg, fixedThetaDeg, surfId, cancelled = null, fine = false) {
  const { m, fx } = lp2DofSetup(lpFaces, zetaFree, fixedZetaDeg, fixedThetaDeg);
  const out = { entry: [], tir: [], seam: [], turn: [], cusp: [], tracks: [] };
  const exact = computeDivergentEdges(lpFaces, m, fx, surfId, lampElev, lampDist, eyeHeight, maxDist, 0.5, true, cancelled, fine);
  const ownName = surfaceFenceName(lpFaces, surfId);
  for (const tk of exact) {
    let dst;
    if (tk.kind === 'cusp') dst = out.cusp;
    else if (tk.kind === 'turn') dst = out.turn;
    else {
      const rest = tk.label.split(ownName).join('');
      dst = rest.includes('entry') ? out.entry : rest.includes('=0') ? out.seam : out.tir;
    }
    if (fine) { for (const pc of tk.pts3) if (!dst.some(q => samePiece(q, pc))) dst.push(pc); }
    else dst.push(...tk.pts3);
    out.tracks.push({ kind: tk.kind, label: tk.label, n: tk.hs.length, h0: tk.hs[0], h1: tk.hs[tk.hs.length - 1], nShadow: tk.nShadow });
  }
  return out;
}

// 修正版去重: 同一条脊线可能被不同种子各追一遍 (参数差 2π 的同一解, 原窗体判重时没折回) → 三维折线逐点重合的只留一条
function samePiece(a, b) {
  const p = a.p, q = b.p;
  if (p.length !== q.length) return false;
  for (let i = 0; i < p.length; i++) if (Math.abs(p[i] - q[i]) > 1e-4) return false;
  return true;
}

export { allFams, famEval3, coneBranchesFull as _coneBranchesFull };

/** 调试: 某层 (太阳 h) 焦散的修正版段 (取向坐标), 以及每段端点的约束余量. */
export function _debugCausticRuns(lpFaces, zetaFree, fixedZetaDeg, fixedThetaDeg, hDeg, nGrid = 101) {
  const { m, fx, P } = lp2DofSetup(lpFaces, zetaFree, fixedZetaDeg, fixedThetaDeg);
  const keep = LP.TH * 180 / Math.PI;
  setSun(hDeg);
  try {
    const raw = sliceRuns2Dof(P, m, fx, nGrid);
    const fin = traceCausticFine(P, m, fx, nGrid), ext = fin;
    const info = x => { const c = Chain3(x[0], x[1], m, fx, P); return { x, g0: c.g0, mE: -SqK - c.qE, qs: Array.from(c.qs || []), J: jana3(x[0], x[1], m, fx, P) }; };
    return { raw, ext, fin, info, m, fx, P };
  } finally { setSun(keep); }
}


/**
 * 二自由度焦散面的外缘 (精确): 外缘点 = 焦散 (Ĵ=0) 与有效区边 (最紧的约束余量 = 0) 的交点. 它在 (h, 取向) 里是一条曲线,
 * 像末端点轨迹一样沿 h 延拓 (每步以上一步为种子做 2×2 Newton, 失败处把 h 步长减半找到轨迹的尽头), 再按锥闭式两支映射成
 * 三维折线 (中点 h 处重解, 加密、端点二分、折叠插点, 与边界线同一套). 种子来自每 2° 一层的修正版焦散段端点.
 * HaloEP 末端点引擎能认出的 (焦散∩某边界) 与这里重合; 它认不出的 (出射临界旁 J 发散、链失效区旁收口) 这里也有.
 * 返回 [{ p, d }].
 */
export function causticRims2Dof(lpFaces, lampElev, lampDist, eyeHeight, maxDist, zetaFree, fixedZetaDeg, fixedThetaDeg, nGrid, cancelled = null) {
  const { m, fx, P } = lp2DofSetup(lpFaces, zetaFree, fixedZetaDeg, fixedThetaDeg);
  const se = lampElev * Math.PI / 180;
  const L = [lampDist * Math.cos(se), 0, lampDist * Math.sin(se)], E = [0, 0, eyeHeight];
  const oldTH = LP.TH * 180 / Math.PI;
  const tracks = [];
  const TWO = 2 * Math.PI;
  const wdist = (a, b) => Math.hypot(a[0] - b[0] - TWO * Math.round((a[0] - b[0]) / TWO), a[1] - b[1] - TWO * Math.round((a[1] - b[1]) / TWO));
  const out = [];
  try {
    const corr = midCorrector(P, m, fx, true);
    const onTrack = (h, x, k) => tracks.some(tk => tk.k === k && tk.hs.some((hh, i) => Math.abs(hh - h) < 1.01 && wdist(tk.xs[i], x) < 0.03));
    // 外缘 = (h, u, v) 里的曲线 { Ĵ = 0, 余量_k = 1e-6 } (余量 1e-6 上稳, 画的时候再用 edgeExact 到精确的边).
    //   伪弧长延拓 (arcTrace): h 也是未知数, 两个外缘点在某个 h 汇合消失 (h 方向的折叠) 时沿曲线绕过去, 另一支接着走.
    const DEG = Math.PI / 180, TAU = 1e-6, Jh = jHat(m, fx, P);
    const rimF = k => y => {
      setSun(y[0] / DEG);
      const j = Jh(y[1], y[2]), mk = marginsAt(y[1], y[2], P, m, fx)[k];
      return Number.isFinite(j) && Number.isFinite(mk) ? [j, mk - TAU] : null;
    };
    // 出射临界边旁 Ĵ 随余量呈 √ 变化: 差分步长要远小于目标余量 (同 edgeNewton)
    const rimOpt = k => ({ eps: k === 1 ? 2e-3 * TAU : 1e-7, tol: [1e-10, 1e-3 * TAU], maxCorr: 0.02, s0: 0.003, sMax: 0.02, sMin: 1e-7, maxN: 3000, cancelled, wrap: [1, 2],
      valid: y => { if (Math.abs(y[0]) > 89.9 * DEG) return false; setSun(y[0] / DEG); return validAt(P, m, fx, y[1], y[2], 1e-7); } });
    for (let hs = -88; hs <= 88; hs += 2) {
      if (cancelled && cancelled()) break;
      setSun(hs);
      const runs = traceCausticFine(P, m, fx, Math.min(nGrid, 71));   // 种子只要大致位置, 粗网格即可
      for (const run of runs) for (const x of [run[0], run[run.length - 1]]) {
        setSun(hs);
        const mg = marginsAt(x[0], x[1], P, m, fx);
        if (!mg.every(Number.isFinite)) continue;
        let k = 0;
        for (let i = 1; i < mg.length; i++) if (mg[i] < mg[k]) k = i;
        if (mg[k] > 1e-2) continue;                                     // 不在有效区边上 (段在别处结束)
        const x0 = edgeNewton(x, k, P, m, fx);
        if (!x0 || onTrack(hs, x0, k)) continue;
        const F = rimF(k), opt = rimOpt(k), r = arcBoth(F, [hs * DEG, x0[0], x0[1]], opt);
        if (!r || r.ys.length < 3) continue;
        tracks.push({ k, F, opt, ys: r.ys, hs: r.ys.map(y => y[0] / DEG), xs: r.ys.map(y => [y[1], y[2]]) });
      }
    }
    if (causticRims2Dof.debug) causticRims2Dof.debug.push(...tracks.map(tk => ({ k: tk.k, h0: tk.hs[0], h1: tk.hs[tk.hs.length - 1], n: tk.hs.length, x0: tk.xs[0], x1: tk.xs[tk.xs.length - 1], ys: tk.ys })));
    // 投影: 两支, 中点在中间 h 重解, 与边界线同一套加密
    const Lden = Math.max(1e-3, Math.min(maxDist, dist3(L, E)) / 150), sag = Math.max(2e-4, 1e-4 * dist3(L, E));
    for (const tk of tracks) {
      // 延拓在余量 1e-6 上做 (稳), 画出来的点用精确交点 (edgeExact)
      //   (外推的几个 Newton 偶尔解不出, 就用余量 1e-6 的点顶上, 只差 0.15°, 不留缺口)
      const exact = (h, x6) => { setSun(h); return edgeExact(x6, tk.k, P, m, fx) || x6; };
      const sample = y => { const h = y[0] / DEG, x6 = [y[1], y[2]]; return { h, x: exact(h, x6), x6, y }; };
      const samples = tk.ys.map(sample);
      // 中点: 弦中点处垂直于弦的超平面与曲线的交点 (折叠两侧的样本之间也对)
      //   (两条边重合的退化段上超平面校正可能失败: 退回定 h 的 edgeNewton, 此时两端 h 不同、中间没有折叠)
      const mid = (a, b) => {
        let y = arcMid(tk.F, a.y, b.y, 0.5, tk.opt);
        if (!y && a.h !== b.h) {
          const h = 0.5 * (a.h + b.h);
          setSun(h);
          const x6 = edgeNewton([0.5 * (a.x6[0] + b.x6[0]), 0.5 * (a.x6[1] + b.x6[1])], tk.k, P, m, fx, 0.05);
          if (x6) y = [h * DEG, x6[0], x6[1]];
        }
        return y ? sample(y) : null;
      };
      const cache = new Map();
      const sky = s0 => {
        setSun(s0.h);
        const ch = Chain3(s0.x[0], s0.x[1], m, fx, P);
        const Dn = Math.sqrt(ch.Dx * ch.Dx + ch.Dy * ch.Dy + ch.Dz * ch.Dz);
        return Dn > 1e-12 ? [ch.Dx / Dn, ch.Dy / Dn, ch.Dz / Dn] : null;
      };
      const branches = s0 => {
        let v = cache.get(s0);
        if (v === undefined) { const D = sky(s0); v = D ? coneBranchesFull(D, s0.h * Math.PI / 180, L, E, maxDist) : [null, null]; cache.set(s0, v); }
        return v;
      };
      const foldAt = s0 => { const D = sky(s0); return D ? coneFold(D, s0.h * Math.PI / 180, L, E) : null; };
      const smp = insertFolds(samples, foldAt, mid, 24);
      for (let b = 0; b < 2; b++) out.push(...refineBranch(smp, s0 => branches(s0)[b], mid, Lden, 9, 24, sag));
    }
  } finally { setSun(oldTH); }
  return out;
}

export const _dbg = { edgeNewton, marginsAt, jHat, edgeExact, arcSolve, arcMid, jacN, nullVec };
