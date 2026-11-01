// 二自由度 (柱晶 Column: α+θ 自由 / Lowitz: α+ζ 自由) 发散光晕面 —
//   (1) Divergentlightform.cs 的 ComputeHaloSurfaceLightCone2Dof / SliceRuns2Dof / RunsFromGrid / ConeMapRuns:
//       每个锥仰角 h 把灯当"太阳", 用光点法在取向平面上追 J=0 焦散 (或选定的一条平行光曲线),
//       焦散点的平行光出射方向 d = −D 代入与一自由度相同的共面闭式 ψ 定位晶体 P = L + t·û;
//   (2) HaloEndpointEngine.cs 的发散光部分 (ComputeDivergentCausticEdges / ComputeDivergentEdges / FenceRuns …):
//       曲面边界线 = 平行光末端点身份的 h-延拓 × 锥闭式.
// 与 C# 逐句对应 (对账见 tools/).

import { LP, OMode, setSun, PathCtx, Chain3, jana3, validAt, axesMeta, traceContours, allFams, famEval3 } from './lp.js';
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
    const c = dot(u, dd), den = 1 - c * c;
    if (den < 1e-12) continue;
    const tR = (c * dot(dd, DDv) - dot(u, DDv)) / den, sR = (c * dot(u, DDv) - dot(dd, DDv)) / den;
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

// 一条参数化曲线 (参数 x: 标量或数组, 任意) 的一支 → 加密好的折线们. 通用于 h-层曲线与 h-延拓轨迹.
//   samples: 参数序列; mapAt(x) → 该支状态或 null; mid(xa, xb) → 中间参数 (可含校正) 或 null.
//   有效 ↔ 无效 的转折处二分 (iters 次) 找精确端点. 返回 [{ p, d }].
function refineBranch(samples, mapAt, mid, Lden, maxDepth = 9, iters = 16) {
  const out = [];
  let P = [], D = [];
  const push = st => { P.push(st.P[0], st.P[1], st.P[2]); D.push(st.s[0], st.s[1], st.s[2]); };
  const flush = () => { if (P.length >= 6) out.push({ p: Float32Array.from(P), d: Float32Array.from(D) }); P = []; D = []; };
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
    if (dep >= maxDepth || !tooFar(sa, sb, Lden)) { push(sb); return; }
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
  const Lden = Math.max(0.1, Math.min(maxDist, dist3(L, E)) / 150);
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
    const samples = tk.hs.map((h, i) => ({ h, sol: tk.sols[i] }));
    const mid = (a, b) => {
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
    for (let b = 0; b < 2; b++) tk.pts3.push(...refineBranch(samples, x => branches(x)[b], mid, Lden));
  }
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
      extendTrack(tk, h0, seedW, slv, dh, cancelled);
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
export function fenceRuns(path, m, fixedAng, surfId) {
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
function coneMapRunsFine(runs, hDeg, m, fx, P, L, E, maxDist, corr, outCurves) {
  const hRad = hDeg * Math.PI / 180;
  const Lden = Math.max(0.1, Math.min(maxDist, dist3(L, E)) / 150);
  for (const run of runs) {
    const cache = new Map();
    const branches = x => {
      let v = cache.get(x);
      if (v === undefined) {
        v = [null, null];
        if (validAt(P, m, fx, x[0], x[1], 1e-4)) {
          const ch = Chain3(x[0], x[1], m, fx, P);
          const Dn = Math.sqrt(ch.Dx * ch.Dx + ch.Dy * ch.Dy + ch.Dz * ch.Dz);
          if (Dn > 1e-12) v = coneBranchesFull([ch.Dx / Dn, ch.Dy / Dn, ch.Dz / Dn], hRad, L, E, maxDist);
        }
        cache.set(x, v);
      }
      return v;
    };
    const mid = (a, b) => corr(a, b);
    for (let b = 0; b < 2; b++) outCurves.push(...refineBranch(run, x => branches(x)[b], mid, Lden));
  }
}

// 中点校正: 焦散 = 沿 ∇J 的 Newton 回到 J=0 (步长受限, 跑远了就用线性中点); 篱笆曲线 = 线性中点 (取样本身很密)
function midCorrector(P, m, fx, caustic) {
  if (!caustic) return (a, b) => [0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])];
  const J = (u, v) => jana3(u, v, m, fx, P);
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
        const fc = [];
        coneMapRunsFine(caus, hDeg, m, fx, P, L, [0, 0, eyeHeight], maxDist, midCorrector(P, m, fx, surfId === 'caustic'), fc);
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
    const runs = surfId === 'caustic' ? sliceRuns2Dof(P, m, fx, nGrid) : fenceRuns(lpFaces, m, fx, surfId);
    if (fine) {
      const corr = midCorrector(P, m, fx, surfId === 'caustic');
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

export { allFams, famEval3 };
