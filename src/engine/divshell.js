// 发散度球壳引擎 — DivergencyShellEngine.cs (static class DivShell) 的 JS 移植.
// 与 C# 逐句对应, 运算次序不变 (与原引擎逐点对账见 tools/cmp-divshell.mjs).
// 工程层面的差别: Parallel.For → 普通循环; CancellationToken → opts.cancelled() (抛 CancelledError);
// double[,] → 一维数组 (行主序); out 参数 → 返回值 (失败返回 null).
//
// 问题: 以眼 E 为心、半径 r 的球壳与一条光路的三维晕面相交成什么曲线? ρ = r/D (D = |L−E|), κ = log10 ρ.
// 一切写在晶体系太阳球 w ∈ S² 上 (理论见 发散度_理论.md §9):
//   链 (闭式):  X=w·u0, S=√(X²+K), P=L·w+S·v, q_E=ℓ_E·w+c_E·S, T_E=√(q_E²−K), D_c=P−(q_E+T_E)·u_E   (K=n²−1)
//   一自由度 R=Rz(α)·R0, R0=Ry(ζ)·Rz(θ): ν=R0ᵀẑ, u_z=−w·ν, e_z=−D_c·ν, p=w·D_c
//   壳方程 Φ(w) = (S_l+ρe_z)² − 2ρ·p·u_z·(S_l+ρe_z) + (ρ²−1)·u_z² = 0
//   篱笆: 入射掠射 X=0 / 接缝 q_k=0 / 蓝边 q_k=−√K (仅纯全反射模式) / 出射临界 q_E=−√K
//   末端点 = 曲线 ∩ 篱笆; 弧 = (w,σ) 伪弧长延拓; 无末端的闭环由网格 Φ 变号补.
// 二自由度 (Column: α+θ / Lowitz: α+ζ): 轮廓 = 焦散 (扩展状态 (w,T,τ,σ) 5 方程延拓) + 各篱笆像 (闭式).

export const DEG = Math.PI / 180.0;
const WALL = 0.02;            // 出射篱笆 margin 小于此值 → 允许括根校正 / 末端检测
const MAX_SKY = 0.5 * DEG;    // 相邻弧点天球弦长上限
const MIN_SKY = 0.12 * DEG;   // 低于此值放大步长
const DSIG = 0.05;            // 种子匹配 / 闭合判定的 σ 容差
const STEP0 = 0.4 * DEG;
const MAX_STEPS = 40000;
const N_FENCE = 2400;
const H6_MAX = 0.007, H6_0 = 0.002;

export class CancelledError extends Error {
  constructor() { super('cancelled'); this.name = 'CancelledError'; }
}
const never = () => false;
function ckCancel(cx) { if (cx.cancelled()) throw new CancelledError(); }

// ================================================================== 向量小件
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const nrm = a => Math.sqrt(dot(a, a));
const scl = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const addSc = (a, b, s) => [a[0] + s * b[0], a[1] + s * b[1], a[2] + s * b[2]];
function dist(a, b) {
  const x = a[0] - b[0], y = a[1] - b[1], z = a[2] - b[2];
  return Math.sqrt(x * x + y * y + z * z);
}
function unit(a) {
  const n = nrm(a);
  return [a[0] / n, a[1] / n, a[2] / n];
}
const unitMid = (a, b) => unit([0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2])]);
function ang(a, b) {
  if (a == null || b == null) return 0.0;
  let c = dot(a, b);
  if (c > 1) c = 1; else if (c < -1) c = -1;
  return Math.acos(c);
}
const fin = x => Number.isFinite(x);
const fin3 = v => v != null && fin(v[0]) && fin(v[1]) && fin(v[2]);
function projT(g, w) {
  const d = dot(g, w);
  g[0] -= d * w[0]; g[1] -= d * w[1]; g[2] -= d * w[2];
}

// .NET Math.IEEERemainder
function ieeeRemainder(x, y) {
  if (Number.isNaN(x)) return x;
  if (Number.isNaN(y)) return y;
  const regularMod = x % y;
  if (Number.isNaN(regularMod)) return NaN;
  if (regularMod === 0 && (x < 0 || Object.is(x, -0))) return -0;
  const alternativeResult = regularMod - (Math.abs(y) * Math.sign(x));
  if (Math.abs(alternativeResult) === Math.abs(regularMod)) {
    const divisionResult = x / y;
    const roundedResult = roundEven(divisionResult);
    if (Math.abs(roundedResult) > Math.abs(divisionResult)) return alternativeResult;
    return regularMod;
  }
  if (Math.abs(alternativeResult) < Math.abs(regularMod)) return alternativeResult;
  return regularMod;
}
// C# Math.Round 默认银行家舍入
function roundEven(x) {
  if (Math.abs(x % 1) === 0.5) return 2 * Math.round(x / 2);
  return Math.round(x);
}

// 3×3 高斯消元 (部分主元); A 行主序 9 元
function solve3(A, b, x) {
  const M = new Array(12);
  for (let i = 0; i < 3; i++) { for (let j = 0; j < 3; j++) M[i * 4 + j] = A[i * 3 + j]; M[i * 4 + 3] = b[i]; }
  for (let c = 0; c < 3; c++) {
    let piv = c, best = Math.abs(M[c * 4 + c]);
    for (let r = c + 1; r < 3; r++) if (Math.abs(M[r * 4 + c]) > best) { best = Math.abs(M[r * 4 + c]); piv = r; }
    if (!(best > 1e-300) || !fin(best)) return false;
    if (piv !== c) for (let j = 0; j < 4; j++) { const t = M[c * 4 + j]; M[c * 4 + j] = M[piv * 4 + j]; M[piv * 4 + j] = t; }
    for (let r = c + 1; r < 3; r++) {
      const f = M[r * 4 + c] / M[c * 4 + c];
      for (let j = c; j < 4; j++) M[r * 4 + j] -= f * M[c * 4 + j];
    }
  }
  for (let i = 2; i >= 0; i--) {
    let s = M[i * 4 + 3];
    for (let j = i + 1; j < 3; j++) s -= M[i * 4 + j] * x[j];
    x[i] = s / M[i * 4 + i];
  }
  return fin(x[0]) && fin(x[1]) && fin(x[2]);
}

// n×n 高斯消元; A 行主序 n·n 元
function solveN(A, b, n, x) {
  const m1 = n + 1, M = new Array(n * m1);
  for (let i = 0; i < n; i++) { for (let j = 0; j < n; j++) M[i * m1 + j] = A[i * n + j]; M[i * m1 + n] = b[i]; }
  for (let c = 0; c < n; c++) {
    let piv = c, best = Math.abs(M[c * m1 + c]);
    for (let r = c + 1; r < n; r++) if (Math.abs(M[r * m1 + c]) > best) { best = Math.abs(M[r * m1 + c]); piv = r; }
    if (!(best > 1e-13) || !fin(best)) return false;
    if (piv !== c) for (let j = 0; j <= n; j++) { const t = M[c * m1 + j]; M[c * m1 + j] = M[piv * m1 + j]; M[piv * m1 + j] = t; }
    for (let r = c + 1; r < n; r++) {
      const f = M[r * m1 + c] / M[c * m1 + c];
      if (f === 0) continue;
      for (let j = c; j <= n; j++) M[r * m1 + j] -= f * M[c * m1 + j];
    }
  }
  for (let i = n - 1; i >= 0; i--) {
    let s = M[i * m1 + n];
    for (let j = i + 1; j < n; j++) s -= M[i * m1 + j] * x[j];
    x[i] = s / M[i * m1 + i];
    if (!fin(x[i])) return false;
  }
  return true;
}

// ================================================================== 篱笆 / 光路常量
// 篱笆 (L-型线 ℓ·w + c·S = k): { name, L, C, Kk, Side }
class PathK {
  constructor(nrmv, refN, partial) {
    this.K = refN * refN - 1.0; this.SqK = Math.sqrt(this.K);
    const U0 = this.U0 = unit(nrmv[0]);
    const UE = this.UE = unit(nrmv[nrmv.length - 1]);
    let Lm = new Array(9);
    for (let i = 0; i < 3; i++)
      for (let j = 0; j < 3; j++)
        Lm[i * 3 + j] = (i === j ? 1.0 : 0.0) - U0[i] * U0[j];
    let Vv = U0.slice();
    this.Fences = [];        // 只含约束型 (side≠0)
    this.AllFences = [];     // 全部几何曲线, 含部分反射模式下不作约束的蓝边 (Side=0); 二自由度显示用
    this.Fences.push({ name: 'entry', L: U0.slice(), C: 0, Kk: 0, Side: +1 });
    for (let m = 1; m < nrmv.length - 1; m++) {
      const u = unit(nrmv[m]);
      const l = [
        Lm[0] * u[0] + Lm[3] * u[1] + Lm[6] * u[2],
        Lm[1] * u[0] + Lm[4] * u[1] + Lm[7] * u[2],
        Lm[2] * u[0] + Lm[5] * u[1] + Lm[8] * u[2]];                    // Lᵀu
      const cc = dot(Vv, u);
      const seam = { name: 'seam' + m, L: l.slice(), C: cc, Kk: 0, Side: -1 };
      const blue = { name: 'blue' + m, L: l.slice(), C: cc, Kk: -this.SqK, Side: partial ? 0.0 : +1.0 };
      this.Fences.push(seam);
      if (!partial) this.Fences.push(blue);                              // 部分反射模式下蓝边只是强度阶跃线
      this.AllFences.push(seam);
      this.AllFences.push(blue);
      const Ln = new Array(9);
      for (let i = 0; i < 3; i++)
        for (let j = 0; j < 3; j++)
          Ln[i * 3 + j] = Lm[i * 3 + j] - 2 * u[i] * (u[0] * Lm[j] + u[1] * Lm[3 + j] + u[2] * Lm[6 + j]);
      const qv = dot(u, Vv);
      Vv = [Vv[0] - 2 * qv * u[0], Vv[1] - 2 * qv * u[1], Vv[2] - 2 * qv * u[2]];
      Lm = Ln;
    }
    this.Lm = Lm; this.Vv = Vv;
    this.LE = [
      Lm[0] * UE[0] + Lm[3] * UE[1] + Lm[6] * UE[2],
      Lm[1] * UE[0] + Lm[4] * UE[1] + Lm[7] * UE[2],
      Lm[2] * UE[0] + Lm[5] * UE[1] + Lm[8] * UE[2]];
    this.CE = dot(Vv, UE);
    this.Exit = { name: 'exit', L: this.LE.slice(), C: this.CE, Kk: -this.SqK, Side: -1 };
    this.Fences.push(this.Exit);
    this.AllFences.unshift(this.Fences[0]);                             // entry
    this.AllFences.push(this.Exit);
  }

  // 扩展链: T (=T_E) 当独立变量, D_c(w,T) = P(w) − (q_E(w)+T)·u_E。M0 = ∂D_c/∂w|_T, g = ∇q_E。返回 q_E
  chainT(w, T, Dc, M0, g) {
    const U0 = this.U0, Lm = this.Lm, Vv = this.Vv, UE = this.UE, LE = this.LE, CE = this.CE;
    const X = dot(U0, w), S = Math.sqrt(X * X + this.K), xs = X / S;
    const qE = dot(LE, w) + CE * S;
    const tEv = -qE - T;
    for (let i = 0; i < 3; i++)
      Dc[i] = Lm[i * 3] * w[0] + Lm[i * 3 + 1] * w[1] + Lm[i * 3 + 2] * w[2] + S * Vv[i] + tEv * UE[i];
    for (let j = 0; j < 3; j++) g[j] = LE[j] + CE * xs * U0[j];
    if (M0 != null)
      for (let i = 0; i < 3; i++)
        for (let j = 0; j < 3; j++)
          M0[i * 3 + j] = Lm[i * 3 + j] + Vv[i] * xs * U0[j] - UE[i] * g[j];
    return qE;
  }

  // 链: false = 出射全反射一侧 (T_E 无定义)。M 非空时同时给 ∂D_c/∂w
  chain(w, Dc, M) {
    const U0 = this.U0, Lm = this.Lm, Vv = this.Vv, UE = this.UE, LE = this.LE, CE = this.CE, K = this.K;
    const X = dot(U0, w), S = Math.sqrt(X * X + K);
    const qE = dot(LE, w) + CE * S;
    const t2 = qE * qE - K;
    let TE;
    if (t2 >= 0) TE = Math.sqrt(t2);
    else if (t2 > -1e-11) TE = 0.0;
    else return false;
    const tEv = -qE - TE;
    for (let i = 0; i < 3; i++)
      Dc[i] = Lm[i * 3] * w[0] + Lm[i * 3 + 1] * w[1] + Lm[i * 3 + 2] * w[2] + S * Vv[i] + tEv * UE[i];
    if (M != null) {
      const xs = X / S;
      const fac = 1.0 + qE / Math.max(TE, 1e-15);           // d(q_E+T_E) = fac·g; 篱笆上发散, 钳到有限大
      for (let j = 0; j < 3; j++) {
        const gj = LE[j] + CE * xs * U0[j];
        for (let i = 0; i < 3; i++)
          M[i * 3 + j] = Lm[i * 3 + j] + Vv[i] * xs * U0[j] - UE[i] * gj * fac;
      }
    }
    return true;
  }

  fval(f, w) {
    const X = dot(this.U0, w);
    return dot(f.L, w) + f.C * Math.sqrt(X * X + this.K) - f.Kk;
  }

  fgrad(f, w) {
    const U0 = this.U0;
    const X = dot(U0, w), S = Math.sqrt(X * X + this.K), xs = f.C * X / S;
    return [f.L[0] + xs * U0[0], f.L[1] + xs * U0[1], f.L[2] + xs * U0[2]];
  }
}

// ================================================================== 场景 (一条光路 × 一层壳, 一自由度)
class Ctx {
  constructor() {
    this.P = null;
    this.R0 = new Array(9).fill(0);
    this.Nu = null;
    this.Sl = 0; this.Rho = 0;
    this.Eye = null;
    this.Radius = 0; this.AzL = 0;
    this.cancelled = never;
  }

  applyR0(v) {
    const R = this.R0;
    return [
      R[0] * v[0] + R[1] * v[1] + R[2] * v[2],
      R[3] * v[0] + R[4] * v[1] + R[5] * v[2],
      R[6] * v[0] + R[7] * v[1] + R[8] * v[2]];
  }

  // → [u_z, e_z, p] 或 null
  scal(w) {
    const Dc = [0, 0, 0];
    if (!this.P.chain(w, Dc, null)) return null;
    const Nu = this.Nu;
    return [-dot(w, Nu), -dot(Dc, Nu), dot(w, Dc)];
  }

  scalG(w, guz, gez, gp) {
    const Dc = [0, 0, 0], M = new Array(9);
    if (!this.P.chain(w, Dc, M)) return null;
    const Nu = this.Nu;
    const uz = -dot(w, Nu), ez = -dot(Dc, Nu), p = dot(w, Dc);
    for (let j = 0; j < 3; j++) {
      guz[j] = -Nu[j];
      gez[j] = -(Nu[0] * M[j] + Nu[1] * M[3 + j] + Nu[2] * M[6 + j]);
      gp[j] = Dc[j] + w[0] * M[j] + w[1] * M[3 + j] + w[2] * M[6 + j];
    }
    return [uz, ez, p];
  }

  phi(w) {
    const s = this.scal(w);
    if (s === null) return NaN;
    const uz = s[0], ez = s[1], p = s[2], Rho = this.Rho;
    const A = this.Sl + Rho * ez;
    return A * A - 2 * Rho * p * uz * A + (Rho * Rho - 1) * uz * uz;
  }

  phiG(w, grad) {
    const guz = [0, 0, 0], gez = [0, 0, 0], gp = [0, 0, 0];
    const s = this.scalG(w, guz, gez, gp);
    if (s === null) return NaN;
    const uz = s[0], ez = s[1], p = s[2];
    const A = this.Sl + this.Rho * ez, r = this.Rho;
    for (let j = 0; j < 3; j++)
      grad[j] = 2 * A * r * gez[j]
              - 2 * r * (gp[j] * uz * A + p * guz[j] * A + p * uz * r * gez[j])
              + 2 * (r * r - 1) * uz * guz[j];
    return A * A - 2 * r * p * uz * A + (r * r - 1) * uz * uz;
  }

  aFun(w) {       // 退化篱笆上的末端函数 A = S_l + ρ e_z
    const s = this.scal(w);
    if (s === null) return NaN;
    return this.Sl + this.Rho * s[1];
  }

  dprime(w) {
    const s = this.scal(w);
    if (s === null) return NaN;
    const uz = s[0], ez = s[1], p = s[2], Sl = this.Sl;
    return uz * uz + ez * ez - 2 * uz * ez * p - Sl * Sl * (1 - p * p);
  }

  dprimeG(w, grad) {
    const guz = [0, 0, 0], gez = [0, 0, 0], gp = [0, 0, 0];
    const s = this.scalG(w, guz, gez, gp);
    if (s === null) return NaN;
    const uz = s[0], ez = s[1], p = s[2], Sl = this.Sl;
    for (let j = 0; j < 3; j++)
      grad[j] = 2 * uz * guz[j] + 2 * ez * gez[j]
              - 2 * (guz[j] * ez * p + uz * gez[j] * p + uz * ez * gp[j])
              + 2 * Sl * Sl * p * gp[j];
    return uz * uz + ez * ez - 2 * uz * ez * p - Sl * Sl * (1 - p * p);
  }

  // → [f1, f2] 或 null
  f12(w, sig) {
    const s = this.scal(w);
    if (s === null) return null;
    const uz = s[0], ez = s[1], p = s[2], Rho = this.Rho;
    const f1 = sig * uz + Rho * ez + this.Sl;
    const f2 = sig * sig + 2 * sig * Rho * p + Rho * Rho - 1;
    return (fin(f1) && fin(f2)) ? [f1, f2] : null;
  }

  // 2×4 雅可比 [∂F/∂w (切平面投影) | ∂F/∂σ], J 行主序 8 元
  jac(w, sig, J) {
    const guz = [0, 0, 0], gez = [0, 0, 0], gp = [0, 0, 0];
    const s = this.scalG(w, guz, gez, gp);
    if (s === null) return false;
    const uz = s[0], p = s[2], Rho = this.Rho;
    const g1 = [0, 0, 0], g2 = [0, 0, 0];
    for (let j = 0; j < 3; j++) { g1[j] = sig * guz[j] + Rho * gez[j]; g2[j] = 2 * sig * Rho * gp[j]; }
    projT(g1, w); projT(g2, w);
    for (let j = 0; j < 3; j++) { J[j] = g1[j]; J[4 + j] = g2[j]; }
    J[3] = uz; J[7] = 2 * sig + 2 * Rho * p;
    for (let i = 0; i < 8; i++) if (!fin(J[i])) return false;
    return true;
  }

  sigmaOf(w) {
    const s = this.scal(w);
    if (s === null) return NaN;
    const uz = s[0], ez = s[1], p = s[2], Rho = this.Rho;
    const lin = -(this.Sl + Rho * ez);
    if (Math.abs(uz) > 1e-7) return lin / uz;
    const d = Rho * Rho * p * p - Rho * Rho + 1;
    if (d < 0) return NaN;
    const s1 = -Rho * p + Math.sqrt(d), s2 = -Rho * p - Math.sqrt(d);
    return Math.abs(s1 * uz - lin) <= Math.abs(s2 * uz - lin) ? s1 : s2;
  }

  margin(f, w) { return f.Side * this.P.fval(f, w); }

  inOmega(w, tol) {
    for (const f of this.P.Fences) if (!(this.margin(f, w) >= -tol)) return false;
    return true;
  }

  minMarginFence(w) {
    let best = null, bm = Infinity;
    for (const f of this.P.Fences) { const m = this.margin(f, w); if (m < bm) { bm = m; best = f; } }
    return best;
  }

  exitMargin(w) { return this.P.Exit == null ? Infinity : this.margin(this.P.Exit, w); }

  phiIn(w) {
    if (!this.inOmega(w, 1e-12)) return NaN;
    return this.phi(w);
  }

  fgradT(f, w) {
    const g = this.P.fgrad(f, w);
    projT(g, w);
    return g;
  }

  // (w,σ) → 世界系: 天球单位方向 ξ (自眼) 与晶体方位 α (度)。失败 → null
  world(w, sig) {
    const Dc = [0, 0, 0];
    if (!this.P.chain(w, Dc, null)) return null;
    const u0 = this.applyR0([-w[0], -w[1], -w[2]]);
    const e0 = this.applyR0([-Dc[0], -Dc[1], -Dc[2]]);
    const Rho = this.Rho;
    const dx = -(Rho * e0[0] + sig * u0[0]), dy = -(Rho * e0[1] + sig * u0[1]);
    const a = this.AzL - Math.atan2(dy, dx);
    const ca = Math.cos(a), sa = Math.sin(a);
    const xi = [-(ca * e0[0] - sa * e0[1]), -(sa * e0[0] + ca * e0[1]), -e0[2]];
    const al = ((a / DEG) % 360.0 + 360.0) % 360.0;
    return fin3(xi) ? { xi, al } : null;
  }

  xi(w, sig) {
    const r = this.world(w, sig);
    return r === null ? null : r.xi;
  }
}

// ================================================================== 篱笆取样 (闭式参数化)
function fenceSamples(p, f, n) {
  const comps = [];
  const nv = f.L, u0 = p.U0;
  const c = f.C, k = f.Kk, K = p.K;
  const dt = dot(nv, u0);
  const nPerp = addSc(nv, u0, -dt);
  const npn = nrm(nPerp);
  if (npn > 1e-12) {
    const nh = scl(nPerp, 1.0 / npn);
    const ph = cross(u0, nh);
    const bof = X => (k - c * Math.sqrt(X * X + K) - dt * X) / npn;
    const d2f = X => { const b = bof(X); return 1.0 - X * X - b * b; };
    const NG = 4001;
    const Xg = new Array(NG), ok = new Array(NG);
    for (let i = 0; i < NG; i++) { Xg[i] = -1.0 + 2.0 * i / (NG - 1); ok[i] = d2f(Xg[i]) >= 0; }
    let ii = 0;
    while (ii < NG) {
      if (!ok[ii]) { ii++; continue; }
      let jj = ii;
      while (jj + 1 < NG && ok[jj + 1]) jj++;
      let lo = Xg[ii], hi = Xg[jj];
      if (ii > 0) {
        let a = Xg[ii - 1], b = Xg[ii];
        for (let t = 0; t < 60; t++) { const m = 0.5 * (a + b); if (d2f(m) >= 0) b = m; else a = m; }
        lo = b;
      }
      if (jj + 1 < NG) {
        let a = Xg[jj], b = Xg[jj + 1];
        for (let t = 0; t < 60; t++) { const m = 0.5 * (a + b); if (d2f(m) >= 0) a = m; else b = m; }
        hi = a;
      }
      const Xs = new Array(n);
      for (let t = 0; t < n; t++) {
        const s = n > 1 ? t / (n - 1) : 0.0;
        Xs[t] = lo + (hi - lo) * 0.5 * (1.0 - Math.cos(Math.PI * s));      // 端点加密
      }
      const pts = new Array(2 * n);
      let q = 0;
      for (let sgi = 0; sgi < 2; sgi++) {
        const sg = sgi === 0 ? +1.0 : -1.0;
        for (let t = 0; t < n; t++) {
          const X = sgi === 0 ? Xs[t] : Xs[n - 1 - t];
          const b_ = bof(X);
          const d = Math.sqrt(Math.max(1.0 - X * X - b_ * b_, 0.0));
          pts[q++] = [
            X * u0[0] + b_ * nh[0] + sg * d * ph[0],
            X * u0[1] + b_ * nh[1] + sg * d * ph[1],
            X * u0[2] + b_ * nh[2] + sg * d * ph[2]];
        }
      }
      comps.push(pts);
      ii = jj + 1;
    }
  } else {
    // ℓ ∥ u0: 篱笆是整圈 w·u0 = X*, X* 由 dot·X + c·S = k 解出
    const Xr = [];
    if (Math.abs(c) < 1e-14) {
      if (Math.abs(dt) > 1e-14) Xr.push(k / dt);
    } else {
      const qa = dt * dt - c * c, qb = -2 * dt * k, qc = k * k - c * c * K;
      const cand = [];
      if (Math.abs(qa) < 1e-14) { if (Math.abs(qb) > 1e-14) cand.push(-qc / qb); }
      else {
        const disc = qb * qb - 4 * qa * qc;
        if (disc >= 0) {
          const sq = Math.sqrt(disc);
          cand.push((-qb + sq) / (2 * qa)); cand.push((-qb - sq) / (2 * qa));
        }
      }
      for (const X of cand)
        if (Math.abs(dt * X + c * Math.sqrt(X * X + K) - k) < 1e-9) Xr.push(X);
    }
    for (const X of Xr) {
      if (Math.abs(X) >= 1.0) continue;
      const r = Math.sqrt(1.0 - X * X);
      let a1 = cross(u0, [1.0, 0.0, 0.0]);
      if (nrm(a1) < 1e-8) a1 = cross(u0, [0.0, 1.0, 0.0]);
      a1 = unit(a1);
      const a2 = cross(u0, a1);
      const m = 2 * n;
      const pts = new Array(m);
      for (let t = 0; t < m; t++) {
        const th = 2.0 * Math.PI * t / m;
        const ct = Math.cos(th), st = Math.sin(th);
        pts[t] = [
          X * u0[0] + r * (ct * a1[0] + st * a2[0]),
          X * u0[1] + r * (ct * a1[1] + st * a2[1]),
          X * u0[2] + r * (ct * a1[2] + st * a2[2])];
      }
      comps.push(pts);
    }
  }
  return comps;
}

// ================================================================== 篱笆上的一维求根
function fenceFootP(P, f, w) {
  w = w.slice();
  for (let it = 0; it < 10; it++) {
    const F = P.fval(f, w);
    if (Math.abs(F) < 1e-15) break;
    const g = P.fgrad(f, w);
    projT(g, w);
    const gg = dot(g, g);
    if (!(gg > 1e-30)) return null;
    w = unit(addSc(w, g, -F / gg));
  }
  return w;
}
const fenceFoot = (cx, f, w) => fenceFootP(cx.P, f, w);   // C# FenceFoot ≡ FenceFootP (FgradT = Fgrad + ProjT)

function fenceBisect(cx, f, fa, fb, fn) {
  let ra = fn(fa), rb = fn(fb);
  if (!fin(ra) || !fin(rb) || ra * rb > 0) return null;
  for (let it = 0; it < 60; it++) {
    const fm = fenceFoot(cx, f, unitMid(fa, fb));
    if (fm === null) return null;
    const rm = fn(fm);
    if (!fin(rm)) return null;
    if (ra * rm <= 0) { fb = fm; rb = rm; } else { fa = fm; ra = rm; }
    if (dist(fa, fb) < 1e-16) break;
  }
  return fenceFoot(cx, f, unitMid(fa, fb));
}

function refineBoundary(cx, f, wOk, wBad, isOk) {
  let lo = wOk, hi = wBad;
  for (let it = 0; it < 50; it++) {
    const mid = fenceFoot(cx, f, unitMid(lo, hi));
    if (mid === null) break;
    if (isOk(mid)) lo = mid; else hi = mid;
  }
  return lo;
}

// 沿篱笆样本找变号区间; 有效/无效样本相邻时先二分出有效域边界 (Ω 的角点) 并当作样本
function fenceRootCandidates(cx, f, comp, vals, isOk, fval) {
  const outL = [];
  for (let i = 0; i + 1 < comp.length; i++) {
    const a = vals[i], b = vals[i + 1];
    const fa = fin(a), fb = fin(b);
    if (fa && fb) {
      if (a * b <= 0) outL.push([comp[i], comp[i + 1]]);
    } else if (fa !== fb) {
      const iOk = fa ? i : i + 1, iBad = fa ? i + 1 : i;
      const wb = refineBoundary(cx, f, comp[iOk], comp[iBad], isOk);
      const vb = fval(wb);
      if (fin(vb) && vb * vals[iOk] <= 0) outL.push([comp[iOk], wb]);
    }
  }
  return outL;
}

// 3 方程 Newton: F1=F2=0, |w|=1。F(w, grad) 返回值并填梯度
function polishPair(F1, F2, w0) {
  const w = w0.slice();
  const g1 = [0, 0, 0], g2 = [0, 0, 0];
  const Jm = new Array(9).fill(0), r = [0, 0, 0], dw = [0, 0, 0];
  for (let it = 0; it < 60; it++) {
    r[0] = F1(w, g1); r[1] = F2(w, g2); r[2] = dot(w, w) - 1.0;
    if (!fin(r[0]) || !fin(r[1]) || !fin3(g1) || !fin3(g2)) return null;
    if (Math.max(Math.abs(r[0]), Math.max(Math.abs(r[1]), Math.abs(r[2]))) < 1e-15) break;
    for (let j = 0; j < 3; j++) { Jm[j] = g1[j]; Jm[3 + j] = g2[j]; Jm[6 + j] = 2 * w[j]; }
    if (!solve3(Jm, [-r[0], -r[1], -r[2]], dw)) return null;
    const n = nrm(dw);
    if (n > 0.3) { const s = 0.3 / n; dw[0] *= s; dw[1] *= s; dw[2] *= s; }
    w[0] += dw[0]; w[1] += dw[1]; w[2] += dw[2];
  }
  return fin3(w) ? w : null;
}

// ================================================================== 末端点
function newEnd(o) {
  return Object.assign({
    w: null, sigma: 0, fence: null, fenceLabel: null, missed: false, degenerate: false, startFail: false,
    used: 0, p: null, dir: null, alphaDeg: 0, inAlphaRange: true,
    fences: new Set(), fen: null, isLoop: false,          // 内部 (C# internal)
  }, o);
}

function findEnds(cx) {
  const outL = [];
  const P = cx.P;
  for (const fen of P.Fences) {
    ckCancel(cx);
    const others = [];
    for (const f of P.Fences) if (f !== fen) others.push(f);
    const othersOk = w => {
      for (const f of others) if (!(cx.margin(f, w) >= -1e-7)) return false;
      return true;
    };
    for (const comp of fenceSamples(P, fen, N_FENCE)) {
      const n = comp.length;
      const okm = new Array(n), ph = new Array(n), sgl = new Array(n), uzA = new Array(n), Aval = new Array(n);
      let allUz0 = true;
      for (let i = 0; i < n; i++) {
        const w = comp[i];
        okm[i] = othersOk(w);
        uzA[i] = -dot(w, cx.Nu);                  // u_z 不依赖链, 出射全反射一侧也有值
        const s = cx.scal(w);
        if (s !== null) {
          const uz = s[0], ez = s[1], p = s[2];
          const A = cx.Sl + cx.Rho * ez;
          ph[i] = A * A - 2 * cx.Rho * p * uz * A + (cx.Rho * cx.Rho - 1) * uz * uz;
          sgl[i] = -A / uz; Aval[i] = A;
        } else { ph[i] = NaN; sgl[i] = NaN; Aval[i] = NaN; }
        if (!(Math.abs(uzA[i]) < 1e-7)) allUz0 = false;
      }
      const vals = new Array(n);
      if (allUz0) {
        // 退化篱笆: 篱笆上 u_z≡0 (如 ζ=−90° 时入射面法向竖直), Φ=(S_l+ρe_z)² 触零不变号 → 末端条件退为 A=0
        for (let i = 0; i < n; i++) vals[i] = (okm[i] && fin(Aval[i])) ? Aval[i] : NaN;
        const fA = x => cx.aFun(x);
        const isOkA = x => othersOk(x) && fin(fA(x));
        for (const pair of fenceRootCandidates(cx, fen, comp, vals, isOkA, fA)) {
          const w = fenceBisect(cx, fen, pair[0], pair[1], fA);
          if (w === null || !(Math.abs(fA(w)) <= 1e-9)) continue;
          if (!othersOkTight(cx, fen, w)) continue;
          const s = cx.scal(w);
          if (s === null) continue;
          const pp = s[2];
          const d = cx.Rho * cx.Rho * pp * pp - cx.Rho * cx.Rho + 1;
          if (d < 0) continue;
          const sq = Math.sqrt(d);
          for (const sig of [-cx.Rho * pp + sq, -cx.Rho * pp - sq]) {
            if (!(sig > 0)) continue;
            if (hasEnd(outL, w, sig)) continue;
            const e = newEnd({ w, sigma: sig, fence: fen.name, fen, degenerate: true });
            e.fences.add(fen.name);
            outL.push(e);
          }
        }
        continue;
      }
      for (let i = 0; i < n; i++)
        vals[i] = (okm[i] && fin(ph[i]) && (sgl[i] > 0 || Math.abs(uzA[i]) < 1e-7)) ? ph[i] : NaN;
      const fPhi = x => cx.phi(x);
      const isOk = x => {
        if (!othersOk(x)) return false;
        if (!fin(cx.phi(x))) return false;
        const sgx = cx.sigmaOf(x);
        return fin(sgx) && sgx > 0;
      };
      for (const pair of fenceRootCandidates(cx, fen, comp, vals, isOk, fPhi)) {
        let w = fenceBisect(cx, fen, pair[0], pair[1], fPhi);
        if (w === null) continue;
        if (fen !== P.Exit) {
          const fenL = fen;
          const w2 = polishPair(
            (x, g) => { const gg = P.fgrad(fenL, x); g[0] = gg[0]; g[1] = gg[1]; g[2] = gg[2]; return P.fval(fenL, x); },
            (x, g) => cx.phiG(x, g), w);
          if (w2 !== null && Math.abs(P.fval(fen, w2)) < 1e-12 && Math.abs(cx.phi(w2)) < 1e-12 && dist(w2, w) < 1e-6)
            w = w2;
        }
        // 出射篱笆上 T_E=√(q_E²−K) 把 1e-16 的舍入放大到 ~1e-8
        const tolPhi = fen === P.Exit ? 1e-6 : 1e-9;
        if (Math.abs(P.fval(fen, w)) > 1e-10 || !(Math.abs(cx.phi(w)) < tolPhi)) continue;
        if (!othersOkTight(cx, fen, w)) continue;
        const sgm = cx.sigmaOf(w);
        if (!(fin(sgm) && sgm > 0)) continue;
        if (hasEnd(outL, w, sgm)) continue;
        const e = newEnd({ w, sigma: sgm, fence: fen.name, fen });
        e.fences.add(fen.name);
        outL.push(e);
      }
    }
  }
  return mergeCornerEnds(outL);
}

function othersOkTight(cx, fen, w) {
  for (const f of cx.P.Fences)
    if (f !== fen && cx.margin(f, w) < -1e-8) return false;
  return true;
}

function hasEnd(L, w, sig) {
  for (const r of L) if (dist(w, r.w) < 1e-8 && Math.abs(sig - r.sigma) < 1e-6) return true;
  return false;
}

// 角点: 同一个 (w,σ) 同时是两条篱笆的末端 → 合并
function mergeCornerEnds(ends) {
  const merged = [];
  for (const e of ends) {
    let hit = null;
    for (const m of merged)
      if (dist(m.w, e.w) < 1e-7 && Math.abs(m.sigma - e.sigma) < 1e-6) { hit = m; break; }
    if (hit === null) { merged.push(e); continue; }
    for (const nm of e.fences) hit.fences.add(nm);
    if (hit.fences.has('exit')) {
      hit.fence = 'exit';
      if (e.fence === 'exit') hit.fen = e.fen;
    }
  }
  for (const m of merged) {
    const names = [...m.fences];
    names.sort();                                // UTF-16 码元序 = StringComparer.Ordinal
    m.fenceLabel = names.join('+');
  }
  return merged;
}

// ================================================================== (w,σ) 切向 / 校正
// 曲线切向 (单位 4-向量) = [J;(w,0)] 的零空间; 两行近乎平行时退到 t_w = w×∇Φ
function tangent(cx, w, sig) {
  const J = new Array(8);
  if (!cx.jac(w, sig, J)) return null;
  const A = new Array(12);
  for (let j = 0; j < 4; j++) { A[j] = J[j]; A[4 + j] = J[4 + j]; }
  A[8] = w[0]; A[9] = w[1]; A[10] = w[2]; A[11] = 0;
  const nvec = [0, 0, 0, 0];
  for (let c = 0; c < 4; c++) {
    const cols = c === 0 ? [1, 2, 3] : c === 1 ? [0, 2, 3] : c === 2 ? [0, 1, 3] : [0, 1, 2];
    const c0 = cols[0], c1 = cols[1], c2 = cols[2];
    const det = A[c0] * (A[4 + c1] * A[8 + c2] - A[4 + c2] * A[8 + c1])
              - A[c1] * (A[4 + c0] * A[8 + c2] - A[4 + c2] * A[8 + c0])
              + A[c2] * (A[4 + c0] * A[8 + c1] - A[4 + c1] * A[8 + c0]);
    nvec[c] = ((c & 1) === 0 ? 1.0 : -1.0) * det;
  }
  const nn = Math.sqrt(nvec[0] * nvec[0] + nvec[1] * nvec[1] + nvec[2] * nvec[2] + nvec[3] * nvec[3]);
  const r1 = Math.sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2] + A[3] * A[3]);
  const r2 = Math.sqrt(A[4] * A[4] + A[5] * A[5] + A[6] * A[6] + A[7] * A[7]);
  if (fin(nn) && r1 > 0 && r2 > 0 && nn > 1e-7 * r1 * r2)
    return [nvec[0] / nn, nvec[1] / nn, nvec[2] / nn, nvec[3] / nn];
  const g = [0, 0, 0];
  if (!fin(cx.phiG(w, g))) return null;
  projT(g, w);
  let tw = cross(w, g);
  const n = nrm(tw);
  if (!fin(n) || n < 1e-14) return null;
  tw = scl(tw, 1.0 / n);
  const uz = J[3];
  let ts = Math.abs(uz) > 1e-9 ? -(J[0] * tw[0] + J[1] * tw[1] + J[2] * tw[2]) / uz : 0.0;
  if (!fin(ts)) ts = 0.0;
  const l = Math.sqrt(1.0 + ts * ts);
  return [tw[0] / l, tw[1] / l, tw[2] / l, ts / l];
}

// Gauss–Newton 最小范数校正到 F1=F2=0, |w|=1。→ { w, sig } 或 null
function gnCorrect(cx, wIn, sigIn) {
  let w = wIn.slice(), sig = sigIn;
  const J = new Array(8), AAt = new Array(9), y = [0, 0, 0];
  for (let it = 0; it < 12; it++) {
    const f = cx.f12(w, sig);
    if (f === null) return null;
    const f1 = f[0], f2 = f[1];
    if (Math.max(Math.abs(f1), Math.abs(f2)) < 1e-13) break;
    if (!cx.jac(w, sig, J)) return null;
    const A = new Array(12);
    for (let j = 0; j < 4; j++) { A[j] = J[j]; A[4 + j] = J[4 + j]; }
    A[8] = w[0]; A[9] = w[1]; A[10] = w[2]; A[11] = 0;
    for (let i = 0; i < 3; i++)
      for (let k = 0; k < 3; k++) {
        let s = 0;
        for (let j = 0; j < 4; j++) s += A[i * 4 + j] * A[k * 4 + j];
        AAt[i * 3 + k] = s;
      }
    if (!solve3(AAt, [-f1, -f2, 0.0], y)) return null;
    const dx = [0, 0, 0, 0];
    for (let j = 0; j < 4; j++) dx[j] = A[j] * y[0] + A[4 + j] * y[1] + A[8 + j] * y[2];
    if (!fin(dx[0]) || !fin(dx[1]) || !fin(dx[2]) || !fin(dx[3])) return null;
    w = unit([w[0] + dx[0], w[1] + dx[1], w[2] + dx[2]]);
    sig += dx[3];
  }
  const g = cx.f12(w, sig);
  if (g === null || Math.max(Math.abs(g[0]), Math.abs(g[1])) > 1e-10) return null;
  return { w, sig };
}

// 括根校正 (出射篱笆附近): 沿 n_c 在 [−M, M] 内找 Φ=0, 只在 Ω 内取值; 取离 0 最近的根。→ { w, m } 或 null
function bracketPhi(cx, w1, nc, Mr) {
  let wRoot = null, mRoot = NaN;
  const ms = [];
  for (let i = 0; i < 33; i++) ms.push(-Mr + 2.0 * Mr * i / 32.0);
  const f = m => cx.phiIn(unit(addSc(w1, nc, m)));
  let vals = [];
  for (const m of ms) vals.push(f(m));
  const extra = [];
  for (let i = 0; i + 1 < ms.length; i++) {
    const aOk = fin(vals[i]), bOk = fin(vals[i + 1]);
    if (aOk === bOk) continue;
    let lo = aOk ? ms[i] : ms[i + 1], hi = aOk ? ms[i + 1] : ms[i];
    for (let t = 0; t < 50; t++) {
      const mid = 0.5 * (lo + hi);
      if (fin(f(mid))) lo = mid; else hi = mid;
    }
    extra.push(lo);
  }
  if (extra.length > 0) {
    for (const e of extra) ms.push(e);
    ms.sort((a, b) => a - b);
    vals = [];
    for (const m of ms) vals.push(f(m));
  }
  let found = false;
  for (let i = 0; i + 1 < ms.length; i++) {
    const a = vals[i], b = vals[i + 1];
    if (!fin(a) || !fin(b) || a * b > 0) continue;
    let lo = ms[i], hi = ms[i + 1], fa = a;
    let ok = true;
    for (let t = 0; t < 55; t++) {
      const mid = 0.5 * (lo + hi), fm = f(mid);
      if (!fin(fm)) { ok = false; break; }
      if (fa * fm <= 0) hi = mid; else { lo = mid; fa = fm; }
    }
    if (!ok) continue;
    const mr = 0.5 * (lo + hi);
    if (!found || Math.abs(mr) < Math.abs(mRoot)) {
      found = true; mRoot = mr; wRoot = unit(addSc(w1, nc, mr));
    }
  }
  return found ? { w: wRoot, m: mRoot } : null;
}

// ================================================================== 追踪
function newRawArc() {
  return {
    W: [], Sig: [], Xi: [], From: -1, To: -1, Closed: false, ViaLoop: false, Stop: 'max_steps',
    FoldIdx: [], FoldW: [], FoldSig: [],
  };
}

class Tracer {
  constructor(c, seeds) { this.cx = c; this.Seeds = seeds; this.fenE = c.P.Exit; }

  nearSeed(w, sig, tol, fence, exclude) {
    let best = tol, bk = -1;
    const Seeds = this.Seeds;
    for (let k = 0; k < Seeds.length; k++) {
      const s = Seeds[k];
      if (k === exclude || s.isLoop) continue;
      if (fence != null && !s.fences.has(fence)) continue;
      if (Math.abs(s.sigma - sig) > DSIG) continue;
      const d = dist(w, s.w);
      if (d < best) { best = d; bk = k; }
    }
    return bk;
  }

  addSeed(w, sig, fen) {
    const e = newEnd({ w, sigma: sig, fence: fen.name, fenceLabel: fen.name + '*', fen, missed: true });
    e.fences.add(fen.name);
    this.Seeds.push(e);
    return this.Seeds.length - 1;
  }

  wallStart(s) {
    const cx = this.cx, fenE = this.fenE;
    let outL = [];
    if (fenE == null) return outL;
    const w = s.w;
    const nIn = unit(scl(cx.fgradT(fenE, w), fenE.Side));
    const tau = unit(cross(w, nIn));
    for (const dl of [+1.0, -1.0])
      for (const ds of [2e-4, 1e-3, 4e-3]) {
        const f1 = fenceFoot(cx, fenE, unit(addSc(w, tau, ds * dl)));
        if (f1 === null) continue;
        const br = bracketPhi(cx, f1, nIn, 0.05);
        if (br === null) continue;
        const wr = br.w, mr = br.m;
        const sig = cx.sigmaOf(wr);
        if (!(fin(sig) && sig > 0)) continue;
        const t4 = [wr[0] - w[0], wr[1] - w[1], wr[2] - w[2], sig - s.sigma];
        const l = Math.sqrt(t4[0] * t4[0] + t4[1] * t4[1] + t4[2] * t4[2] + t4[3] * t4[3]);
        if (!(l > 0)) continue;
        for (let j = 0; j < 4; j++) t4[j] /= l;
        outL.push({ W: wr, Sig: sig, T: t4, M: Math.abs(mr) });
        break;
      }
    if (outL.length === 2 && !(outL[0].M < 1e-6 && outL[1].M < 1e-6)) {
      const keep = outL[0].M <= outL[1].M ? outL[0] : outL[1];
      outL = [keep];
    }
    return outL;
  }

  // 先切向 GN 起步 (首段天球弦长 ≤ MAX_SKY); 出射端 (相切) GN 常因 √ 病态失败, 再退到贴壁起步
  startDirs(s) {
    const cx = this.cx, fenE = this.fenE;
    let outL = [];
    const w = s.w, sig = s.sigma;
    const xi0 = cx.xi(w, sig);
    const t = tangent(cx, w, sig);
    if (t !== null) {
      for (const sgn of [+1.0, -1.0]) {
        const dirn = [sgn * t[0], sgn * t[1], sgn * t[2], sgn * t[3]];
        for (const st of [STEP0 / 4, STEP0 / 32, STEP0 / 256, STEP0 / 2048]) {
          const w1 = unit([w[0] + st * dirn[0], w[1] + st * dirn[1], w[2] + st * dirn[2]]);
          if (!cx.inOmega(w1, 1e-9)) continue;
          const g = gnCorrect(cx, w1, sig + st * dirn[3]);
          if (g === null) continue;
          const wc = g.w, sc = g.sig;
          if (!cx.inOmega(wc, 1e-9) || !(sc > 0)) continue;
          if (dist(wc, w) < 1e-12 && Math.abs(sc - sig) < 1e-12) continue;
          if (ang(cx.xi(wc, sc), xi0) > MAX_SKY && st > STEP0 / 2048) continue;
          outL.push({ W: wc, Sig: sc, T: dirn, M: 0 });
          break;
        }
      }
    }
    if (outL.length === 0 && fenE != null && (s.fences.has('exit') || cx.exitMargin(w) < 3 * WALL))
      outL = this.wallStart(s);
    if (!s.isLoop && outL.length === 2) {
      const d0 = ang(cx.xi(outL[0].W, outL[0].Sig), xi0), d1 = ang(cx.xi(outL[1].W, outL[1].Sig), xi0);
      let dd = 0; for (let j = 0; j < 4; j++) dd += outL[0].T[j] * outL[1].T[j];
      if (dd < -0.5) outL = [d0 <= d1 ? outL[0] : outL[1]];
    }
    return outL;
  }

  trace(si, sd) {
    const cx = this.cx, fenE = this.fenE, Seeds = this.Seeds;
    const s0 = Seeds[si];
    const arc = newRawArc(); arc.From = si;
    arc.W.push(s0.w.slice()); arc.Sig.push(s0.sigma); arc.Xi.push(cx.xi(s0.w, s0.sigma));
    arc.W.push(sd.W.slice()); arc.Sig.push(sd.Sig); arc.Xi.push(cx.xi(sd.W, sd.Sig));
    let w = sd.W.slice(), sig = sd.Sig, tang = sd.T;
    let step = STEP0 / 4;
    let armed = false;
    for (let it = 0; it < MAX_STEPS; it++) {
      if ((it & 63) === 0) ckCancel(cx);
      const t = tangent(cx, w, sig);
      if (t === null) { arc.Stop = 'tangent_fail'; break; }
      let dtt = 0; for (let j = 0; j < 4; j++) dtt += t[j] * tang[j];
      if (dtt < 0) for (let j = 0; j < 4; j++) t[j] = -t[j];

      let accepted = false, ev = 0;         // ev: 1 = leave, 2 = exit_end
      let wc = null, evA = null, evB = null, xi1 = null, sc = NaN, evSig = NaN, ch = 0;
      let failed = false;
      while (true) {
        const w1 = unit([w[0] + step * t[0], w[1] + step * t[1], w[2] + step * t[2]]);
        const s1 = sig + step * t[3];
        const gc = gnCorrect(cx, w1, s1);
        let good = gc !== null && fin(gc.sig)
                   && dist(gc.w, w1) < 2.0 * step + 1e-9 && Math.abs(gc.sig - s1) < 0.5;
        if (good) {
          wc = gc.w; sc = gc.sig;
          if (!cx.inOmega(wc, 1e-9)) {
            const mf = cx.minMarginFence(wc);
            if (mf === fenE) good = false;      // 出射篱笆只能相切地到达, 交给括根分支
            else { ev = 1; evA = wc; evSig = sc; break; }
          }
        }
        if (!good) {
          const nearExit = fenE != null && (cx.exitMargin(w) < 3 * WALL || cx.exitMargin(w1) < 3 * WALL);
          let haveR = false;
          if (nearExit) {
            const nc = cross(w1, [t[0], t[1], t[2]]);
            if (nrm(nc) > 1e-12) {
              const br = bracketPhi(cx, w1, unit(nc), Math.max(2.0 * step, 1e-4));
              if (br !== null) {
                const sr = cx.sigmaOf(br.w);
                if (fin(sr) && Math.abs(sr - s1) < 0.5) { haveR = true; wc = br.w; sc = sr; }
              }
            }
          }
          if (!haveR) {
            if (nearExit && cx.exitMargin(w) < 1e-7) {
              // 当前点已经贴在篱笆上: 这就是末端
              const f0 = fenceFoot(cx, fenE, w);
              if (f0 !== null) { ev = 2; evA = f0; evB = null; break; }
            }
            if (nearExit && step < STEP0 / 2048) {
              // 已缩到最小步长仍括不到根: 篱笆上 Φ 在前后足点间变号 → 相切末端
              const f0 = fenceFoot(cx, fenE, w);
              const f1 = fenceFoot(cx, fenE, w1);
              if (f0 !== null && f1 !== null) {
                const pa = cx.phi(f0), pb = cx.phi(f1);
                if (fin(pa) && fin(pb) && pa * pb <= 0) { ev = 2; evA = f0; evB = f1; break; }
              }
            }
            step *= 0.5;
            if (step < STEP0 / 4096) { arc.Stop = 'correct_fail'; failed = true; break; }
            continue;
          }
        }
        xi1 = cx.xi(wc, sc);
        ch = ang(xi1, arc.Xi[arc.Xi.length - 1]);
        if (ch > MAX_SKY && step > STEP0 / 4096) { step *= 0.5; continue; }
        accepted = true;
        break;
      }
      if (failed) break;

      if (ev === 2) {
        const wE = evB === null ? evA : fenceBisect(cx, fenE, evA, evB, x => cx.phi(x));
        if (wE === null) { arc.Stop = 'exit_end_fail'; break; }
        const sE = cx.sigmaOf(wE);
        let k = this.nearSeed(wE, sE, 3e-3, 'exit', arc.W.length < 6 ? si : -1);
        if (k < 0) { k = this.addSeed(wE, sE, fenE); arc.Stop = 'missed_end'; } else arc.Stop = 'end';
        arc.W.push(Seeds[k].w.slice()); arc.Sig.push(Seeds[k].sigma); arc.Xi.push(cx.xi(Seeds[k].w, Seeds[k].sigma));
        arc.To = k;
        break;
      }
      if (ev === 1) {
        let loW = w.slice(), hiW = evA, loS = sig, hiS = evSig;
        for (let q = 0; q < 60; q++) {
          const wm = unitMid(loW, hiW), sm = 0.5 * (loS + hiS);
          const g = gnCorrect(cx, wm, sm);
          if (g === null || !cx.inOmega(g.w, 1e-9)) { hiW = wm; hiS = sm; }
          else { loW = g.w; loS = g.sig; }
          if (dist(loW, hiW) < 1e-15) break;
        }
        const mf = cx.minMarginFence(loW);
        let k = this.nearSeed(loW, loS, 3e-3, mf.name, arc.W.length < 6 ? si : -1);
        if (k < 0) { k = this.addSeed(loW, loS, mf); arc.Stop = 'missed_end'; } else arc.Stop = 'end';
        arc.W.push(Seeds[k].w.slice()); arc.Sig.push(Seeds[k].sigma); arc.Xi.push(cx.xi(Seeds[k].w, Seeds[k].sigma));
        arc.To = k;
        break;
      }
      if (!accepted) break;

      w = wc; sig = sc; tang = t;
      arc.W.push(w.slice()); arc.Sig.push(sig); arc.Xi.push(xi1);
      if (ch < MIN_SKY) step = Math.min(step * 1.5, STEP0);
      if (dist(w, arc.W[1]) > 4.0 * STEP0) armed = true;
      if (arc.W.length > 4) {
        let k = this.nearSeed(w, sig, Math.max(1.5 * step, 1e-4), null, si);
        if (k >= 0 && Seeds[k].fence === 'exit' && dist(w, Seeds[k].w) > 1e-6) k = -1;
        if (k >= 0 && ang(xi1, cx.xi(Seeds[k].w, Seeds[k].sigma)) < MAX_SKY) {
          arc.W.push(Seeds[k].w.slice()); arc.Sig.push(Seeds[k].sigma); arc.Xi.push(cx.xi(Seeds[k].w, Seeds[k].sigma));
          arc.To = k; arc.Stop = 'end';
          break;
        }
        if (armed && dist(w, arc.W[1]) < 1.5 * step && Math.abs(sig - arc.Sig[1]) < DSIG) {
          // 收口那一段也要守天球弦长上限
          if (ang(xi1, arc.Xi[1]) > MAX_SKY && step > STEP0 / 4096) step *= 0.5;
          else {
            arc.W.push(arc.W[1].slice()); arc.Sig.push(arc.Sig[1]); arc.Xi.push(arc.Xi[1]);
            arc.Closed = true; arc.To = si; arc.Stop = 'closed';
            break;
          }
        }
      }
    }
    this.markFolds(arc);
    return arc;
  }

  // 折叠点 (两个落位重合, D'=0): D' 的局部极小 (<1e-3) → 三方程 Newton (D', Φ) 抛光
  markFolds(arc) {
    const cx = this.cx;
    const n = arc.W.length;
    if (n < 5) return;
    const Dv = new Array(n);
    for (let i = 0; i < n; i++) Dv[i] = cx.dprime(arc.W[i]);
    for (let i = 1; i + 1 < n; i++) {
      if (!(fin(Dv[i]) && Dv[i] < 1e-3 && Dv[i] <= Dv[i - 1] && Dv[i] <= Dv[i + 1])) continue;
      const w = polishPair((x, g) => cx.dprimeG(x, g), (x, g) => cx.phiG(x, g), arc.W[i]);
      if (w === null) continue;
      if (!(Math.abs(cx.dprime(w)) <= 1e-9) || !(Math.abs(cx.phi(w)) <= 1e-9)) continue;
      if (dist(w, arc.W[i]) > 2.0 * STEP0) continue;
      arc.FoldIdx.push(i); arc.FoldW.push(w); arc.FoldSig.push(cx.sigmaOf(w));
    }
  }
}

function traceArcs(cx, seeds) {
  const tr = new Tracer(cx, seeds);
  const arcs = [];
  for (let si = 0; si < seeds.length; si++) {          // seeds 可能在追踪中增长 (漏端), 同 C#
    const s = seeds[si];
    if (s.isLoop || s.used !== 0) continue;
    const dirs = tr.startDirs(s);
    if (dirs.length === 0) s.startFail = true;
    for (const d of dirs) {
      if (s.used > 0) break;                          // 一个末端只发一条弧
      const a = tr.trace(si, d);
      s.used++;
      if (a.To >= 0 && a.To !== si) seeds[a.To].used++;
      arcs.push(a);
    }
  }
  return arcs;
}

// 没有末端点的闭环: 网格 Φ 变号 → 远离已追踪弧的候选 → 双向追踪并合并
function findLoops(cx, seeds, arcs) {
  const d = 0.5;
  const nLat = roundEven(180.0 / d), nLon = roundEven(360.0 / d) + 1;
  const W = new Array(nLat * nLon), F = new Float64Array(nLat * nLon);
  const cLon = new Array(nLon), sLon = new Array(nLon);
  for (let j = 0; j < nLon; j++) { const lon = (-180.0 + j * d) * DEG; cLon[j] = Math.cos(lon); sLon[j] = Math.sin(lon); }
  for (let i = 0; i < nLat; i++) {
    if ((i & 15) === 0) ckCancel(cx);
    const lat = (-90.0 + d / 2 + i * d) * DEG, cl = Math.cos(lat), sl = Math.sin(lat);
    for (let j = 0; j < nLon; j++) {
      const w = [cl * cLon[j], cl * sLon[j], sl];
      W[i * nLon + j] = w;
      F[i * nLon + j] = cx.inOmega(w, 1e-6) ? cx.phi(w) : NaN;
    }
  }
  const cand = [];
  for (let i = 0; i + 1 < nLat; i++)
    for (let j = 0; j < nLon; j++) {
      const a = F[i * nLon + j], b = F[(i + 1) * nLon + j];
      if (fin(a) && fin(b) && a * b < 0) cand.push(unitMid(W[i * nLon + j], W[(i + 1) * nLon + j]));
    }
  for (let i = 0; i < nLat; i++)
    for (let j = 0; j + 1 < nLon; j++) {
      const a = F[i * nLon + j], b = F[i * nLon + j + 1];
      if (fin(a) && fin(b) && a * b < 0) cand.push(unitMid(W[i * nLon + j], W[i * nLon + j + 1]));
    }

  const tr = new Tracer(cx, seeds);
  const found = [];
  const haveW = [], haveS = [];
  for (const a of arcs) { for (const x of a.W) haveW.push(x); for (const x of a.Sig) haveS.push(x); }
  const tolW = 1.0 * DEG;
  for (const w0 of cand) {
    ckCancel(cx);
    const sig0 = cx.sigmaOf(w0);
    if (!fin(sig0)) continue;
    const g = gnCorrect(cx, w0, sig0);
    if (g === null) continue;
    const cw = g.w, cs = g.sig;
    if (!cx.inOmega(cw, 1e-8) || !(cs > 0)) continue;
    let close = false;
    for (let q = 0; q < haveW.length; q++)
      if (Math.abs(haveS[q] - cs) < 0.1 && dist(haveW[q], cw) < tolW) { close = true; break; }
    if (close) continue;
    const seed = newEnd({ w: cw, sigma: cs, fence: 'loop', fenceLabel: 'loop', isLoop: true });
    seeds.push(seed);
    const si = seeds.length - 1;
    const dirs = tr.startDirs(seed);
    if (dirs.length === 0) continue;
    const a1 = tr.trace(si, dirs[0]);
    if (a1.To >= 0 && a1.To !== si) seeds[a1.To].used++;
    let res;
    if (a1.Closed || dirs.length === 1) {
      a1.ViaLoop = true;
      a1.From = -1;                                // 起点是网格种子, 不是末端
      if (a1.Closed) a1.To = -1;
      res = a1;
    } else {
      const a2 = tr.trace(si, dirs[1]);
      if (a2.To >= 0 && a2.To !== si) seeds[a2.To].used++;
      res = newRawArc();
      res.From = a1.To; res.To = a2.To; res.ViaLoop = true; res.Stop = a1.Stop + '|' + a2.Stop;
      const n1 = a1.W.length;
      for (let q = n1 - 1; q >= 0; q--) { res.W.push(a1.W[q]); res.Sig.push(a1.Sig[q]); res.Xi.push(a1.Xi[q]); }
      for (let q = 2; q < a2.W.length; q++) { res.W.push(a2.W[q]); res.Sig.push(a2.Sig[q]); res.Xi.push(a2.Xi[q]); }
      for (const fi of a1.FoldIdx) res.FoldIdx.push(n1 - 1 - fi);
      for (const fi of a2.FoldIdx) res.FoldIdx.push(n1 - 2 + fi);
      for (const x of a1.FoldW) res.FoldW.push(x); for (const x of a2.FoldW) res.FoldW.push(x);
      for (const x of a1.FoldSig) res.FoldSig.push(x); for (const x of a2.FoldSig) res.FoldSig.push(x);
    }
    found.push(res);
    for (const x of res.W) haveW.push(x); for (const x of res.Sig) haveS.push(x);
  }
  return found;
}

// ================================================================== 独立核对 (与闭式链无关的逐面 snell)
function snellTrace(N, refN, u, dOut) {
  const eta = 1.0 / refN;
  const nIn = N[0];
  let cosI = -dot(u, nIn);
  if (cosI < -1e-9) return false;
  if (cosI < 0) cosI = 0;
  const s2 = eta * eta * (1 - cosI * cosI);
  if (s2 > 1) return false;
  const c1 = eta * cosI - Math.sqrt(1 - s2);
  let v = [eta * u[0] + c1 * nIn[0], eta * u[1] + c1 * nIn[1], eta * u[2] + c1 * nIn[2]];
  for (let k = 1; k + 1 < N.length; k++) {
    const dn = dot(v, N[k]);
    if (dn <= 0) return false;
    v = addSc(v, N[k], -2 * dn);
  }
  const nO = N[N.length - 1];
  const cosI2 = dot(v, nO);
  if (cosI2 <= 0) return false;
  const s22 = refN * refN * (1 - cosI2 * cosI2);
  if (s22 > 1 + 1e-9) return false;
  const c2 = Math.sqrt(Math.max(0.0, 1 - s22)) - refN * cosI2;
  const o = [refN * v[0] + c2 * nO[0], refN * v[1] + c2 * nO[1], refN * v[2] + c2 * nO[2]];
  const l = nrm(o);
  if (!(l > 1e-12)) return false;
  dOut[0] = o[0] / l; dOut[1] = o[1] / l; dOut[2] = o[2] / l;
  return true;
}

// ==================================================================
// 二自由度族 (Column: α+θ 自由, ζ 定 / Lowitz: α+ζ 自由, θ 定)
// 落位: b(w) = σ·w + ρ·D_c(w) 是晶体系"眼→灯"单位方向, 条件 ν·b(w) = sin Σ;
//   二自由度 ν(τ) = c0 + c1 cosτ + c2 sinτ ⇒ A cosτ + B sinτ = S_l − C 闭式解 τ: 壳流形是 w-球上的 4 叶图。
// 轮廓 = 焦散 (Ĵ ≡ T·J = 0, 扩展状态延拓) + 各条篱笆的像 (闭式)。
// ==================================================================
export function fenceKind(name) {
  if (name === 'entry') return 1;
  if (name === 'exit') return 2;
  if (name.startsWith('seam')) return 3;
  if (name.startsWith('blue')) return 4;
  return 0;
}

class Ctx2 {
  constructor() {
    this.P = null;
    this.column = true;                      // Column: ζ0 固定;  Lowitz: θ0 固定
    this.Fixed = 0;                          // 弧度
    this.C0 = null; this.C1 = null; this.C2 = null; this.Gen = null;
    this.Sl = 0; this.Rho = 0;
    this.Eye = null;
    this.Radius = 0; this.AzL = 0;
    this.cancelled = never;
  }

  nu(tau) {
    const c = Math.cos(tau), s = Math.sin(tau), C0 = this.C0, C1 = this.C1, C2 = this.C2;
    return [C0[0] + C1[0] * c + C2[0] * s, C0[1] + C1[1] * c + C2[1] * s, C0[2] + C1[2] * c + C2[2] * s];
  }
  dnu(tau) {
    const c = Math.cos(tau), s = Math.sin(tau), C1 = this.C1, C2 = this.C2;
    return [-C1[0] * s + C2[0] * c, -C1[1] * s + C2[1] * c, -C1[2] * s + C2[2] * c];
  }

  // R1(τ) = Ry(ζ)·Rz(θ), 第二自由度取 τ (R 行主序 9 元)
  r1(tau, R) {
    const z = this.column ? this.Fixed : tau, t = this.column ? tau : this.Fixed;
    const cz = Math.cos(z), sz = Math.sin(z), ct = Math.cos(t), st = Math.sin(t);
    R[0] = cz * ct; R[1] = -cz * st; R[2] = sz;
    R[3] = st; R[4] = ct; R[5] = 0;
    R[6] = -sz * ct; R[7] = sz * st; R[8] = cz;
  }

  // Ĵ = T·J。w 不必单位, T 独立 —— 供差分与扩展追踪用
  jhat(w, T, tau) {
    const Dc = [0, 0, 0], M0 = new Array(9), g = [0, 0, 0];
    const qE = this.P.chainT(w, T, Dc, M0, g);
    const a1 = this.nu(tau), a2 = this.Gen;
    const wa1 = cross(w, a1), wa2 = cross(w, a2);
    const p1 = cross(a1, Dc), p2 = cross(a2, Dc);
    for (let i = 0; i < 3; i++) {
      p1[i] += M0[i * 3] * wa1[0] + M0[i * 3 + 1] * wa1[1] + M0[i * 3 + 2] * wa1[2];
      p2[i] += M0[i * 3] * wa2[0] + M0[i * 3 + 1] * wa2[1] + M0[i * 3 + 2] * wa2[2];
    }
    const c1 = dot(g, wa1), c2 = dot(g, wa2);
    const Ja = dot(Dc, cross(p1, p2));
    const Jb = c2 * dot(Dc, cross(p1, this.P.UE)) - c1 * dot(Dc, cross(p2, this.P.UE));
    return T * Ja - qE * Jb;
  }

  // 第 (a,c) 叶上的闭式落位: a=±1 选 σ 根, c=±1 选 τ 根。onExit: w 在出射篱笆上 (T 取 0)。→ [T, τ, σ] 或 null
  sheet(w, a, c, onExit) {
    const P = this.P, Rho = this.Rho;
    const Dc = [0, 0, 0], g = [0, 0, 0];
    const qE = P.chainT(w, 0.0, Dc, null, g);
    if (qE > 0) return null;
    const t2 = qE * qE - P.K;
    let T;
    if (onExit || (t2 < 0 && t2 > -1e-9)) T = 0.0;
    else if (t2 < 0) return null;
    else T = Math.sqrt(t2);
    Dc[0] -= T * P.UE[0]; Dc[1] -= T * P.UE[1]; Dc[2] -= T * P.UE[2];
    const p = dot(w, Dc);
    const disc = Rho * Rho * p * p - Rho * Rho + 1;
    if (disc < 0) return null;
    const sig = -Rho * p + a * Math.sqrt(disc);
    if (!(sig > 0)) return null;
    const b = [sig * w[0] + Rho * Dc[0], sig * w[1] + Rho * Dc[1], sig * w[2] + Rho * Dc[2]];
    const A = dot(b, this.C1), B = dot(b, this.C2), Cc = dot(b, this.C0);
    const Rr = Math.sqrt(A * A + B * B);
    if (!(Rr > 1e-13)) return null;
    const kap = (this.Sl - Cc) / Rr;
    if (!(Math.abs(kap) <= 1.0)) return null;
    const tau = Math.atan2(B, A) + c * Math.acos(kap);
    return [T, tau, sig];
  }

  othersOk(w, skip, tol) {
    const P = this.P;
    for (const f of P.Fences)
      if (f !== skip && !(f.Side * P.fval(f, w) >= -tol)) return false;
    return true;
  }

  // (w,T,τ,σ) → 天球单位方向 ξ (总是返回数组; 调用方用 fin3 判有效, 同 C# World 的返回值)
  world(x) {
    const w = unit([x[0], x[1], x[2]]);
    const Dc = [0, 0, 0], g = [0, 0, 0];
    this.P.chainT(w, x[3], Dc, null, g);
    const R = new Array(9); this.r1(x[4], R);
    const u0 = [0, 0, 0], e0 = [0, 0, 0];
    for (let i = 0; i < 3; i++) {
      u0[i] = -(R[i * 3] * w[0] + R[i * 3 + 1] * w[1] + R[i * 3 + 2] * w[2]);
      e0[i] = -(R[i * 3] * Dc[0] + R[i * 3 + 1] * Dc[1] + R[i * 3 + 2] * Dc[2]);
    }
    const Rho = this.Rho;
    const dx = -(Rho * e0[0] + x[5] * u0[0]), dy = -(Rho * e0[1] + x[5] * u0[1]);
    const al = this.AzL - Math.atan2(dy, dx), ca = Math.cos(al), sa = Math.sin(al);
    return [-(ca * e0[0] - sa * e0[1]), -(sa * e0[0] + ca * e0[1]), -e0[2]];
  }

  alphaOf(x) {
    const w = unit([x[0], x[1], x[2]]);
    const Dc = [0, 0, 0], g = [0, 0, 0];
    this.P.chainT(w, x[3], Dc, null, g);
    const R = new Array(9); this.r1(x[4], R);
    let ux = 0, uy = 0, ex = 0, ey = 0;
    for (let j = 0; j < 3; j++) {
      ux -= R[j] * w[j]; uy -= R[3 + j] * w[j];
      ex -= R[j] * Dc[j]; ey -= R[3 + j] * Dc[j];
    }
    const Rho = this.Rho;
    return this.AzL - Math.atan2(-(Rho * ey + x[5] * uy), -(Rho * ex + x[5] * ux));
  }
}

// ------------------------------------------------------------------ 篱笆像 (显式)
function fenceCuts(cx, fen, outRuns) {
  const onExit = fen === cx.P.Exit;
  for (const comp of fenceSamples(cx.P, fen, 1200)) {
    ckCancel(cx);
    const n = comp.length;
    const okW = new Array(n);
    for (let i = 0; i < n; i++) okW[i] = cx.othersOk(comp[i], fen, 1e-9);
    for (let a = +1; a >= -1; a -= 2)
      for (let c = +1; c >= -1; c -= 2) {
        const aa = a, cc = c;
        const stateAt = w => {
          if (!cx.othersOk(w, fen, 1e-9)) return null;
          const s = cx.sheet(w, aa, cc, onExit);
          if (s === null) return null;
          return [w[0], w[1], w[2], s[0], s[1], s[2]];
        };
        const st = new Array(n);
        for (let i = 0; i < n; i++) st[i] = okW[i] ? stateAt(comp[i]) : null;
        // 有效样本的连续段; 段两端向无效侧二分出精确端点
        let i0 = 0;
        while (i0 < n) {
          if (st[i0] === null) { i0++; continue; }
          let j0 = i0;
          while (j0 + 1 < n && st[j0 + 1] !== null) j0++;
          const run = [];
          if (i0 > 0) { const e = refineRunEnd(cx, fen, comp[i0], comp[i0 - 1], stateAt); if (e !== null) run.push(e); }
          for (let k = i0; k <= j0; k++) {
            if (run.length > 0) subdivide(cx, fen, run[run.length - 1], st[k], stateAt, run, 0);
            run.push(st[k]);
          }
          if (j0 + 1 < n) {
            const e = refineRunEnd(cx, fen, comp[j0], comp[j0 + 1], stateAt);
            if (e !== null) { subdivide(cx, fen, run[run.length - 1], e, stateAt, run, 0); run.push(e); }
          }
          if (run.length >= 2) outRuns.push(run);
          i0 = j0 + 1;
        }
      }
  }
}

function refineRunEnd(cx, fen, wOk, wBad, stateAt) {
  let lo = wOk, hi = wBad, best = null;
  for (let it = 0; it < 48; it++) {
    const mid = fenceFootP(cx.P, fen, unitMid(lo, hi));
    if (mid === null) break;
    const s = stateAt(mid);
    if (s !== null) { lo = mid; best = s; } else hi = mid;
    if (dist(lo, hi) < 1e-15) break;
  }
  return best;
}

// 相邻两点天球弦长超限就沿篱笆取中点插入 (递归)
function subdivide(cx, fen, sa, sb, stateAt, run, depth) {
  if (depth > 10) return;
  const xa = cx.world(sa);
  if (!fin3(xa)) return;
  const xb = cx.world(sb);
  if (!fin3(xb)) return;
  if (ang(xa, xb) <= MAX_SKY) return;
  const wm = fenceFootP(cx.P, fen, unitMid([sa[0], sa[1], sa[2]], [sb[0], sb[1], sb[2]]));
  if (wm === null) return;
  const sm = stateAt(wm);
  if (sm === null) return;
  subdivide(cx, fen, sa, sm, stateAt, run, depth + 1);
  run.push(sm);
  subdivide(cx, fen, sm, sb, stateAt, run, depth + 1);
}

// ------------------------------------------------------------------ 焦散: 扩展状态 x = (w1,w2,w3,T,τ,σ)
// E 5 元, J 行主序 5×6 (可为 null)
function evalE(cx, x, E, J) {
  const w = [x[0], x[1], x[2]];
  const T = x[3], tau = x[4], sig = x[5], rho = cx.Rho;
  const Dc = [0, 0, 0], M0 = new Array(9), g = [0, 0, 0];
  const qE = cx.P.chainT(w, T, Dc, M0, g);
  const nu = cx.nu(tau);
  const p = dot(w, Dc);
  E[0] = dot(w, w) - 1.0;
  E[1] = T * T - qE * qE + cx.P.K;
  E[2] = -(sig * dot(w, nu) + rho * dot(Dc, nu)) + cx.Sl;
  E[3] = sig * sig + 2 * sig * rho * p + rho * rho - 1.0;
  E[4] = cx.jhat(w, T, tau);
  for (let i = 0; i < 5; i++) if (!fin(E[i])) return false;
  if (J == null) return true;
  const dnu = cx.dnu(tau);
  const UE = cx.P.UE;
  for (let j = 0; j < 3; j++) {
    const mtNu = M0[j] * nu[0] + M0[3 + j] * nu[1] + M0[6 + j] * nu[2];
    const mtW = M0[j] * w[0] + M0[3 + j] * w[1] + M0[6 + j] * w[2];
    J[j] = 2 * w[j];
    J[6 + j] = -2 * qE * g[j];
    J[12 + j] = -(sig * nu[j] + rho * mtNu);
    J[18 + j] = 2 * sig * rho * (Dc[j] + mtW);
  }
  J[3] = 0; J[4] = 0; J[5] = 0;
  J[9] = 2 * T; J[10] = 0; J[11] = 0;
  J[15] = rho * dot(UE, nu); J[16] = -(sig * dot(w, dnu) + rho * dot(Dc, dnu)); J[17] = -dot(w, nu);
  J[21] = -2 * sig * rho * dot(w, UE); J[22] = 0; J[23] = 2 * sig + 2 * rho * p;
  const h = 1e-6;
  for (let k = 0; k < 5; k++) {
    const wp = [x[0], x[1], x[2]], wm = [x[0], x[1], x[2]];
    let Tp = T, Tm = T, tp = tau, tm = tau;
    if (k < 3) { wp[k] += h; wm[k] -= h; } else if (k === 3) { Tp += h; Tm -= h; } else { tp += h; tm -= h; }
    J[24 + k] = (cx.jhat(wp, Tp, tp) - cx.jhat(wm, Tm, tm)) / (2 * h);
  }
  J[29] = 0;
  for (let i = 0; i < 30; i++) if (!fin(J[i])) return false;
  return true;
}

// 切向: 解 [J; t_prevᵀ]·t = e₆ 再归一
function tangent6(cx, x, tPrev) {
  const E = new Array(5), J = new Array(30);
  if (!evalE(cx, x, E, J)) return null;
  const A = new Array(36), b = [0, 0, 0, 0, 0, 0], t = [0, 0, 0, 0, 0, 0];
  for (let i = 0; i < 30; i++) A[i] = J[i];
  for (let j = 0; j < 6; j++) A[30 + j] = tPrev[j];
  b[5] = 1.0;
  if (!solveN(A, b, 6, t)) return null;
  let n = 0; for (let j = 0; j < 6; j++) n += t[j] * t[j];
  n = Math.sqrt(n);
  if (!(n > 0) || !fin(n)) return null;
  for (let j = 0; j < 6; j++) t[j] /= n;
  return t;
}

// 校正: Newton 解 [J; tᵀ]·Δ = [−E; 0] (沿切向不动), 到 |E|∞ < 1e-11
function correct6(cx, x0, t) {
  const x = x0.slice();
  const E = new Array(5), J = new Array(30);
  const A = new Array(36), b = [0, 0, 0, 0, 0, 0], d = [0, 0, 0, 0, 0, 0];
  for (let it = 0; it < 12; it++) {
    if (!evalE(cx, x, E, J)) return null;
    let mx = 0; for (let i = 0; i < 5; i++) mx = Math.max(mx, Math.abs(E[i]));
    if (mx < 1e-11) return x;
    for (let i = 0; i < 5; i++) { for (let j = 0; j < 6; j++) A[i * 6 + j] = J[i * 6 + j]; b[i] = -E[i]; }
    for (let j = 0; j < 6; j++) A[30 + j] = t[j];
    b[5] = 0;
    if (!solveN(A, b, 6, d)) return null;
    let dn = 0; for (let j = 0; j < 6; j++) dn += d[j] * d[j];
    dn = Math.sqrt(dn);
    if (dn > 0.2) { const s = 0.2 / dn; for (let j = 0; j < 6; j++) d[j] *= s; }
    for (let j = 0; j < 6; j++) x[j] += d[j];
  }
  if (!evalE(cx, x, E, null)) return null;
  let m2 = 0; for (let i = 0; i < 5; i++) m2 = Math.max(m2, Math.abs(E[i]));
  return m2 < 1e-9 ? x : null;
}

// 有效裕度: min(T, σ, 出射以外各约束篱笆的 margin)。which: 0=T 1=σ 2+=篱笆下标。→ { m, which }
function minMargin6(cx, x) {
  let m = x[3], which = 0;
  if (x[5] < m) { m = x[5]; which = 1; }
  const w = unit([x[0], x[1], x[2]]);
  const P = cx.P;
  for (let k = 0; k < P.Fences.length; k++) {
    const f = P.Fences[k];
    if (f === P.Exit) continue;
    const v = f.Side * P.fval(f, w);
    if (v < m) { m = v; which = 2 + k; }
  }
  return { m, which };
}

function dist6(a, b) {
  let s = 0;
  for (let j = 0; j < 6; j++) {
    let d = a[j] - b[j];
    if (j === 4) d = ieeeRemainder(d, 2 * Math.PI);
    s += d * d;
  }
  return Math.sqrt(s);
}

// 从 x0 沿 t0 方向追一支。→ { pts (含 x0), endWhich (≥0 撞边界, −1 闭合, −2 失败/步数用尽), stop }
function traceDir6(cx, x0, t0) {
  const pts = [x0.slice()];
  let x = x0.slice(), t = t0, dirPrev;
  const out = { pts, endWhich: -2, stop: 'max_steps' };
  dirPrev = cx.world(x);
  if (!fin3(dirPrev)) { out.stop = 'world_fail'; return out; }
  let h = H6_0, armed = false;
  for (let it = 0; it < 30000; it++) {
    if ((it & 63) === 0) ckCancel(cx);
    const tn = tangent6(cx, x, t);
    if (tn === null) { out.stop = 'tangent_fail'; return out; }
    let xc = null, dir = null, ch = 0;
    while (true) {
      const xp = new Array(6);
      for (let j = 0; j < 6; j++) xp[j] = x[j] + h * tn[j];
      xc = correct6(cx, xp, tn);
      if (xc === null || dist6(xc, xp) > 2 * h) {
        h *= 0.5;
        if (h < 1e-8) { out.stop = 'correct_fail'; return out; }
        continue;
      }
      const mm = minMargin6(cx, xc);
      if (mm.m < 0) {
        // 越界: 在步长上二分出边界点 (亮缘的末端)
        let lo = 0, hi = h, xlo = x;
        for (let q = 0; q < 44; q++) {
          const mid = 0.5 * (lo + hi);
          const xm = new Array(6);
          for (let j = 0; j < 6; j++) xm[j] = x[j] + mid * tn[j];
          const xmc = correct6(cx, xm, tn);
          if (xmc === null || minMargin6(cx, xmc).m < 0) hi = mid; else { lo = mid; xlo = xmc; }
        }
        pts.push(xlo);
        out.endWhich = mm.which; out.stop = 'end';
        return out;
      }
      dir = cx.world(xc);
      if (!fin3(dir)) { h *= 0.5; if (h < 1e-8) { out.stop = 'world_fail'; return out; } continue; }
      ch = ang(dir, dirPrev);
      if (ch > MAX_SKY && h > 1e-7) { h *= 0.5; continue; }
      break;
    }
    x = xc; t = tn; dirPrev = dir;
    pts.push(x);
    if (ch < MIN_SKY) h = Math.min(1.5 * h, H6_MAX);
    if (dist6(x, x0) > 4 * H6_MAX) armed = true;
    if (armed && dist6(x, x0) < 1.5 * h) {
      const d0 = cx.world(x0);
      if (ang(dir, d0) > MAX_SKY && h > 1e-7) h *= 0.5;
      else { pts.push(x0.slice()); out.endWhich = -1; out.stop = 'closed'; return out; }
    }
  }
  return out;
}

// 种子: 4 叶上 1° 经纬网格, Ĵ 在相邻有效节点间变号 → 沿边二分。→ { seeds, degenerate }
function causticSeeds(cx) {
  let nNormal = 0, nFlat = 0;
  const d = 1.0;
  const nLat = roundEven(180.0 / d), nLon = roundEven(360.0 / d) + 1;
  const W = new Array(nLat * nLon), ok = new Array(nLat * nLon);
  const cLon = new Array(nLon), sLon = new Array(nLon);
  for (let j = 0; j < nLon; j++) { const lon = (-180.0 + j * d) * DEG; cLon[j] = Math.cos(lon); sLon[j] = Math.sin(lon); }
  for (let i = 0; i < nLat; i++) {
    const lat = (-90.0 + d / 2 + i * d) * DEG, cl = Math.cos(lat), sl = Math.sin(lat);
    for (let j = 0; j < nLon; j++) {
      const w = [cl * cLon[j], cl * sLon[j], sl];
      W[i * nLon + j] = w;
      ok[i * nLon + j] = cx.othersOk(w, null, 1e-9);
    }
  }
  const seeds = [];
  const F = new Float64Array(nLat * nLon);
  for (let a = +1; a >= -1; a -= 2)
    for (let c = +1; c >= -1; c -= 2) {
      ckCancel(cx);
      const aa = a, cc = c;
      const jAt = w => {
        if (!cx.othersOk(w, null, 1e-9)) return NaN;
        const s = cx.sheet(w, aa, cc, false);
        if (s === null) return NaN;
        return cx.jhat(w, s[0], s[1]);
      };
      let fMax = 0;
      let nValid = 0, nBig = 0;
      for (let i = 0; i < nLat * nLon; i++) {
        F[i] = ok[i] ? jAt(W[i]) : NaN;
        if (!fin(F[i])) continue;
        nValid++;
        const af = Math.abs(F[i]);
        if (af > fMax) fMax = af;
        if (af > 1e-6) nBig++;
      }
      if (nValid === 0) continue;                                 // 这一叶上没有有效节点
      // J≡0 的光路 (平行面 / 镜面型且第二转轴是镜面法线): 网格上只剩舍入噪声 —— 按"绝大多数节点"判
      if (nBig <= 0.02 * nValid) { nFlat++; continue; }
      nNormal++;
      const noise = 1e-9 * fMax;
      const edge = (wa, wb, fa, fb) => {
        let lo = wa, hi = wb, flo = fa;
        for (let q = 0; q < 30; q++) {
          const mid = unitMid(lo, hi);
          const fm = jAt(mid);
          if (!fin(fm)) return;                              // 叶在边内失效: 不当种子
          if (flo * fm <= 0) hi = mid; else { lo = mid; flo = fm; }
        }
        const w = unitMid(lo, hi);
        const s = cx.sheet(w, aa, cc, false);
        if (s === null) return;
        seeds.push([w[0], w[1], w[2], s[0], s[1], s[2]]);
      };
      for (let i = 0; i + 1 < nLat; i++)
        for (let j = 0; j < nLon; j++) {
          const f0 = F[i * nLon + j], f1 = F[(i + 1) * nLon + j];
          if (fin(f0) && fin(f1) && f0 * f1 < 0 && Math.max(Math.abs(f0), Math.abs(f1)) > noise)
            edge(W[i * nLon + j], W[(i + 1) * nLon + j], f0, f1);
        }
      for (let i = 0; i < nLat; i++)
        for (let j = 0; j + 1 < nLon; j++) {
          const f0 = F[i * nLon + j], f1 = F[i * nLon + j + 1];
          if (fin(f0) && fin(f1) && f0 * f1 < 0 && Math.max(Math.abs(f0), Math.abs(f1)) > noise)
            edge(W[i * nLon + j], W[i * nLon + j + 1], f0, f1);
        }
    }
  return { seeds, degenerate: nNormal === 0 && nFlat > 0 };
}

function traceCaustic(cx, res, outCurves, outEnds) {
  const cs = causticSeeds(cx);
  if (cs.degenerate) { res.note = '该光路在这个族里 J≡0 (出射方向恒等于入射方向, 或不随第二个转角变): 二自由度退化, 没有焦散'; return; }
  const have = [];
  const covered = s => {
    for (const y of have)
      if (Math.abs(y[5] - s[5]) < 0.05 && Math.abs(y[3] - s[3]) < 0.05 && dist6(y, s) < 0.02) return true;
    return false;
  };
  const axes = [];
  for (let k = 0; k < 6; k++) { const ax = [0, 0, 0, 0, 0, 0]; ax[k] = 1.0; axes.push(ax); }
  for (const s0 of cs.seeds) {
    ckCancel(cx);
    if (covered(s0)) continue;
    // 种子先校正到曲线上 (闭式叶已满足 E1–E4, 只差 Ĵ 的二分残差)
    let t0 = null;
    for (const ax of axes) { t0 = tangent6(cx, s0, ax); if (t0 !== null) break; }
    if (t0 === null) continue;
    const x0 = correct6(cx, s0, t0);
    if (x0 === null || minMargin6(cx, x0).m < 0 || covered(x0)) continue;
    t0 = tangent6(cx, x0, t0);
    if (t0 === null) continue;
    const a1 = traceDir6(cx, x0, t0);
    let curve;
    if (a1.endWhich === -1) curve = a1.pts;
    else {
      const tneg = new Array(6); for (let j = 0; j < 6; j++) tneg[j] = -t0[j];
      const a2 = traceDir6(cx, x0, tneg);
      curve = [];
      for (let q = a1.pts.length - 1; q >= 0; q--) curve.push(a1.pts[q]);
      for (let q = 1; q < a2.pts.length; q++) curve.push(a2.pts[q]);
      if (a1.endWhich >= 0) outEnds.push(a1.pts[a1.pts.length - 1]);
      if (a2.endWhich >= 0) outEnds.push(a2.pts[a2.pts.length - 1]);
      if (a1.stop !== 'end' && a1.stop !== 'closed') res.badArcs.push('caustic:' + a1.stop);
      if (a2.stop !== 'end' && a2.stop !== 'closed') res.badArcs.push('caustic:' + a2.stop);
    }
    if (curve.length >= 2) outCurves.push(curve);
    for (const y of curve) have.push(y);
  }
}

/**
 * 二自由度族在一层球壳上的轮廓: 焦散 (亮缘) 与各条篱笆的像 (= DivShell.Compute2Dof)。
 * fam='column': R = Rz(α)·Ry(fixedDeg)·Rz(τ), τ=θ 自由;  fam='lowitz': R = Rz(α)·Ry(τ)·Rz(fixedDeg), τ=ζ 自由。
 * wantFence(name): 要不要画名为 name 的篱笆 ('entry'/'exit'/'seam1'/'blue1'…), 可为 null。
 * opts.cancelled(): 返回 true 则抛 CancelledError;  opts.verify: 独立 snell 核对 (填 maxVerifyErr/maxJfdErr/worstJInfo)。
 * 返回 { curves:[{kind,name,p,dir,closed}], endP, endDir, cuspP, cuspDir, rho, kappa, baseDist, lampElevDeg, radius,
 *        nPoints, nCaustic, nFence, badArcs, note, maxVerifyErr, maxJfdErr, worstJInfo }  (点 = [x,y,z])
 */
export function computeShell2Dof(faceNormals, fam, fixedDeg, refN, lamp, eye, radius, partial,
  wantCaustic, wantFence, opts = {}) {
  const verify = !!opts.verify;
  const res = {
    curves: [], endP: [], endDir: [], cuspP: [], cuspDir: [],
    rho: 0, kappa: 0, baseDist: 0, lampElevDeg: 0, radius,
    nPoints: 0, nCaustic: 0, nFence: 0,
    badArcs: [], note: '', maxVerifyErr: NaN, maxJfdErr: NaN, worstJInfo: '',
  };
  const bx = lamp[0] - eye[0], by = lamp[1] - eye[1], bz = lamp[2] - eye[2];
  const D = Math.sqrt(bx * bx + by * by + bz * bz);
  res.baseDist = D;
  if (!(D > 1e-9) || !(radius > 0) || faceNormals == null || faceNormals.length < 2 || !(refN > 1.0)) {
    res.note = '参数无效';
    return res;
  }
  let rho = radius / D;
  if (Math.abs(rho - 1.0) < 1e-9) rho = 1.0 + 1e-9;
  res.rho = rho; res.kappa = Math.log10(rho);
  res.lampElevDeg = Math.asin(Math.max(-1.0, Math.min(1.0, bz / D))) / DEG;

  const fx = fixedDeg * DEG;
  const cx = new Ctx2();
  cx.P = new PathK(faceNormals, refN, partial);
  cx.column = !(fam === 'lowitz' || fam === 'Lowitz' || fam === 1);
  cx.Fixed = fx;
  cx.Sl = bz / D; cx.Rho = rho; cx.Eye = eye.slice(); cx.Radius = radius;
  cx.AzL = Math.atan2(by, bx);
  if (typeof opts.cancelled === 'function') cx.cancelled = opts.cancelled;
  if (cx.column) {
    // ν(θ) = (−sinζ0 cosθ, sinζ0 sinθ, cosζ0);  θ 的转轴 = 晶体 c 轴
    cx.C0 = [0.0, 0.0, Math.cos(fx)];
    cx.C1 = [-Math.sin(fx), 0.0, 0.0];
    cx.C2 = [0.0, Math.sin(fx), 0.0];
    cx.Gen = [0.0, 0.0, 1.0];
    if (Math.abs(Math.sin(fx)) < 1e-9) { res.note = 'ζ=0/180°: θ 与 α 退化成同一个转动, 这不是二自由度族'; return res; }
  } else {
    // ν(ζ) = (−sinζ cosθ0, sinζ sinθ0, cosζ);  ζ 的转轴在晶体系是 (sinθ0, cosθ0, 0)
    cx.C0 = [0.0, 0.0, 0.0];
    cx.C1 = [0.0, 0.0, 1.0];
    cx.C2 = [-Math.cos(fx), Math.sin(fx), 0.0];
    cx.Gen = [Math.sin(fx), Math.cos(fx), 0.0];
  }

  const emit = (run, kind, name, closed) => {
    const Pw = [], Dw = [];
    for (const x of run) {
      const xi = cx.world(x);
      if (!fin3(xi)) continue;
      Dw.push(xi);
      Pw.push([eye[0] + radius * xi[0], eye[1] + radius * xi[1], eye[2] + radius * xi[2]]);
    }
    if (Pw.length < 2) return;
    res.curves.push({ kind, name, p: Pw, dir: Dw, closed });
    res.nPoints += Pw.length;
    if (kind === 0) res.nCaustic++; else res.nFence++;
  };

  if (wantFence != null)
    for (const f of cx.P.AllFences) {
      if (!wantFence(f.name)) continue;
      const runs = [];
      fenceCuts(cx, f, runs);
      for (const r of runs) emit(r, fenceKind(f.name), f.name, false);
    }

  const cRuns = [], cEnds = [];
  if (wantCaustic) {
    traceCaustic(cx, res, cRuns, cEnds);
    for (const r of cRuns) {
      const closed = r.length > 2 && dist6(r[0], r[r.length - 1]) < 1e-9;
      emit(r, 0, 'caustic', closed);
      // 尖点: 天球像弦方向反转的孤立点
      const cu = res.curves[res.curves.length - 1];
      if (cu === undefined) throw new RangeError('Index was out of range');   // 同 C# (空表取 [−1])
      if (cu.kind === 0 && cu.dir.length >= 7) {
        const chd = new Array(cu.dir.length - 1);
        for (let i = 0; i + 1 < cu.dir.length; i++) chd[i] = sub(cu.dir[i + 1], cu.dir[i]);
        for (let i = 1; i < chd.length; i++) {
          if (!(dot(chd[i], chd[i - 1]) < 0)) continue;
          const l0 = Math.min(nrm(chd[i]), nrm(chd[i - 1]));
          let wide = 0;
          for (let j = Math.max(0, i - 4); j < Math.min(chd.length, i + 4); j++) wide = Math.max(wide, nrm(chd[j]));
          if (wide > 20 * Math.max(l0, 1e-12)) { res.cuspP.push(cu.p[i]); res.cuspDir.push(cu.dir[i]); }
        }
      }
    }
    for (const x of cEnds) {
      const xi = cx.world(x);
      if (!fin3(xi)) continue;
      res.endDir.push(xi);
      res.endP.push([eye[0] + radius * xi[0], eye[1] + radius * xi[1], eye[2] + radius * xi[2]]);
    }
  }
  if (res.badArcs.length > 0) res.note = '焦散追踪异常: ' + res.badArcs.join(',');

  if (verify) verify2(cx, faceNormals, refN, lamp, eye, cRuns, res);
  return res;
}

// .NET "E2" / "F3" (InvariantCulture): 1.23E-005 / 0.123
function fmtE(x, d) {
  if (Number.isNaN(x)) return 'NaN';
  if (x === Infinity) return 'Infinity';
  if (x === -Infinity) return '-Infinity';
  const s = x.toExponential(d), k = s.indexOf('e');
  const ex = parseInt(s.slice(k + 1), 10);
  return (Object.is(x, -0) ? '-' : '') + s.slice(0, k) + 'E' + (ex < 0 ? '-' : '+') + String(Math.abs(ex)).padStart(3, '0');
}
function fmtF(x, d) {
  if (Number.isNaN(x)) return 'NaN';
  if (x === Infinity) return 'Infinity';
  if (x === -Infinity) return '-Infinity';
  return (Object.is(x, -0) ? '-' : '') + x.toFixed(d);
}

// 独立核对 (逐面 snell, 世界系 R = Rz(α)·R1(τ)): (1) 焦散抽样点出射方向指向眼; (2) 对 (α,τ) 差分的 J_fd 在焦散上 ≈0
function verify2(cx, faceNormals, refN, lamp, eye, cRuns, res) {
  let worst = 0, worstJ = 0;
  const trace = (al, tau, u) => {
    const R = new Array(9); cx.r1(tau, R);
    const ca = Math.cos(al), sa = Math.sin(al);
    const N = new Array(faceNormals.length);
    for (let f = 0; f < faceNormals.length; f++) {
      const nf = unit(faceNormals[f]);
      const m = [
        R[0] * nf[0] + R[1] * nf[1] + R[2] * nf[2],
        R[3] * nf[0] + R[4] * nf[1] + R[5] * nf[2],
        R[6] * nf[0] + R[7] * nf[1] + R[8] * nf[2]];
      N[f] = [ca * m[0] - sa * m[1], sa * m[0] + ca * m[1], m[2]];
    }
    const o = [0, 0, 0];
    return snellTrace(N, refN, u, o) ? o : null;
  };
  for (const run of cRuns) {
    const n = run.length, stride = Math.max(1, Math.trunc(n / 40));
    for (let i = 0; i < n; i += stride) {
      const x = run[i];
      if (minMargin6(cx, x).m < 1e-4) continue;          // 贴着篱笆: 逐面追踪的判据不稳, 不比
      const xi = cx.world(x);
      if (!fin3(xi)) continue;
      const Pw = [eye[0] + cx.Radius * xi[0], eye[1] + cx.Radius * xi[1], eye[2] + cx.Radius * xi[2]];
      const u = unit(sub(Pw, lamp));
      const want = unit(sub(eye, Pw));
      const al = cx.alphaOf(x), tau = x[4];
      const got = trace(al, tau, u);
      if (got === null) { worst = Math.max(worst, 1.0); continue; }
      worst = Math.max(worst, dist(got, want));
      const dd = 1e-5;
      const ap = trace(al + dd, tau, u), am = trace(al - dd, tau, u);
      const tp = trace(al, tau + dd, u), tm = trace(al, tau - dd, u);
      if (ap === null || am === null || tp === null || tm === null) continue;
      const dA = scl(sub(ap, am), 1 / (2 * dd)), dT = scl(sub(tp, tm), 1 / (2 * dd));
      if (x[3] < 0.05) continue;                                 // 出射临界旁差分误差 ∝1/T³: 不比
      const jfd = Math.abs(dot(got, cross(dA, dT)));
      if (jfd > worstJ) {
        worstJ = jfd;
        const mm = minMargin6(cx, x);
        res.worstJInfo = 'Jfd=' + fmtE(jfd, 2) + ' T=' + fmtE(x[3], 2) + ' sig=' + fmtF(x[5], 3)
          + ' minMargin=' + fmtE(mm.m, 2) + '(#' + mm.which + ') |dA|=' + fmtE(nrm(dA), 2)
          + ' |dT|=' + fmtE(nrm(dT), 2) + ' Jhat=' + fmtE(cx.jhat([x[0], x[1], x[2]], x[3], x[4]), 2);
      }
    }
  }
  res.maxVerifyErr = worst; res.maxJfdErr = worstJ;
}

// ================================================================== 主入口 (一自由度)
/**
 * 一条光路在一层球壳上的晕线 (= DivShell.Compute)。
 * faceNormals: 晶体系面法向 (入射面, 各反射面…, 出射面)。取向族 R = Rz(α)·Ry(ζ)·Rz(θ), α ∈ [alphaMin, alphaMax] (度)。
 * lamp/eye: 世界坐标 (米, Z-up)。radius: 球壳半径 (米, 以眼为心)。partial: 中间反射允许部分反射。
 * opts.cancelled(): 返回 true 则抛 CancelledError;  opts.verify: 独立 snell 核对 (填 maxVerifyErr)。
 * 返回 { ends:[{w,sigma,fence,fenceLabel,missed,degenerate,startFail,used,p,dir,alphaDeg,inAlphaRange}],
 *        arcs:[{w,sigma,dir,p,alphaDeg,from,to,closed,viaLoop,stop,foldIdx,foldP,foldDir}],
 *        polylines, polylineDirs, cuspP, rho, kappa, baseDist, lampElevDeg, radius,
 *        ok, nMissed, nUnused, nOverused, nPoints, badArcs, note, maxVerifyErr }
 */
export function computeShell(faceNormals, zetaDeg, thetaDeg, refN, lamp, eye, radius, partial,
  alphaMinDeg, alphaMaxDeg, opts = {}) {
  const verify = !!opts.verify;
  const res = {
    ends: [], arcs: [], polylines: [], polylineDirs: [], cuspP: [],
    rho: 0, kappa: 0, baseDist: 0, lampElevDeg: 0, radius,
    ok: false, nMissed: 0, nUnused: 0, nOverused: 0, nPoints: 0,
    badArcs: [], note: '', maxVerifyErr: NaN,
  };
  const bx = lamp[0] - eye[0], by = lamp[1] - eye[1], bz = lamp[2] - eye[2];
  const D = Math.sqrt(bx * bx + by * by + bz * bz);
  res.baseDist = D;
  if (!(D > 1e-9) || !(radius > 0) || faceNormals == null || faceNormals.length < 2 || !(refN > 1.0)) {
    res.note = '参数无效';
    return res;
  }
  let rho = radius / D;
  if (Math.abs(rho - 1.0) < 1e-9) rho = 1.0 + 1e-9;      // ρ=1 是生灭层 (σ 可到 0), 稍偏一点避开退化
  res.rho = rho; res.kappa = Math.log10(rho);
  res.lampElevDeg = Math.asin(Math.max(-1.0, Math.min(1.0, bz / D))) / DEG;

  const cx = new Ctx();
  cx.P = new PathK(faceNormals, refN, partial);
  cx.Sl = bz / D; cx.Rho = rho;
  cx.Eye = eye.slice(); cx.Radius = radius;
  cx.AzL = Math.atan2(by, bx);
  if (typeof opts.cancelled === 'function') cx.cancelled = opts.cancelled;
  const z = zetaDeg * DEG, t = thetaDeg * DEG;
  const cz = Math.cos(z), sz = Math.sin(z), cth = Math.cos(t), sth = Math.sin(t);
  // R0 = Ry(ζ)·Rz(θ)
  const R0 = cx.R0;
  R0[0] = cz * cth; R0[1] = -cz * sth; R0[2] = sz;
  R0[3] = sth; R0[4] = cth; R0[5] = 0;
  R0[6] = -sz * cth; R0[7] = sz * sth; R0[8] = cz;
  cx.Nu = [R0[6], R0[7], R0[8]];

  const seeds = findEnds(cx);
  const raw = traceArcs(cx, seeds);
  for (const a of findLoops(cx, seeds, raw)) raw.push(a);

  // ---------------- 输出: 世界坐标
  for (const s of seeds) {
    if (s.isLoop) continue;
    const wr = cx.world(s.w, s.sigma);
    if (wr !== null) {
      const xi = wr.xi;
      s.dir = xi; s.alphaDeg = wr.al;
      s.p = [eye[0] + radius * xi[0], eye[1] + radius * xi[1], eye[2] + radius * xi[2]];
      s.inAlphaRange = inRange(wr.al, alphaMinDeg, alphaMaxDeg);
    }
    if (s.fenceLabel == null) s.fenceLabel = s.fence;
    res.ends.push(s);
  }
  const endIndex = new Map();
  for (let i = 0; i < res.ends.length; i++) endIndex.set(res.ends[i], i);

  for (const a of raw) {
    const o = {
      w: [], sigma: [], dir: [], p: [], alphaDeg: [], from: -1, to: -1,
      closed: a.Closed, viaLoop: a.ViaLoop, stop: a.Stop, foldIdx: [], foldP: [], foldDir: [],
    };
    o.from = (a.From >= 0 && a.From < seeds.length && !seeds[a.From].isLoop) ? endIndex.get(seeds[a.From]) : -1;
    o.to = (a.To >= 0 && a.To < seeds.length && !seeds[a.To].isLoop) ? endIndex.get(seeds[a.To]) : -1;
    for (let i = 0; i < a.W.length; i++) {
      const wr = cx.world(a.W[i], a.Sig[i]);
      if (wr === null) continue;
      const xi = wr.xi;
      o.w.push(a.W[i]); o.sigma.push(a.Sig[i]); o.dir.push(xi); o.alphaDeg.push(wr.al);
      o.p.push([eye[0] + radius * xi[0], eye[1] + radius * xi[1], eye[2] + radius * xi[2]]);
    }
    for (let q = 0; q < a.FoldW.length; q++) {
      const wr = cx.world(a.FoldW[q], a.FoldSig[q]);
      if (wr === null) continue;
      if (!inRange(wr.al, alphaMinDeg, alphaMaxDeg)) continue;
      const xi = wr.xi;
      o.foldIdx.push(a.FoldIdx[q]);
      o.foldDir.push(xi);
      o.foldP.push([eye[0] + radius * xi[0], eye[1] + radius * xi[1], eye[2] + radius * xi[2]]);
    }
    res.arcs.push(o);
    res.nPoints += o.p.length;
    clipByAlpha(o, alphaMinDeg, alphaMaxDeg, eye, radius, res);
    findCusps(o, res);
  }

  // ---------------- 记账
  for (const s of res.ends) {
    if (s.missed) res.nMissed++;
    if (s.used === 0) res.nUnused++;
    if (s.used > 1) res.nOverused++;
  }
  for (const a of res.arcs)
    for (const part of a.stop.split('|'))
      if (part !== 'end' && part !== 'closed') { res.badArcs.push(a.stop); break; }
  res.ok = res.nMissed === 0 && res.nUnused === 0 && res.nOverused === 0 && res.badArcs.length === 0 && (res.ends.length % 2 === 0);
  if (!res.ok)
    res.note = '记账异常: 漏端' + res.nMissed + ' 未用' + res.nUnused + ' 重用' + res.nOverused
             + (res.badArcs.length > 0 ? ' 坏弧[' + res.badArcs.join(',') + ']' : '')
             + (res.ends.length % 2 !== 0 ? ' 末端数为奇' : '');

  if (verify) res.maxVerifyErr = verify1(cx, faceNormals, refN, lamp, eye, res);
  // 去掉内部字段 (C# internal: Fences / Fen / IsLoop)
  for (const s of res.ends) { delete s.fences; delete s.fen; delete s.isLoop; }
  return res;
}

function inRange(alDeg, aMin, aMax) {
  const span = aMax - aMin;
  if (!(span < 360.0 - 1e-9)) return true;
  if (span < 0) return false;
  const rel = ((alDeg - aMin) % 360.0 + 360.0) % 360.0;
  return rel <= span + 1e-9;
}

// 按 α 区间把弧裁成折线 (区间边界处在相邻两点间线性插值)
function clipByAlpha(a, aMin, aMax, eye, radius, res) {
  const n = a.p.length;
  if (n < 2) return;
  let curP = [], curD = [];
  const flush = () => {
    if (curP.length >= 2) { res.polylines.push(curP); res.polylineDirs.push(curD); }
    curP = []; curD = [];
  };
  const full = !(aMax - aMin < 360.0 - 1e-9);
  for (let i = 0; i < n; i++) {
    const inI = full || inRange(a.alphaDeg[i], aMin, aMax);
    if (i > 0 && !full) {
      const inPrev = inRange(a.alphaDeg[i - 1], aMin, aMax);
      if (inPrev !== inI) {
        // 在 (i−1, i) 之间找 α 穿越区间端点的位置: 按 α 差线性内插
        const a0 = a.alphaDeg[i - 1], a1 = a.alphaDeg[i];
        const da = ((a1 - a0 + 540.0) % 360.0) - 180.0;
        const bound = nearestBound(a0, da, aMin, aMax);
        let fr = Math.abs(da) > 1e-12 ? (((bound - a0 + 540.0) % 360.0) - 180.0) / da : 0.5;
        if (fr < 0) fr = 0; else if (fr > 1) fr = 1;
        const D0 = a.dir[i - 1], D1 = a.dir[i];
        const dmid = unit([
          D0[0] + fr * (D1[0] - D0[0]),
          D0[1] + fr * (D1[1] - D0[1]),
          D0[2] + fr * (D1[2] - D0[2])]);
        const pmid = [eye[0] + radius * dmid[0], eye[1] + radius * dmid[1], eye[2] + radius * dmid[2]];
        if (inPrev) { curP.push(pmid); curD.push(dmid); flush(); }
        else { curP.push(pmid); curD.push(dmid); }
      }
    }
    if (inI) { curP.push(a.p[i]); curD.push(a.dir[i]); }
    else if (curP.length > 0) flush();
  }
  flush();
}

function nearestBound(a0, da, aMin, aMax) {
  // 沿 α 的行进方向, 先碰到哪个端点
  const toMin = ((aMin - a0 + 540.0) % 360.0) - 180.0, toMax = ((aMax - a0 + 540.0) % 360.0) - 180.0;
  const minAhead = toMin * da >= 0 && Math.abs(toMin) <= Math.abs(da) + 1e-9;
  const maxAhead = toMax * da >= 0 && Math.abs(toMax) <= Math.abs(da) + 1e-9;
  if (minAhead && maxAhead) return Math.abs(toMin) <= Math.abs(toMax) ? aMin : aMax;
  if (minAhead) return aMin;
  if (maxAhead) return aMax;
  return Math.abs(toMin) <= Math.abs(toMax) ? aMin : aMax;
}

// 尖点: 弧上天球像弦方向反转的孤立点 (两侧 |弦| 需显著回升)
function findCusps(a, res) {
  const n = a.dir.length;
  if (n < 7) return;
  const ch = new Array(n - 1);
  for (let i = 0; i + 1 < n; i++) ch[i] = sub(a.dir[i + 1], a.dir[i]);
  for (let i = 1; i < ch.length; i++) {
    if (!(dot(ch[i], ch[i - 1]) < 0)) continue;
    const l0 = Math.min(nrm(ch[i]), nrm(ch[i - 1]));
    let wide = 0;
    for (let j = Math.max(0, i - 4); j < Math.min(ch.length, i + 4); j++) wide = Math.max(wide, nrm(ch[j]));
    if (wide > 20 * Math.max(l0, 1e-12)) res.cuspP.push(a.p[i]);
  }
}

// 抽样核对: 每条弧约 40 点, 逐面 snell (世界系 R=Rz(α)·R0) 复算出射方向
function verify1(cx, faceNormals, refN, lamp, eye, res) {
  let worst = 0;
  for (const a of res.arcs) {
    const n = a.p.length, stride = Math.max(1, Math.trunc(n / 40));
    for (let i = 0; i < n; i += stride) {
      if (cx.P.Fences.length > 0) {
        let nearFence = false;
        for (const f of cx.P.Fences) if (cx.margin(f, a.w[i]) < 1e-6) { nearFence = true; break; }
        if (nearFence) continue;                  // 篱笆上逐面追踪的判据不稳, 不比
      }
      const al = a.alphaDeg[i] * DEG, ca = Math.cos(al), sa = Math.sin(al);
      const N = new Array(faceNormals.length);
      for (let f = 0; f < faceNormals.length; f++) {
        const m = cx.applyR0(unit(faceNormals[f]));
        N[f] = [ca * m[0] - sa * m[1], sa * m[0] + ca * m[1], m[2]];
      }
      const u = unit(sub(a.p[i], lamp));
      const want = unit(sub(eye, a.p[i]));
      const got = [0, 0, 0];
      if (!snellTrace(N, refN, u, got)) { worst = Math.max(worst, 1.0); continue; }
      worst = Math.max(worst, dist(got, want));
    }
  }
  return worst;
}
