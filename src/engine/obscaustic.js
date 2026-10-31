// 观测者视角发散焦散引擎 — DivergentObserverCaustic.cs (static class ObsCaustic) 的 JS 移植.
// 与 C# 逐句对应, 运算次序保持不变 (tools/cmp-obs.mjs 与 C# 参考 tools/ref-obs 逐点对账).
//
// 理论 (halotheories/divergent/observer caustic 1dof):
//   取向族 R(s) = Rot_W(s)·R0; ξ = R^T·û 晶体系入射传播; F(ξ) = Snell 链; b(s) = R(s)^T·𝔻, 𝔻 = L − E.
//   折叠焦散 χ = det(ξ, F, V) = 0, V = DF⁻¹(Ŵ×F) − Ŵ×ξ (伴随式) —— 不含 s、不含 𝔻, 灯只管配对与可见性.
//   二自由度: χ_内(ξ)=0 曲线 → 每点解内层角 → 闭式外层角.
//
// 移植约定 (与 C# 的差别只在表示层):
//   · 3×3 矩阵 = 嵌套数组 M[i][j] (C# double[,] 的 M[i,j]);
//   · C# 的 `out bool ok` → 传入对象 o, 函数里置 o.ok (Field 签名 = (xi, o) => value);
//   · List<T> → 数组; Tuple → 小对象 (RetroPoints: {xi, elev}; VerifySegs: {count, worst});
//   · (float) 输出 (Seg.pts3 / retroCircles) → Float32Array (存入即 round-to-nearest-even, 等同 Math.fround);
//   · Lang.T / Lang.English → 本模块的 setLang('zh'|'en') + 内部 T(zh, en);
//   · AppendLine 的换行用 "\n" (C# 在 Windows 上是 "\r\n");
//   · Dictionary<string,Seg> 的枚举次序 (Fold2Dof 段的输出顺序) 用 NetDict 精确仿真 (槽位序, 删除槽 LIFO 复用).

// ================================================== 语言开关 (Lang.T / Lang.English)
let _english = false;
/** 'zh' (默认) | 'en' — 只影响 Diagnose / SelfTest 等人话文本. */
export function setLang(code) { _english = code === 'en'; }
export function isEnglish() { return _english; }
const T = (zh, en) => (_english ? en : zh);

// ================================================== 上下文与输出

/** 一自由度场景 (二自由度时内层角固定后也用它). */
export class Ctx {
  constructor(o = {}) {
    this.faces = o.faces ?? null;        // 紧凑面索引 (ParsePath 输出)
    this.baseN = o.baseN ?? null;        // 基准法向表 (未施加 R0), 按紧凑索引
    this.n = o.n ?? 1.31;                // 折射率
    this.r0 = o.r0 ?? null;              // 3×3
    this.w = o.w ?? null;                // 实验室旋转轴 (单位)
    this.what = o.what ?? null;          // R0^T·W (晶体系轴, 恒定)
    this.dvec = o.dvec ?? null;          // 𝔻 = L − E
    this.eyePos = o.eyePos ?? null;      // E (世界)
    this.maxDist = o.maxDist ?? 80.0;    // 距眼截断
  }
  finish() {
    this.w = nz(this.w);
    this.what = nz(mtv(this.r0, this.w));   // R0^T·W
  }
}

/** 一段可见焦散 (天球方向 + 空间脊线). */
export class Seg {
  constructor(kind = 'K1', branch = 0) {
    this.kind = kind;          // K1 折叠 / K5:entry 等缝像 / K1-2dof
    this.branch = branch;      // 配对分支 0/1
    this.truncated = false;    // 顶到 MaxDist (无界分支)
    this.xi = [];
    this.d = [];               // 天球方向
    this.p = [];               // 空间脊线
    this.alpha = [];
    this.r = [];
    this.t = [];
    this.inner = [];           // 二自由度: 逐点内层角 b
  }
  /** C# Pts3(): (float) 打包 — Float32Array 存入即 fround. */
  pts3() {
    const P = this.p;
    const a = new Float32Array(P.length * 3);
    for (let i = 0; i < P.length; i++) {
      a[3 * i] = P[i][0]; a[3 * i + 1] = P[i][1]; a[3 * i + 2] = P[i][2];
    }
    return a;
  }
}

// ================================================== 向量工具
const cr = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const dt = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const ad = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
const sb = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const sc = (a, k) => [a[0] * k, a[1] * k, a[2] * k];
const nrm = a => Math.sqrt(dt(a, a));
function nz(a) { const m = nrm(a); return m < 1e-300 ? a : sc(a, 1.0 / m); }
const det3 = (a, b, c) => dt(cr(a, b), c);
const mv = (M, v) => [
  M[0][0] * v[0] + M[0][1] * v[1] + M[0][2] * v[2],
  M[1][0] * v[0] + M[1][1] * v[1] + M[1][2] * v[2],
  M[2][0] * v[0] + M[2][1] * v[1] + M[2][2] * v[2]];
const mtv = (M, v) => [
  M[0][0] * v[0] + M[1][0] * v[1] + M[2][0] * v[2],
  M[0][1] * v[0] + M[1][1] * v[1] + M[2][1] * v[2],
  M[0][2] * v[0] + M[1][2] * v[1] + M[2][2] * v[2]];
const isNaN_ = Number.isNaN;

/** 绕单位轴 w 转 a 作用于 v (Rodrigues). */
export function rotV(w, a, v) {
  const c = Math.cos(a), s = Math.sin(a);
  const wxv = cr(w, v);
  const wd = dt(w, v);
  return [
    v[0] * c + wxv[0] * s + w[0] * wd * (1 - c),
    v[1] * c + wxv[1] * s + w[1] * wd * (1 - c),
    v[2] * c + wxv[2] * s + w[2] * wd * (1 - c)];
}

export function rotZ(a) {
  const c = Math.cos(a), s = Math.sin(a);
  return [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]];
}
export function rotY(a) {
  const c = Math.cos(a), s = Math.sin(a);
  return [[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]];
}
export function mm(A, B) {
  const C = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let i = 0; i < 3; i++)
    for (let j = 0; j < 3; j++)
      C[i][j] = A[i][0] * B[0][j] + A[i][1] * B[1][j] + A[i][2] * B[2][j];
  return C;
}

/** 右手切基: e1×e2 = v (v 单位). 返回 [e1, e2]. */
function perpBasis(v) {
  let k = 0;
  let am = Math.abs(v[0]);
  if (Math.abs(v[1]) < am) { k = 1; am = Math.abs(v[1]); }
  if (Math.abs(v[2]) < am) k = 2;
  const a = [0, 0, 0]; a[k] = 1.0;
  const e1 = nz(cr(a, v));
  const e2 = cr(v, e1);
  return [e1, e2];
}

// ================================================== Snell 链 (解析雅可比)

export class ChainOut {
  constructor() {
    this.f = null;          // 晶体系出射传播
    this.df = null;         // 3×3 雅可比 (wantJac 时)
    this.valid = false;
    this.entry = 0; this.exitApp = 0; this.exitTir = 0;
    this.mirApp = null;     // 每面镜的逼近量 (Float64Array, 长 = 面数 − 2)
  }
}

// 雅可比工作区 (行主序 3×3, 模块级复用: 引擎同步执行, 每个 Worker 一份). 只换存储, 不换算术:
// 系数矩阵先算好 (C# 在三重循环里逐项重算, 值相同), 求和仍是 s = 0; s += …(p = 0,1,2).
let _dfA = new Float64Array(9), _dfB = new Float64Array(9);
const _co = new Float64Array(9);
function rank1Apply(diag, k, m) {          // _dfA ← [(i==p ? diag : 0) − k·m_i·m_p]·_dfA
  for (let i = 0; i < 3; i++)
    for (let p = 0; p < 3; p++) _co[3 * i + p] = (i === p ? diag : 0.0) - k * m[i] * m[p];
  const A = _dfA, Tm = _dfB;
  for (let i = 0; i < 3; i++)
    for (let j = 0; j < 3; j++) {
      let s = 0;
      s += _co[3 * i] * A[j];
      s += _co[3 * i + 1] * A[3 + j];
      s += _co[3 * i + 2] * A[6 + j];
      Tm[3 * i + j] = s;
    }
  _dfA = Tm; _dfB = A;
}

/** 晶体系 Snell 链 ξ → F, 附解析 DF (界面事件 = 沿法向的秩一修正之积). */
export function chain(xi, c, wantJac) {
  const o = new ChainOut();
  const n = c.n, mu = 1.0 / n;
  const faces = c.faces, baseN = c.baseN;
  const last = faces.length - 1;
  const m0 = baseN[faces[0]];
  const c1 = -dt(xi, m0);
  o.entry = c1;
  o.valid = c1 > 0;
  const c1c = Math.abs(c1) < 1e-15 ? 1e-15 : c1;
  let disc0 = 1.0 - mu * mu * (1.0 - c1c * c1c);
  if (disc0 <= 0) { o.valid = false; disc0 = 1e-300; }
  const c2 = Math.sqrt(disc0);

  // zeta 用三个标量 (C# 是数组; (mu·c1c − c2) 每分量重算, 值相同)
  const kz = mu * c1c - c2;
  let z0 = mu * xi[0] + kz * m0[0], z1 = mu * xi[1] + kz * m0[1], z2 = mu * xi[2] + kz * m0[2];

  if (wantJac) {
    const k = mu * (1.0 - mu * c1c / c2);
    const A = _dfA;
    for (let i = 0; i < 3; i++)
      for (let j = 0; j < 3; j++)
        A[3 * i + j] = (i === j ? mu : 0.0) - k * m0[i] * m0[j];
  }

  o.mirApp = new Float64Array(Math.max(0, last - 1));
  for (let q = 1; q < last; q++) {
    const m = baseN[faces[q]];
    const a = z0 * m[0] + z1 * m[1] + z2 * m[2];
    o.mirApp[q - 1] = a;
    if (a <= 0) o.valid = false;
    z0 -= 2.0 * a * m[0]; z1 -= 2.0 * a * m[1]; z2 -= 2.0 * a * m[2];
    if (wantJac) rank1Apply(1.0, 2.0, m);
  }

  const me = baseN[faces[last]];
  const ce = z0 * me[0] + z1 * me[1] + z2 * me[2];
  o.exitApp = ce;
  const discE = 1.0 - n * n * (1.0 - ce * ce);
  o.exitTir = discE;
  if (ce <= 0 || discE <= 0) o.valid = false;
  const cec = Math.abs(ce) < 1e-15 ? 1e-15 : ce;
  const c2e = Math.sqrt(discE > 0 ? discE : 1e-300);
  const kf = n * cec - c2e;
  o.f = [n * z0 - kf * me[0], n * z1 - kf * me[1], n * z2 - kf * me[2]];
  if (wantJac) {
    rank1Apply(n, n * (1.0 - n * cec / c2e), me);
    const D = _dfA;
    o.df = [[D[0], D[1], D[2]], [D[3], D[4], D[5]], [D[6], D[7], D[8]]];
  }
  return o;
}

// ================================================== χ 场 (伴随式连续延拓)

/** χ' = det(M)·χ, 用伴随式避开 DF⁻¹ 的奇异; 右手切基下光滑、符号稳定. out.ok ← 链有效. */
export function chi(xi, c, what, out = null) {
  const ch = chain(xi, c, true);
  if (out) out.ok = ch.valid;
  if (ch.df === null) return NaN;
  const [e1, e2] = perpBasis(xi);
  const [f1, f2] = perpBasis(nz(ch.f));
  const DFe1 = mv(ch.df, e1);
  const DFe2 = mv(ch.df, e2);
  const M11 = dt(f1, DFe1), M12 = dt(f1, DFe2);
  const M21 = dt(f2, DFe1), M22 = dt(f2, DFe2);
  const wxF = cr(what, ch.f);
  const g1 = dt(f1, wxF), g2 = dt(f2, wxF);
  const detM = M11 * M22 - M12 * M21;
  const x1 = M22 * g1 - M12 * g2;
  const x2 = -M21 * g1 + M11 * g2;
  const wxXi = cr(what, xi);
  const Vs = [0, 0, 0];
  for (let i = 0; i < 3; i++) Vs[i] = x1 * e1[i] + x2 * e2[i] - detM * wxXi[i];
  return det3(xi, ch.f, Vs);
}

/** Field = (xi, o) => value, 并置 o.ok. */
export function chiFieldOf(c, what) {
  return (xi, o) => chi(xi, c, what, o);
}

/** 缝场: which ∈ entry / exitTir / exitApp / mir{i}. 零线 = 掠射/TIR 缝. */
export function seamFieldOf(c, which) {
  const isMir = which.startsWith('mir');
  const k = isMir ? parseIntNet(which.substring(3)) : -1;
  return (xi, o) => {
    const ch = chain(xi, c, false);
    o.ok = true;
    if (which === 'entry') return ch.entry;
    if (which === 'exitTir') return ch.exitTir;
    if (which === 'exitApp') return ch.exitApp;
    if (isMir) {
      if (ch.mirApp === null || k >= ch.mirApp.length) { o.ok = false; return NaN; }
      return ch.mirApp[k];
    }
    o.ok = false; return NaN;
  };
}
// int.Parse: 非整数串在 C# 抛异常
function parseIntNet(s) {
  if (!/^\s*[+-]?\d+\s*$/.test(s)) throw new Error(`FormatException: '${s}'`);
  return parseInt(s, 10);
}

/** 缝曲线上的物理有效性 (其余条件 ≥ −tol, 允许相邻缝同时为零). */
export function seamValid(xi, c, which, tol = 1e-7) {
  const ch = chain(xi, c, false);
  if (which !== 'entry' && ch.entry <= -tol) return false;
  if (which !== 'exitApp' && ch.exitApp <= -tol) return false;
  if (which !== 'exitTir' && ch.exitTir <= -tol) return false;
  if (ch.mirApp !== null)
    for (let i = 0; i < ch.mirApp.length; i++)
      if (which !== 'mir' + i && ch.mirApp[i] <= -tol) return false;
  return true;
}

// ================================================== 球面零线追踪 (1-D 圆步进 + 二分)

function bisectArc(f, p1, p2, iters = 45) {
  const o = { ok: false };
  let v1 = f(p1, o);
  for (let i = 0; i < iters; i++) {
    const pm = nz(ad(p1, p2));
    const vm = f(pm, o);
    if (!o.ok || isNaN_(vm)) break;
    if (Math.sign(vm) === Math.sign(v1)) { p1 = pm; v1 = vm; }
    else p2 = pm;
  }
  return nz(ad(p1, p2));
}

function findSeeds(f, nScan, nSamp) {
  const seeds = [];
  const o = { ok: false };
  for (let a = 0; a < nScan; a++) {
    const phi = 2 * Math.PI * a / nScan;
    const cp = Math.cos(phi), sp = Math.sin(phi);
    let prev = null, pv = 0, pok = false;
    for (let i = 0; i < nSamp; i++) {
      const th = 0.015 + (Math.PI - 0.03) * i / (nSamp - 1.0);
      const p = [Math.sin(th) * cp, Math.sin(th) * sp, Math.cos(th)];
      const v = f(p, o); const ok = o.ok;
      if (pok && ok && !isNaN_(v) && !isNaN_(pv) && Math.sign(v) !== Math.sign(pv))
        seeds.push(bisectArc(f, prev, p));
      prev = p; pv = v; pok = ok && !isNaN_(v);
    }
  }
  const nPar = Math.max(8, Math.trunc(nScan / 2));
  for (let a = 0; a < nPar; a++) {
    const tv = 0.08 + (Math.PI - 0.16) * a / (nPar - 1.0);
    const st = Math.sin(tv), ct = Math.cos(tv);
    let first = null, prev = null, fv = 0, pv = 0, pok = false, fok = false;
    for (let i = 0; i < nSamp; i++) {
      const ph = 2 * Math.PI * i / nSamp;
      const p = [st * Math.cos(ph), st * Math.sin(ph), ct];
      const v = f(p, o); const ok = o.ok;
      if (i === 0) { first = p; fv = v; fok = ok && !isNaN_(v); }
      else if (pok && ok && !isNaN_(v) && Math.sign(v) !== Math.sign(pv))
        seeds.push(bisectArc(f, prev, p));
      prev = p; pv = v; pok = ok && !isNaN_(v);
    }
    if (pok && fok && Math.sign(pv) !== Math.sign(fv)) seeds.push(bisectArc(f, prev, first));
  }
  return seeds;
}

function marchDirs(f, p, h) {
  const [e1, e2] = perpBasis(p);
  const res = [];
  const NP = 72;
  const cand = new Array(NP), vals = new Array(NP), oks = new Array(NP);
  const o = { ok: false };
  for (let i = 0; i < NP; i++) {
    const ph = 2 * Math.PI * i / NP;
    const q = [0, 0, 0];
    for (let k = 0; k < 3; k++)
      q[k] = Math.cos(h) * p[k] + Math.sin(h) * (Math.cos(ph) * e1[k] + Math.sin(ph) * e2[k]);
    cand[i] = nz(q);
    vals[i] = f(cand[i], o); oks[i] = o.ok && !isNaN_(vals[i]);
  }
  for (let i = 0; i < NP; i++) {
    const j = (i + 1) % NP;
    if (oks[i] && oks[j] && Math.sign(vals[i]) !== Math.sign(vals[j])) {
      const q = bisectArc(f, cand[i], cand[j]);
      const tau = sb(q, sc(p, dt(q, p)));
      if (nrm(tau) > 1e-12) res.push(nz(tau));
    }
  }
  return res;
}

/** 追踪 field=0 的全部曲线 (纯 1-D 圆步进 + 二分, 无网格). 返回 曲线[] (每条 = 单位向量[]). */
export function traceZeros(f, h = 0.008, maxSteps = 6000, maxCurves = 40, nScan = 96, nSamp = 1440,
  cancelled = null) {
  const seeds = findSeeds(f, nScan, nSamp);
  const curves = [];
  const visited = [];
  const seen = (p, tol) => {
    const ct = Math.cos(tol);
    for (const q of visited) if (dt(q, p) > ct) return true;
    return false;
  };
  const o = { ok: false };
  const NW = 27;
  const cand = new Array(NW), vals = new Array(NW), oks = new Array(NW);

  for (const seed of seeds) {
    if (curves.length >= maxCurves) break;
    if (cancelled && cancelled()) break;
    if (seen(seed, 1.6 * h)) continue;
    const dirs = marchDirs(f, seed, h);
    if (dirs.length === 0) continue;
    const halves = [];
    for (let d = 0; d < Math.min(2, dirs.length); d++) {
      const pts = [seed];
      let p = seed, tau = dirs[d], hh = h;
      for (let step = 0; step < maxSteps; step++) {
        if (cancelled && cancelled()) break;
        const enu = cr(p, tau);
        let best = NaN, bi = -1;
        for (let i = 0; i < NW; i++) {
          const w = -1.25 + 2.5 * i / (NW - 1.0);
          const q = [0, 0, 0];
          for (let k = 0; k < 3; k++)
            q[k] = Math.cos(hh) * p[k] + Math.sin(hh) * (Math.cos(w) * tau[k] + Math.sin(w) * enu[k]);
          cand[i] = nz(q);
          vals[i] = f(cand[i], o); oks[i] = o.ok && !isNaN_(vals[i]);
        }
        for (let i = 0; i < NW - 1; i++) {
          if (!oks[i] || !oks[i + 1] || Math.sign(vals[i]) === Math.sign(vals[i + 1])) continue;
          const wm = (-1.25 + 2.5 * i / (NW - 1.0)) + (-1.25 + 2.5 * (i + 1) / (NW - 1.0));
          if (bi < 0 || Math.abs(wm) < Math.abs(best)) { best = wm; bi = i; }
        }
        if (bi < 0) {
          if (hh > 0.0015) { hh *= 0.5; continue; }
          break;
        }
        const qq = bisectArc(f, cand[bi], cand[bi + 1]);
        const turn = Math.abs(best) / 2;
        if (turn > 0.5 && hh > 0.0015) { hh *= 0.5; continue; }
        if (turn < 0.15) hh = Math.min(hh * 1.25, h * 2.0);
        const nt = sb(qq, sc(p, dt(qq, p)));
        if (nrm(nt) < 1e-14) break;
        tau = nz(nt); p = qq; pts.push(p);
        if (pts.length > 12 && dt(p, seed) > Math.cos(2.0 * hh) && dt(tau, dirs[d]) > 0) break;
      }
      halves.push(pts);
    }
    let curve;
    if (halves.length === 2 && halves[1].length > 1) {
      curve = [];
      for (let i = halves[1].length - 1; i >= 0; i--) curve.push(halves[1][i]);
      for (let i = 1; i < halves[0].length; i++) curve.push(halves[0][i]);
    } else curve = halves[0];
    if (curve.length < 6) continue;
    curves.push(curve);
    for (const q of curve) visited.push(q);
  }
  return curves;
}

/** 把追踪折线在 field=0 上加密到步长 ≤ maxStep (弧度), 横向二分投回零集. */
export function densify(curve, f, maxStep) {
  if (curve === null || curve === undefined || curve.length < 2) return curve;
  const outp = [curve[0]];
  for (let i = 0; i + 1 < curve.length; i++) {
    const p0 = curve[i], p1 = curve[i + 1];
    const dot = Math.max(-1, Math.min(1, dt(p0, p1)));
    const ang = Math.acos(dot);
    if (ang < 1e-12) continue;
    const k = Math.max(1, Math.ceil(ang / maxStep));
    let nu = cr(p0, p1);
    const nn = nrm(nu);
    nu = nn < 1e-15 ? [0, 0, 1] : sc(nu, 1.0 / nn);
    const sa = Math.sin(ang);
    for (let j = 1; j <= k; j++) {
      const tt = j / k;
      const s0 = Math.sin((1 - tt) * ang) / sa, s1 = Math.sin(tt * ang) / sa;
      const q = nz([s0 * p0[0] + s1 * p1[0], s0 * p0[1] + s1 * p1[1], s0 * p0[2] + s1 * p1[2]]);
      outp.push(projectTransverse(q, nu, f, 3e-3));
    }
  }
  return outp;
}

function projectTransverse(p, nu, f, delta) {
  const o = { ok: false };
  const f0 = f(p, o);
  if (!o.ok || isNaN_(f0)) return p;
  for (const dd of [delta, 4 * delta]) {
    const pa = nz(sb(p, sc(nu, dd)));
    const pb = nz(ad(p, sc(nu, dd)));
    const fa = f(pa, o); const oka = o.ok;
    const fb = f(pb, o); const okb = o.ok;
    let lo = null, hi = null;
    if (oka && !isNaN_(fa) && Math.sign(fa) !== Math.sign(f0)) { lo = pa; hi = p; }
    else if (okb && !isNaN_(fb) && Math.sign(fb) !== Math.sign(f0)) { lo = p; hi = pb; }
    if (lo !== null) return bisectArc(f, lo, hi, 42);
  }
  return p;
}

// ================================================== α 配对 (锥上闭式) 与可见性

export class Branch {
  constructor() {
    this.alpha = 0;
    this.d = null; this.u = null;
    this.r = 0; this.t = 0; this.du2 = 0;
    this.exists = false; this.vis = false;   // vis = 链有效 且 几何可见
    this.visGeom = false;                    // 只看几何 (存在 α 根 且 r>0 且 t>0)
  }
}

/** 解 C=0: A·cos s + B·sin s + C0 = 0 的两支, 并给出天球方向与 r/t. */
export function pairAlpha(xi, c) {
  const ch = chain(xi, c, false);
  const q = mv(c.r0, cr(xi, ch.f));
  const w = c.w, Dv = c.dvec;
  const A = dt(q, Dv) - dt(w, q) * dt(w, Dv);
  const B = dt(cr(w, q), Dv);
  const C0 = dt(w, q) * dt(w, Dv);
  const RAB = Math.sqrt(A * A + B * B);
  const exists = RAB >= Math.abs(C0) - 1e-14;
  const phi0 = Math.atan2(B, A);
  const arg = RAB < 1e-300 ? 2.0 : -C0 / RAB;
  const dlt = Math.acos(Math.max(-1, Math.min(1, arg)));
  const FLab0 = mv(c.r0, ch.f);
  const XiLab0 = mv(c.r0, xi);

  const outb = new Array(2);
  for (let k = 0; k < 2; k++) {
    const alpha = phi0 + (k === 0 ? dlt : -dlt);
    const b = new Branch();
    b.alpha = alpha; b.exists = exists;
    b.d = sc(rotV(w, alpha, FLab0), -1.0);
    b.u = rotV(w, alpha, XiLab0);
    const DxU = cr(b.d, b.u);
    const du2 = dt(DxU, DxU);
    b.du2 = du2;
    if (du2 < 1e-30) { b.r = NaN; b.t = NaN; }
    else {
      b.r = dt(cr(Dv, b.u), DxU) / du2;
      b.t = dt(cr(Dv, b.d), DxU) / du2;
    }
    b.visGeom = exists && !isNaN_(b.r) && b.r > 0 && b.t > 0;
    b.vis = ch.valid && b.visGeom;
    outb[k] = b;
  }
  return outb;
}

// ================================================== 主管线: 场 → 可见天球/空间弧段

/** 通用管线: 追踪零线 → 致密化 → 配对 → 可见窗口细化 → 段列表. */
export function curvesFromField(c, f, kind, extraValid, traceH = 0.008, coarse = 3e-4,
  windowPts = 400, cancelled = null, maxSegs = 4000) {
  const segs = [];
  const curves = traceZeros(f, traceH, 6000, 40, 96, 1440, cancelled);
  for (const cur of curves) {
    if (cancelled && cancelled()) break;
    if (segs.length >= maxSegs) break;
    const cd = densify(cur, f, coarse);
    const vis = new Array(cd.length);
    const brs = new Array(cd.length);
    for (let i = 0; i < cd.length; i++) {
      brs[i] = pairAlpha(cd[i], c);
      const ex = extraValid === null || extraValid === undefined || extraValid(cd[i]);
      vis[i] = [brs[i][0].visGeom && ex, brs[i][1].visGeom && ex];
    }
    for (let bi = 0; bi < 2; bi++) {
      let i0 = -1;
      for (let i = 0; i <= cd.length; i++) {
        const on = i < cd.length && vis[i][bi];
        if (on && i0 < 0) i0 = i;
        if (!on && i0 >= 0) {
          const a = Math.max(0, i0 - 2), b = Math.min(cd.length - 1, i - 1 + 2);
          const sub = cd.slice(a, b + 1);
          let arc = 0;
          for (let q = 0; q + 1 < sub.length; q++)
            arc += Math.acos(Math.max(-1, Math.min(1, dt(sub[q], sub[q + 1]))));
          const fine = densify(sub, f, Math.max(arc / windowPts, 1e-7));
          emitRuns(c, fine, bi, kind, extraValid, segs);
          i0 = -1;
        }
      }
    }
  }
  return segs;
}

function emitRuns(c, pts, bi, kind, extraValid, segs) {
  let cur = null;
  for (const xi of pts) {
    const br = pairAlpha(xi, c)[bi];
    const ok = br.visGeom && (extraValid === null || extraValid === undefined || extraValid(xi)) && br.r <= c.maxDist;
    if (!ok) { if (cur !== null && cur.p.length >= 2) segs.push(cur); cur = null; continue; }
    if (cur === null) cur = new Seg(kind, bi);
    cur.xi.push(xi); cur.d.push(br.d); cur.alpha.push(br.alpha);
    cur.r.push(br.r); cur.t.push(br.t);
    cur.p.push(sc(br.d, br.r));
    if (br.r > 0.98 * c.maxDist) cur.truncated = true;
  }
  if (cur !== null && cur.p.length >= 2) segs.push(cur);
}

/** K1 折叠焦散 (一自由度). */
export function fold1Dof(c, traceH = 0.008, cancelled = null) {
  const f = chiFieldOf(c, c.what);
  return curvesFromField(c, f, 'K1', xi => chain(xi, c, false).valid, traceH, 3e-4, 400, cancelled);
}

/** K5 缝像 (掠射 / TIR). */
export function seam1Dof(c, which, traceH = 0.008, cancelled = null) {
  const f = seamFieldOf(c, which);
  // 缝上链的严格有效性必然踩边界 -> 用容差版 SeamValid 作为唯一可见性来源
  return curvesFromField(c, f, 'K5:' + which, xi => seamValid(xi, c, which, 1e-6),
    traceH, 3e-4, 400, cancelled, 600);
}

/** 该光路的全部缝名 (入射掠射 + 各镜面掠射 + 出射 TIR). */
export function seamNames(c) {
  const l = ['entry', 'exitTir'];
  for (let i = 0; i < Math.max(0, c.faces.length - 2); i++) l.push('mir' + i);
  return l;
}

// ================================================== K4 渐近圆 (逆反射点)

/** 解 F(ξ) = −ξ (逆反射). 返回 [{xi, elev}] (ξ*, 圆仰角 rad); 无界族专属. */
export function retroPoints(c, nSeed = 4000) {
  const sols = [];
  const gold = Math.PI * (3 - Math.sqrt(5));
  const cands = [];
  for (let i = 0; i < nSeed; i++) {
    const z = 1 - 2 * (i + 0.5) / nSeed;
    const rr = Math.sqrt(Math.max(0, 1 - z * z));
    const p = [rr * Math.cos(gold * i), rr * Math.sin(gold * i), z];
    const ch = chain(p, c, false);
    if (!ch.valid) continue;
    if (nrm(ad(ch.f, p)) < 0.25) cands.push(p);
  }
  for (const seed of cands) {
    let p = seed.slice();
    let good = false;
    for (let it = 0; it < 60; it++) {
      const ch = chain(p, c, true);
      if (!ch.valid || ch.df === null) break;
      const G = ad(ch.f, p);
      if (nrm(G) < 1e-13) { good = true; break; }
      const [e1, e2] = perpBasis(p);
      const J1 = ad(mv(ch.df, e1), e1);
      const J2 = ad(mv(ch.df, e2), e2);
      // 最小二乘 2×2 正规方程
      const a11 = dt(J1, J1), a12 = dt(J1, J2), a22 = dt(J2, J2);
      const b1 = -dt(J1, G), b2 = -dt(J2, G);
      const det = a11 * a22 - a12 * a12;
      if (Math.abs(det) < 1e-30) break;
      const s1 = (b1 * a22 - a12 * b2) / det, s2 = (a11 * b2 - a12 * b1) / det;
      let step = [s1 * e1[0] + s2 * e2[0], s1 * e1[1] + s2 * e2[1], s1 * e1[2] + s2 * e2[2]];
      const sn = nrm(step);
      if (sn > 0.2) step = sc(step, 0.2 / sn);
      p = nz(ad(p, step));
    }
    if (!good) continue;
    const chk = chain(p, c, false);
    if (!chk.valid || nrm(ad(chk.f, p)) > 1e-10) continue;
    let dup = false;
    for (const s of sols) if (Math.abs(dt(s.xi, p)) > 0.99999) { dup = true; break; }
    if (dup) continue;
    const elev = Math.asin(Math.max(-1, Math.min(1, dt(mv(c.r0, p), c.w))));
    sols.push({ xi: p, elev });
  }
  return sols;
}

/** 渐近圆画在远天球上 (半径 skyR, 以眼为心) — 观测者视角的 K4. 返回 Float32Array[] (xyz 交错). */
export function retroCircles(c, skyR = 200.0, nPts = 721) {
  const outl = [];
  for (const s of retroPoints(c)) {
    const v0 = mv(c.r0, s.xi);
    const arr = new Float32Array(nPts * 3);
    for (let i = 0; i < nPts; i++) {
      const d = rotV(c.w, 2 * Math.PI * i / (nPts - 1.0), v0);
      arr[3 * i] = c.eyePos[0] + skyR * d[0];
      arr[3 * i + 1] = c.eyePos[1] + skyR * d[1];
      arr[3 * i + 2] = c.eyePos[2] + skyR * d[2];
    }
    outl.push(arr);
  }
  return outl;
}

// ================================================== 二自由度 (归约: χ_内=0 曲线 → 解内层角 → 闭式外层角)

/** Column = (α 绕天顶, θ 绕晶体 c 轴) 定 ζ;  Lowitz = (α 绕天顶, ζ 绕 Rz(α)ŷ) 定 θ. */
export const Mode2 = Object.freeze({ Column: 0, Lowitz: 1 });

export class Ctx2 {
  constructor(o = {}) {
    this.faces = o.faces ?? null;
    this.baseN = o.baseN ?? null;
    this.n = o.n ?? 1.31;
    this.mode = o.mode ?? Mode2.Column;
    this.fixedRad = o.fixedRad ?? 0;    // Column: ζ; Lowitz: θ
    this.dvec = o.dvec ?? null;
    this.eyePos = o.eyePos ?? null;
    this.maxDist = o.maxDist ?? 80.0;
  }

  /** 内层角 b 固定后的一自由度场景 (外层角 = α, 轴 = 实验室 ẑ). */
  inner(b) {
    const zeta = this.mode === Mode2.Column ? this.fixedRad : b;
    const theta = this.mode === Mode2.Column ? b : this.fixedRad;
    const c = new Ctx({
      faces: this.faces, baseN: this.baseN, n: this.n,
      r0: mm(rotY(zeta), rotZ(theta)),
      w: [0.0, 0.0, 1.0],
      dvec: this.dvec, eyePos: this.eyePos, maxDist: this.maxDist,
    });
    c.finish();
    return c;
  }

  /** 内层生成元 (体系恒定) —— χ_内 是纯 ξ 函数. */
  omegaInner() {
    if (this.mode === Mode2.Column) return [0.0, 0.0, 1.0];          // 绕晶体 c 轴
    const t = this.fixedRad;
    return [Math.sin(t), Math.cos(t), 0.0];                          // Rz(−θ)·ŷ
  }

  /** 外层生成元 Ω_α = R^T·ẑ = Rz(−θ)Ry(−ζ)ẑ (只依赖内层角, 与 α 无关). */
  omegaOuter(b) {
    const zeta = this.mode === Mode2.Column ? this.fixedRad : b;
    const theta = this.mode === Mode2.Column ? b : this.fixedRad;
    const sz = Math.sin(zeta), cz = Math.cos(zeta);
    const st = Math.sin(theta), ct = Math.cos(theta);
    return [-sz * ct, sz * st, cz];
  }
}

/** 沿 χ_内=0 曲线上的点 ξ, 解 χ_外(ξ; b)=0 求内层角 b (可能多解). */
export function solveInner(xi, c2, nScan = 360) {
  const res = [];
  const probe = c2.inner(0.0);
  const o = { ok: false };
  const g = b => {
    const v = chi(xi, probe, c2.omegaOuter(b), o);   // chi 只用 faces/baseN/n
    return o.ok ? v : NaN;
  };
  let prev = NaN;
  for (let i = 0; i <= nScan; i++) {
    const b = 2 * Math.PI * i / nScan;
    const v = g(b);
    if (i > 0 && !isNaN_(prev) && !isNaN_(v) && Math.sign(prev) !== Math.sign(v)) {
      let lo = 2 * Math.PI * (i - 1) / nScan, hi = b, vlo = prev;
      for (let k = 0; k < 50; k++) {
        const mid = 0.5 * (lo + hi);
        const vm = g(mid);
        if (isNaN_(vm)) break;
        if (Math.sign(vm) === Math.sign(vlo)) { lo = mid; vlo = vm; }
        else hi = mid;
      }
      res.push(0.5 * (lo + hi));
    }
    prev = v;
  }
  return res;
}

// .NET Dictionary<string,T> 的枚举次序仿真: 按槽位序; Remove 的槽进 LIFO 空闲链, 下次 Add 先复用.
class NetDict {
  constructor() { this.slotKey = []; this.slotVal = []; this.slotUsed = []; this.free = []; this.idx = new Map(); }
  get(k) { const i = this.idx.get(k); return i === undefined ? undefined : this.slotVal[i]; }
  set(k, v) {
    let i = this.idx.get(k);
    if (i === undefined) {
      if (this.free.length > 0) i = this.free.pop();
      else { i = this.slotKey.length; this.slotKey.push(null); this.slotVal.push(null); this.slotUsed.push(false); }
      this.idx.set(k, i); this.slotKey[i] = k; this.slotUsed[i] = true;
    }
    this.slotVal[i] = v;
  }
  remove(k) {
    const i = this.idx.get(k);
    if (i === undefined) return false;
    this.idx.delete(k); this.slotUsed[i] = false; this.slotKey[i] = null; this.slotVal[i] = null;
    this.free.push(i);
    return true;
  }
  *entries() {
    for (let i = 0; i < this.slotKey.length; i++) if (this.slotUsed[i]) yield [this.slotKey[i], this.slotVal[i]];
  }
}

/** 二自由度观测者焦散: 追 χ_内=0 (灯无关) → 每点解内层角 → 闭式外层角.
 *  返回散点式弧段 (按曲线顺序串起来). */
export function fold2Dof(c2, traceH = 0.01, maxPts = 400, cancelled = null) {
  const probe = c2.inner(0.0);
  const omIn = c2.omegaInner();
  const f = chiFieldOf(probe, omIn);
  const curves = traceZeros(f, traceH, 6000, 40, 96, 1440, cancelled);
  const segs = [];
  for (const cur of curves) {
    if (cancelled && cancelled()) break;
    const step = Math.max(1, Math.trunc(cur.length / Math.max(1, maxPts)));
    // 每个 (分支, 内层解序号) 各串一条
    const open = new NetDict();
    for (let i = 0; i < cur.length; i += step) {
      if (cancelled && cancelled()) break;
      const xi = cur[i];
      const bs = solveInner(xi, c2);
      const live = new Set();
      for (let bi = 0; bi < bs.length; bi++) {
        const cc = c2.inner(bs[bi]);
        const brs = pairAlpha(xi, cc);
        for (let k = 0; k < 2; k++) {
          if (!brs[k].vis || brs[k].r > cc.maxDist) continue;
          const key = bi + '/' + k;
          live.add(key);
          let sg = open.get(key);
          if (sg === undefined) {
            sg = new Seg('K1-2dof', k);
            open.set(key, sg);
          }
          sg.xi.push(xi); sg.d.push(brs[k].d); sg.alpha.push(brs[k].alpha);
          sg.r.push(brs[k].r); sg.t.push(brs[k].t); sg.inner.push(bs[bi]);
          sg.p.push(sc(brs[k].d, brs[k].r));
          if (brs[k].r > 0.98 * cc.maxDist) sg.truncated = true;
        }
      }
      // 断开的 key 收尾
      const drop = [];
      for (const [k] of open.entries()) if (!live.has(k)) drop.push(k);
      for (const k of drop) {
        if (open.get(k).p.length >= 2) segs.push(open.get(k));
        open.remove(k);
      }
    }
    for (const [, v] of open.entries()) if (v.p.length >= 2) segs.push(v);
  }
  return segs;
}

// ================================================== 便捷构造 (与 Divergentlightform 的 (ζ,θ,α) 约定对齐)

/** 一自由度场景: sweep ∈ "alpha"/"zeta"/"theta". R(α,ζ,θ) = Rz(α)·Ry(ζ)·Rz(θ). */
export function make1Dof(faces, baseN, n, sweep, alphaRad, zetaRad, thetaRad, dvec, eyePos, maxDist) {
  let R0, W;
  if (sweep === 'zeta') {
    R0 = rotZ(alphaRad + thetaRad);
    W = mv(rotZ(alphaRad), [0.0, 1.0, 0.0]);
  } else if (sweep === 'theta') {
    R0 = mm(rotZ(alphaRad), rotY(zetaRad));
    W = mv(R0, [0.0, 0.0, 1.0]);
  } else {
    R0 = mm(rotY(zetaRad), rotZ(thetaRad));
    W = [0.0, 0.0, 1.0];
  }
  const c = new Ctx({ faces, baseN, n, r0: R0, w: W, dvec, eyePos, maxDist });
  c.finish();
  return c;
}

// ================================================== .NET 数字格式 (自检文本与 C# 逐字一致)
const INF = '∞';   // .NET 5+ (ICU, en-US / zh-CN) 的 PositiveInfinitySymbol
function fmtF(x, d) {               // "F{d}"
  if (isNaN_(x)) return 'NaN';
  if (x === Infinity) return INF;
  if (x === -Infinity) return '-' + INF;
  if (Object.is(x, -0)) return '-' + (0).toFixed(d);
  return x.toFixed(d);
}
function fmtE(x, d) {               // "E{d}": 尾数 d 位小数, 指数带符号且至少三位
  if (isNaN_(x)) return 'NaN';
  if (x === Infinity) return INF;
  if (x === -Infinity) return '-' + INF;
  const neg = x < 0 || Object.is(x, -0);
  const s = Math.abs(x).toExponential(d);
  const m = /^([\d.]+)e([+-])(\d+)$/.exec(s);
  return (neg ? '-' : '') + m[1] + 'E' + m[2] + m[3].padStart(3, '0');
}
function fmtP1(x) {                 // "P1" (en-US: "n%"): 精确十进制在第 3 位小数处舍入再移位
  if (isNaN_(x)) return 'NaN';
  const neg = x < 0;
  const s = Math.abs(x).toFixed(3);
  const [ip, fp] = s.split('.');
  const pct = (BigInt(ip) * 100n + BigInt(fp.slice(0, 2))).toString() + '.' + fp[2];
  return (neg ? '-' : '') + pct + '%';
}
const padL = (s, w) => String(s).padStart(w);    // {x,w}
const padR = (s, w) => String(s).padEnd(w);      // {x,-w}
const boolNet = b => (b ? 'True' : 'False');

// ================================================== 自检

/** 自检: 板晶 3-5 的 χ=0 应等价于 Bravais 对称通过 (s=1); 逆反射点上 χ=0. */
export function selfTest(baseN) {
  const lines = [];
  const se = 10.0 * Math.PI / 180;
  const dv = [15 * Math.cos(se), 0.0, 15 * Math.sin(se)];
  const eye = [0, 0, 0];

  // (1) 板晶 3-5: 水平最小偏向 21.84°
  const cPlate = make1Dof([2, 4], baseN, 1.31, 'alpha', 0, 0, 0, dv, eye, 80);
  let devMin = 999;
  for (let i = 0; i < 4000; i++) {
    const ang = (15.0 + 74.5 * i / 3999.0) * Math.PI / 180;
    const xi = [-Math.cos(ang), Math.sin(ang), 0.0];
    const ch = chain(xi, cPlate, false);
    if (!ch.valid) continue;
    const d = Math.acos(Math.max(-1, Math.min(1, dt(xi, ch.f)))) * 180 / Math.PI;
    if (d < devMin) devMin = d;
  }
  const dTheory = 2 * (Math.asin(1.31 * 0.5) - Math.PI / 6) * 180 / Math.PI;
  lines.push(`[1] 板晶3-5 最小偏向 ${fmtF(devMin, 6)}° vs 理论 ${fmtF(dTheory, 6)}°  Δ=${fmtE(Math.abs(devMin - dTheory), 2)}`);

  // (2) Parry 6-3-5 / 3-5-8 的逆反射点上 χ 应为 0
  const cPar = make1Dof([5, 2, 4], baseN, 1.31, 'alpha', 0, -Math.PI / 2, 0, dv, eye, 80);
  for (const rp of retroPoints(cPar)) {
    const o = { ok: false };
    const chiV = chi(rp.xi, cPar, cPar.what, o);
    const ch = chain(rp.xi, cPar, false);
    lines.push(`[2] 6-3-5 逆反射 仰角 ${fmtF(rp.elev * 180 / Math.PI, 3)}°  |F+ξ|=${fmtE(nrm(ad(ch.f, rp.xi)), 2)}  χ=${fmtE(chiV, 2)}`);
  }

  // (3) 配对一致性: 位置重构 |r·D − (𝔻 + t·û)|
  let worst = 0, cnt = 0;
  for (const sg of fold1Dof(cPar))
    for (let i = 0; i < sg.p.length; i++) {
      const br = pairAlpha(sg.xi[i], cPar);
      const k = sg.branch;
      const P2 = ad(dv, sc(br[k].u, br[k].t));
      worst = Math.max(worst, nrm(sb(sg.p[i], P2)));
      cnt++;
    }
  lines.push(`[3] Parry 6-3-5 折叠段位置重构最大残差 ${fmtE(worst, 2)} m  (N=${cnt})`);
  return lines.map(l => l + '\n').join('');
}

/** 世界系独立 Snell 链 (不复用晶体系 Chain) — 运行时校验脊线是真实晕点. */
function worldChain(u, wn, faces, n, tol = 0.0) {
  const eta = 1.0 / n;
  const last = faces.length - 1;
  const a = wn[faces[0]];
  let ci = -dt(u, a);
  if (ci <= -tol) return null;
  if (ci < 0) ci = 0;
  const d0 = 1 - eta * eta * (1 - ci * ci);
  if (d0 <= 0) return null;
  const v = [0, 0, 0];
  for (let i = 0; i < 3; i++) v[i] = eta * u[i] + (eta * ci - Math.sqrt(d0)) * a[i];
  for (let q = 1; q < last; q++) {
    const m = wn[faces[q]];
    const dd = dt(v, m);
    if (dd <= -tol) return null;
    for (let i = 0; i < 3; i++) v[i] -= 2 * dd * m[i];
  }
  const me = wn[faces[last]];
  const ce = dt(v, me);
  let de = 1 - n * n * (1 - ce * ce);
  if (ce <= -tol || de <= -tol) return null;
  if (de < 0) de = 0;
  const F = [0, 0, 0];
  for (let i = 0; i < 3; i++) F[i] = n * v[i] - (n * ce - Math.sqrt(de)) * me[i];
  return F;
}

function rotAxisMat(w, a) {
  const c = Math.cos(a), s = Math.sin(a);
  const M = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let i = 0; i < 3; i++)
    for (let j = 0; j < 3; j++)
      M[i][j] = (i === j ? c : 0.0) + (1 - c) * w[i] * w[j];
  M[0][1] += -s * w[2]; M[0][2] += s * w[1];
  M[1][0] += s * w[2]; M[1][2] += -s * w[0];
  M[2][0] += -s * w[1]; M[2][1] += s * w[0];
  return M;
}

/** 校验一组脊线: 用世界系独立链验证 灯→晶体→眼. r0Of(seg, i) 给每点的 R0.
 *  返回 {count, worst} (点数, 最大残差). */
export function verifySegs(segs, c, r0Of, tol = 0.0) {
  let cnt = 0, worst = 0;
  for (const sg of segs)
    for (let i = 0; i < sg.p.length; i++) {
      const Rw = mm(rotAxisMat(c.w, sg.alpha[i]), r0Of(sg, i));
      const wn = new Array(c.baseN.length);
      for (let f = 0; f < c.baseN.length; f++) wn[f] = mv(Rw, c.baseN[f]);
      const u = nz(sb(sg.p[i], c.dvec));
      const e = worldChain(u, wn, c.faces, c.n, tol);
      cnt++;
      if (e === null) { worst = Infinity; continue; }
      worst = Math.max(worst, nrm(ad(e, sg.d[i])));   // 出射传播应 = −D
    }
  return { count: cnt, worst };
}

/** 运行时自检: 一自由度 (含窄窗口/缝像) + 二自由度, 全部世界系独立校验. */
export function selfTest2(baseN) {
  const lines = [];
  const se = 10.0 * Math.PI / 180;
  const dv = [15 * Math.cos(se), 0.0, 15 * Math.sin(se)];
  const eye = [0, 0, 0];

  for (const [name, faces, zeta] of [
    ['Parry 6-3-5', [5, 2, 4], -Math.PI / 2],
    ['Parry 3-5-8', [2, 4, 7], -Math.PI / 2],
    ['Parry 3-5', [2, 4], -Math.PI / 2],
    ['板晶 3-5', [2, 4], 0.0]]) {
    const c = make1Dof(faces, baseN, 1.31, 'alpha', 0, zeta, 0, dv, eye, 80);
    const R0 = mm(rotY(zeta), rotZ(0));
    const fold = fold1Dof(c);
    const v = verifySegs(fold, c, () => R0);
    let bd = 0, inf = 0;
    for (const sg of fold) { if (sg.truncated) inf++; else bd++; }
    lines.push(`[1dof] ${padR(name, 12)} 段 ${padL(fold.length, 2)} (有界${bd}/无界${inf})  ` +
      `点 ${padL(v.count, 5)}  世界链残差 ${fmtE(v.worst, 2)}`);
  }

  {
    const zeta = -Math.PI / 2;
    const c = make1Dof([5, 2, 4], baseN, 1.31, 'alpha', 0, zeta, 0, dv, eye, 80);
    const R0 = mm(rotY(zeta), rotZ(0));
    let tot = 0, w = 0, nseg = 0;
    for (const nm of seamNames(c)) {
      const sm = seam1Dof(c, nm);
      nseg += sm.length;
      const v = verifySegs(sm, c, () => R0, 1e-6);   // 缝恰在边界上, 需容差
      tot += v.count; w = Math.max(w, v.worst);
    }
    lines.push(`[缝像] Parry 6-3-5 全缝   段 ${padL(nseg, 2)}              点 ${padL(tot, 5)}  世界链残差 ${fmtE(w, 2)}`);
  }

  for (const [name, faces] of [['3-5', [2, 4]], ['6-3-5', [5, 2, 4]]]) {
    const c2 = new Ctx2({
      faces, baseN, n: 1.31, mode: Mode2.Column,
      fixedRad: -Math.PI / 2, dvec: dv, eyePos: eye, maxDist: 80,
    });
    const fold = fold2Dof(c2, 0.01, 200);
    const cInner = c2.inner(0.0);
    // 二自由度: 每点的 R0 = Ry(ζ)·Rz(θ_i), θ_i 为该点的内层解
    const v = verifySegs(fold, cInner, (sg, i) => mm(rotY(-Math.PI / 2), rotZ(sg.inner[i])));
    lines.push(`[2dof] Column ${padR(name, 6)}      段 ${padL(fold.length, 2)}              ` +
      `点 ${padL(v.count, 5)}  世界链残差 ${fmtE(v.worst, 2)}`);
  }
  return lines.map(l => l + '\n').join('');
}

// ================================================== 体检: "这条光路到底有没有焦散, 没有是为什么"

/** 焦散存在性的逐级体检结果. */
export class Diag {
  constructor() {
    this.validFrac = 0;        // 链有效域占入射球面的比例
    this.posFrac = 0;          // 有效域内 χ>0 的比例 (0 或 1 ⟹ 不变号)
    this.chiSignChange = false;// χ 是否变号 ⟹ 上游是否存在焦散曲线
    this.curves = 0; this.curvePts = 0;   // 追到的 χ=0 曲线
    this.noAlphaRoot = 0;      // 致密化后: 锥不可解 (|C0|>√(A²+B²)) 的点数
    this.behindEye = 0;        // r ≤ 0 (晶体在眼后方)
    this.behindLamp = 0;       // t ≤ 0 (晶体在灯后方)
    this.beyondMaxDist = 0;    // r > 截断距离
    this.visible = 0;          // 最终可见点数
    this.segs = 0;             // 最终弧段数
    this.verdict = '';         // 人话结论
  }
  toString() { return this.verdict; }
}

/** 逐级体检: 无焦散时给出确切原因 (真的没有 / 有但不可见 / 被截断). */
export function diagnose(c, nSphere = 60000, traceH = 0.008, cancelled = null) {
  const d = new Diag();
  // (1) 球面均匀采样: 有效域大小 + χ 是否变号
  const gold = Math.PI * (3 - Math.sqrt(5));
  let nv = 0, npos = 0;
  const o = { ok: false };
  for (let i = 0; i < nSphere; i++) {
    const z = 1 - 2 * (i + 0.5) / nSphere;
    const rr = Math.sqrt(Math.max(0, 1 - z * z));
    const p = [rr * Math.cos(gold * i), rr * Math.sin(gold * i), z];
    const v = chi(p, c, c.what, o);
    if (!o.ok || isNaN_(v)) continue;
    nv++; if (v > 0) npos++;
  }
  d.validFrac = nSphere > 0 ? nv / nSphere : 0;
  d.posFrac = nv > 0 ? npos / nv : 0;
  d.chiSignChange = nv > 0 && npos > 0 && npos < nv;

  if (nv === 0) {
    d.verdict = T('该光路在此取向下不通 (链有效域为空) — 没有晕, 自然没有焦散。',
      'This path is blocked at this orientation (empty valid domain) — no halo, hence no caustic.');
    return d;
  }
  if (!d.chiSignChange) {
    d.verdict = _english
      ? 'There is genuinely no caustic (not a computation failure): χ keeps a constant ' +
        (npos === nv ? 'positive' : 'negative') + ' sign over ' +
        'the whole valid domain, so its zero set is empty ⟹ the sky map is nowhere singular and this ' +
        'path produces no fold line at this orientation. Its sharp edge comes from the grazing/TIR ' +
        'seam instead — tick "With seams".'
      : '本来就没有焦散 (不是没算出来): χ 在整个有效域内恒' + (npos === nv ? '正' : '负') + '号, 零集为空 ⟹ 天球映射处处非奇异, ' +
        '该光路在此取向下不产生折叠亮线。它的锐利边界来自掠射/TIR 焊缝, 请勾选"含缝像"。';
    return d;
  }

  // (2) χ=0 曲线确实存在 → 逐级看可见性在哪一关被卡掉
  const f = chiFieldOf(c, c.what);
  const curves = traceZeros(f, traceH, 6000, 40, 96, 1440, cancelled);
  d.curves = curves.length;
  for (const cur of curves) {
    const cd = densify(cur, f, 3e-4);
    d.curvePts += cd.length;
    for (const xi of cd) {
      if (!chain(xi, c, false).valid) continue;
      for (const br of pairAlpha(xi, c)) {
        if (!br.exists) { d.noAlphaRoot++; continue; }
        if (isNaN_(br.r) || br.r <= 0) { d.behindEye++; continue; }
        if (br.t <= 0) { d.behindLamp++; continue; }
        if (br.r > c.maxDist) { d.beyondMaxDist++; continue; }
        d.visible++;
      }
    }
  }
  // 权威结果来自 Fold1Dof (含可见窗口细化); 上面的粗计数只用于"没有时归因"
  const segs = fold1Dof(c, traceH, cancelled);
  d.segs = segs.length;
  let finePts = 0;
  for (const sg of segs) finePts += sg.p.length;

  if (d.segs > 0) {
    d.visible = finePts;
    d.verdict = _english
      ? `Caustic found: ${d.segs} arc(s) / ${finePts} points.`
      : `有焦散: ${d.segs} 段 / ${finePts} 个点。`;
  } else if (d.visible > 0)
    d.verdict = _english
      ? `The caustic curve exists and ${d.visible} coarse points fall in the visible region, but no ` +
        'continuous arc survives refinement (visible window too narrow) — try nudging lamp elevation/distance.'
      : `焦散曲线存在, 粗扫有 ${d.visible} 个点落在可见区, 但细化后未形成连续弧段 ` +
        '(可见窗口过窄) —— 可尝试微调灯仰角/距离。';
  else if (d.behindLamp > 0 && d.behindLamp >= d.behindEye && d.behindLamp >= d.beyondMaxDist)
    d.verdict = T('焦散曲线存在, 但全部落在**灯的后方** (t<0) —— 光没法从灯出发打到那里, 物理上不可见。',
      'The caustic curve exists but lies entirely BEHIND the lamp (t<0) — light cannot reach it, so it is physically invisible.');
  else if (d.beyondMaxDist > 0 && d.beyondMaxDist >= d.behindEye)
    d.verdict = _english
      ? `The caustic curve exists but lies entirely beyond the cutoff distance (${fmtF(c.maxDist, 0)} m) — raise "Max distance" to see it.`
      : `焦散曲线存在, 但全部超出截断距离 (${fmtF(c.maxDist, 0)} m) —— 调大"最大距离"即可看到。`;
  else if (d.behindEye > 0)
    d.verdict = T('焦散曲线存在, 但全部落在**眼的后方** (r<0) —— 不在视野内。',
      'The caustic curve exists but lies entirely BEHIND the eye (r<0) — outside the field of view.');
  else if (d.noAlphaRoot > 0)
    d.verdict = T('焦散曲线存在, 但锥上配对无解 (基线圆与散射面不相交) —— 该灯位下取不到取向。',
      'The caustic curve exists but the cone pairing has no root (baseline circle misses the scattering plane) — no orientation is reachable at this lamp position.');
  else
    d.verdict = T('焦散曲线存在但无可见点 (原因未归类)。',
      'The caustic curve exists but has no visible point (cause unclassified).');
  return d;
}

/** 批量体检 (供自检/排查用). paths: "3-5" 这类串 (注意: 同 C#, 面号直接 −1, 不走紧凑映射). */
export function diagnoseReport(baseN, zetaDeg, thetaDeg, lampElevDeg, lampDist, eyeHeight, maxDist, ...paths) {
  let s = '';
  const se = lampElevDeg * Math.PI / 180;
  const dv = [lampDist * Math.cos(se), 0.0, lampDist * Math.sin(se) - eyeHeight];
  const eye = [0.0, 0.0, eyeHeight];
  for (const ps of paths) {
    const parts = ps.split('-');
    const faces = new Array(parts.length);
    for (let i = 0; i < parts.length; i++) faces[i] = parseIntNet(parts[i]) - 1;
    const c = make1Dof(faces, baseN, 1.31, 'alpha', 0,
      zetaDeg * Math.PI / 180, thetaDeg * Math.PI / 180, dv, eye, maxDist);
    const d = diagnose(c);
    s += `${padR(ps, 10)} 有效域${padL(fmtP1(d.validFrac), 6)}  χ>0占${padL(fmtP1(d.posFrac), 6)}  ` +
      `变号=${padR(boolNet(d.chiSignChange), 5)}  曲线${padL(d.curves, 2)}/${padL(d.curvePts, 5)}点  ` +
      `无α根${padL(d.noAlphaRoot, 5)} 眼后${padL(d.behindEye, 5)} 灯后${padL(d.behindLamp, 5)} ` +
      `超距${padL(d.beyondMaxDist, 5)} 可见${padL(d.visible, 5)}` + '\n';
    s += '           → ' + d.verdict + '\n';
  }
  return s;
}

// ================================================== 窗体侧辅助 (Divergentlightform.cs 的 ParsePath / ObsBaseNormals)

const CRYSTAL_NORMALS = [
  [0, 0, 1],                           // 面1: Basal 顶
  [0, 0, -1],                          // 面2: Basal 底
  [1, 0, 0],                           // 面3: Prism
  [0.5, Math.sqrt(3) / 2.0, 0],        // 面4
  [-0.5, Math.sqrt(3) / 2.0, 0],       // 面5
  [-1, 0, 0],                          // 面6
  [-0.5, -Math.sqrt(3) / 2.0, 0],      // 面7
  [0.5, -Math.sqrt(3) / 2.0, 0]];      // 面8

/** 紧凑索引 ci 的晶体系基准法向 (锥面: 方位同柱面, 抬起 δ=顶角/2, 上锥 +z / 下锥 −z). */
export function baseNormalFromCompact(ci, upApexDeg = 56.0, loApexDeg = 56.0) {
  if (ci < 8) return [CRYSTAL_NORMALS[ci][0], CRYSTAL_NORMALS[ci][1], CRYSTAL_NORMALS[ci][2]];
  const p = (ci < 14 ? ci - 8 : ci - 14) + 2;                  // 柱面 3..8 的行号 2..7
  const delta = (ci < 14 ? upApexDeg : loApexDeg) * Math.PI / 360.0;
  const cd = Math.cos(delta), sd = Math.sin(delta);
  return [cd * CRYSTAL_NORMALS[p][0], cd * CRYSTAL_NORMALS[p][1], (ci < 14 ? 1.0 : -1.0) * sd];
}

/** 基准法向表 (紧凑索引 0..19), 同 ObsBaseNormals(). */
export function obsBaseNormals(upApexDeg = 56.0, loApexDeg = 56.0) {
  const b = new Array(20);
  for (let f = 0; f < 20; f++) b[f] = baseNormalFromCompact(f, upApexDeg, loApexDeg);
  return b;
}

/** 光路串 → 紧凑面索引 (1-8→0-7, 13-18→8-13, 23-28→14-19); 非法返回 null. 同 ParsePath. */
export function parsePath(code) {
  if (code === null || code === undefined || code.trim() === '') return null;
  const parts = code.split(/[- ]/).filter(s => s.length > 0);
  const faces = [];
  for (const p of parts) {
    const t = p.trim();
    if (!/^[+-]?\d+$/.test(t)) return null;
    const f = parseInt(t, 10);
    const ci = (f >= 1 && f <= 8) ? f - 1
      : (f >= 13 && f <= 18) ? 8 + (f - 13)
        : (f >= 23 && f <= 28) ? 14 + (f - 23) : -1;   // 紧凑索引 (含锥面)
    if (ci < 0) return null;
    faces.push(ci);
  }
  return faces.length >= 2 ? faces : null;
}
