// 球壳法晕面网格 —— 骨架 = 特殊点 / 特殊线, 三角形按参数直接铺.
//
// 解决什么: 把一张晕面 (二自由度: 焦散或某条边界的像; 一自由度: 晕面) 铺成三角形, 边界与折痕精确地是网格的边,
//   相邻两层之间的对应由参数决定, 不靠 "天空里找最近点" 去猜.
// 工具: 以眼为心、半径 r 的球壳截晕面得到一层轮廓 (引擎闭式 + Newton, 每层直接算). 轮廓上的特殊点 —— 末端点
//   (碰到另一条边界)、折点 (天空里掉头)、焦点 (整段映到一点)、截断点 (截断球面切口) —— 跨层连成特殊线
//   (边界线、折叠线、焦线). 引擎另外直接算出的三维边界线 (opt.lines) 也是骨架的一部分: 给末端点 / 折点贴上
//   "落在哪条线、线上哪个位置" 的标签 (跨层只连同一条线上相邻的交点), 给出临界半径 (拓扑在这里变), 两层之间
//   沿它插入它的点当顶点 (网格的边就是精确边界线).
// 步骤:
//   1) 层: ln r 均匀 opt.rows 层; 临界半径 (三维边界线离眼距离的极值与端点) 两侧、灯所在球壳两侧各加一层.
//   2) 每层: 轮廓 → 特殊点, 末端点 / 折点按三维边界线贴标签.
//   3) 相邻两层的特殊点按类型、标签连起来 → 特殊线 (链号).
//   4) 每层轮廓在 "两层都还在的特殊线" 处切段; 两端连着同一对特殊线的段互相对应, 按段内弧长比例直接铺成三角形带.
//      每层的顶点固定 (特殊点 + 段内按点距均匀取), 上下两条带共用; 有标签的末端之间沿三维边界线插入它的点.
//   5) 加密 (算新层): 有段或特殊点对不上 (这两层之间有特殊线出生 / 消失 = 拓扑在变) → 二分夹到相对 EPSR;
//      对应点在天空里隔得 > GAPSKY、边界线在两层之间弯出弦高、三层对应点的径向弯曲超差、灯旁半径步太大 → 中间加层.
//   6) 拓扑变化处 (相邻两层几乎重合, 或临界半径两侧那一对): 对不上的段, 小的收口, 其余逐点就近缝;
//      灯两侧那一对层收到灯上 (晕面在灯上是锥顶).
// 输入: compute(r, { seeds }) → { arcs: [{ p, d, closed }], focal: [[x, y, z] …], states } (p / d: 一层弧上的点与单位视方向,
//   [x, y, z] 数组或长 3n 的数组; 开弧两端是引擎的精确末端点), opt = { maxD, rMin, rMax, rows, skyStep (弧度),
//   lines: [{ p (Float32Array 3n), kind ('cusp' / 其它) }], critical (半径数组, 可无: 由 lines 算), cancelled, log }
// 输出: { pos, dir, h, nrm (每个三角形三个顶点依次排), nTri, lines: [{ r, gap, curves: [{ p, d }] }] (各层的线),
//   special: [{ kind: 'fold' | 'focal', p, d }] (折叠线 / 焦线), stat }

import { criticalRadii } from './mesh2s.js';

const DEG = Math.PI / 180;

export function meshShellSurface(compute, L, E, opt = {}) {
  const t0 = Date.now();
  const Dx = L[0] - E[0], Dy = L[1] - E[1], Dz = L[2] - E[2], D = Math.hypot(Dx, Dy, Dz);
  const maxD = opt.maxD || 80, rMin = opt.rMin || 1e-3 * D, rMax = opt.rMax || maxD + D;
  const N0 = opt.rows || 72;
  let SKY = opt.skyStep > 0 ? opt.skyStep : 0.5 * DEG;
  let GAPSKY = Math.max(1.2 * DEG, 3 * SKY), MAXEDGE = Math.max(4 * GAPSKY, 5 * DEG);
  const CAPMAX = 1 * DEG, VBUDGET = opt.budget || 90000;
  const EPSR = opt.epsR || 5e-4, CRIT = 2e-4, PINL = opt.pinLamp || 1e-3;
  const MAXCALLS = opt.maxCalls || 900, MAXROWS = opt.maxRows || 700;
  const FOLDW = Math.max(0.05 * DEG, 0.25 * SKY);          // 折点检测: 前后各沿线走这么长 (天空弧长) 比方向
  const tol3 = r => 0.01 + 1e-3 * r;                       // 三维容差 (弦高、径向弯曲)
  const cancelled = opt.cancelled || null;
  const log = opt.log || null;
  const stat = { rows: 0, calls: 0, events: 0, caps: 0, loose: 0, looseDeg: 0, folds: 0, focal: 0, labelled: 0, unlabelled: 0, lineVerts: 0, ms: 0 };

  const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const angU = (a, b) => { const c = dot(a, b); return c >= 1 ? 0 : c <= -1 ? Math.PI : Math.acos(c); };
  const unit = a => { const l = Math.hypot(a[0], a[1], a[2]) || 1; return [a[0] / l, a[1] / l, a[2] / l]; };
  const dist3 = (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
  const at = (r, u) => [E[0] + r * u[0], E[1] + r * u[1], E[2] + r * u[2]];

  // ---- 三维边界线 (引擎直接算的): 弧长参数、包围盒
  const LN = [];
  for (const l of opt.lines || []) {
    const p = l && (l.p || l);
    if (!p || p.length < 6) continue;
    const n = p.length / 3, s = new Float64Array(n), lo = [Infinity, Infinity, Infinity], hi = [-Infinity, -Infinity, -Infinity];
    for (let i = 0; i < n; i++) {
      if (i) s[i] = s[i - 1] + Math.hypot(p[3 * i] - p[3 * i - 3], p[3 * i + 1] - p[3 * i - 2], p[3 * i + 2] - p[3 * i - 1]);
      for (let k = 0; k < 3; k++) { lo[k] = Math.min(lo[k], p[3 * i + k]); hi[k] = Math.max(hi[k], p[3 * i + k]); }
    }
    LN.push({ id: LN.length, kind: l.kind || 'edge', p, n, s, lo, hi });
  }
  // q 附近 (tol 之内) 的所有线 (kindOk 过滤), 每条线取最近点: → [{ id, s, d }] (可空).
  //   同一条边界线常有几份重合的副本 (不同分组 / 镜像各给一次), 只取最近的一条会在相邻层间换来换去 —— 全留下, 有共同的线就算同一条
  const nearLines = (q, kindOk, tol) => {
    const out = [];
    for (const ln of LN) {
      if (kindOk && !kindOk(ln.kind)) continue;
      const bx = Math.max(ln.lo[0] - q[0], 0, q[0] - ln.hi[0]), by = Math.max(ln.lo[1] - q[1], 0, q[1] - ln.hi[1]), bz = Math.max(ln.lo[2] - q[2], 0, q[2] - ln.hi[2]);
      if (Math.hypot(bx, by, bz) > tol) continue;
      const p = ln.p;
      let best = null, bd = tol;
      for (let i = 0; i + 1 < ln.n; i++) {
        const ax = p[3 * i], ay = p[3 * i + 1], az = p[3 * i + 2], ex = p[3 * i + 3] - ax, ey = p[3 * i + 4] - ay, ez = p[3 * i + 5] - az;
        const qx = q[0] - ax, qy = q[1] - ay, qz = q[2] - az, l2 = ex * ex + ey * ey + ez * ez;
        let u = l2 > 0 ? (qx * ex + qy * ey + qz * ez) / l2 : 0; u = u < 0 ? 0 : u > 1 ? 1 : u;
        const d = Math.hypot(qx - u * ex, qy - u * ey, qz - u * ez);
        if (d < bd) { bd = d; best = { id: ln.id, s: ln.s[i] + u * (ln.s[i + 1] - ln.s[i]), d }; }
      }
      if (best) out.push(best);
    }
    if (!out.length) return null;
    // 只留 "差不多最近" 的: 重合的副本 (距离相同) 都留; 两条线交叉处, 另一支的线 (远好几倍) 不留 —— 否则两个末端标签相同、交叉而过时分不清
    const dmin = Math.min(...out.map(o => o.d)), qr = Math.hypot(q[0] - E[0], q[1] - E[1], q[2] - E[2]);
    return out.filter(o => o.d <= 3 * dmin + 2e-6 * (1 + qr));
  };
  // 两个标签的共同线里, 线上弧长差最小的那条: → { id, sa, sb, ds } 或 null
  const sharedLine = (la, lb) => {
    if (!la || !lb) return null;
    let best = null;
    for (const a of la) for (const b of lb) if (a.id === b.id) { const ds = Math.abs(a.s - b.s); if (!best || ds < best.ds) best = { id: a.id, sa: a.s, sb: b.s, ds }; }
    return best;
  };
  // 线 id 上弧长 s0 → s1 之间 (不含两端) 的点, 按 s0 → s1 的顺序
  const linePts = (id, s0, s1) => {
    const ln = LN[id], out = [], a = Math.min(s0, s1), b = Math.max(s0, s1);
    for (let i = 0; i < ln.n; i++) if (ln.s[i] > a + 1e-9 && ln.s[i] < b - 1e-9) out.push([ln.p[3 * i], ln.p[3 * i + 1], ln.p[3 * i + 2]]);
    if (s0 > s1) out.reverse();
    return out;
  };
  // 两层之间这段线可信吗: 线上的点离眼距离落在两层半径之间 (略放宽), 且离两端连线不超过两端距离的一半 ——
  //   否则 (共同的线在两端之间绕了远路、或是另一处交点) 不用它
  const lineOk = (Q, P0, P1, r0, r1) => {
    const lo = Math.min(r0, r1) * (1 - 1e-3) - 1e-4, hi = Math.max(r0, r1) * (1 + 1e-3) + 1e-4, d01 = dist3(P0, P1);
    for (const q of Q) { const rq = Math.hypot(q[0] - E[0], q[1] - E[1], q[2] - E[2]); if (rq < lo || rq > hi) return false; }
    return chordDev(P0, P1, Q) <= 0.5 * d01 + 1e-6;
  };
  const chordDev = (P0, P1, Q) => {
    let m = 0;
    const ex = P1[0] - P0[0], ey = P1[1] - P0[1], ez = P1[2] - P0[2], l2 = ex * ex + ey * ey + ez * ez;
    for (const q of Q) {
      const qx = q[0] - P0[0], qy = q[1] - P0[1], qz = q[2] - P0[2];
      let u = l2 > 0 ? (qx * ex + qy * ey + qz * ez) / l2 : 0; u = u < 0 ? 0 : u > 1 ? 1 : u;
      m = Math.max(m, Math.hypot(qx - u * ex, qy - u * ey, qz - u * ez));
    }
    return m;
  };

  // ---- 顶点表
  let VP = new Float64Array(3 * 65536), VU = new Float64Array(3 * 65536), nv = 0;
  const addV = (p, u) => {
    if (3 * nv + 3 > VP.length) { const a = new Float64Array(2 * VP.length); a.set(VP); VP = a; const b = new Float64Array(2 * VU.length); b.set(VU); VU = b; }
    VP[3 * nv] = p[0]; VP[3 * nv + 1] = p[1]; VP[3 * nv + 2] = p[2]; VU[3 * nv] = u[0]; VU[3 * nv + 1] = u[1]; VU[3 * nv + 2] = u[2];
    return nv++;
  };
  const vu = v => [VU[3 * v], VU[3 * v + 1], VU[3 * v + 2]];
  const vang = (a, b) => angU(vu(a), vu(b));
  const tris = [];
  let nSkip = 0;
  const emit = (a, b, c) => {
    if (a === b || b === c || c === a) return;
    if (vang(a, b) > MAXEDGE || vang(b, c) > MAXEDGE || vang(c, a) > MAXEDGE) { nSkip++; return; }
    tris.push(a, b, c);
  };

  // ---- 一层: 引擎 → 弧 (天空累计弧长) → 特殊点
  let calls = 0;
  const seedsOf = (...rows) => {
    const runs = [];
    for (const rw of rows) if (rw && rw.states) for (const run of rw.states) if (run && run.length) runs.push(run);
    if (!runs.length) return null;
    const per = Math.max(3, Math.min(16, Math.floor(80 / runs.length))), out = [];
    for (const run of runs) { const st = Math.max(1, Math.floor(run.length / per)); for (let i = st >> 1; i < run.length; i += st) out.push(run[i]); out.push(run[0], run[run.length - 1]); }
    return out;
  };
  // 弧: U = 单位视方向 (3n), S = 累计天空弧长; 折点检测
  const makeArc = (row, U, n, closed, headType, tailType, joints = []) => {
    const S = new Float64Array(n);
    for (let i = 1; i < n; i++) S[i] = S[i - 1] + angU([U[3 * i - 3], U[3 * i - 2], U[3 * i - 1]], [U[3 * i], U[3 * i + 1], U[3 * i + 2]]);
    const len = S[n - 1] + (closed && n > 2 ? angU([U[3 * n - 3], U[3 * n - 2], U[3 * n - 1]], [U[0], U[1], U[2]]) : 0);
    const arc = { U, n, S, len, closed, sp: [], ai: row.arcs.length };
    const uAt = i => [U[3 * i], U[3 * i + 1], U[3 * i + 2]];
    const mk = (t, type, u) => ({ t, type, u, p: at(row.r, u), arc, lab: null, chain: -1 });
    if (!closed) { arc.sp.push(mk(0, headType, uAt(0))); arc.sp.push(mk(S[n - 1], tailType, uAt(n - 1))); }
    // 接头 (引擎分段给的同一条轮廓并接处; 退化弧两侧的两支在焦点处相接 —— 两张面交于焦线): 也是特殊点, 段在这里切开
    for (const j of joints) if ((j > 0 && j < n - 1) || closed) arc.sp.push(mk(S[Math.min(j, n - 1)], 'joint', uAt(Math.min(j, n - 1))));
    // 折点: 前后各沿线走 FOLDW 的两段弦反向 (夹角 > 120°); 连着几个点都满足的取最尖的那个 (沿来向投影最远)
    if (n >= 5 && len > 2 * FOLDW) {
      const Sx = k => closed ? (k < 0 ? S[k + n] - len : k >= n ? S[k - n] + len : S[k]) : S[k];
      const Ux = k => uAt(closed ? ((k % n) + n) % n : k);
      const lo = closed ? 0 : 1, hi = closed ? n : n - 1;
      let h = lo - 1, j = lo + 1, cl = null;
      const flush = () => {
        if (!cl) return;
        let best = cl.i0, bv = -Infinity;
        for (let k = cl.i0; k <= cl.i1; k++) { const v = dot(Ux(k), cl.a); if (v > bv) { bv = v; best = k; } }
        const bi = closed ? ((best % n) + n) % n : best;
        if (closed || (S[bi] > FOLDW && S[n - 1] - S[bi] > FOLDW)) { const f = mk(S[bi], 'fold', uAt(bi)); f.sharp = -cl.cmin; arc.sp.push(f); }
        cl = null;
      };
      for (let i = lo; i < hi; i++) {
        if (!closed && (S[i] < FOLDW || S[n - 1] - S[i] < FOLDW)) { flush(); continue; }
        while (Sx(i) - Sx(h) > FOLDW && h < i - 1) h++;
        while (Sx(i) - Sx(h) < FOLDW && h > i - n) h--;
        if (j <= i) j = i + 1;
        while (Sx(j) - Sx(i) < FOLDW && j < i + n - 1 && (closed || j < n - 1)) j++;
        const ui = Ux(i), uh = Ux(h), uj = Ux(j);
        const a = [ui[0] - uh[0], ui[1] - uh[1], ui[2] - uh[2]], b = [uj[0] - ui[0], uj[1] - ui[1], uj[2] - ui[2]];
        const la = Math.hypot(a[0], a[1], a[2]), lb = Math.hypot(b[0], b[1], b[2]);
        if (la > 1e-12 && lb > 1e-12 && dot(a, b) < -0.5 * la * lb) {
          const cs = dot(a, b) / (la * lb);
          if (!cl) cl = { i0: i, i1: i, a: [a[0] / la, a[1] / la, a[2] / la], cmin: cs }; else { cl.i1 = i; cl.cmin = Math.min(cl.cmin, cs); }
        } else flush();
      }
      flush();
      const fs = arc.sp.filter(q => q.type === 'fold').sort((x, y) => x.t - y.t);
      for (let k = 1; k < fs.length; k++) if (fs[k].t - fs[k - 1].t < FOLDW) fs[k].dup = true;
      if (closed && fs.length > 1 && fs[0].t + len - fs[fs.length - 1].t < FOLDW) fs[fs.length - 1].dup = true;
      arc.sp = arc.sp.filter(q => !q.dup);
    }
    arc.sp.sort((x, y) => x.t - y.t);
    // 标签: 末端点 → 最近的非尖点脊线的三维边界线; 折点 → 最近的尖点脊线 (焦散的折点就在尖点脊线上)
    for (const q of arc.sp) {
      if (q.type === 'cut' || q.type === 'joint') continue;
      // 容差只按 r 的比例 (所有边界线都在眼睛处汇合: 带常数项的话眼旁每个末端点都挨着所有线, 跨层连错)
      const lab = nearLines(q.p, q.type === 'fold' ? (k => k === 'cusp') : (k => k !== 'cusp'), Math.max(1e-6, 2e-3 * row.r));
      if (lab) { q.lab = lab; stat.labelled++; } else if (q.type === 'end') stat.unlabelled++;
    }
    row.arcs.push(arc);
    for (const q of arc.sp) row.sp.push(q);
    return arc;
  };
  const buildRow = (r, nbA, nbB, kind) => {
    if (cancelled && cancelled()) throw new Error('cancelled');
    calls++;
    const s = compute(r, { seeds: seedsOf(nbA, nbB) }) || { arcs: [] };
    const row = { r, kind, arcs: [], sp: [], raw: [], states: s.states || null, lmin: Infinity, lmax: -Infinity };
    for (const q of s.focal || []) {
      if (Math.hypot(q[0] - L[0], q[1] - L[1], q[2] - L[2]) > maxD) continue;
      const u = unit([q[0] - E[0], q[1] - E[1], q[2] - E[2]]);
      row.sp.push({ t: 0, type: 'focal', u, p: at(r, u), arc: null, lab: null, chain: -1 });
    }
    const pre = [];                                   // 截断后的弧 (还没并接): { V, m, closed, h, t }
    for (const arc of s.arcs || []) {
      const flat = !Array.isArray(arc.d[0]);
      let n = flat ? arc.d.length / 3 : arc.d.length;
      if (n < 1) continue;
      const U = new Float64Array(3 * n);
      for (let i = 0; i < n; i++) {
        const x = flat ? arc.d[3 * i] : arc.d[i][0], y = flat ? arc.d[3 * i + 1] : arc.d[i][1], z = flat ? arc.d[3 * i + 2] : arc.d[i][2];
        const l = Math.hypot(x, y, z) || 1; U[3 * i] = x / l; U[3 * i + 1] = y / l; U[3 * i + 2] = z / l;
      }
      // 绕了好几圈的弧 (刚出生的闭环比追踪步长还小, 引擎认不出回到起点): 在第一次回到起点处截断, 当闭环
      let closed = !!arc.closed;
      if (!closed && n > 8) {
        const d0 = [U[0], U[1], U[2]];
        let len = 0, mxStep = 0, mxFar = 0;
        for (let i = 1; i < n; i++) { const st = angU([U[3 * i - 3], U[3 * i - 2], U[3 * i - 1]], [U[3 * i], U[3 * i + 1], U[3 * i + 2]]); len += st; mxStep = Math.max(mxStep, st); }
        for (let i = 0; i < n; i++) mxFar = Math.max(mxFar, angU(d0, [U[3 * i], U[3 * i + 1], U[3 * i + 2]]));
        if (mxFar > 0 && len > 4 * Math.PI * mxFar) {
          let away = false;
          for (let i = 1; i < n; i++) {
            const f = angU(d0, [U[3 * i], U[3 * i + 1], U[3 * i + 2]]);
            if (!away) { if (f > 0.5 * mxFar) away = true; continue; }
            if (f <= Math.max(1.5 * mxStep, 1e-9)) { n = i + 1; closed = true; break; }
          }
        }
      }
      if (closed && n > 2 && angU([U[0], U[1], U[2]], [U[3 * n - 3], U[3 * n - 2], U[3 * n - 1]]) < 1e-12) n--;   // 首尾重复点
      // 截断球面 (离灯 maxD) 之内的连续段; 截口处线性插值出精确的截断点
      const ok = new Uint8Array(n);
      for (let i = 0; i < n; i++) {
        const p = at(r, [U[3 * i], U[3 * i + 1], U[3 * i + 2]]), dl = Math.hypot(p[0] - L[0], p[1] - L[1], p[2] - L[2]);
        ok[i] = dl <= maxD ? 1 : 0;
        if (dl < row.lmin) row.lmin = dl; if (dl > row.lmax) row.lmax = dl;
      }
      const cut = (a, b) => {           // a 在内 b 在外
        let lo = 0, hi = 1;
        const ua = [U[3 * a], U[3 * a + 1], U[3 * a + 2]], ub = [U[3 * b], U[3 * b + 1], U[3 * b + 2]];
        const uAt = t => unit([ua[0] + t * (ub[0] - ua[0]), ua[1] + t * (ub[1] - ua[1]), ua[2] + t * (ub[2] - ua[2])]);
        for (let it = 0; it < 40; it++) { const m = 0.5 * (lo + hi), p = at(r, uAt(m)); if (Math.hypot(p[0] - L[0], p[1] - L[1], p[2] - L[2]) <= maxD) lo = m; else hi = m; }
        return uAt(lo);
      };
      const all = ok.every(Boolean);
      if (all) { pre.push({ V: n === U.length / 3 ? U : U.slice(0, 3 * n), m: n, closed: closed && n > 2, h: 'end', t: 'end', J: [] }); continue; }
      // 有截断: 闭环从一个出界点起算, 当开弧切
      let order = [...Array(n).keys()];
      if (closed) { const k0 = ok.indexOf(0); order = order.map(i => (i + k0) % n); }
      const runs = [];
      for (let q = 0; q < order.length;) {
        if (!ok[order[q]]) { q++; continue; }
        let e = q; while (e + 1 < order.length && ok[order[e + 1]]) e++;
        runs.push([q, e]); q = e + 1;
      }
      for (const [q0, q1] of runs) {
        const idx = order.slice(q0, q1 + 1), head = q0 > 0 || closed ? cut(order[q0], order[(q0 - 1 + order.length) % order.length]) : null;
        const tail = q1 < order.length - 1 || closed ? cut(order[q1], order[(q1 + 1) % order.length]) : null;
        const m = idx.length + (head ? 1 : 0) + (tail ? 1 : 0), V = new Float64Array(3 * m);
        let k = 0;
        if (head) { V.set(head, 0); k = 1; }
        for (const i of idx) { V[3 * k] = U[3 * i]; V[3 * k + 1] = U[3 * i + 1]; V[3 * k + 2] = U[3 * i + 2]; k++; }
        if (tail) V.set(tail, 3 * k);
        if (m < 2) continue;
        pre.push({ V, m, closed: false, h: head ? 'cut' : 'end', t: tail ? 'cut' : 'end', J: [] });
      }
    }
    // 退化的弧 (整条映到天上同一点, 天空直径 < 1e-5 rad; 如一自由度 1-3-2 ζ90 每层都是): 焦点, 不当弧 (二自由度引擎已先做了这一步)
    for (let i = pre.length - 1; i >= 0; i--) {
      const a = pre[i], u0 = [a.V[0], a.V[1], a.V[2]];
      let mx = 0;
      for (let k = 1; k < a.m && mx < 1e-5; k++) mx = Math.max(mx, angU(u0, [a.V[3 * k], a.V[3 * k + 1], a.V[3 * k + 2]]));
      if (mx >= 1e-5) continue;
      pre.splice(i, 1);
      if (Math.hypot(...at(r, u0).map((x, k) => x - L[k])) <= maxD && !row.sp.some(q => q.type === 'focal' && angU(q.u, u0) < 1e-5)) row.sp.push({ t: 0, type: 'focal', u: u0, p: at(r, u0), arc: null, lab: null, chain: -1 });
    }
    // 重复的弧: 整条落在另一条弧上 (天空里 1e-4 rad 以内; 两条折线取样不同, 弯处弦高本身就到 1e-5 量级) 的去掉. 焦点附近引擎有时给出两条状态空间里不同、天空像完全
    //   重合的弧 (3-5 column ζ−90 Σ10 r≈0.65), 不去掉的话下面会把它们首尾并成闭环, 相邻层间拓扑乱跳
    {
      const DUPTOL = 1e-4;
      const near = (a, u) => {                       // u 到弧 a 折线的最小夹角 (近似: 到各段的弦)
        let bd = Infinity;
        for (let i = 0; i + 1 < a.m; i++) {
          const ax = a.V[3 * i], ay = a.V[3 * i + 1], az = a.V[3 * i + 2], ex = a.V[3 * i + 3] - ax, ey = a.V[3 * i + 4] - ay, ez = a.V[3 * i + 5] - az;
          const qx = u[0] - ax, qy = u[1] - ay, qz = u[2] - az, l2 = ex * ex + ey * ey + ez * ez;
          let t = l2 > 0 ? (qx * ex + qy * ey + qz * ez) / l2 : 0; t = t < 0 ? 0 : t > 1 ? 1 : t;
          const d = Math.hypot(qx - t * ex, qy - t * ey, qz - t * ez);
          if (d < bd) bd = d;
        }
        if (a.m === 1) bd = Math.hypot(u[0] - a.V[0], u[1] - a.V[1], u[2] - a.V[2]);
        return bd;
      };
      const lenA = a => { let l = 0; for (let i = 1; i < a.m; i++) l += angU([a.V[3 * i - 3], a.V[3 * i - 2], a.V[3 * i - 1]], [a.V[3 * i], a.V[3 * i + 1], a.V[3 * i + 2]]); return l; };
      const order = pre.map((a, i) => [lenA(a), i]).sort((x, y) => x[0] - y[0]);
      const gone = new Set();
      for (const [, i] of order) {
        const a = pre[i];
        for (let j = 0; j < pre.length; j++) {
          if (j === i || gone.has(j) || pre[j].m < 2) continue;
          let all = true;
          const st = Math.max(1, Math.floor(a.m / 24));
          for (let k = 0; k < a.m && all; k += st) if (near(pre[j], [a.V[3 * k], a.V[3 * k + 1], a.V[3 * k + 2]]) > DUPTOL) all = false;
          if (all && near(pre[j], [a.V[3 * a.m - 3], a.V[3 * a.m - 2], a.V[3 * a.m - 1]]) > DUPTOL) all = false;
          if (all) { gone.add(i); break; }
        }
      }
      if (gone.size) { const kept = pre.filter((_, i) => !gone.has(i)); pre.length = 0; pre.push(...kept); stat.dups = (stat.dups || 0) + gone.size; }
    }
    // 首尾相接的弧并成一条: 引擎把一条轮廓分成几段给 (边界像在 σ 的一支没了、另一支接着走的地方断开; 两条段在那里
    //   端点重合), 这种接头不是末端点 (不落在任何边界线上). 并接后若真是尖角, 折点检测会认出来; 首尾相接的成闭环
    const JOIN = 1e-6;
    const endU = (a, tail) => { const i = tail ? a.m - 1 : 0; return [a.V[3 * i], a.V[3 * i + 1], a.V[3 * i + 2]]; };
    const flip = a => { const V = new Float64Array(3 * a.m); for (let i = 0; i < a.m; i++) { V[3 * i] = a.V[3 * (a.m - 1 - i)]; V[3 * i + 1] = a.V[3 * (a.m - 1 - i) + 1]; V[3 * i + 2] = a.V[3 * (a.m - 1 - i) + 2]; } return { V, m: a.m, closed: false, h: a.t, t: a.h, J: a.J.map(j => a.m - 1 - j) }; };
    const cat = (a, b) => { const V = new Float64Array(3 * (a.m + b.m - 1)); V.set(a.V.subarray(0, 3 * a.m)); V.set(b.V.subarray(3, 3 * b.m), 3 * a.m); return { V, m: a.m + b.m - 1, closed: false, h: a.h, t: b.t, J: [...a.J, a.m - 1, ...b.J.map(j => j + a.m - 1)] }; };
    for (let again = true; again;) {
      again = false;
      for (let i = 0; i < pre.length && !again; i++) {
        const a = pre[i];
        if (a.closed) continue;
        if (a.m > 2 && a.h === 'end' && a.t === 'end' && angU(endU(a, false), endU(a, true)) < JOIN) { a.m--; a.closed = true; a.J.push(0); again = true; break; }
        for (let j = i + 1; j < pre.length && !again; j++) {
          const b = pre[j];
          if (b.closed) continue;
          let c = null;
          if (a.t === 'end' && b.h === 'end' && angU(endU(a, true), endU(b, false)) < JOIN) c = cat(a, b);
          else if (a.t === 'end' && b.t === 'end' && angU(endU(a, true), endU(b, true)) < JOIN) c = cat(a, flip(b));
          else if (a.h === 'end' && b.t === 'end' && angU(endU(a, false), endU(b, true)) < JOIN) c = cat(b, a);
          else if (a.h === 'end' && b.h === 'end' && angU(endU(a, false), endU(b, false)) < JOIN) c = cat(flip(a), b);
          if (c) { pre[i] = c; pre.splice(j, 1); again = true; }
        }
      }
    }
    for (const a of pre) {
      makeArc(row, a.V, a.m, a.closed, a.h, a.t, a.J);
      row.raw.push({ p: Float32Array.from({ length: 3 * a.m }, (_, k) => E[k % 3] + r * a.V[k]), d: Float32Array.from(a.V.subarray(0, 3 * a.m)) });
    }
    return row;
  };

  // ---- 弧上任意位置 (天空弧长 t) 的视方向
  const uOn = (arc, t) => {
    const S = arc.S, n = arc.n, U = arc.U;
    if (arc.closed) { t = ((t % arc.len) + arc.len) % arc.len; if (t > S[n - 1]) { const f = (t - S[n - 1]) / Math.max(1e-15, arc.len - S[n - 1]); return unit([U[3 * n - 3] + f * (U[0] - U[3 * n - 3]), U[3 * n - 2] + f * (U[1] - U[3 * n - 2]), U[3 * n - 1] + f * (U[2] - U[3 * n - 1])]); } }
    if (t <= 0) return [U[0], U[1], U[2]];
    if (t >= S[n - 1]) return [U[3 * n - 3], U[3 * n - 2], U[3 * n - 1]];
    let lo = 0, hi = n - 1;
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (S[m] <= t) lo = m; else hi = m; }
    const f = (t - S[lo]) / Math.max(1e-15, S[hi] - S[lo]);
    return unit([U[3 * lo] + f * (U[3 * hi] - U[3 * lo]), U[3 * lo + 1] + f * (U[3 * hi + 1] - U[3 * lo + 1]), U[3 * lo + 2] + f * (U[3 * hi + 2] - U[3 * lo + 2])]);
  };

  // ---- 相邻两层的特殊点连起来: 同类型; 两个都有标签的只连同一条线 (按线上弧长最近), 否则按三维距离;
  //   互为最近才连, 没标签的还要明显近于次近 (无歧义)
  const link = (A, B) => {
    const pairs = [];
    // 代价 = 视方向的夹角 (球壳上一切按角度量: 近眼处各层是同一条线按比例缩放, 三维距离里半径差占了大头, 挨着的
    //   特殊点分不清); 两个都有标签的还要有共同的线
    const cost = (a, b) => {
      if (a.type !== b.type) return Infinity;
      if (a.lab && b.lab && !sharedLine(a.lab, b.lab)) return Infinity;
      const d = angU(a.u, b.u);
      return d <= 0.35 ? d : Infinity;
    };
    const nA = A.sp.length, nB = B.sp.length;
    if (!nA || !nB) return pairs;
    const C = new Float64Array(nA * nB);
    for (let i = 0; i < nA; i++) for (let j = 0; j < nB; j++) C[i * nB + j] = cost(A.sp[i], B.sp[j]);
    const bestB = new Int32Array(nA).fill(-1), bestA = new Int32Array(nB).fill(-1), sec = new Float64Array(nA).fill(Infinity);
    for (let i = 0; i < nA; i++) { let bv = Infinity; for (let j = 0; j < nB; j++) { const c = C[i * nB + j]; if (c < bv) { sec[i] = bv; bv = c; bestB[i] = j; } else if (c < sec[i]) sec[i] = c; } if (bv === Infinity) bestB[i] = -1; }
    for (let j = 0; j < nB; j++) { let bv = Infinity; for (let i = 0; i < nA; i++) { const c = C[i * nB + j]; if (c < bv) { bv = c; bestA[j] = i; } } if (bv === Infinity) bestA[j] = -1; }
    // 第一遍: 互为最近且明显近于次近的 (无歧义) 先连
    const usedB = new Uint8Array(nB), done = new Uint8Array(nA), defer = [];
    for (let i = 0; i < nA; i++) {
      const j = bestB[i];
      if (j < 0 || bestA[j] !== i) { if (j >= 0) defer.push(i); continue; }
      if (sec[i] < Infinity && C[i * nB + j] > 0.6 * sec[i]) { defer.push(i); continue; }
      pairs.push([i, j]); usedB[j] = 1; done[i] = 1;
    }
    // 第二遍: 有歧义的 (两条弧的末端在同一条线上挨着 —— 两弧在临界半径附近相接 / 分开) 按弧的对应定: 同一条弧上已连上的
    //   特殊点连到了对面哪条弧, 就在那条弧上找; 唯一才连
    const arcMap = new Map();                 // A 的弧 → B 的弧 (按已连上的特殊点投票)
    for (const [i, j] of pairs) {
      const aa = A.sp[i].arc, bb = B.sp[j].arc;
      if (!aa || !bb) continue;
      const m = arcMap.get(aa) || new Map(); m.set(bb, (m.get(bb) || 0) + 1); arcMap.set(aa, m);
    }
    //   一条弧上没定的特殊点一起指派 (对称的 V 形小弧: 两端挨着、标签相同, 单个看分不清; 整体取总夹角最小的配法).
    //   弧的对应: 已连上的特殊点投票; 没有就取弧中点视方向最近的
    const mid = arc => uOn(arc, 0.5 * arc.len);
    const groups = new Map();
    for (const i of defer) { const a = A.sp[i]; if (!a.arc) continue; if (!groups.has(a.arc)) groups.set(a.arc, []); groups.get(a.arc).push(i); }
    for (const [aa, is] of groups) {
      let bb = null;
      const m = arcMap.get(aa);
      if (m) { let bv = 0, tie = false; for (const [k, v] of m) { if (v > bv) { bv = v; bb = k; tie = false; } else if (v === bv) tie = true; } if (tie) bb = null; }
      if (!bb) { const ma = mid(aa); let bd = 0.35; for (const arc of B.arcs) { const d = angU(ma, mid(arc)); if (d < bd) { bd = d; bb = arc; } } }
      if (!bb) continue;
      const js = [];
      for (let j = 0; j < nB; j++) if (!usedB[j] && B.sp[j].arc === bb) js.push(j);
      if (!js.length) continue;
      // 小集合穷举 (单射), 大的贪心
      let best = null, bestC = Infinity;
      if (is.length <= 5 && js.length <= 6) {
        const rec = (k, used, acc, cur) => {
          if (acc >= bestC) return;
          if (k === is.length) { if (cur.length) { bestC = acc; best = cur.slice(); } return; }
          rec(k + 1, used, acc + 1e3, cur);                     // 这一个不连 (罚分)
          for (const j of js) { if (used.has(j)) continue; const c = C[is[k] * nB + j]; if (c === Infinity) continue; used.add(j); cur.push([is[k], j]); rec(k + 1, used, acc + c, cur); cur.pop(); used.delete(j); }
        };
        rec(0, new Set(), 0, []);
      } else {
        best = [];
        const cand = [];
        for (const i of is) for (const j of js) if (C[i * nB + j] < Infinity) cand.push([C[i * nB + j], i, j]);
        cand.sort((x, y) => x[0] - y[0]);
        const ui = new Set(), uj = new Set();
        for (const [, i, j] of cand) if (!ui.has(i) && !uj.has(j)) { ui.add(i); uj.add(j); best.push([i, j]); }
      }
      for (const [i, j] of best || []) { pairs.push([i, j]); usedB[j] = 1; done[i] = 1; }
    }
    return pairs;
  };
  let nextChain = 0;
  const relink = rows => {
    nextChain = 0;
    for (const q of rows[0].sp) q.chain = nextChain++;
    for (let k = 0; k + 1 < rows.length; k++) {
      const A = rows[k], B = rows[k + 1];
      for (const q of B.sp) q.chain = -1;
      for (const [i, j] of link(A, B)) B.sp[j].chain = A.sp[i].chain;
      for (const q of B.sp) if (q.chain < 0) q.chain = nextChain++;
    }
  };

  // ---- 两层之间: 只在两层都还在的特殊线处切段 → 段的对应
  //   段 = { arc, t0, t1 (闭环可越过 len), k0, k1 (两端特殊线的链号; 不延续的端点给唯一的负号), loop }
  let uniq = -1;
  const piecesOf = (row, keep) => {
    const out = [];
    for (const arc of row.arcs) {
      const br = arc.sp.filter(q => keep.has(q.chain));
      if (!arc.closed) {
        const pts = [{ t: 0, q: arc.sp[0] }];
        for (const q of br) if (q !== arc.sp[0] && q !== arc.sp[arc.sp.length - 1]) pts.push({ t: q.t, q });
        pts.push({ t: arc.sp[arc.sp.length - 1].t, q: arc.sp[arc.sp.length - 1] });
        const key = q => keep.has(q.chain) ? q.chain : uniq--;
        for (let i = 0; i + 1 < pts.length; i++) out.push({ arc, t0: pts[i].t, t1: pts[i + 1].t, k0: key(pts[i].q), k1: key(pts[i + 1].q), loop: false });
      } else if (!br.length) {
        out.push({ arc, t0: 0, t1: arc.len, k0: 'L', k1: 'L', loop: true });
      } else {
        for (let i = 0; i < br.length; i++) {
          const a = br[i], b = br[(i + 1) % br.length];
          out.push({ arc, t0: a.t, t1: i + 1 < br.length ? b.t : b.t + arc.len, k0: a.chain, k1: b.chain, loop: false });
        }
      }
    }
    for (const pc of out) { pc.len = pc.t1 - pc.t0; pc.mid = uOn(pc.arc, 0.5 * (pc.t0 + pc.t1)); }
    return out;
  };
  // 相邻两层 A < B 的分析: 对应段 (pairs: [pa, pb, rev]), 对不上的段, 对不上的特殊点, 天空间距, 边界线弦高
  const analyze = (A, B) => {
    const chA = new Set(A.sp.map(q => q.chain)), chB = new Set(B.sp.map(q => q.chain));
    const keep = new Set([...chA].filter(c => chB.has(c)));
    const PA = piecesOf(A, keep), PB = piecesOf(B, keep);
    const used = new Set(), pairs = [], unA = [];
    for (const pa of PA) {
      let best = null, bd = Infinity, rev = false;
      for (const pb of PB) {
        if (used.has(pb)) continue;
        let ok = false, rv = false;
        if (pa.loop && pb.loop) { const q = pb.len / Math.max(pa.len, 1e-15); ok = q > 0.4 && q < 2.5 && angU(pa.mid, pb.mid) < 0.5 * Math.max(pa.len, pb.len) + GAPSKY; }
        else if (!pa.loop && !pb.loop && typeof pa.k0 === 'number' && pa.k0 >= 0 && pa.k1 >= 0) {
          if (pb.k0 === pa.k0 && pb.k1 === pa.k1) ok = true;
          else if (pb.k0 === pa.k1 && pb.k1 === pa.k0) { ok = true; rv = true; }
        }
        if (!ok) continue;
        const d = angU(pa.mid, pb.mid);
        if (d < bd) { bd = d; best = pb; rev = rv; }
      }
      // 走向再按形状核一遍: 两端 + 四分点的夹角之和, 正反取小 (挨得很近的两个折点在相邻层间可能互换了链号 ——
      //   两种连法总夹角一样; 只看链号会把整段反着拉链)
      if (best && !pa.loop) {
        const cst = r2 => { let c = 0; for (const f of [0, 0.25, 0.5, 0.75, 1]) c += angU(uOn(pa.arc, pa.t0 + f * pa.len), uOn(best.arc, best.t0 + (r2 ? 1 - f : f) * best.len)); return c; };
        const cf = cst(false), cr = cst(true);
        if (rev ? cf < 0.5 * cr : cr < 0.5 * cf) rev = !rev;
      }
      if (best) { used.add(best); pairs.push([pa, best, rev]); } else unA.push(pa);
    }
    const unB = PB.filter(pb => !used.has(pb));
    const spA = A.sp.filter(q => !keep.has(q.chain)), spB = B.sp.filter(q => !keep.has(q.chain));
    // 天空间距: 对应段上 9 个比例点的最大夹角
    let gap = 0;
    for (const [pa, pb, rev] of pairs) for (let k = 0; k <= 8; k++) {
      const f = k / 8, fb = rev ? 1 - f : f;
      gap = Math.max(gap, angU(uOn(pa.arc, pa.t0 + f * pa.len), uOn(pb.arc, pb.t0 + fb * pb.len)));
    }
    // 边界线在两层之间的弦高 (同一条线上相连的两个末端点之间)
    let lineDev = 0;
    for (const qa of A.sp) {
      if (!qa.lab || !keep.has(qa.chain)) continue;
      const qb = B.sp.find(q => q.chain === qa.chain), sh = qb && sharedLine(qa.lab, qb.lab);
      if (!sh) continue;
      const Q = linePts(sh.id, sh.sa, sh.sb);
      if (lineOk(Q, qa.p, qb.p, A.r, B.r)) lineDev = Math.max(lineDev, chordDev(qa.p, qb.p, Q));
    }
    const event = unA.length > 0 || unB.length > 0 || spA.length > 0 || spB.length > 0;
    return { pairs, unA, unB, spA, spB, gap, lineDev, event, keep };
  };

  // ---- 1) 初始层
  const crit = [];
  {
    const cr = opt.critical || (LN.length ? criticalRadii(LN.map(l => l.p), L, E) : []);
    for (const rc of cr.slice().sort((p, q) => p - q)) {
      if (!(rc * (1 - 2 * CRIT) > rMin && rc * (1 + 2 * CRIT) < rMax) || Math.abs(rc / D - 1) < 2 * PINL) continue;
      if (crit.length && rc / crit[crit.length - 1] - 1 < 4 * CRIT) continue;
      crit.push(rc);
    }
  }
  const rs = [];
  for (let k = 0; k < N0; k++) {
    const r = rMin * Math.pow(rMax / rMin, k / Math.max(1, N0 - 1));
    if (Math.abs(r / D - 1) > PINL && !crit.some(rc => Math.abs(r / rc - 1) < 2 * CRIT)) rs.push([r, 'init']);
  }
  if (rMin < D * (1 - PINL) && rMax > D * (1 + PINL)) rs.push([D * (1 - PINL), 'lampLo'], [D * (1 + PINL), 'lampHi']);
  for (const rc of crit) rs.push([rc * (1 - CRIT), 'critLo'], [rc * (1 + CRIT), 'critHi']);
  rs.sort((p, q) => p[0] - q[0]);
  let rows = [];
  for (const [r, kind] of rs) rows.push(buildRow(r, rows[rows.length - 1], null, kind));
  if (log) log(`initial ${rows.length} shells, ${crit.length} critical radii, ${LN.length} lines`);
  // 顶点预算 (opt.budget): 铺满大半个天空的晕面每层的线有几百度长. 按初始层估一下顶点数 (加密大约再多一半),
  //   超了就把点距按比例放大 (至多 2°); 加密过程中超了预算就不再为天空间距 / 弦高加层, 只为拓扑变化二分
  const vEst = rows => { let v = 0; for (const rw of rows) for (const a of rw.arcs) v += Math.max(a.closed ? 8 : 4, Math.ceil(a.len / SKY)) + a.sp.length; return v; };
  {
    const est = 1.5 * vEst(rows) * Math.max(1, N0 / Math.max(1, rows.length));
    if (est > VBUDGET) { SKY = Math.min(2 * DEG, SKY * est / VBUDGET); GAPSKY = Math.max(1.2 * DEG, 3 * SKY); MAXEDGE = Math.max(4 * GAPSKY, 5 * DEG); }
  }

  // ---- 5) 加密
  const pairSkip = (A, B) => A.kind === 'lampLo' && B.kind === 'lampHi' || A.kind === 'critLo' && B.kind === 'critHi' && B.r / A.r - 1 < 3 * CRIT;
  const nearLamp = (A, B) => {
    const dl = row => { if (row.rl === undefined) { let m = Infinity; for (const arc of row.arcs) for (let i = 0; i < arc.n; i += Math.max(1, arc.n >> 6)) { const p = at(row.r, [arc.U[3 * i], arc.U[3 * i + 1], arc.U[3 * i + 2]]); m = Math.min(m, Math.hypot(p[0] - L[0], p[1] - L[1], p[2] - L[2])); } row.rl = m; } return row.rl; };
    return B.r - A.r > 0.35 * Math.max(Math.min(dl(A), dl(B)), 2e-3 * D);
  };
  for (let pass = 0; pass < 28; pass++) {
    relink(rows);
    const ins = new Map(), over = vEst(rows) > 2.5 * VBUDGET;
    const AN = [];
    for (let k = 0; k + 1 < rows.length; k++) {
      const A = rows[k], B = rows[k + 1];
      if (!A.arcs.length && !B.arcs.length && !A.sp.length && !B.sp.length) { AN.push(null); continue; }
      const an = analyze(A, B);
      AN.push(an);
      if (pairSkip(A, B)) continue;
      const ratio = B.r / A.r - 1, mid = Math.sqrt(A.r * B.r);
      if (an.event) { if (ratio > EPSR) ins.set(mid, [A, B]); continue; }
      if (ratio <= 4 * EPSR) continue;
      if (over) continue;
      if (an.gap > GAPSKY || nearLamp(A, B) || an.lineDev > tol3(B.r)) ins.set(mid, [A, B]);
    }
    // 径向弯曲: 三层 (A, B, C) 上对应段的同一比例点不在一条直线上 (超了预算不做) (按 ln r 插值) → 两侧各加一层
    for (let k = 1; k + 1 < rows.length && !over; k++) {
      const a1 = AN[k - 1], a2 = AN[k];
      if (!a1 || !a2 || a1.event || a2.event) continue;
      const A = rows[k - 1], B = rows[k], C = rows[k + 1];
      if (pairSkip(A, B) || pairSkip(B, C) || B.r / A.r - 1 <= 4 * EPSR || C.r / B.r - 1 <= 4 * EPSR) continue;
      const lam = Math.log(B.r / A.r) / Math.log(C.r / A.r);
      let dev = 0;
      for (const [pa, pb, rv1] of a1.pairs) {
        if (pb.loop) continue;
        const nx = a2.pairs.find(q => q[0].arc === pb.arc && q[0].k0 === pb.k0 && q[0].k1 === pb.k1);
        if (!nx) continue;
        const [, pc, rv2] = nx;
        for (const f of [0.25, 0.5, 0.75]) {
          const fa = rv1 ? 1 - f : f, fc = rv2 ? 1 - f : f;
          const PA = at(A.r, uOn(pa.arc, pa.t0 + fa * pa.len)), PB = at(B.r, uOn(pb.arc, pb.t0 + f * pb.len)), PC = at(C.r, uOn(pc.arc, pc.t0 + fc * pc.len));
          dev = Math.max(dev, dist3(PB, [PA[0] + lam * (PC[0] - PA[0]), PA[1] + lam * (PC[1] - PA[1]), PA[2] + lam * (PC[2] - PA[2])]));
        }
      }
      if (dev > 2 * tol3(B.r)) { ins.set(Math.sqrt(A.r * B.r), [A, B]); ins.set(Math.sqrt(B.r * C.r), [B, C]); }
    }
    if (!ins.size || calls >= MAXCALLS || rows.length >= MAXROWS) break;
    const add = [];
    for (const [r, [A, B]] of [...ins.entries()].sort((p, q) => p[0] - q[0])) {
      if (calls >= MAXCALLS || rows.length + add.length >= MAXROWS) break;
      add.push(buildRow(r, A, B, 'ref'));
    }
    rows = rows.concat(add).sort((p, q) => p.r - q.r);
    if (log) log(`pass ${pass}: +${add.length} shells (${rows.length})`);
  }
  relink(rows);
  if (opt.debugFinal) opt.debugFinal(rows, analyze);      // 测试 / 排查用: 看最终各层的特殊点与对应

  // ---- 4) / 6) 三角形
  // 每层的顶点: 每条弧上 = 全部特殊点 + 相邻两个特殊点之间按点距均匀取; 上下两条带共用
  const vertsOf = (row, arc) => {
    if (arc.vt) return arc.vt;
    const ts = arc.sp.filter(q => q.type !== 'focal').map(q => q.t);
    if (!arc.closed) { ts.push(0, arc.S[arc.n - 1]); }
    const T = [...new Set(ts.map(t => +t.toFixed(12)))].sort((a, b) => a - b);
    const all = [];
    const spans = arc.closed ? (T.length ? T.map((t, i) => [t, i + 1 < T.length ? T[i + 1] : T[0] + arc.len]) : [[0, arc.len]]) : T.slice(0, -1).map((t, i) => [t, T[i + 1]]);
    // 每段按点距取, 但小弧也要保住形状: 闭环全周至少 8 段, 开弧至少 4 段 (比点距还短的小透镜原来只剩一个点, 铺不出三角形)
    const minSeg = arc.closed ? 8 : 4;
    for (const [a, b] of spans) {
      const m = Math.max(1, Math.ceil((b - a) / SKY), arc.len > 0 ? Math.ceil(minSeg * (b - a) / arc.len) : 1);
      for (let k = 0; k < m; k++) all.push(a + (b - a) * k / m);
    }
    if (!arc.closed) all.push(T[T.length - 1]);
    const vt = { t: Float64Array.from(all), id: new Int32Array(all.length) };
    for (let i = 0; i < all.length; i++) { const u = uOn(arc, all[i]); vt.id[i] = addV(at(row.r, u), u); }
    arc.vt = vt;
    return vt;
  };
  // 段 → 顶点序号串 (按 t0 → t1) 与各点在段内的比例
  const pieceVerts = (row, pc) => {
    const vt = vertsOf(row, pc.arc), n = vt.t.length, ids = [], fr = [];
    const L0 = pc.arc.len, span = Math.max(1e-15, pc.t1 - pc.t0);
    if (pc.loop) {
      for (let i = 0; i < n; i++) { ids.push(vt.id[i]); fr.push((vt.t[i] - vt.t[0]) / L0); }
      ids.push(vt.id[0]); fr.push(1);
      return { ids, fr };
    }
    const eps = 1e-10;
    const push = (i, t) => { ids.push(vt.id[i]); fr.push((t - pc.t0) / span); };
    if (!pc.arc.closed) {
      for (let i = 0; i < n; i++) if (vt.t[i] >= pc.t0 - eps && vt.t[i] <= pc.t1 + eps) push(i, vt.t[i]);
    } else {                                           // 闭环: 从 t0 起绕着走到 t1 (可越过起点, 整周时回到起点)
      let i0 = -1;
      for (let i = 0; i < n; i++) if (vt.t[i] >= pc.t0 - eps) { i0 = i; break; }
      let off = 0;
      if (i0 < 0) { i0 = 0; off = L0; }
      for (let k = 0; k <= n; k++) {
        const i = (i0 + k) % n, tt = vt.t[i] + off + (i0 + k >= n ? L0 : 0);
        if (tt > pc.t1 + eps) break;
        push(i, tt);
      }
    }
    return { ids, fr };
  };
  // 两串顶点按比例拉链
  const zipFr = (A, B, i0 = 0, i1 = A.ids.length - 1, j0 = 0, j1 = B.ids.length - 1) => {
    let i = i0, j = j0;
    while (i < i1 || j < j1) {
      if (j >= j1 || (i < i1 && A.fr[i + 1] <= B.fr[j + 1])) { emit(A.ids[i], A.ids[i + 1], B.ids[j]); i++; }
      else { emit(A.ids[i], B.ids[j + 1], B.ids[j]); j++; }
    }
  };
  const revV = V => ({ ids: V.ids.slice().reverse(), fr: V.fr.slice().reverse().map(f => 1 - f) });
  // 段端的特殊点在两层上落在同一条三维边界线上: 两层之间沿这条线插入它的点 (边界线就是网格的边)
  const sideLine = (A, B, qa, qb) => {
    const sh = qa && qb ? sharedLine(qa.lab, qb.lab) : null;
    if (!sh) return null;
    const Q = linePts(sh.id, sh.sa, sh.sb);
    if (!Q.length || !lineOk(Q, qa.p, qb.p, A.r, B.r)) return null;
    return Q.map(p => { stat.lineVerts++; return addV(p, unit([p[0] - E[0], p[1] - E[1], p[2] - E[2]])); });
  };
  const spAt = (row, pc, end) => {                      // 段一端的特殊点 (按链号)
    const k = end ? pc.k1 : pc.k0;
    return typeof k === 'number' && k >= 0 ? row.sp.find(q => q.chain === k && q.arc === pc.arc) : null;
  };
  const loopAlign = (VA, VB) => {             // 闭环: VB 转到离 VA 起点最近的点起, 两个走向里取四分点对得近的; 比例按累计天空弧长
    const ring = VB.ids.slice(0, -1), m = ring.length;
    let j0 = 0, bd = Infinity;
    for (let j = 0; j < m; j++) { const d = vang(VA.ids[0], ring[j]); if (d < bd) { bd = d; j0 = j; } }
    const mk = sgn => { const ids = []; for (let k = 0; k <= m; k++) ids.push(ring[((j0 + sgn * k) % m + m) % m]); return ids; };
    const frOf = ids => { const c = [0]; for (let k = 1; k < ids.length; k++) c.push(c[k - 1] + vang(ids[k - 1], ids[k])); const t = c[c.length - 1] || 1; return c.map(x => x / t); };
    const fa = frOf(VA.ids), A2 = { ids: VA.ids, fr: fa };
    const at4 = (V, f) => { let k = 0; while (k + 1 < V.fr.length && V.fr[k + 1] < f) k++; return V.ids[k]; };
    const cand = [1, -1].map(sg => { const ids = mk(sg); const V = { ids, fr: frOf(ids) }; let c = 0; for (const f of [0.25, 0.5, 0.75]) c += vang(at4(A2, f), at4(V, f)); return [c, V]; });
    return [A2, cand[0][0] <= cand[1][0] ? cand[0][1] : cand[1][1]];
  };
  const strip = (A, B, pa, pb, rev) => {
    let VA = pieceVerts(A, pa), VB = pieceVerts(B, pb);
    if (pa.loop && pb.loop) { if (VA.ids.length > 2 && VB.ids.length > 2) { [VA, VB] = loopAlign(VA, VB); zipFr(VA, VB); } return; }
    if (rev) VB = revV(VB);
    if (VA.ids.length < 2 || VB.ids.length < 2) {                    // 缩成一点的段: 扇形
      const P = VA.ids.length < 2 ? VB : VA, c = VA.ids.length < 2 ? VA.ids[0] : VB.ids[0];
      if (c === undefined) return;
      for (let i = 0; i + 1 < P.ids.length; i++) emit(P.ids[i], P.ids[i + 1], c);
      return;
    }
    let i0 = 0, j0 = 0, i1 = VA.ids.length - 1, j1 = VB.ids.length - 1;
    // 起点一侧
    const s0 = VA.ids.length > 2 && VB.ids.length > 2 ? sideLine(A, B, spAt(A, pa, false), spAt(B, pb, rev)) : null;
    if (s0) {
      const S = [VA.ids[0], ...s0, VB.ids[0]], m = S.length >> 1;
      for (let k = 0; k + 1 < S.length; k++) emit(S[k], S[k + 1], k < m ? VA.ids[1] : VB.ids[1]);
      emit(VA.ids[1], S[m], VB.ids[1]);
      i0 = 1; j0 = 1;
    }
    const s1 = VA.ids.length > 2 && VB.ids.length > 2 && i1 - i0 >= 1 && j1 - j0 >= 1 ? sideLine(A, B, spAt(A, pa, true), spAt(B, pb, !rev)) : null;
    if (s1) {
      const S = [VA.ids[i1], ...s1, VB.ids[j1]], m = S.length >> 1;
      for (let k = 0; k + 1 < S.length; k++) emit(S[k], S[k + 1], k < m ? VA.ids[i1 - 1] : VB.ids[j1 - 1]);
      emit(VA.ids[i1 - 1], S[m], VB.ids[j1 - 1]);
      i1--; j1--;
    }
    zipFr(VA, VB, i0, i1, j0, j1);
  };
  let lampV = -1;
  const capTo = (row, pc, toLamp) => {
    const V = pieceVerts(row, pc);
    if (!V.ids.length) return;
    let c;
    if (toLamp) { if (lampV < 0) lampV = addV(L, [Dx / D, Dy / D, Dz / D]); c = lampV; }
    else {
      let x = 0, y = 0, z = 0;
      for (const v of V.ids) { x += VU[3 * v]; y += VU[3 * v + 1]; z += VU[3 * v + 2]; }
      const u = unit([x, y, z]); c = addV(at(row.r, u), u);
    }
    for (let i = 0; i + 1 < V.ids.length; i++) emit(V.ids[i], V.ids[i + 1], c);
    stat.caps++;
  };
  // 对不上的段 (相邻两层几乎重合时): 每个点找对面对不上的段上视方向最近的点 (LIM 之内), 同一段上脚点单调的一串一串拉链
  const LIM = 2.5 * GAPSKY;
  const looseSew = (A, B, UA, UB) => {
    const VA = UA.map(pc => pieceVerts(A, pc)), VB = UB.map(pc => pieceVerts(B, pc));
    let unc = 0;
    const pass = (S, T, flip) => {
      S.forEach(w => {
        const ft = w.ids.map(v => { let bq = -1, bj = -1, bd = LIM; T.forEach((x, q) => { for (let j = 0; j < x.ids.length; j++) { const d = vang(v, x.ids[j]); if (d < bd) { bd = d; bq = q; bj = j; } } }); return bq < 0 ? null : [bq, bj]; });
        for (let s0 = 0; s0 < ft.length;) {
          if (!ft[s0]) { if (s0 + 1 < w.ids.length) unc += vang(w.ids[s0], w.ids[s0 + 1]); s0++; continue; }
          let e = s0, dir = 0;
          while (e + 1 < ft.length && ft[e + 1] && ft[e + 1][0] === ft[s0][0]) {
            const sg = Math.sign(ft[e + 1][1] - ft[e][1]);
            if ((sg && dir && sg !== dir) || Math.abs(ft[e + 1][1] - ft[e][1]) > 12) break;
            if (sg) dir = sg; e++;
          }
          if (e > s0) {
            const x = T[ft[s0][0]], j0 = ft[s0][1], j1 = ft[e][1], st = j1 >= j0 ? 1 : -1, Q = [];
            for (let j = j0; j !== j1 + st; j += st) Q.push(x.ids[j]);
            const P = w.ids.slice(s0, e + 1);
            const a = { ids: flip ? Q : P, fr: (flip ? Q : P).map((_, k, arr) => k / Math.max(1, arr.length - 1)) };
            const b = { ids: flip ? P : Q, fr: (flip ? P : Q).map((_, k, arr) => k / Math.max(1, arr.length - 1)) };
            zipFr(a, b);
          }
          s0 = Math.max(e, s0 + 1);
        }
      });
    };
    if (VA.length && VB.length) { pass(VA, VB, false); pass(VB, VA, true); }
    else for (const w of VA.length ? VA : VB) for (let i = 0; i + 1 < w.ids.length; i++) unc += vang(w.ids[i], w.ids[i + 1]);
    return unc;
  };
  for (let k = 0; k + 1 < rows.length; k++) {
    const A = rows[k], B = rows[k + 1];
    if (!A.arcs.length && !B.arcs.length) continue;
    const an = analyze(A, B);
    for (const [pa, pb, rev] of an.pairs) strip(A, B, pa, pb, rev);
    if (!an.unA.length && !an.unB.length) continue;
    // 灯两侧那一对层: 晕面若过灯, 灯旁的段缩向灯的方向 (锥顶) —— 对不上的、靠近灯方向 (20° 内) 的段收到灯上
    if (A.kind === 'lampLo' && B.kind === 'lampHi') {
      const ul = [Dx / D, Dy / D, Dz / D], UA = [], UB = [];
      for (const [row, list, out] of [[A, an.unA, UA], [B, an.unB, UB]]) for (const pc of list) {
        const V = pieceVerts(row, pc);
        if (V.ids.every(v => angU(vu(v), ul) < 20 * DEG)) capTo(row, pc, true); else out.push(pc);
      }
      if (UA.length || UB.length) { const unc = looseSew(A, B, UA, UB); if (unc > 1.5 * SKY) { stat.loose++; stat.looseDeg += unc / DEG; } }
      continue;
    }
    stat.events++;
    if (log) log(`event r ${A.r.toFixed(5)} .. ${B.r.toFixed(5)}: pieces ${an.unA.length} | ${an.unB.length} unmatched`);
    // 小段 (正在缩没 / 刚出生) 收口, 其余就近缝
    const tiny = pc => pc.len < Math.min(CAPMAX, 30 * Math.sqrt(2 * Math.max(EPSR, B.r / A.r - 1)));
    const UA = [], UB = [];
    for (const pc of an.unA) { if (tiny(pc)) capTo(A, pc, false); else UA.push(pc); }
    for (const pc of an.unB) { if (tiny(pc)) capTo(B, pc, false); else UB.push(pc); }
    if (UA.length || UB.length) {
      const unc = looseSew(A, B, UA, UB);
      if (unc > 1.5 * SKY) { stat.loose++; stat.looseDeg += unc / DEG; }
    }
  }

  // ---- 特殊线输出 (折叠线、焦线): 各链按半径排的点
  const chains = new Map();
  for (const row of rows) for (const q of row.sp) {
    if (q.type !== 'fold' && q.type !== 'focal') continue;
    let c = chains.get(q.chain);
    if (!c) { c = { kind: q.type, pts: [], dirs: [] }; chains.set(q.chain, c); }
    c.pts.push(...q.p); c.dirs.push(...q.u);
  }
  const special = [];
  for (const c of chains.values()) if (c.pts.length >= 6) { special.push({ kind: c.kind, p: Float32Array.from(c.pts), d: Float32Array.from(c.dirs) }); if (c.kind === 'fold') stat.folds++; else stat.focal++; }

  // ---- 输出: 每个三角形三个顶点依次排; 法向 = 共顶点且朝向相近 (夹角 < 60°) 的面法向之和 (至多看 48 个)
  const nT = tris.length / 3, pos = new Float32Array(9 * nT), dir = new Float32Array(9 * nT), hh = new Float32Array(3 * nT), nrm = new Float32Array(9 * nT);
  const TN = new Float64Array(3 * nT), cnt = new Uint32Array(nv + 1);
  for (let t = 0; t < nT; t++) {
    const a = tris[3 * t], b = tris[3 * t + 1], c = tris[3 * t + 2];
    const ux = VP[3 * b] - VP[3 * a], uy = VP[3 * b + 1] - VP[3 * a + 1], uz = VP[3 * b + 2] - VP[3 * a + 2], vx = VP[3 * c] - VP[3 * a], vy = VP[3 * c + 1] - VP[3 * a + 1], vz = VP[3 * c + 2] - VP[3 * a + 2];
    TN[3 * t] = uy * vz - uz * vy; TN[3 * t + 1] = uz * vx - ux * vz; TN[3 * t + 2] = ux * vy - uy * vx;
    cnt[a + 1]++; cnt[b + 1]++; cnt[c + 1]++;
  }
  for (let i = 0; i < nv; i++) cnt[i + 1] += cnt[i];
  const off = cnt.slice(), inc = new Uint32Array(3 * nT);
  for (let i = 0; i < tris.length; i++) inc[off[tris[i]]++] = (i / 3) | 0;
  for (let t = 0; t < nT; t++) {
    const tl = Math.hypot(TN[3 * t], TN[3 * t + 1], TN[3 * t + 2]);
    for (let k = 0; k < 3; k++) {
      const v = tris[3 * t + k], o = 9 * t + 3 * k;
      let x = 0, y = 0, z = 0;
      const nInc = cnt[v + 1] - cnt[v], stp = nInc > 48 ? nInc / 48 : 1;
      if (stp > 1) { x = TN[3 * t]; y = TN[3 * t + 1]; z = TN[3 * t + 2]; }
      for (let q = 0; q * stp < nInc; q++) {
        const j = inc[cnt[v] + Math.floor(q * stp)], c = TN[3 * j] * TN[3 * t] + TN[3 * j + 1] * TN[3 * t + 1] + TN[3 * j + 2] * TN[3 * t + 2];
        if (j !== t && Math.abs(c) < 0.5 * Math.hypot(TN[3 * j], TN[3 * j + 1], TN[3 * j + 2]) * tl) continue;
        const sg = c < 0 ? -1 : 1;
        x += sg * TN[3 * j]; y += sg * TN[3 * j + 1]; z += sg * TN[3 * j + 2];
      }
      const nl = Math.hypot(x, y, z) || 1;
      pos[o] = VP[3 * v]; pos[o + 1] = VP[3 * v + 1]; pos[o + 2] = VP[3 * v + 2];
      dir[o] = VU[3 * v]; dir[o + 1] = VU[3 * v + 1]; dir[o + 2] = VU[3 * v + 2];
      nrm[o] = x / nl; nrm[o + 1] = y / nl; nrm[o + 2] = z / nl;
      hh[3 * t + k] = Math.asin((L[2] - VP[3 * v + 2]) / (Math.hypot(VP[3 * v] - L[0], VP[3 * v + 1] - L[1], VP[3 * v + 2] - L[2]) || 1)) / DEG;
    }
  }
  // 各层的线 (画层线用). gap = 与里面一层在天空里隔多远 (对应段的最大夹角), 画的时候按累计的 gap 挑着画
  const lines = [];
  let prev = null;
  for (const row of rows) {
    if (!row.raw.length) continue;
    let gap = Infinity;
    if (prev) { const an = analyze(prev, row); gap = an.pairs.length ? an.gap : GAPSKY; }
    lines.push({ r: row.r, gap, curves: row.raw });
    prev = row;
  }
  stat.rows = rows.filter(r => r.arcs.length).length; stat.calls = calls; stat.ms = Date.now() - t0; stat.skipped = nSkip;
  return { pos, dir, h: hh, nrm, nTri: nT, lines, special, stat, rowR: rows.filter(r => r.arcs.length).map(r => r.r) };
}
