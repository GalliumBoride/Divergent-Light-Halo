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
  // 交界层 (这一段的最里、最外一层, 相邻段各算一次): 顶点按传入的点距 SKY0 取 —— 各段传进来的一样; 段内的点距可以按
  //   这一段的预算放大, 交界层不跟着变, 两段在交界层上的顶点才逐位相同 (否则两侧的三角形边对不上, 屏幕上是一圈细缝)
  const SKY0 = SKY;
  let GAPSKY = Math.max(1.2 * DEG, 3 * SKY), MAXEDGE = Math.max(4 * GAPSKY, 5 * DEG);
  const CAPMAX = 1 * DEG, VBUDGET = opt.budget || 90000;
  const EPSR = opt.epsR || 5e-4, CRIT = 2e-4, PINL = opt.pinLamp || 1e-3;
  // 层数上限随初始层数 (一段只有七八层初始层时不许加到几百层: 那是相邻层反复对不上在空转, 每层都是整圈轮廓,
  //   三角形与用时一起爆掉); softMs: 加密用掉这么久就收手, 用已有的层出网格 (不整段作废)
  const MAXROWS = opt.maxRows || Math.max(60, Math.min(700, 12 * N0 + 80)), MAXCALLS = opt.maxCalls || Math.round(1.3 * MAXROWS);
  const SOFTMS = opt.softMs || Infinity;
  const FOLDW = Math.max(0.05 * DEG, 0.25 * SKY);          // 折点检测: 前后各沿线走这么长 (天空弧长) 比方向 (最小尺度)
  const CORNW = Math.max(FOLDW, SKY), FOLDC = Math.cos(100 * DEG);   // 最大尺度; 两段弦夹角超过 120° (最小尺度) / 100° (更大尺度) 算折回
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
    // 闭合的线 (首尾同一点; 如极点线): 线上两点之间走短的那一边, 可以越过起止点
    const closed = n > 3 && Math.hypot(p[0] - p[3 * n - 3], p[1] - p[3 * n - 2], p[2] - p[3 * n - 1]) < 1e-5;
    LN.push({ id: LN.length, kind: l.kind || 'edge', p, n, s, lo, hi, tot: s[n - 1], closed });
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
    for (const a of la) for (const b of lb) if (a.id === b.id) {
      let ds = Math.abs(a.s - b.s);
      if (LN[a.id].closed && ds > 0.5 * LN[a.id].tot) ds = LN[a.id].tot - ds;
      if (!best || ds < best.ds) best = { id: a.id, sa: a.s, sb: b.s, ds };
    }
    return best;
  };
  // 线 id 上弧长 s0 → s1 之间 (不含两端) 的点, 按 s0 → s1 的顺序 (闭合的线走短的那一边)
  const linePts = (id, s0, s1) => {
    const ln = LN[id], out = [], a = Math.min(s0, s1), b = Math.max(s0, s1);
    const pt = i => [ln.p[3 * i], ln.p[3 * i + 1], ln.p[3 * i + 2]];
    if (ln.closed && b - a > 0.5 * ln.tot) {
      // 越过起止点 (第 0 点与最后一点是同一个点, 只取一次)
      if (s0 < s1) { for (let i = ln.n - 2; i >= 0; i--) if (ln.s[i] < s0 - 1e-9) out.push(pt(i)); for (let i = ln.n - 2; i >= 1; i--) if (ln.s[i] > s1 + 1e-9) out.push(pt(i)); }
      else { for (let i = 1; i < ln.n; i++) if (ln.s[i] > s0 + 1e-9) out.push(pt(i)); for (let i = 1; i < ln.n - 1; i++) if (ln.s[i] < s1 - 1e-9) out.push(pt(i)); }
      return out;
    }
    for (let i = 0; i < ln.n; i++) if (ln.s[i] > a + 1e-9 && ln.s[i] < b - 1e-9) out.push(pt(i));
    if (s0 > s1) out.reverse();
    return out;
  };
  // 这些点离眼的距离都夹在两层半径之间吗 (tol: 相对放宽量). 两层各有一个末端点落在同一条边界线上、两点之间的那段线又没出过
  //   这两层的半径范围 ⟹ 那段线就是晕面在这两层之间的边
  const inBand = (Q, r0, r1, tol) => {
    const lo = Math.min(r0, r1) * (1 - tol), hi = Math.max(r0, r1) * (1 + tol);
    for (const q of Q) { const rq = Math.hypot(q[0] - E[0], q[1] - E[1], q[2] - E[2]); if (rq < lo || rq > hi) return false; }
    return true;
  };
  const BAND = 2e-5;                       // 远小于临界半径两侧那对层的间距 (2·CRIT): 线的极值点不算 "夹在" 相邻的下一对层之间
  // q0 → q1 沿边界线走的点 (不含两端; 可空): 同一条线; 或两条在端点相接的线 (晕面的角点). 走不通 (不在这两层之间) 返回 null
  const linePath = (q0, q1, r0, r1) => {
    if (!q0.lab || !q1.lab) return null;
    const sh = sharedLine(q0.lab, q1.lab);
    if (sh) { const Q = linePts(sh.id, sh.sa, sh.sb); if (inBand(Q, r0, r1, BAND)) return Q; }
    for (const a of q0.lab) for (const b of q1.lab) {
      if (a.id === b.id) continue;
      const la = LN[a.id], lb = LN[b.id];
      if (la.closed || lb.closed) continue;
      for (const ea of [0, la.n - 1]) for (const eb of [0, lb.n - 1]) {
        const ca = [la.p[3 * ea], la.p[3 * ea + 1], la.p[3 * ea + 2]], cb = [lb.p[3 * eb], lb.p[3 * eb + 1], lb.p[3 * eb + 2]];
        if (dist3(ca, cb) > 1e-3 * Math.max(r0, r1)) continue;
        const Q = [...linePts(a.id, a.s, la.s[ea]), ca, ...linePts(b.id, lb.s[eb], b.s)];
        if (inBand(Q, r0, r1, BAND)) return Q;
      }
    }
    return null;
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
  // VF: 顶点的标记 (自检用): 1 = 可以在晕面边界上 (弧的末端点 / 截断点, 边界线上的点, 灯), 2 / 4 = 在这一段最里 / 最外一层上
  let VP = new Float64Array(3 * 65536), VU = new Float64Array(3 * 65536), VF = new Uint8Array(65536), nv = 0;
  const addV = (p, u, f = 0) => {
    if (3 * nv + 3 > VP.length) { const a = new Float64Array(2 * VP.length); a.set(VP); VP = a; const b = new Float64Array(2 * VU.length); b.set(VU); VU = b; const c = new Uint8Array(2 * VF.length); c.set(VF); VF = c; }
    VP[3 * nv] = p[0]; VP[3 * nv + 1] = p[1]; VP[3 * nv + 2] = p[2]; VU[3 * nv] = u[0]; VU[3 * nv + 1] = u[1]; VU[3 * nv + 2] = u[2];
    VF[nv] = f;
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
  let calls = 0, cutSeq = 0;
  const seedsOf = (dense, ...rows) => {
    const runs = [];
    for (const rw of rows) if (rw && rw.states) for (const run of rw.states) if (run && run.length) runs.push(run);
    if (!runs.length) return null;
    const per = dense ? Math.max(8, Math.min(120, Math.floor(600 / runs.length))) : Math.max(3, Math.min(16, Math.floor(80 / runs.length))), out = [];
    for (const run of runs) { const st = Math.max(1, Math.floor(run.length / per)); for (let i = st >> 1; i < run.length; i += st) out.push(run[i]); out.push(run[0], run[run.length - 1]); }
    return out;
  };
  // 弧: U = 单位视方向 (3n), S = 累计天空弧长; 折点检测
  const makeArc = (row, U, n, closed, headType, tailType, joints = [], cuts = null) => {
    const S = new Float64Array(n);
    for (let i = 1; i < n; i++) S[i] = S[i - 1] + angU([U[3 * i - 3], U[3 * i - 2], U[3 * i - 1]], [U[3 * i], U[3 * i + 1], U[3 * i + 2]]);
    const len = S[n - 1] + (closed && n > 2 ? angU([U[3 * n - 3], U[3 * n - 2], U[3 * n - 1]], [U[0], U[1], U[2]]) : 0);
    const arc = { U, n, S, len, closed, sp: [], ai: row.arcs.length };
    const uAt = i => [U[3 * i], U[3 * i + 1], U[3 * i + 2]];
    const mk = (t, type, u) => ({ t, type, u, p: at(row.r, u), arc, lab: null, chain: -1 });
    if (!closed) {
      const h = mk(0, headType, uAt(0)), t = mk(S[n - 1], tailType, uAt(n - 1));
      // 截断点的配对号: 同一条轮廓上被截断球面切掉的一段, 它两头的两个截断点同号 (joinPieces 只把同号的接上)
      if (cuts) { if (headType === 'cut') h.cutId = cuts[0]; if (tailType === 'cut') t.cutId = cuts[1]; }
      arc.sp.push(h); arc.sp.push(t);
    }
    // 接头 (引擎分段给的同一条轮廓并接处; 退化弧两侧的两支在焦点处相接 —— 两张面交于焦线): 也是特殊点, 段在这里切开
    for (const j of joints) if ((j > 0 && j < n - 1) || closed) arc.sp.push(mk(S[Math.min(j, n - 1)], 'joint', uAt(Math.min(j, n - 1))));
    // 折点 (轮廓在天空里折回 —— 曲面的折痕): 前后各沿线走 W 的两段弦夹角 > 120° (W = FOLDW) 或 > 100° (W = 2·FOLDW … 到一个天空步长)
    //   (多尺度: 眼旁折痕尖到点距以下, 离眼远了变成圆角 —— 圆角半径还小于天空步长时网格照样需要一个顶点在尖上);
    //   连着几个点都满足的取最尖的那个 (沿 "来向 − 去向" 投影最远)
    if (n >= 5 && len > 2 * FOLDW) {
      const Sx = k => closed ? (k < 0 ? S[k + n] - len : k >= n ? S[k - n] + len : S[k]) : S[k];
      const Ux = k => uAt(closed ? ((k % n) + n) % n : k);
      const lo = closed ? 0 : 1, hi = closed ? n : n - 1;
      const WS = [FOLDW];
      // 大尺度只用在比窗口长得多的弧上 (闭环 16 倍、开弧 8 倍): 灯旁一两度的小圈上, 大窗口处处 "折回", 认出的折点没有意义
      while (WS.length < 4 && 2 * WS[WS.length - 1] <= CORNW * 1.0001 && 2 * WS[WS.length - 1] <= (closed ? 1 / 16 : 1 / 8) * len) WS.push(2 * WS[WS.length - 1]);
      const hs = WS.map(() => lo - 1), js = WS.map(() => lo + 1);
      // 先逐点算: 触发的最小尺度 KT (−1 = 不触发)、该尺度的余弦与尖的方向
      const KT = new Int8Array(n).fill(-1), CS = new Float64Array(n), DIR = new Array(n);
      for (let i = lo; i < hi; i++) {
        if (!closed && (S[i] < FOLDW || S[n - 1] - S[i] < FOLDW)) continue;
        const ui = Ux(i);
        for (let w = 0; w < WS.length; w++) {
          const W = WS[w];
          if (!closed && (S[i] < W || S[n - 1] - S[i] < W)) break;
          let h = hs[w], j = js[w];
          while (Sx(i) - Sx(h) > W && h < i - 1) h++;
          while (Sx(i) - Sx(h) < W && h > i - n) h--;
          if (j <= i) j = i + 1;
          while (Sx(j) - Sx(i) < W && j < i + n - 1 && (closed || j < n - 1)) j++;
          hs[w] = h; js[w] = j;
          const uh = Ux(h), uj = Ux(j);
          const a = [ui[0] - uh[0], ui[1] - uh[1], ui[2] - uh[2]], b = [uj[0] - ui[0], uj[1] - ui[1], uj[2] - ui[2]];
          const la = Math.hypot(a[0], a[1], a[2]), lb = Math.hypot(b[0], b[1], b[2]);
          if (la < 1e-12 || lb < 1e-12) continue;
          const cs = dot(a, b) / (la * lb);
          // 触发的最小尺度定尖的方向: 真尖点 (小窗口就折回) 用小窗口的弦, 尖取得准 (大窗口的弦带着两支的弯曲, 尖会偏);
          //   只有大窗口才折回的圆角才用大窗口
          if (cs < (w ? FOLDC : -0.5)) { KT[i] = w; CS[i] = cs; DIR[i] = unit([a[0] / la - b[0] / lb, a[1] / la - b[1] / lb, a[2] / la - b[2] / lb]); break; }
        }
      }
      // 再把连着触发的点聚成簇, 每簇一个折点. 闭环从一个不触发的点开始绕一圈 (跨过起点的簇不被拆成两半 ——
      //   拆开就在同一个折角两侧各出一个折点, 两者间距在去重阈值上下摆动, 相邻层时有时无); 整圈都触发 (比窗口还小的圈) 不算
      let cl = null;
      // 一簇 → 折点. 簇里在最小触发尺度上触发的点若分成几串 (相隔不到一个大窗口的两个发夹弯, 大窗口下连成一簇), 每串一个折点
      //   —— 否则灯旁挨得近的两条折叠线在某些层上并成一个折点, 线就断成几截
      const flush = () => {
        if (!cl) return;
        const subs = [];
        let cur = null;
        for (let k = cl.i0; k <= cl.i1; k++) {
          const i = closed ? ((k % n) + n) % n : k;
          if (KT[i] === cl.k) {
            if (!cur) subs.push(cur = { i0: k, i1: k, a: DIR[i], cmin: CS[i] });
            else { cur.i1 = k; if (CS[i] < cur.cmin) { cur.cmin = CS[i]; cur.a = DIR[i]; } }
          } else cur = null;
        }
        for (const sb of subs) {
          // 尖: 沿该串最尖处的 "来向 − 去向" 投影最远的点, 在这一串两侧各放宽到簇的边界之内找
          const lo2 = Math.max(cl.i0, sb.i0 - (sb.i1 - sb.i0) - 2), hi2 = Math.min(cl.i1, sb.i1 + (sb.i1 - sb.i0) + 2);
          let best = sb.i0, bv = -Infinity;
          for (let k = subs.length > 1 ? lo2 : cl.i0; k <= (subs.length > 1 ? hi2 : cl.i1); k++) { const v = dot(Ux(k), sb.a); if (v > bv) { bv = v; best = k; } }
          const bi = closed ? ((best % n) + n) % n : best;
          if (closed || (S[bi] > FOLDW && S[n - 1] - S[bi] > FOLDW)) arc.sp.push(mk(S[bi], 'fold', uAt(bi)));
        }
        cl = null;
      };
      let s0 = lo;
      if (closed) { s0 = -1; for (let i = 0; i < n; i++) if (KT[i] < 0) { s0 = i; break; } }
      if (s0 >= 0) for (let q = closed ? s0 + 1 : lo; q < (closed ? s0 + n : hi); q++) {
        const i = closed ? q % n : q, kt = KT[i];
        if (kt >= 0) {
          if (!cl) cl = { i0: q, i1: q, a: DIR[i], cmin: CS[i], k: kt };
          else { cl.i1 = q; if (kt < cl.k || (kt === cl.k && CS[i] < cl.cmin)) { cl.cmin = CS[i]; cl.a = DIR[i]; cl.k = kt; } }
        } else flush();
      }
      flush();
      const fs = arc.sp.filter(q => q.type === 'fold').sort((x, y) => x.t - y.t);
      for (let k = 1; k < fs.length; k++) if (fs[k].t - fs[k - 1].t < FOLDW) fs[k].dup = true;
      if (closed && fs.length > 1 && fs[0].t + len - fs[fs.length - 1].t < FOLDW) fs[fs.length - 1].dup = true;
      arc.sp = arc.sp.filter(q => !q.dup);
    }
    // 折点与接头落在同一处 (引擎分段的并接处常常就是轮廓在天空里折回的地方): 两个特殊点挤在一起, 相邻层间次序随意对调,
    //   段就对不上 (当成拓扑在变). 留接头 (每层都有, 身份稳定), 折点并进去 (foldAt: 画折叠线时照样算折点)
    {
      const js = arc.sp.filter(q => q.type === 'joint');
      if (js.length) for (const f of arc.sp) {
        if (f.type !== 'fold') continue;
        const j = js.find(q => { let d = Math.abs(q.t - f.t); if (closed) d = Math.min(d, len - d); return d < FOLDW; });
        if (j) { j.foldAt = true; f.dup = true; }
      }
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
  const buildRow = (r, nbA, nbB, kind, dense = false) => {
    if (cancelled && cancelled()) throw new Error('cancelled');
    calls++;
    const s = compute(r, { seeds: seedsOf(dense, nbA, nbB) }) || { arcs: [] };
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
      const cut0 = cutSeq; cutSeq += runs.length + 1;
      for (let ri = 0; ri < runs.length; ri++) {
        const [q0, q1] = runs[ri];
        const idx = order.slice(q0, q1 + 1), head = q0 > 0 || closed ? cut(order[q0], order[(q0 - 1 + order.length) % order.length]) : null;
        const tail = q1 < order.length - 1 || closed ? cut(order[q1], order[(q1 + 1) % order.length]) : null;
        const m = idx.length + (head ? 1 : 0) + (tail ? 1 : 0), V = new Float64Array(3 * m);
        let k = 0;
        if (head) { V.set(head, 0); k = 1; }
        for (const i of idx) { V[3 * k] = U[3 * i]; V[3 * k + 1] = U[3 * i + 1]; V[3 * k + 2] = U[3 * i + 2]; k++; }
        if (tail) V.set(tail, 3 * k);
        if (m < 2) continue;
        // 第 ri 段的头与上一段的尾之间是切掉的一段 (闭环: 最后一段的尾接回第一段的头)
        pre.push({ V, m, closed: false, h: head ? 'cut' : 'end', t: tail ? 'cut' : 'end', J: [], hc: cut0 + ri, tc: closed && ri === runs.length - 1 ? cut0 : cut0 + ri + 1 });
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
    // 并集 (opt.union; 降自由度的光路): 各条边界曲线的像互相重叠着铺满同一条轮廓 (取向的每条纤维两头都落在边界上, 轮廓上
    //   每一点至少被盖两次), 每条像的末端多半落在别的像中间 —— 不是轮廓的真末端. 不并的话: 重叠的几张面各铺各的, 假末端
    //   在相邻层间时有时无 (层层当成拓扑在变, 二分加层到上限), 接不上的地方留细缝.
    //   并法: 长的先留; 后面的每条只留没被已留的盖住的部分, 而且从已留弧的端点接着走 (接头逐位相同, 下面的首尾并接
    //   把它们并成一条). 已留的线段登记在方向空间的格子里, 查一个点只看周围 27 格
    if (opt.union && pre.length > 1) {
      const TOLU = 1e-4, G = 96, cells = new Map(), segs = [], ends = [];
      const cellOf = (x, y, z) => (Math.floor((x + 1) * G) * 256 + Math.floor((y + 1) * G)) * 256 + Math.floor((z + 1) * G);
      const reg = a => {
        for (let i = 0; i + 1 < a.m + (a.closed ? 1 : 0); i++) {
          const j = (i + 1) % a.m, si = segs.length / 6;
          const ax = a.V[3 * i], ay = a.V[3 * i + 1], az = a.V[3 * i + 2], bx = a.V[3 * j], by = a.V[3 * j + 1], bz = a.V[3 * j + 2];
          segs.push(ax, ay, az, bx, by, bz);
          const n = Math.max(1, Math.ceil(Math.hypot(bx - ax, by - ay, bz - az) * 2 * G));
          let last = -1;
          for (let k = 0; k <= n; k++) {
            const f = k / n, c = cellOf(ax + f * (bx - ax), ay + f * (by - ay), az + f * (bz - az));
            if (c === last) continue;
            last = c;
            const l = cells.get(c);
            if (l) { if (l[l.length - 1] !== si) l.push(si); } else cells.set(c, [si]);
          }
        }
        if (!a.closed) ends.push([a.V[0], a.V[1], a.V[2]], [a.V[3 * a.m - 3], a.V[3 * a.m - 2], a.V[3 * a.m - 1]]);
      };
      const covered = (x, y, z) => {                  // 这个方向离已留的线不到 TOLU
        const cx = Math.floor((x + 1) * G), cy = Math.floor((y + 1) * G), cz = Math.floor((z + 1) * G);
        for (let dx = -1; dx <= 1; dx++) for (let dy = -1; dy <= 1; dy++) for (let dz = -1; dz <= 1; dz++) {
          const l = cells.get(((cx + dx) * 256 + cy + dy) * 256 + cz + dz);
          if (!l) continue;
          for (const si of l) {
            const ax = segs[6 * si], ay = segs[6 * si + 1], az = segs[6 * si + 2], ex = segs[6 * si + 3] - ax, ey = segs[6 * si + 4] - ay, ez = segs[6 * si + 5] - az;
            const qx = x - ax, qy = y - ay, qz = z - az, l2 = ex * ex + ey * ey + ez * ez;
            let t = l2 > 0 ? (qx * ex + qy * ey + qz * ez) / l2 : 0; t = t < 0 ? 0 : t > 1 ? 1 : t;
            if (Math.hypot(qx - t * ex, qy - t * ey, qz - t * ez) <= TOLU) return true;
          }
        }
        return false;
      };
      // 已留弧的端点里, 落在线段 p → q 上 (2·TOLU 之内) 的那个: 新的一截从它起
      const endOn = (p, q) => {
        const ex = q[0] - p[0], ey = q[1] - p[1], ez = q[2] - p[2], l2 = ex * ex + ey * ey + ez * ez;
        let best = null, bd = 2 * TOLU;
        for (const e of ends) {
          const qx = e[0] - p[0], qy = e[1] - p[1], qz = e[2] - p[2];
          let t = l2 > 0 ? (qx * ex + qy * ey + qz * ez) / l2 : 0; t = t < 0 ? 0 : t > 1 ? 1 : t;
          const d = Math.hypot(qx - t * ex, qy - t * ey, qz - t * ez);
          if (d < bd) { bd = d; best = e; }
        }
        return best;
      };
      const lenOfPre = a => { let l = 0; for (let i = 1; i < a.m; i++) l += angU([a.V[3 * i - 3], a.V[3 * i - 2], a.V[3 * i - 1]], [a.V[3 * i], a.V[3 * i + 1], a.V[3 * i + 2]]); return l; };
      const order = pre.map(a => [lenOfPre(a), a]).sort((x, y) => y[0] - x[0]).map(x => x[1]);
      const out = [];
      let nCut = 0;
      for (const a of order) {
        const P = i => [a.V[3 * i], a.V[3 * i + 1], a.V[3 * i + 2]];
        const on = new Uint8Array(a.m);
        let nOn = 0;
        if (segs.length) for (let i = 0; i < a.m; i++) if (covered(a.V[3 * i], a.V[3 * i + 1], a.V[3 * i + 2])) { on[i] = 1; nOn++; }
        if (nOn === a.m) { nCut++; continue; }                               // 整条被盖住
        if (nOn === 0) { out.push(a); reg(a); continue; }
        nCut++;
        // 没被盖住的各截 (闭环从一个被盖住的点起算, 当开弧切)
        const idx = [];
        if (a.closed) { const k0 = on.indexOf(1); for (let k = 0; k <= a.m; k++) idx.push((k0 + k) % a.m); } else for (let k = 0; k < a.m; k++) idx.push(k);
        const fresh = [];
        for (let q = 0; q < idx.length;) {
          if (on[idx[q]]) { q++; continue; }
          let e = q; while (e + 1 < idx.length && !on[idx[e + 1]]) e++;
          const pts = [];
          if (q > 0) { const pv = P(idx[q - 1]), cu = P(idx[q]); pts.push(endOn(pv, cu) || pv); }
          for (let k = q; k <= e; k++) pts.push(P(idx[k]));
          if (e + 1 < idx.length) { const cu = P(idx[e]), nx = P(idx[e + 1]); pts.push(endOn(cu, nx) || nx); }
          if (pts.length >= 2) {
            const V = new Float64Array(3 * pts.length);
            pts.forEach((p, k) => { V[3 * k] = p[0]; V[3 * k + 1] = p[1]; V[3 * k + 2] = p[2]; });
            fresh.push({ V, m: pts.length, closed: false, h: q === 0 && !a.closed ? a.h : 'end', t: e === idx.length - 1 && !a.closed ? a.t : 'end', J: [], hc: q === 0 && !a.closed ? a.hc : undefined, tc: e === idx.length - 1 && !a.closed ? a.tc : undefined });
          }
          q = e + 1;
        }
        for (const f of fresh) { out.push(f); reg(f); }
      }
      if (nCut) { pre.length = 0; pre.push(...out); stat.unionCuts = (stat.unionCuts || 0) + nCut; }
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
    // 两个末端算不算同一点: 1e-6 rad 以内; 都不在边界线上的放宽到 1e-4 rad (对称面上左右两半的轮廓各自停在离对称面
    //   千分之几度处, 不并的话晕面沿对称面留一道零宽的缝 —— 两条挨着的 "边")
    const near2 = (u, v) => {
      const g = angU(u, v);
      if (g < JOIN) return true;
      if (g >= 1e-4) return false;
      const tolL = Math.max(1e-6, 2e-3 * r);
      return !(LN.length && (nearLines(at(r, u), k => k !== 'cusp', tolL) || nearLines(at(r, v), k => k !== 'cusp', tolL)));
    };
    const endU = (a, tail) => { const i = tail ? a.m - 1 : 0; return [a.V[3 * i], a.V[3 * i + 1], a.V[3 * i + 2]]; };
    const flip = a => { const V = new Float64Array(3 * a.m); for (let i = 0; i < a.m; i++) { V[3 * i] = a.V[3 * (a.m - 1 - i)]; V[3 * i + 1] = a.V[3 * (a.m - 1 - i) + 1]; V[3 * i + 2] = a.V[3 * (a.m - 1 - i) + 2]; } return { V, m: a.m, closed: false, h: a.t, t: a.h, J: a.J.map(j => a.m - 1 - j), hc: a.tc, tc: a.hc }; };
    const cat = (a, b) => { const V = new Float64Array(3 * (a.m + b.m - 1)); V.set(a.V.subarray(0, 3 * a.m)); V.set(b.V.subarray(3, 3 * b.m), 3 * a.m); return { V, m: a.m + b.m - 1, closed: false, h: a.h, t: b.t, J: [...a.J, a.m - 1, ...b.J.map(j => j + a.m - 1)], hc: a.hc, tc: b.tc }; };
    for (let again = true; again;) {
      again = false;
      for (let i = 0; i < pre.length && !again; i++) {
        const a = pre[i];
        if (a.closed) continue;
        if (a.m > 2 && a.h === 'end' && a.t === 'end' && angU(endU(a, false), endU(a, true)) < JOIN) { a.m--; a.closed = true; again = true; break; }
        for (let j = i + 1; j < pre.length && !again; j++) {
          const b = pre[j];
          if (b.closed) continue;
          let c = null;
          if (a.t === 'end' && b.h === 'end' && near2(endU(a, true), endU(b, false))) c = cat(a, b);
          else if (a.t === 'end' && b.t === 'end' && near2(endU(a, true), endU(b, true))) c = cat(a, flip(b));
          else if (a.h === 'end' && b.t === 'end' && near2(endU(a, false), endU(b, true))) c = cat(b, a);
          else if (a.h === 'end' && b.h === 'end' && near2(endU(a, false), endU(b, false))) c = cat(flip(a), b);
          if (c) { pre[i] = c; pre.splice(j, 1); again = true; }
        }
      }
    }
    // 首尾差一点没合上的弧 (两端相距不到半个点距、也不到弧长的 2%, 而且都不在边界线上): 当闭环. 引擎的轮廓在对称面上、
    //   灯旁、几条边界像换手处常差这么一点 (千分之几度到零点几度); 留着当开弧的话, 两个末端挤在一起, 相邻层间谁连谁
    //   是随意的, 8 字形的闭环还可能两层从不同的那一次过叉点处断开 —— 按弧长比例对应就错开半圈
    for (const a of pre) {
      if (a.closed || a.m < 8 || a.h !== 'end' || a.t !== 'end') continue;
      const u0 = endU(a, false), u1 = endU(a, true), g = angU(u0, u1);
      if (g >= 0.5 * SKY0) continue;
      let len = 0;
      for (let i = 1; i < a.m; i++) len += angU([a.V[3 * i - 3], a.V[3 * i - 2], a.V[3 * i - 1]], [a.V[3 * i], a.V[3 * i + 1], a.V[3 * i + 2]]);
      if (g >= 0.02 * len) continue;
      const tolL = Math.max(1e-6, 2e-3 * r);
      if (LN.length && (nearLines(at(r, u0), k => k !== 'cusp', tolL) || nearLines(at(r, u1), k => k !== 'cusp', tolL))) continue;
      a.closed = true; a.J = a.J.filter(j => j > 0 && j < a.m - 1);
      stat.closedArcs = (stat.closedArcs || 0) + 1;
    }
    for (const a of pre) {
      // 并集模式下接头不当特殊点: 那里只是两条边界像换手的地方 (留哪条、从哪接, 随各条的长短变), 位置在相邻层间乱跳
      //   (对称的闭环上, 接头这一层在左、下一层在右 → 按链号对上的段差着十几度, 层层当成对不上)
      makeArc(row, a.V, a.m, a.closed, a.h, a.t, opt.union ? [] : a.J, [a.hc, a.tc]);
      row.raw.push({ p: Float32Array.from({ length: 3 * a.m }, (_, k) => E[k % 3] + r * a.V[k]), d: Float32Array.from(a.V.subarray(0, 3 * a.m)) });
    }
    // 转折缝穿过这一层的地方: 轮廓上多取一个顶点 (arc.extraT). 转折缝是晕面的折痕; 那里没有特殊点的弧 (引擎一口气走过去,
    //   没在那里分段) 顶点跨在折痕两边, 网格削掉折痕的尖, 离线一二十厘米. 只是多一个顶点, 不是特殊点 —— 不改段的划分与对应
    for (const ln of LN) {
      if (ln.kind !== 'turn') continue;
      const rr = i => Math.hypot(ln.p[3 * i] - E[0], ln.p[3 * i + 1] - E[1], ln.p[3 * i + 2] - E[2]);
      let ra = rr(0);
      for (let i = 0; i + 1 < ln.n; i++) {
        const rb = rr(i + 1), a = ra - r, b = rb - r;
        ra = rb;
        if (!(a === 0 || a * b < 0)) continue;
        const f = a / (a - b), u = unit([ln.p[3 * i] + f * (ln.p[3 * i + 3] - ln.p[3 * i]) - E[0], ln.p[3 * i + 1] + f * (ln.p[3 * i + 4] - ln.p[3 * i + 1]) - E[1], ln.p[3 * i + 2] + f * (ln.p[3 * i + 5] - ln.p[3 * i + 2]) - E[2]]);
        let best = null, bd = 1e-3;
        for (const arc of row.arcs) {
          const m = arc.n + (arc.closed ? 1 : 0);
          for (let k = 0; k + 1 < m; k++) {
            const k1 = (k + 1) % arc.n, px = arc.U[3 * k], py = arc.U[3 * k + 1], pz = arc.U[3 * k + 2], ex = arc.U[3 * k1] - px, ey = arc.U[3 * k1 + 1] - py, ez = arc.U[3 * k1 + 2] - pz, l2 = ex * ex + ey * ey + ez * ez;
            let t = l2 > 0 ? ((u[0] - px) * ex + (u[1] - py) * ey + (u[2] - pz) * ez) / l2 : 0; t = t < 0 ? 0 : t > 1 ? 1 : t;
            const d = Math.hypot(u[0] - px - t * ex, u[1] - py - t * ey, u[2] - pz - t * ez);
            if (d < bd) { bd = d; best = { arc, t: arc.S[k] + t * ((k + 1 < arc.n ? arc.S[k + 1] : arc.len) - arc.S[k]) }; }
          }
        }
        if (!best) continue;
        const arc = best.arc, gapTo = t => { let d = Math.abs(t - best.t); if (arc.closed) d = Math.min(d, arc.len - d); return d; }, lim = 0.25 * SKY0;
        if (arc.sp.some(q => gapTo(q.t) < lim) || (arc.extraT || []).some(t => gapTo(t) < lim)) continue;
        if (!arc.closed && (best.t < lim || arc.S[arc.n - 1] - best.t < lim)) continue;
        (arc.extraT || (arc.extraT = [])).push(best.t);
      }
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
      const sh = a.lab && b.lab ? sharedLine(a.lab, b.lab) : null;
      if (a.lab && b.lab && !sh) return Infinity;
      const d = angU(a.u, b.u);
      if (d <= 0.35) return d;
      // 隔得远也连: 两个末端点在同一条边界线上, 两点之间的那段线没出过这两层的半径范围 —— 就是同一条边 (晕面在这里几乎与
      //   球壳相切时, 末端沿线跑得很快: 3-1-5-8 Lowitz Σ87.5 入射面的极点线整圈都在 r 15.797–15.812 之内, 层距 0.004 m 末端挪 17°)
      if (sh && a.type === 'end') { const Q = linePts(sh.id, sh.sa, sh.sb); if (Q.length && inBand(Q, A.r, B.r, BAND)) return d; }
      return Infinity;
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
        for (let i = 0; i + 1 < pts.length; i++) out.push({ arc, t0: pts[i].t, t1: pts[i + 1].t, k0: key(pts[i].q), k1: key(pts[i + 1].q), q0: pts[i].q, q1: pts[i + 1].q, loop: false });
      } else if (!br.length) {
        out.push({ arc, t0: 0, t1: arc.len, k0: 'L', k1: 'L', loop: true });
      } else {
        for (let i = 0; i < br.length; i++) {
          const a = br[i], b = br[(i + 1) % br.length];
          out.push({ arc, t0: a.t, t1: i + 1 < br.length ? b.t : b.t + arc.len, k0: a.chain, k1: b.chain, q0: a, q1: b, loop: false });
        }
      }
    }
    for (const pc of out) { pc.len = pc.t1 - pc.t0; pc.mid = uOn(pc.arc, 0.5 * (pc.t0 + pc.t1)); }
    return out;
  };
  // 对应段的一端: 两层的末端点落在同一条边界线上, 而末端沿线跑得比层距快 (晕面在这里几乎与球壳相切 —— 临界半径旁, 或
  //   边界线整段几乎躺在一层球壳上). 这时按整段弧长比例对应会把这一头拉歪 (对应点隔出十几度, 整条带被当成对不上).
  //   正确的对应: 长的那一层多出来的一头对着两层之间的那段边界线, 剩下的才对着短的那一层.
  //   → 每端 null (不用管) 或 { who: 'A' | 'B' (哪一层多出来), x (多出来的弧长), Q (那段线上的点, A 的末端 → B 的末端) }
  const LONGSIDE = () => Math.max(2 * SKY, 1 * DEG);
  const sideTrim = (A, B, pa, pb, rev) => {
    const out = [null, null];
    if (pa.loop || pb.loop) return out;
    for (let e = 0; e < 2; e++) {
      const eb = rev ? 1 - e : e, qa = e ? pa.q1 : pa.q0, qb = eb ? pb.q1 : pb.q0;
      if (!qa || !qb) continue;
      // 截断点也一样: 截断球面的切口就是这里的边界 (远处截断球面几乎与球壳相切, 切口沿轮廓一层挪十几度); 没有线上的点可插,
      //   多出来的一头直接对着两个截断点的连线
      const isCut = qa.type === 'cut' && qb.type === 'cut';
      if (!isCut && (qa.type !== 'end' || qb.type !== 'end' || !qa.lab || !qb.lab)) continue;
      const d = angU(qa.u, qb.u);
      if (d <= LONGSIDE()) continue;
      let Q = [];
      if (!isCut) {
        const sh = sharedLine(qa.lab, qb.lab);
        if (!sh) continue;
        Q = linePts(sh.id, sh.sa, sh.sb);
        if (!Q.length || !inBand(Q, A.r, B.r, BAND)) continue;
      }
      // 对面那一层的弧上离本层末端最近的位置 (从它自己的末端往里找, 至多到 60%)
      const near = (pc, fromEnd, u) => {
        const lim = 0.6 * pc.len, st = Math.max(Math.min(0.5 * SKY, lim / 48), lim / 400);
        let bx = 0, bd = Infinity;
        for (let x = 0; x <= lim; x += st) { const dd = angU(uOn(pc.arc, fromEnd ? pc.t1 - x : pc.t0 + x), u); if (dd < bd) { bd = dd; bx = x; } }
        return [bx, bd, st];
      };
      const [xb, db, sb] = near(pb, eb === 1, qa.u), [xa, da, sa] = near(pa, e === 1, qb.u);
      if (xb > 1.5 * sb && db < 0.7 * d && db <= da) out[e] = { who: 'B', x: xb, Q };
      else if (xa > 1.5 * sa && da < 0.7 * d) out[e] = { who: 'A', x: xa, Q };
    }
    // 两头都从同一层上截, 截完得还有剩
    for (const [w, pc] of [['A', pa], ['B', pb]]) {
      const x = (out[0] && out[0].who === w ? out[0].x : 0) + (out[1] && out[1].who === w ? out[1].x : 0);
      if (x > 0.9 * pc.len) { out[0] = out[1] = null; break; }
    }
    return out;
  };
  // ---- 两层之间, 边界线也是轮廓的一部分 ("接合"): 同一层上两个没延续到对面的末端点, 若都在边界线上、沿线从一个走到另一个的
  //   那一段整个夹在这两层之间, 那段线就是晕面在这两层之间的边 —— 把它两头的弧接成一条 (拼接段). 几种拓扑变化是同一个道理:
  //     一条弧的两端互相接上: 弧从线上出生 / 缩回线上消失 (月牙), 或闭环碰到线后张开成弧;
  //     两条弧的末端接上: 两条弧在线离眼最远 (最近) 处并成一条 / 一条分成两条.
  //   拼接段照常与对面一层的段对应 (两端的链号, 或闭环对闭环); 顶点 = 各弧自己的顶点 + 线上的点 (网格的边 = 精确边界线).
  //   拼接段带一条合成的 "弧" (各段的采样点 + 线上的点; 只用来按比例取方向、量间距), 这样后面的代码不用区分
  const mkVirtual = (seq, closed) => {
    const us = [];
    for (const it of seq) {
      if (it.Q) { for (const p of it.Q) us.push(unit([p[0] - E[0], p[1] - E[1], p[2] - E[2]])); continue; }
      const n = Math.max(2, Math.ceil(it.pc.len / (0.5 * SKY)));
      for (let k = 0; k <= n; k++) us.push(uOn(it.pc.arc, it.pc.t0 + (it.rev ? 1 - k / n : k / n) * it.pc.len));
    }
    const n = us.length, U = new Float64Array(3 * n), S = new Float64Array(n);
    for (let i = 0; i < n; i++) { U[3 * i] = us[i][0]; U[3 * i + 1] = us[i][1]; U[3 * i + 2] = us[i][2]; if (i) S[i] = S[i - 1] + angU(us[i - 1], us[i]); }
    const len = S[n - 1] + (closed ? angU(us[n - 1], us[0]) : 0), sp = [];
    for (const it of seq) if (it.pc) { const c = it.pc.arc.sp.find(q => q.type === 'cut'); if (c) sp.push(c); }
    const first = seq[0], last = seq[seq.length - 1];
    const pc = {
      arc: { U, n, S, len, closed, sp, virt: true }, t0: 0, t1: len, len, loop: closed, virt: seq,
      k0: closed ? 'L' : first.rev ? first.pc.k1 : first.pc.k0, k1: closed ? 'L' : last.rev ? last.pc.k0 : last.pc.k1,
      q0: closed ? null : first.rev ? first.pc.q1 : first.pc.q0, q1: closed ? null : last.rev ? last.pc.q0 : last.pc.q1,
    };
    pc.mid = uOn(pc.arc, 0.5 * len);
    return pc;
  };
  const joinPieces = (A, B, row, P, joined) => {
    // 候选: 段的自由端 (链不延续到对面), 末端点或截断点
    const ends = [], port = P.map(() => [-1, -1]);
    P.forEach((pc, i) => {
      if (pc.loop) return;
      for (const e of [0, 1]) {
        const q = e ? pc.q1 : pc.q0, k = e ? pc.k1 : pc.k0;
        if (q && (q.type === 'end' || q.type === 'cut') && !(typeof k === 'number' && k >= 0)) { port[i][e] = ends.length; ends.push({ i, e, q }); }
      }
    });
    if (ends.length < 2) return P;
    const cand = [];
    for (let a = 0; a < ends.length; a++) for (let b = a + 1; b < ends.length; b++) {
      const qa = ends[a].q, qb = ends[b].q;
      if (qa.type !== qb.type) continue;
      if (qa.type === 'cut' && (qa.cutId === undefined || qa.cutId !== qb.cutId)) continue;      // 截断点只接同一段切口的两头
      let Q = qa.type === 'end' ? linePath(qa, qb, A.r, B.r) : null;
      // 一条弧的两端 (月牙): 那段线可以比对面那一层还往外鼓出去一点 (至多 0.5%) —— 弧在对面那层上还没出现 (引擎在
      //   出生点旁找不到那条很短的弧), 月牙就是这条弧到边界线之间的全部
      if (Q === null && qa.type === 'end' && ends[a].i === ends[b].i) Q = row === B ? linePath(qa, qb, B.r * (1 - 5e-3), B.r) : linePath(qa, qb, A.r, A.r * (1 + 5e-3));
      const exact = Q !== null;
      if (!exact) {
        // 没有边界线可走 (一自由度没有传线; 或截断球面的切口): 挨得近的直接连 —— 两个末端在这两层之间相遇消失 (两条弧并成
        //   一条 / 一条小弧缩没 / 截断球面刚切进闭环), 它们之间的边界只有这么一小段. 都有标签却走不通的 (两条不相干的边界线
        //   挨着) 只在很近时连
        const lim = qa.type === 'cut' ? MAXEDGE : qa.lab && qb.lab ? GAPSKY : 2.5 * GAPSKY;
        if (angU(qa.u, qb.u) > lim) continue;
        Q = [];
      }
      let l = 0, u = qa.u;
      for (const p of Q) { const v = unit([p[0] - E[0], p[1] - E[1], p[2] - E[2]]); l += angU(u, v); u = v; }
      cand.push({ a, b, Q, l: l + angU(u, qb.u), exact });
    }
    if (!cand.length) return P;
    cand.sort((x, y) => x.l - y.l);                       // 短的先接: 每个末端只接一次
    const jn = new Array(ends.length).fill(null);
    for (const c of cand) if (!jn[c.a] && !jn[c.b]) {
      jn[c.a] = { to: c.b, Q: c.Q }; jn[c.b] = { to: c.a, Q: c.Q.slice().reverse() };
      // 沿精确边界线接上的不再算 "对不上的特殊点"; 直接连的仍算 (照常二分到底, 连线只管把最后剩的小口补上)
      if (c.exact) { joined.add(ends[c.a].q); joined.add(ends[c.b].q); }
    }
    const out = [], done = new Uint8Array(P.length);
    for (let i = 0; i < P.length; i++) {
      if (done[i]) continue;
      const j0 = port[i][0] >= 0 ? jn[port[i][0]] : null, j1 = port[i][1] >= 0 ? jn[port[i][1]] : null;
      if (!j0 && !j1) { out.push(P[i]); done[i] = 1; continue; }
      // 找这一串的起点: 从 i 的 0 端往回走, 走到没接的端 (开链) 或回到 i 的 1 端 (闭环)
      let s = i, ps = 0, closed = false;
      for (let guard = 0; guard <= P.length; guard++) {
        const j = port[s][ps] >= 0 ? jn[port[s][ps]] : null;
        if (!j) break;
        const nx = ends[j.to];
        if (nx.i === i && nx.e === 1) { closed = true; s = i; ps = 0; break; }
        s = nx.i; ps = 1 - nx.e;
      }
      // 从 (s, ps) 起: 进 ps 端、出另一端, 沿接合到下一段 …
      const seq = [];
      let cur = s, pin = ps;
      for (let guard = 0; guard <= P.length; guard++) {
        seq.push({ pc: P[cur], rev: pin === 1 }); done[cur] = 1;
        const pout = 1 - pin, j = port[cur][pout] >= 0 ? jn[port[cur][pout]] : null;
        if (!j) break;
        seq.push({ Q: j.Q });
        const nx = ends[j.to];
        if (nx.i === s && nx.e === ps) { closed = true; break; }
        cur = nx.i; pin = nx.e;
      }
      if (seq.length === 1) { out.push(P[i]); continue; }
      // 一整条开弧, 两端直接连上 (没有边界线): 就是一个闭环 (引擎没把最后一小截走完 —— 灯旁的水滴形轮廓, 尖上差一点;
      //   或闭环刚在没有线的边界上张开). 当真闭环用 (克隆一条闭合的弧, 顶点仍是原弧的): 鞍点、穹顶收口都照常
      //   (两端之间有一小段边界线、比半个点距还短的也一样: 灯旁各条边界线都汇到灯上)
      if (closed && seq.length === 2) {
        const pc = seq[0].pc, a = pc.arc, u0 = [a.U[0], a.U[1], a.U[2]];
        if (!a.closed && pc.t0 === 0 && pc.t1 === a.S[a.n - 1] && seq[1].Q.every(q => angU(u0, unit([q[0] - E[0], q[1] - E[1], q[2] - E[2]])) < 0.5 * SKY) && angU(u0, [a.U[3 * a.n - 3], a.U[3 * a.n - 2], a.U[3 * a.n - 1]]) < 0.5 * SKY || !seq[1].Q.length && !a.closed && pc.t0 === 0 && pc.t1 === a.S[a.n - 1]) {
          const len = a.S[a.n - 1] + angU([a.U[0], a.U[1], a.U[2]], [a.U[3 * a.n - 3], a.U[3 * a.n - 2], a.U[3 * a.n - 1]]);
          const arc = { U: a.U, n: a.n, S: a.S, len, closed: true, sp: a.sp, ai: a.ai, base: a };
          const lp = { arc, t0: 0, t1: len, k0: 'L', k1: 'L', loop: true, len };
          lp.mid = uOn(arc, 0.5 * len);
          out.push(lp);
          continue;
        }
      }
      const v = mkVirtual(seq, closed);
      out.push(v);
    }
    return out;
  };
  // 这一圈上 8 个点到那一圈的最大夹角 (那一圈按采样点算)
  const loopDist = (pa, pb) => {
    const b = pb.arc, st = Math.max(1, b.n >> 7);
    let mx = 0;
    for (let k = 0; k < 8; k++) {
      const u = uOn(pa.arc, pa.t0 + (k + 0.5) / 8 * pa.len);
      let best = -Infinity;
      for (let i = 0; i < b.n; i += st) { const c = u[0] * b.U[3 * i] + u[1] * b.U[3 * i + 1] + u[2] * b.U[3 * i + 2]; if (c > best) best = c; }
      mx = Math.max(mx, best >= 1 ? 0 : Math.acos(Math.max(-1, best)));
    }
    return mx;
  };
  // 月牙: 拼接成的闭环只有一条弧 + 一段线. 弧与线按比例对应的点得都挨着 (闭环碰到边界线后张开成的弧不是月牙: 两端在线上
  //   挨得很近, 弧却有一百多度长 —— 它对着的是对面一层的闭环). 返回弧与线的最大间距; 不是月牙返回 −1
  const lensGap = pc => {
    if (!pc.virt || !pc.loop || pc.virt.length !== 2) return -1;
    const sub = pc.virt[0], a = sub.pc, ua = f => uOn(a.arc, a.t0 + (sub.rev ? 1 - f : f) * a.len);
    const W = [ua(1), ...pc.virt[1].Q.map(p => unit([p[0] - E[0], p[1] - E[1], p[2] - E[2]])), ua(0)], cw = [0];   // 线: 弧的出口端 → 入口端
    for (let k = 1; k < W.length; k++) cw.push(cw[k - 1] + angU(W[k - 1], W[k]));
    const lw = cw[cw.length - 1];
    let g = 0;
    for (let k = 1; k < 16; k++) {
      const s = k / 16, tw = s * lw;
      let j = 1; while (j + 1 < W.length && cw[j] < tw) j++;
      const f = (tw - cw[j - 1]) / Math.max(1e-15, cw[j] - cw[j - 1]), p = W[j - 1], q = W[j];
      g = Math.max(g, angU(unit([p[0] + f * (q[0] - p[0]), p[1] + f * (q[1] - p[1]), p[2] + f * (q[2] - p[2])]), ua(1 - s)));
    }
    return g > 0.5 * MAXEDGE ? -1 : g;
  };
  // 相邻两层 A < B 的分析: 对应段 (pairs: [pa, pb, rev, trim]), 对不上的段, 对不上的特殊点, 天空间距, 边界线弦高
  const analyze = (A, B) => {
    const chA = new Set(A.sp.map(q => q.chain)), chB = new Set(B.sp.map(q => q.chain));
    const keep = new Set([...chA].filter(c => chB.has(c)));
    const joined = new Set();
    const PA = joinPieces(A, B, A, piecesOf(A, keep), joined), PB = joinPieces(A, B, B, piecesOf(B, keep), joined);
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
        // 闭环: 中点落在哪由起点定 (随意), 改按形状比 —— 这一圈上 8 个点到那一圈的最大距离
        const d = pa.loop ? loopDist(pa, pb) : angU(pa.mid, pb.mid);
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
    let unB = PB.filter(pb => !used.has(pb));
    // 沿边界线接上的末端点不算 "对不上的特殊点" (那段线就是它们在这两层之间的去处)
    const spA = A.sp.filter(q => !keep.has(q.chain) && !joined.has(q)), spB = B.sp.filter(q => !keep.has(q.chain) && !joined.has(q));
    // A closed trace can start anywhere and run in either direction. Measure
    // separation with the same cyclic correspondence used for triangulation.
    let gap = 0;
    const matched = [];
    for (const pair of pairs) {
      const [pa, pb, rev] = pair;
      let start = rev ? pb.t1 : pb.t0, sign = rev ? -1 : 1;
      if (pa.loop && pb.loop) {
        const u = uOn(pa.arc, pa.t0), b = pb.arc;
        let best = -Infinity;
        for (let i = 0; i < b.n; i++) {
          const c = u[0] * b.U[3*i] + u[1] * b.U[3*i+1] + u[2] * b.U[3*i+2];
          if (c > best) { best = c; start = b.S[i]; }
        }
        // 闭环自己打叉 (过灯的 8 字形: 起点常常就在叉点上) 时, 离起点最近的地方有两处 —— 两次经过叉点. 差不多近的都试,
        //   两个走向也都试, 取四分点对得最近的 (只看最近的那个点会把整圈错开半圈)
        const thr = Math.cos(Math.min(Math.PI, (best >= 1 ? 0 : Math.acos(best)) + Math.max(SKY, 0.5 * DEG)));
        const qs = [0.25, 0.5, 0.75].map(f => uOn(pa.arc, pa.t0 + f * pa.len));
        let bc = Infinity, lastI = -99;
        for (let i = 0; i < b.n; i++) {
          if (u[0] * b.U[3*i] + u[1] * b.U[3*i+1] + u[2] * b.U[3*i+2] < thr) continue;
          if (i - lastI < 3 && i - lastI > 0 && bc < Infinity) { lastI = i; continue; }       // 同一处挨着的点不必都试
          lastI = i;
          for (const sg of [1, -1]) {
            let cs = 0;
            for (let k = 0; k < 3; k++) cs += angU(qs[k], uOn(b, b.S[i] + sg * (k + 1) / 4 * pb.len));
            if (cs < bc) { bc = cs; start = b.S[i]; sign = sg; }
          }
        }
      }
      let g = 0;
      // 末端沿边界线滑得快的一头先截掉 (那一头对着边界线, 见 sideTrim), 间距只量剩下互相对着的部分
      const trim = sideTrim(A, B, pa, pb, rev);
      if (trim[0] || trim[1]) {
        const cutOf = (w, e) => trim[e] && trim[e].who === w ? trim[e].x : 0;
        const a0 = pa.t0 + cutOf('A', 0), a1 = pa.t1 - cutOf('A', 1);
        // pb 在自己的参数里: pa 的 0 端对着 pb 的 (rev ? 1 : 0) 端
        const b0 = pb.t0 + cutOf('B', rev ? 1 : 0), b1 = pb.t1 - cutOf('B', rev ? 0 : 1);
        for (let k = 0; k <= 8; k++) {
          const f = k / 8;
          g = Math.max(g, angU(uOn(pa.arc, a0 + f * (a1 - a0)), uOn(pb.arc, rev ? b1 - f * (b1 - b0) : b0 + f * (b1 - b0))));
        }
        pair.push(trim);
      } else
      for (let k = 0; k <= 8; k++) {
        const f = k / 8;
        g = Math.max(g, angU(uOn(pa.arc, pa.t0 + f * pa.len), uOn(pb.arc, start + sign * f * pb.len)));
      }
      // 对应点隔得比三角形边长上限还远的一对段, 铺出来的三角形会被整条丢掉 (留一条空带): 这不是 "对应", 是这两层之间
      //   拓扑变了 (开弧两端还在、中段掐成了闭环 …) 或者有一层缺了东西 —— 当对不上的段: 先二分, 分到底就近缝
      if (g > MAXEDGE) { if (log && opt.debugGap) log(`gap reject r ${A.r.toFixed(5)}..${B.r.toFixed(5)}: g ${(g / DEG).toFixed(2)}° > ${(MAXEDGE / DEG).toFixed(2)}°, rev ${rev}, loop ${pa.loop}/${pb.loop}, len ${(pa.len / DEG).toFixed(2)}/${(pb.len / DEG).toFixed(2)}, k ${pa.k0}|${pa.k1} ${pb.k0}|${pb.k1}, trim ${pair.length > 3}`); unA.push(pa); unB.push(pb); }
      else { matched.push(pair); gap = Math.max(gap, g); }
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
    // 月牙 (一条弧 + 一段边界线围成, 对面没有对应的闭环): 不算对不上, 间距 = 弧与线的最大间距
    const lens = [];
    const takeLens = (row, list) => list.filter(pc => { const g = lensGap(pc); if (g < 0) return true; lens.push({ row, pc }); gap = Math.max(gap, g); return false; });
    const unA2 = takeLens(A, unA); unB = takeLens(B, unB);
    const event = unA2.length > 0 || unB.length > 0 || spA.length > 0 || spB.length > 0;
    return { pairs: matched, unA: unA2, unB, spA, spB, gap, lineDev, event, keep, lens };
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
    // 首末两层 = 这一段的边界, 半径取准 (pow 的舍入会让它与相邻段差一个末位), 一定算
    const seam = k === 0 || k === N0 - 1;
    const r = k === 0 ? rMin : k === N0 - 1 ? rMax : rMin * Math.pow(rMax / rMin, k / Math.max(1, N0 - 1));
    if (seam || (Math.abs(r / D - 1) > PINL && !crit.some(rc => Math.abs(r / rc - 1) < 2 * CRIT))) rs.push([r, 'init']);
  }
  if (rMin < D * (1 - PINL) && rMax > D * (1 + PINL)) rs.push([D * (1 - PINL), 'lampLo'], [D * (1 + PINL), 'lampHi']);
  for (const rc of crit) rs.push([rc * (1 - CRIT), 'critLo'], [rc * (1 + CRIT), 'critHi']);
  rs.sort((p, q) => p[0] - q[0]);
  let rows = [];
  for (const [r, kind] of rs) {
    // 交界层不带种子算: 相邻两段各算一次, 都是 "引擎全搜" 这同一个确定的计算, 轮廓逐点相同 (带相邻层的种子的话, 两段的
    //   相邻层不同, 闭环的起点、取样都不同)
    const seam = r === rMin || r === rMax;
    const row = buildRow(r, seam ? null : rows[rows.length - 1], null, kind);
    row.seam = seam;
    rows.push(row);
  }
  if (log) log(`initial ${rows.length} shells, ${crit.length} critical radii, ${LN.length} lines`);
  // 顶点预算 (opt.budget): 铺满大半个天空的晕面每层的线有几百度长. 按初始层估一下顶点数 (加密大约再多一半),
  //   超了就把点距按比例放大 (至多 2°); 加密过程中超了预算就不再为天空间距 / 弦高加层, 只为拓扑变化二分
  const vEst = rows => { let v = 0; for (const rw of rows) for (const a of rw.arcs) v += Math.max(a.closed ? 8 : 4, Math.ceil(a.len / SKY)) + a.sp.length; return v; };
  {
    const est = 1.5 * vEst(rows) * Math.max(1, N0 / Math.max(1, rows.length));
    if (est > VBUDGET) { SKY = Math.min(2 * DEG, SKY * est / VBUDGET); GAPSKY = Math.max(1.2 * DEG, 3 * SKY); MAXEDGE = Math.max(4 * GAPSKY, 5 * DEG); }
  }

  // ---- 层的一致性校验: 不盲信单独一层.
  //   轮廓总长随半径连续 (掐断、合并、出生、消失、被截断球面切短时都连续). 夹在 A、B 之间的一层 M 若总长跳出了两邻层的范围
  //   (odd), 而且放进去以后两侧都对不上 (useless: 只连这三层试一下) —— 这一层多半是引擎漏了 (或重了) 一条轮廓: 它帮不上
  //   忙, 反而把两侧的带子都毁掉 (3-1-5 Parry Σ87.5: 掐断点旁一层丢了 211° 的闭环, 两侧各留一圈空带).
  //   先加密种子重算, 再把半径挪一挪 (离开掐断点) 重算; 都不行就不要这一层, 这两层之间不再二分 (就近缝).
  const lenOf = row => { if (row.slen === undefined) { let l = 0; for (const a of row.arcs) l += a.len; row.slen = l; } return row.slen; };
  const odd = (M, A, B) => {
    const la = lenOf(A), lb = lenOf(B), lm = lenOf(M), tol = 0.5 * Math.abs(la - lb) + 0.02 * Math.max(la, lb) + 0.5 * DEG;
    return lm < Math.min(la, lb) - tol || lm > Math.max(la, lb) + tol;
  };
  const useless = (M, A, B) => { relink([A, M, B]); return analyze(A, M).event && analyze(M, B).event; };
  //   只在两邻层挨得很近 (半径差 < 2%) 时查: 隔得远的两层之间, 轮廓总长本来就可以不单调 (灯旁、铺满全天的晕面), 这时
  //   拒掉中间的层反而留下一条几米宽、没法再分的带 (一自由度 1-3-2 Σ12: r 24.5–30 之间整段只能就近缝, 21.9 m² 假面)
  const bad = (M, A, B) => B.r / A.r - 1 < 0.02 && odd(M, A, B) && useless(M, A, B);
  const checkedRow = (r, A, B, kind) => {
    let M = buildRow(r, A, B, kind);
    if (!bad(M, A, B)) return M;
    stat.retry = (stat.retry || 0) + 1;
    M = buildRow(r, A, B, kind, true);
    if (!bad(M, A, B)) return M;
    for (const f of [0.35, 0.65]) {
      const r2 = A.r * Math.pow(B.r / A.r, f);
      if (!(r2 > A.r * (1 + 0.2 * EPSR) && r2 < B.r * (1 - 0.2 * EPSR))) continue;
      M = buildRow(r2, A, B, kind, true);
      if (!bad(M, A, B)) return M;
    }
    stat.rejected = (stat.rejected || 0) + 1;
    if (log) log(`shell r ${r.toFixed(5)} rejected: contour length ${(lenOf(M) / DEG).toFixed(1)}° between ${(lenOf(A) / DEG).toFixed(1)}° and ${(lenOf(B) / DEG).toFixed(1)}°, unmatched on both sides`);
    return null;
  };
  // 初始层也查 (交界层、临界半径 / 灯两侧的层不动): 可疑的重算, 还不行就去掉 (留给两邻层之间的二分)
  for (let k = 1; k + 1 < rows.length; k++) {
    const A = rows[k - 1], M = rows[k], B = rows[k + 1];
    if (M.seam || M.kind !== 'init' || !bad(M, A, B)) continue;
    const M2 = checkedRow(M.r, A, B, 'init');
    if (M2) rows[k] = M2; else { rows.splice(k, 1); k--; }
  }

  // ---- 5) 加密
  const noSplit = new Set(), pkey = (A, B) => A.r + ':' + B.r;
  const pairSkip = (A, B) => A.kind === 'lampLo' && B.kind === 'lampHi' || A.kind === 'critLo' && B.kind === 'critHi' && B.r / A.r - 1 < 3 * CRIT || noSplit.has(pkey(A, B));
  const nearLamp = (A, B) => {
    const dl = row => { if (row.rl === undefined) { let m = Infinity; for (const arc of row.arcs) for (let i = 0; i < arc.n; i += Math.max(1, arc.n >> 6)) { const p = at(row.r, [arc.U[3 * i], arc.U[3 * i + 1], arc.U[3 * i + 2]]); m = Math.min(m, Math.hypot(p[0] - L[0], p[1] - L[1], p[2] - L[2])); } row.rl = m; } return row.rl; };
    return B.r - A.r > 0.35 * Math.max(Math.min(dl(A), dl(B)), 2e-3 * D);
  };
  const refine = () => {
  for (let pass = 0; pass < 28; pass++) {
    relink(rows);
    const ins = new Map(), over = vEst(rows) > 2.5 * VBUDGET;
    const AN = [], why = { ev: 0, gap: 0, lamp: 0, line: 0, rad: 0, radR: [] };
    for (let k = 0; k + 1 < rows.length; k++) {
      const A = rows[k], B = rows[k + 1];
      if (!A.arcs.length && !B.arcs.length && !A.sp.length && !B.sp.length) { AN.push(null); continue; }
      const an = analyze(A, B);
      AN.push(an);
      if (pairSkip(A, B)) continue;
      const ratio = B.r / A.r - 1, mid = Math.sqrt(A.r * B.r);
      if (an.event) { if (ratio > EPSR) { ins.set(mid, [A, B]); why.ev++; } continue; }
      // A budget limits visual density, but a band wider than the triangle
      // acceptance limit would be dropped altogether. Resolve it first.
      if (an.gap > 0.7 * MAXEDGE && ratio > EPSR) { ins.set(mid, [A, B]); continue; }
      if (ratio <= 4 * EPSR) continue;
      if (over) continue;
      if (an.gap > GAPSKY || nearLamp(A, B) || an.lineDev > tol3(B.r)) { ins.set(mid, [A, B]); if (log) { if (an.gap > GAPSKY) why.gap++; else if (an.lineDev > tol3(B.r)) why.line++; else why.lamp++; } }
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
      if (dev > 2 * tol3(B.r)) { ins.set(Math.sqrt(A.r * B.r), [A, B]); ins.set(Math.sqrt(B.r * C.r), [B, C]); if (log) { why.rad++; if (why.radR.length < 6) why.radR.push(`${B.r.toFixed(3)}:${dev.toFixed(3)}`); } }
    }
    if (!ins.size) break;
    if (calls >= MAXCALLS || rows.length >= MAXROWS || Date.now() - t0 > SOFTMS) { stat.capped = Date.now() - t0 > SOFTMS ? 'time' : 'rows'; break; }
    const add = [];
    for (const [r, [A, B]] of [...ins.entries()].sort((p, q) => p[0] - q[0])) {
      if (calls >= MAXCALLS || rows.length + add.length >= MAXROWS || Date.now() - t0 > SOFTMS) break;
      const M = checkedRow(r, A, B, 'ref');
      if (M) add.push(M); else noSplit.add(pkey(A, B));
    }
    rows = rows.concat(add).sort((p, q) => p.r - q.r);
    if (log) log(`pass ${pass}: +${add.length} shells (${rows.length}) [事件 ${why.ev}, 天空间距 ${why.gap}, 边界线弦高 ${why.line}, 近灯 ${why.lamp}, 径向弯曲 ${why.rad} ${why.radR.join(' ')}]`);
  }
  };
  // 带路重算: 相邻两层对不上, 一层有、另一层没有的一段轮廓 —— 可能是真的在这里出生 / 消失, 也可能是引擎在没有的那一层漏了
  //   (全局搜索找种子有死角: 3-5-6-8 Lowitz Σ10.5 焦散在灯距附近 r 14.94–15.05 找不到上面那条三十几度的弧, 再往外就找得到;
  //   初始层由里往外算, 只有里面一层的种子, 外面先找到的弧带不进来). 用有的那一层的状态点当种子 (加密) 把没有的那一层重算
  //   一遍: 找回来的长度超过缺口的一半就换上. 由外往里、由里往外各扫一遍 (换上的层接着给下一层带路). 交界层不动 (相邻两段
  //   要算出同一条轮廓)
  const reseed = () => {
    let changed = 0;
    const unLen = list => { let l = 0; for (const pc of list) l += pc.len; return l; };
    for (const dir of [-1, 1]) {
      for (let k = dir < 0 ? rows.length - 2 : 0; k >= 0 && k + 1 < rows.length; k += dir) {
        if (calls >= MAXCALLS + 40 || Date.now() - t0 > SOFTMS) return changed;
        const A = rows[k], B = rows[k + 1];
        if (!A.arcs.length && !B.arcs.length) continue;
        relink([A, B]);
        const an = analyze(A, B);
        if (!an.unA.length && !an.unB.length) continue;
        const la = unLen(an.unA), lb = unLen(an.unB);
        // 穷的一侧: 由外往里扫时修 A (用 B 带路), 由里往外扫时修 B
        if (dir < 0 ? !(lb > la + 2 * SKY0) : !(la > lb + 2 * SKY0)) continue;
        const poor = dir < 0 ? A : B, rich = dir < 0 ? B : A, kp = dir < 0 ? k : k + 1;
        if (poor.seam || poor.reseeded) continue;
        const other = dir < 0 ? rows[k - 1] || null : rows[k + 2] || null;
        const R2 = buildRow(poor.r, dir < 0 ? other : rich, dir < 0 ? rich : other, poor.kind, true);
        R2.reseeded = true; poor.reseeded = true;
        if (lenOf(R2) > lenOf(poor) + 0.5 * Math.abs(lb - la)) { rows[kp] = R2; changed++; }
      }
    }
    return changed;
  };
  refine();
  for (let it = 0; it < 2; it++) {
    const n = reseed();
    if (!n) break;
    stat.reseeded = (stat.reseeded || 0) + n;
    if (log) log(`reseed: ${n} shells recomputed with seeds from the neighbour that has more`);
    refine();
  }
  // 预算收尾: 拓扑变化多的一段, 二分加的层不受预算管 (每层都是整圈轮廓), 顶点数会超出几倍. 加密完再按实际层数把段内点距
  //   放大到预算以内 (至多 3°; 交界层照旧按 SKY0, 与相邻段对得上)
  {
    const v = vEst(rows);
    if (v > 1.25 * VBUDGET && SKY < 3 * DEG) {
      // 只放大点距. 层间距 (GAPSKY) 不跟着放大: 层已经按原来的 GAPSKY 加好了, 再放大它, 就近缝的 "够得着" 范围跟着变成
      //   二十几度, 会把天上隔着 20° 的两个小闭环缝到一起 (3-5-6-8 Lowitz Σ34.5 蓝边 2 的两个尖顶之间多出一张假面)
      SKY = Math.min(3 * DEG, SKY * v / VBUDGET); MAXEDGE = Math.max(MAXEDGE, 2.5 * SKY);
      stat.coarsened = +(SKY / DEG).toFixed(2);
    }
  }
  relink(rows);
  if (opt.debugFinal) opt.debugFinal(rows, analyze);      // 测试 / 排查用: 看最终各层的特殊点与对应

  // ---- 4) / 6) 三角形
  // 每层的顶点: 每条弧上 = 全部特殊点 + 相邻两个特殊点之间按点距均匀取; 上下两条带共用
  const vertsOf = (row, arc) => {
    if (arc.base) return vertsOf(row, arc.base);       // 克隆的闭合弧 (joinPieces): 顶点是原弧的
    if (arc.vt) return arc.vt;
    const ts = arc.sp.filter(q => q.type !== 'focal').map(q => q.t);
    if (arc.extraT) ts.push(...arc.extraT);
    if (!arc.closed) { ts.push(0, arc.S[arc.n - 1]); }
    const T = [...new Set(ts.map(t => +t.toFixed(12)))].sort((a, b) => a - b);
    const all = [];
    const spans = arc.closed ? (T.length ? T.map((t, i) => [t, i + 1 < T.length ? T[i + 1] : T[0] + arc.len]) : [[0, arc.len]]) : T.slice(0, -1).map((t, i) => [t, T[i + 1]]);
    // 每段按点距取, 但小弧也要保住形状: 闭环全周至少 8 段, 开弧至少 4 段 (比点距还短的小透镜原来只剩一个点, 铺不出三角形)
    const minSeg = arc.closed ? 8 : 4, sky = row.seam ? SKY0 : SKY;
    for (const [a, b] of spans) {
      const m = Math.max(1, Math.ceil((b - a) / sky), arc.len > 0 ? Math.ceil(minSeg * (b - a) / arc.len) : 1);
      for (let k = 0; k < m; k++) all.push(a + (b - a) * k / m);
    }
    if (!arc.closed) all.push(T[T.length - 1]);
    const vt = { t: Float64Array.from(all), id: new Int32Array(all.length) };
    const fs = row.seam ? (row.r === rMin ? 2 : 4) : 0;
    for (let i = 0; i < all.length; i++) { const u = uOn(arc, all[i]); vt.id[i] = addV(at(row.r, u), u, fs | (!arc.closed && (i === 0 || i === all.length - 1) ? 1 : 0)); }
    arc.vt = vt;
    return vt;
  };
  const cumFr = ids => { const c = [0]; for (let k = 1; k < ids.length; k++) c.push(c[k - 1] + vang(ids[k - 1], ids[k])); const t = c[c.length - 1] || 1; return c.map(x => x / t); };
  const lineV = Q => Q.map(p => { stat.lineVerts++; return addV(p, unit([p[0] - E[0], p[1] - E[1], p[2] - E[2]]), 1); });
  // 段 → 顶点序号串 (按 t0 → t1) 与各点在段内的比例
  const pieceVerts = (row, pc) => {
    if (pc.virt) {                                     // 拼接段: 各弧自己的顶点 (与上下两条带共用) + 线上的点
      const ids = [];
      for (const it of pc.virt) {
        if (it.Q) { if (!it.V) it.V = lineV(it.Q); ids.push(...it.V); continue; }
        const V = pieceVerts(row, it.pc);
        if (it.rev) for (let k = V.ids.length - 1; k >= 0; k--) ids.push(V.ids[k]); else ids.push(...V.ids);
      }
      if (pc.loop) ids.push(ids[0]);
      return { ids, fr: cumFr(ids) };
    }
    const vt = vertsOf(row, pc.arc), n = vt.t.length, ids = [], fr = [], ts = [];       // ts: 各顶点在弧上的参数 (鞍点处切段用)
    const L0 = pc.arc.len, span = Math.max(1e-15, pc.t1 - pc.t0);
    if (pc.loop) {
      for (let i = 0; i < n; i++) { ids.push(vt.id[i]); fr.push((vt.t[i] - vt.t[0]) / L0); ts.push(vt.t[i]); }
      ids.push(vt.id[0]); fr.push(1); ts.push(vt.t[0] + L0);
      return { ids, fr, ts };
    }
    const eps = 1e-10;
    const push = (i, t) => { ids.push(vt.id[i]); fr.push((t - pc.t0) / span); ts.push(t); };
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
    return { ids, fr, ts };
  };
  // 两串顶点按比例拉链
  const zipFr = (A, B, i0 = 0, i1 = A.ids.length - 1, j0 = 0, j1 = B.ids.length - 1, em = emit) => {
    let i = i0, j = j0;
    while (i < i1 || j < j1) {
      if (j >= j1 || (i < i1 && A.fr[i + 1] <= B.fr[j + 1])) { em(A.ids[i], A.ids[i + 1], B.ids[j]); i++; }
      else { em(A.ids[i], B.ids[j + 1], B.ids[j]); j++; }
    }
  };
  const emitAny = (a, b, c) => { if (a !== b && b !== c && c !== a) tris.push(a, b, c); };      // 不查边长 (对应关系可信的地方)
  const revV = V => ({ ids: V.ids.slice().reverse(), fr: V.fr.slice().reverse().map(f => 1 - f) });
  // 段端的特殊点在两层上落在同一条三维边界线上: 两层之间沿这条线插入它的点 (边界线就是网格的边)
  const sideLine = (A, B, qa, qb) => {
    const sh = qa && qb ? sharedLine(qa.lab, qb.lab) : null;
    if (!sh) return null;
    const Q = linePts(sh.id, sh.sa, sh.sb);
    if (!Q.length || !lineOk(Q, qa.p, qb.p, A.r, B.r)) return null;
    return Q.map(p => { stat.lineVerts++; return addV(p, unit([p[0] - E[0], p[1] - E[1], p[2] - E[2]]), 1); });
  };
  const spAt = (row, pc, end) => {                      // 段一端的特殊点 (链延续到对面的才算)
    const k = end ? pc.k1 : pc.k0;
    return typeof k === 'number' && k >= 0 ? (end ? pc.q1 : pc.q0) || null : null;
  };
  const loopAlign = (VA, VB) => {             // 闭环: VB 转到离 VA 起点最近的点起, 两个走向里取四分点对得近的; 比例按累计天空弧长
    const ring = VB.ids.slice(0, -1), m = ring.length;
    let bd = Infinity;
    for (let j = 0; j < m; j++) bd = Math.min(bd, vang(VA.ids[0], ring[j]));
    const mk = (j0, sgn) => { const ids = []; for (let k = 0; k <= m; k++) ids.push(ring[((j0 + sgn * k) % m + m) % m]); return ids; };
    const frOf = ids => { const c = [0]; for (let k = 1; k < ids.length; k++) c.push(c[k - 1] + vang(ids[k - 1], ids[k])); const t = c[c.length - 1] || 1; return c.map(x => x / t); };
    const fa = frOf(VA.ids), A2 = { ids: VA.ids, fr: fa };
    const at4 = (V, f) => { let k = 0; while (k + 1 < V.fr.length && V.fr[k + 1] < f) k++; return V.ids[k]; };
    // 起点: 离 VA 起点差不多最近的几个点都试 (闭环自己打叉时最近的地方有两处, 见 analyze), 两个走向都试, 取 7 个等分点对得最近的
    let best = null, bc = Infinity;
    const lim = bd + Math.max(SKY, 0.5 * DEG);
    for (let j0 = 0; j0 < m; j0++) {
      if (vang(VA.ids[0], ring[j0]) > lim) continue;
      for (const sg of [1, -1]) {
        const ids = mk(j0, sg), V = { ids, fr: frOf(ids) };
        let c = 0;
        for (let k = 1; k < 8; k++) c += vang(at4(A2, k / 8), at4(V, k / 8));
        if (c < bc) { bc = c; best = V; }
      }
    }
    return [A2, best];
  };
  const strip = (A, B, pa, pb, rev, trim = null) => {
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
    // 末端沿边界线滑得快的一头 (analyze 里的 sideTrim): 多出来的那一层的这一头与两层之间的那段边界线拉链
    //   (线: 对面那层的末端 → 这层的末端; 弧: 这层上从截口走到末端), 剩下的部分再与对面那层按比例拉链
    const longSide = (tr, e) => {
      const VX = tr.who === 'B' ? VB : VA, pcX = tr.who === 'B' ? pb : pa, n = VX.ids.length;
      const fx = tr.x / Math.max(1e-15, pcX.len);
      let m;
      if (e === 0) { m = 0; while (m + 1 < n - 1 && VX.fr[m + 1] <= fx + 1e-9) m++; if (m < 1) return false; }
      else { m = n - 1; while (m - 1 > 0 && VX.fr[m - 1] >= 1 - fx - 1e-9) m--; if (m > n - 2) return false; }
      const Qv = lineV(tr.Q), endA = e === 0 ? VA.ids[0] : VA.ids[VA.ids.length - 1], endB = e === 0 ? VB.ids[0] : VB.ids[VB.ids.length - 1];
      const X = tr.who === 'B' ? [endA, ...Qv, endB] : [endB, ...Qv.reverse(), endA], Y = [];
      if (e === 0) for (let k = m; k >= 0; k--) Y.push(VX.ids[k]); else for (let k = m; k <= n - 1; k++) Y.push(VX.ids[k]);
      // 不查边长: 这一头夹在边界 (线 / 截断球面的切口) 与弧之间, 没有线上的点可插时 (截断点) 三角形又细又长, 但位置是对的
      zipFr({ ids: X, fr: cumFr(X) }, { ids: Y, fr: cumFr(Y) }, 0, X.length - 1, 0, Y.length - 1, emitAny);
      if (tr.who === 'B') { if (e === 0) j0 = m; else j1 = m; } else if (e === 0) i0 = m; else i1 = m;
      stat.longSides = (stat.longSides || 0) + 1;
      return true;
    };
    const t0 = !!(trim && trim[0]) && longSide(trim[0], 0), t1 = !!(trim && trim[1]) && longSide(trim[1], 1);
    // 起点一侧
    const s0 = !t0 && VA.ids.length > 2 && VB.ids.length > 2 && i1 - i0 >= 1 && j1 - j0 >= 1 ? sideLine(A, B, spAt(A, pa, false), spAt(B, pb, rev)) : null;
    if (s0) {
      const S = [VA.ids[0], ...s0, VB.ids[0]], m = S.length >> 1;
      for (let k = 0; k + 1 < S.length; k++) emit(S[k], S[k + 1], k < m ? VA.ids[1] : VB.ids[1]);
      emit(VA.ids[1], S[m], VB.ids[1]);
      i0 = 1; j0 = 1;
    }
    const s1 = !t1 && VA.ids.length > 2 && VB.ids.length > 2 && i1 - i0 >= 1 && j1 - j0 >= 1 ? sideLine(A, B, spAt(A, pa, true), spAt(B, pb, !rev)) : null;
    if (s1) {
      const S = [VA.ids[i1], ...s1, VB.ids[j1]], m = S.length >> 1;
      for (let k = 0; k + 1 < S.length; k++) emit(S[k], S[k + 1], k < m ? VA.ids[i1 - 1] : VB.ids[j1 - 1]);
      emit(VA.ids[i1 - 1], S[m], VB.ids[j1 - 1]);
      i1--; j1--;
    }
    if (t0 || t1) {
      const sub = (V, a, b) => { const f0 = V.fr[a], sp = Math.max(1e-15, V.fr[b] - f0); return { ids: V.ids.slice(a, b + 1), fr: V.fr.slice(a, b + 1).map(f => (f - f0) / sp) }; };
      zipFr(sub(VA, i0, i1), sub(VB, j0, j1));
    } else zipFr(VA, VB, i0, i1, j0, j1);
  };
  // 月牙 (analyze 里的 lens): 弧与两端之间的那段边界线拉链 (网格的边 = 精确边界线)
  const lensZip = (row, pc) => {
    const sub = pc.virt[0], jn = pc.virt[1];
    let V = pieceVerts(row, sub.pc);
    if (sub.rev) V = revV(V);
    if (V.ids.length < 2) return;
    if (!jn.V) jn.V = lineV(jn.Q);
    const X = [V.ids[0], ...jn.V.slice().reverse(), V.ids[V.ids.length - 1]];       // 线存的是出口端 → 入口端
    zipFr({ ids: V.ids, fr: cumFr(V.ids) }, { ids: X, fr: cumFr(X) });
    stat.lineCaps = (stat.lineCaps || 0) + 1;
  };
  let lampV = -1;
  const capTo = (row, pc, toLamp) => {
    const V = pieceVerts(row, pc);
    if (!V.ids.length) return;
    let c;
    if (toLamp) { if (lampV < 0) lampV = addV(L, [Dx / D, Dy / D, Dz / D], 1); c = lampV; }
    else {
      let x = 0, y = 0, z = 0;
      for (const v of V.ids) { x += VU[3 * v]; y += VU[3 * v + 1]; z += VU[3 * v + 2]; }
      const u = unit([x, y, z]); c = addV(at(row.r, u), u);
    }
    for (let i = 0; i + 1 < V.ids.length; i++) emit(V.ids[i], V.ids[i + 1], c);
    stat.caps++;
  };
  // 闭环收口 (晕面在这里到了头: 离眼最远 / 最近处与球壳相切). 小的、圆的: 扇形收到中心; 大的 (晕面的顶很平, 最后一层的
  //   闭环还有几十度长 —— 多半是又扁又长的一条): 从相距最远的两个顶点处分成两半, 两半互相拉链 (三角形横跨窄的方向).
  //   不查三角形边长: 这一圈里面本来就该铺满
  const capLoop = (row, pc) => {
    const V = pieceVerts(row, pc), ring = V.ids[0] === V.ids[V.ids.length - 1] ? V.ids.slice(0, -1) : V.ids.slice(), n = ring.length;
    if (n < 3) return;
    let x = 0, y = 0, z = 0;
    for (const v of ring) { x += VU[3 * v]; y += VU[3 * v + 1]; z += VU[3 * v + 2]; }
    const c = unit([x, y, z]);
    let rad = 0;
    for (const v of ring) rad = Math.max(rad, angU(c, vu(v)));
    if (rad <= 0.5 * MAXEDGE || n < 6) { capTo(row, pc, false); return; }
    let bi = 0, bj = 1, bd = -1;
    const st = Math.max(1, Math.floor(n / 400));
    for (let i = 0; i < n; i += st) for (let j = i + 1; j < n; j += st) { const d = vang(ring[i], ring[j]); if (d > bd) { bd = d; bi = i; bj = j; } }
    const C1 = ring.slice(bi, bj + 1), C2 = [];
    for (let k = bi; k !== bj; k = (k - 1 + n) % n) C2.push(ring[k]);
    C2.push(ring[bj]);
    const f1 = cumFr(C1), f2 = cumFr(C2);
    for (let i = 0, j = 0; i < C1.length - 1 || j < C2.length - 1;) {
      if (j >= C2.length - 1 || (i < C1.length - 1 && f1[i + 1] <= f2[j + 1])) { if (C1[i] !== C2[j] && C1[i + 1] !== C2[j]) tris.push(C1[i], C1[i + 1], C2[j]); i++; }
      else { if (C1[i] !== C2[j] && C1[i] !== C2[j + 1]) tris.push(C1[i], C2[j + 1], C2[j]); j++; }
    }
    stat.caps++;
  };
  // 对不上的段 (相邻两层几乎重合时): 每个点找对面对不上的段上视方向最近的点 (LIM 之内), 同一段上脚点单调的一串一串拉链.
  //   ok(w, i, q, j) 可再筛候选 (同一层自缝时排除沿线的近邻). 对面一个候选也没有的点连成串, 放进 lone (调用方再处理)
  const LIM = 2.5 * GAPSKY + (stat.coarsened ? SKY : 0);
  const sewPass = (S, T, flip, ok, lone) => {
    S.forEach((w, wi) => {
      const ft = w.ids.map((v, i) => { let bq = -1, bj = -1, bd = LIM; T.forEach((x, q) => { for (let j = 0; j < x.ids.length; j++) { if (ok && !ok(wi, i, q, j)) continue; const d = vang(v, x.ids[j]); if (d < bd) { bd = d; bq = q; bj = j; } } }); return bq < 0 ? null : [bq, bj]; });
      // 截断出来的弧 (带截断切口) 上的孤点串标上 fromCut, 不自缝: 那是截断球面把轮廓一点点切短, 不是曲面折回
      //   (外界 r = 截断距离 + 灯距处曾把左右两段缝到一起)
      if (lone) { let run = null; const push = rn => { if (w.noSelf) rn.fromCut = true; lone.push(rn); }; for (let i = 0; i < ft.length; i++) { if (!ft[i]) (run || (run = [])).push(w.ids[i]); else if (run) { push(run); run = null; } } if (run) push(run); }
      for (let s0 = 0; s0 < ft.length;) {
        if (!ft[s0]) { s0++; continue; }
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
  const runLen = runs => { let s = 0; for (const w of runs) for (let i = 0; i + 1 < w.length; i++) s += vang(w[i], w[i + 1]); return s; };
  // 同一层里对不上、对面也没有可缝的点 (lone 串): 曲面在这一层与相邻一层之间到了离眼最远 (或最近) 处折回 —— 这一层的两条轮廓
  //   (或一条闭环的两侧) 在两层之间汇合消失, 中间那条带子是折回处. 把两侧直接缝上 (带子在这一层的球面上, 离真曲面不超过
  //   两层的半径差). 候选: 别的串上, 或同一串上沿线隔得远 (沿线弧长 > 3 倍直线角距, 且隔 3 个点以上) 的点. 返回仍没缝上的长度
  const selfSew = runs0 => {
    const runs = runs0.filter(w => !w.fromCut);
    let cutLen = runLen(runs0.filter(w => w.fromCut));
    if (!runs.length) return cutLen;
    const V = runs.map(ids => {
      const c = [0];
      for (let k = 1; k < ids.length; k++) c.push(c[k - 1] + vang(ids[k - 1], ids[k]));
      return { ids, c };
    });
    const ok = (a, i, b, j) => {
      if (a !== b) return true;
      if (Math.abs(i - j) < 3) return false;
      return Math.abs(V[a].c[i] - V[a].c[j]) > 3 * vang(V[a].ids[i], V[a].ids[j]);
    };
    const lone = [];
    sewPass(V, V, false, ok, lone);
    return runLen(lone) + cutLen;
  };
  const looseSew = (A, B, UA, UB) => {
    const cutArc = pc => pc.arc.sp.some(q => q.type === 'cut');
    const VA = UA.map(pc => Object.assign(pieceVerts(A, pc), { noSelf: cutArc(pc) })), VB = UB.map(pc => Object.assign(pieceVerts(B, pc), { noSelf: cutArc(pc) }));
    const loneA = [], loneB = [];
    sewPass(VA, VB, false, null, loneA); sewPass(VB, VA, true, null, loneB);
    // 只在对面一层这里什么都没有 (对不上的段全在这一侧) 时自缝: 两侧都有对不上的段是别的拓扑变化, 孤点只是离得远, 不是折回.
    //   两侧都有时, 孤点串扇形收到对面最近的那个点上: 典型是闭环在两层之间张开成弧 (晕面几乎与球壳相切, r 差 0.004 m 轮廓在天上
    //   挪好几度), 闭环上消失的那一段缩向边界线, 对面那层的弧端就在边界线上 (3-5-6-8 Lowitz Σ34.5 蓝边 2 r≈13.44:
    //   闭环上 15° 的一段消失, 原来天上留一道 1°×4° 的缝). 截断出来的弧不收 (同自缝)
    const fanTo = (runs, V) => {
      let unc = 0;
      for (const w of runs) {
        if (w.fromCut || w.length < 2) { unc += runLen([w]); continue; }
        let bv = -1, bd = 3 * LIM;
        for (const x of V) for (const v of x.ids) for (const q of w) { const d = vang(q, v); if (d < bd) { bd = d; bv = v; } }
        if (bv < 0) { unc += runLen([w]); continue; }
        for (let i = 0; i + 1 < w.length; i++) emit(w[i], w[i + 1], bv);
        stat.caps++;
      }
      return unc;
    };
    return (VB.length ? fanTo(loneA, VB) : selfSew(loneA)) + (VA.length ? fanTo(loneB, VA) : selfSew(loneB));
  };
  // ---- 鞍点 (晕面内部的拓扑变化: 两条轮廓在两层之间碰到一起、换了搭档 —— 闭环掐成两个 / 两条弧交换端点 / 弧上掐下
  //   一个闭环; 反过来也一样). 每层上挨得最近的那一对 "不相邻" 的点是双曲线的两个顶点, 两层的四个点围着鞍点.
  //   各段在这些点处切开; 切出来的小段两两对应成 "臂" (一端是同一条特殊线、另一端是鞍点, 或两端都是鞍点), 每条臂照常
  //   拉链; 四条臂在中间围出的四边形 (两层各两个点, 交替) 补两个三角形. 这样铺出来的面不重叠、不留洞
  //   (就近缝 looseSew: 两个方向各缝一遍, 互相重叠, 中间还剩一个洞). 认不出来 (不是鞍点) → false, 调用方退回就近缝
  const saddleSew = (A, B, UA, UB) => {
    const why = c => { if (log) log(`saddle? r ${A.r.toFixed(5)}..${B.r.toFixed(5)}: no (${c}), pieces ${UA.length}|${UB.length}`); return false; };
    if (!UA.length || !UB.length || UA.length + UB.length > 8 || UA.some(pc => pc.virt) || UB.some(pc => pc.virt)) return why(1);
    const LIMS = 2.5 * GAPSKY, cl = Math.cos(LIMS);
    const mk = (row, U) => U.map(pc => { const V = pieceVerts(row, pc); return { pc, V, cyc: V.ids.length > 2 && V.ids[0] === V.ids[V.ids.length - 1], cuts: [] }; });
    const SA = mk(A, UA), SB = mk(B, UB);
    if (SA.some(s => !s.V.ts) || SB.some(s => !s.V.ts)) return why(2);
    // 一层上最近的一对不相邻的点 (不同的弧上; 或同一条弧上沿线隔得比直线距离远 3 倍以上). near: 只在这些方向的 LIMS 之内找;
    //   used: 已经定下的鞍点附近不再找
    const closest = (S, near, used) => {
      const C = [];
      S.forEach((s, p) => {
        const n = s.V.ids.length - (s.cyc ? 1 : 0), al = s.pc.arc.len;
        for (let i = 0; i < n; i++) {
          const u = vu(s.V.ids[i]);
          if (near && !near.some(q => dot(u, q) > cl)) continue;
          if (used.some(q => dot(u, q) > cl)) continue;
          C.push({ p, i, id: s.V.ids[i], u, arc: s.pc.arc, t: s.pc.arc.closed ? ((s.V.ts[i] % al) + al) % al : s.V.ts[i] });
        }
      });
      let best = null, bd = LIMS, cb = cl;
      for (let x = 0; x < C.length; x++) for (let y = x + 1; y < C.length; y++) {
        const a = C[x], b = C[y], c = dot(a.u, b.u);
        if (c <= cb) continue;
        const d = c >= 1 ? 0 : Math.acos(c);
        if (a.arc === b.arc) {
          let al = Math.abs(a.t - b.t);
          if (a.arc.closed) al = Math.min(al, a.arc.len - al);
          if (al <= 3 * d || al < 2 * SKY) continue;
        }
        bd = d; cb = c; best = { a, b, d };
      }
      return best;
    };
    const sites = [], usedA = [], usedB = [];
    for (let it = 0; it < 4; it++) {
      const a = closest(SA, null, usedA);
      if (!a) break;
      const b = closest(SB, [a.a.u, a.b.u], usedB);
      usedA.push(a.a.u, a.b.u);
      if (!b) continue;
      if (Math.max(angU(a.a.u, b.a.u), angU(a.a.u, b.b.u), angU(a.b.u, b.a.u), angU(a.b.u, b.b.u)) > LIMS) continue;
      usedB.push(b.a.u, b.b.u);
      const k = sites.length;
      sites.push({ m: [a.a.id, a.b.id], n: [b.a.id, b.b.id] });
      SA[a.a.p].cuts.push([a.a.i, k]); SA[a.b.p].cuts.push([a.b.i, k]);
      SB[b.a.p].cuts.push([b.a.i, k]); SB[b.b.p].cuts.push([b.b.i, k]);
    }
    if (!sites.length) return why(3);
    // 切段: 每段 { pc (给 strip 用的小段), d0, d1 (两端的身份: 'c' + 链号 / 's' + 鞍点号), v0, v1 (两端的顶点) }
    const subsOf = S => {
      const out = [];
      for (const s of S) {
        const pc = s.pc, V = s.V, n = V.ids.length;
        if (!s.cuts.length) return null;                                   // 有一段不沾鞍点: 不是单纯的鞍点
        const cs = s.cuts.slice().sort((x, y) => x[0] - y[0]);
        for (let k = 1; k < cs.length; k++) if (cs[k][0] === cs[k - 1][0]) return null;
        const sub = (i0, i1, t0, t1, e0, e1) => {                          // e = [身份, 链号 (给 strip 找边界线) 或 null, 特殊点]
          const q = { arc: pc.arc, t0, t1, len: t1 - t0, loop: false, k0: e0[1] === null ? uniq-- : e0[1], k1: e1[1] === null ? uniq-- : e1[1], q0: e0[2], q1: e1[2] };
          q.mid = uOn(pc.arc, 0.5 * (t0 + t1));
          out.push({ pc: q, d0: e0[0], d1: e1[0], v0: V.ids[i0], v1: V.ids[i1] });
        };
        const cutE = c => ['s' + c[1], null, null];
        if (s.cyc && pc.loop) {                                            // 闭环: k 个切点 → k 段 (一个切点: 整圈, 从切点回到切点)
          for (let k = 0; k < cs.length; k++) {
            const a = cs[k], b = cs[(k + 1) % cs.length], wrap = k + 1 >= cs.length;
            sub(a[0], b[0], V.ts[a[0]], V.ts[b[0]] + (wrap ? pc.arc.len : 0), cutE(a), cutE(b));
          }
          continue;
        }
        if (cs[0][0] <= 0 || cs[cs.length - 1][0] >= n - 1) return null;   // 切点落在段的端点上
        const head = [typeof pc.k0 === 'number' && pc.k0 >= 0 ? 'c' + pc.k0 : 'u' + (uniq--), typeof pc.k0 === 'number' && pc.k0 >= 0 ? pc.k0 : null, pc.q0 || null];
        const tail = [typeof pc.k1 === 'number' && pc.k1 >= 0 ? 'c' + pc.k1 : 'u' + (uniq--), typeof pc.k1 === 'number' && pc.k1 >= 0 ? pc.k1 : null, pc.q1 || null];
        sub(0, cs[0][0], pc.t0, V.ts[cs[0][0]], head, cutE(cs[0]));
        for (let k = 0; k + 1 < cs.length; k++) sub(cs[k][0], cs[k + 1][0], V.ts[cs[k][0]], V.ts[cs[k + 1][0]], cutE(cs[k]), cutE(cs[k + 1]));
        sub(cs[cs.length - 1][0], n - 1, V.ts[cs[cs.length - 1][0]], pc.t1, cutE(cs[cs.length - 1]), tail);
      }
      return out;
    };
    const XA = subsOf(SA), XB = subsOf(SB);
    if (!XA || !XB) return why(4);
    // 臂: 两端身份相同的一对 (多个候选取中点最近的); 走向: 身份不同的两端直接定, 相同的按四分点
    const gapOf = (pa, pb, rev) => { let g = 0; for (let k = 0; k <= 8; k++) g = Math.max(g, angU(uOn(pa.arc, pa.t0 + k / 8 * pa.len), uOn(pb.arc, pb.t0 + (rev ? 1 - k / 8 : k / 8) * pb.len))); return g; };
    const same = (a, b) => (a.d0 === b.d0 && a.d1 === b.d1) || (a.d0 === b.d1 && a.d1 === b.d0);
    const revOf = (a, b) => {
      if (a.d0 !== a.d1) return a.d0 !== b.d0;
      const cst = r2 => { let c = 0; for (const f of [0.25, 0.5, 0.75]) c += angU(uOn(a.pc.arc, a.pc.t0 + f * a.pc.len), uOn(b.pc.arc, b.pc.t0 + (r2 ? 1 - f : f) * b.pc.len)); return c; };
      return cst(true) < cst(false);
    };
    const arms = [], usedX = new Set();
    for (const a of XA) {
      let best = null, bd = Infinity;
      for (const b of XB) {
        if (usedX.has(b) || !same(a, b)) continue;
        const d = angU(a.pc.mid, b.pc.mid);
        if (d < bd) { bd = d; best = b; }
      }
      if (!best) continue;
      const rev = revOf(a, best);
      if (gapOf(a.pc, best.pc, rev) > MAXEDGE) return why(6);
      usedX.add(best); usedX.add(a);
      arms.push({ a, b: best, rev });
    }
    if (!arms.length) return why(5);
    // 没配上的小段. 同一层上两条、两端是同样的两个鞍点 (s_k, s_l, k ≠ l): 晕面的一道折回脊斜穿过这两层 —— 这一层上两条轮廓
    //   中间那一截面对面挨着, 到对面那层已经汇合消失; 两截互相拉链. 两端是同一个鞍点的小段 (对面那层在鞍点处的发夹弯):
    //   留着, 它与各臂在鞍点处围成的小多边形由最后的补洞补上. 别的 (带着特殊线的端却没有对应) → 认不出来
    const ridges = [];
    for (const [X, row] of [[XA, A], [XB, B]]) {
      const left = X.filter(x => !usedX.has(x));
      for (let i = 0; i < left.length; i++) {
        const a = left[i];
        if (usedX.has(a)) continue;
        if (a.d0[0] !== 's' || a.d1[0] !== 's') return why(4);
        if (a.d0 === a.d1) { if (a.pc.len > 2 * LIMS) return why(7); continue; }
        const b = left.find((x, j) => j > i && !usedX.has(x) && same(a, x));
        if (!b) return why(7);
        const rev = a.d0 !== b.d0;
        if (gapOf(a.pc, b.pc, rev) > MAXEDGE) return why(6);
        usedX.add(a); usedX.add(b);
        ridges.push({ row, a, b, rev });
      }
    }
    for (const { a, b, rev } of arms) strip(A, B, a.pc, b.pc, rev);
    for (const { row, a, b, rev } of ridges) {
      const VX = pieceVerts(row, a.pc);
      let VY = pieceVerts(row, b.pc);
      if (rev) VY = revV(VY);
      if (VX.ids.length >= 2 && VY.ids.length >= 2) zipFr({ ids: VX.ids, fr: cumFr(VX.ids) }, { ids: VY.ids, fr: cumFr(VY.ids) });
    }
    // 每个鞍点: 四条臂的鞍点端连出的四条边正好是 m0-n0, m0-n1, m1-n0, m1-n1 时, 中间的四边形补两个三角形
    //   (不是这种的 —— 有折回脊 / 发夹弯 —— 留给最后的补洞)
    for (let k = 0; k < sites.length; k++) {
      const { m, n } = sites[k], got = new Set();
      let ok = true;
      for (const { a, b, rev } of arms) for (const e of [0, 1]) {
        if ((e ? a.d1 : a.d0) !== 's' + k) continue;
        const va = e ? a.v1 : a.v0, vb = (e ? !rev : rev) ? b.v1 : b.v0;
        const i = m.indexOf(va), j = n.indexOf(vb);
        if (i < 0 || j < 0) ok = false; else got.add(i * 2 + j);
      }
      if (!ok || got.size !== 4) continue;
      if (vang(m[0], m[1]) <= vang(n[0], n[1])) { emit(m[0], n[0], m[1]); emit(m[0], m[1], n[1]); }
      else { emit(n[0], m[0], n[1]); emit(n[0], n[1], m[1]); }
    }
    stat.saddles = (stat.saddles || 0) + sites.length;
    return true;
  };
  for (let k = 0; k + 1 < rows.length; k++) {
    const A = rows[k], B = rows[k + 1];
    if (!A.arcs.length && !B.arcs.length) continue;
    const an = analyze(A, B);
    for (const [pa, pb, rev, trim] of an.pairs) strip(A, B, pa, pb, rev, trim);
    for (const { row, pc } of an.lens) lensZip(row, pc);
    if (!an.unA.length && !an.unB.length) continue;
    // 灯两侧那一对层: 晕面若过灯, 灯旁的段缩向灯的方向 (锥顶) —— 对不上的、靠近灯方向 (20° 内) 的段收到灯上
    if (A.kind === 'lampLo' && B.kind === 'lampHi') {
      const ul = [Dx / D, Dy / D, Dz / D], UA = [], UB = [];
      for (const [row, list, out] of [[A, an.unA, UA], [B, an.unB, UB]]) for (const pc of list) {
        const V = pieceVerts(row, pc);
        if (V.ids.every(v => angU(vu(v), ul) < 20 * DEG)) capTo(row, pc, true); else out.push(pc);
      }
      // 剩下的: 多半是两个闭环在灯这个点上碰着 (另一侧是一个绕过灯的闭环) —— 就是掐断点在灯上的鞍点
      if ((UA.length || UB.length) && !saddleSew(A, B, UA, UB)) { const unc = looseSew(A, B, UA, UB); if (unc > 1.5 * SKY) { stat.loose++; stat.looseDeg += unc / DEG; } }
      continue;
    }
    stat.events++;
    if (log) log(`event r ${A.r.toFixed(5)} .. ${B.r.toFixed(5)}: pieces ${an.unA.length} | ${an.unB.length} unmatched`);
    // 小段 (正在缩没 / 刚出生) 收口, 其余就近缝
    const tiny = pc => pc.len < Math.min(CAPMAX, 30 * Math.sqrt(2 * Math.max(EPSR, B.r / A.r - 1)));
    const UA = [], UB = [];
    for (const pc of an.unA) { if (tiny(pc)) capTo(A, pc, false); else UA.push(pc); }
    for (const pc of an.unB) { if (tiny(pc)) capTo(B, pc, false); else UB.push(pc); }
    // 一侧只剩一个闭环、对面什么也没有: 晕面像穹顶一样收成一点 (离眼最远 / 最近处与球壳在一点相切), 闭环收成扇形.
    //   (两个闭环同时消失 —— 一里一外 —— 是一圈折回, 留给下面的同层自缝; 单个闭环对边的沿线距离只有直线的 π/2 倍, 自缝不缝它)
    //   3-5-6-8 Parry Σ34.5: r≈14.43 处最后一层还剩 8.8° 的闭环, 原来中间留一个五边形的洞
    //   几个闭环同时消失: 互不相套的 (几个互不相干的尖顶 / 几条并排的扁长条) 各收各的; 套在一起的才是一圈折回
    const domes = (row, U) => {
      // 点 c 在不在闭环 ring 里面: 绕 c 的方位角转过的圈数
      const inside = (c, ring) => {
        const e1 = unit(Math.abs(c[2]) < 0.9 ? [-c[1], c[0], 0] : [0, -c[2], c[1]]), e2 = [c[1] * e1[2] - c[2] * e1[1], c[2] * e1[0] - c[0] * e1[2], c[0] * e1[1] - c[1] * e1[0]];
        let w = 0, prev = null;
        for (let k = 0; k <= ring.length; k++) {
          const u = vu(ring[k % ring.length]), x = dot(u, e1), y = dot(u, e2);
          if (x === 0 && y === 0) return true;
          const b = Math.atan2(y, x);
          if (prev !== null) { let d = b - prev; if (d > Math.PI) d -= 2 * Math.PI; else if (d < -Math.PI) d += 2 * Math.PI; w += d; }
          prev = b;
        }
        return Math.abs(w) > Math.PI;
      };
      const info = U.map(pc => {
        if (!pc.loop) return null;
        const V = pieceVerts(row, pc), ring = V.ids[0] === V.ids[V.ids.length - 1] ? V.ids.slice(0, -1) : V.ids;
        let x = 0, y = 0, z = 0;
        for (const v of ring) { x += VU[3 * v]; y += VU[3 * v + 1]; z += VU[3 * v + 2]; }
        return { c: unit([x, y, z]), ring };
      });
      return U.filter((pc, i) => {
        const a = info[i];
        if (!a) return true;
        // 套在一起的 (一个的中心落在另一个里面): 一圈折回, 不单独收
        for (let j = 0; j < U.length; j++) if (j !== i && info[j] && (inside(a.c, info[j].ring) || inside(info[j].c, a.ring))) return true;
        capLoop(row, pc);
        return false;
      });
    };
    if (!UB.length && UA.length) { const U2 = domes(A, UA); UA.length = 0; UA.push(...U2); }
    if (!UA.length && UB.length) { const U2 = domes(B, UB); UB.length = 0; UB.push(...U2); }
    if ((UA.length || UB.length) && !saddleSew(A, B, UA, UB)) {
      const unc = looseSew(A, B, UA, UB);
      if (unc > 1.5 * SKY) { stat.loose++; stat.looseDeg += unc / DEG; }
    }
  }

  // ---- 特殊线输出 (折叠线、焦线): 各链按半径排的点
  const chains = new Map();
  rows.forEach((row, ri) => {
    for (const q of row.sp) {
      if (q.type !== 'fold' && q.type !== 'focal' && !q.foldAt) continue;
      let c = chains.get(q.chain);
      if (!c) { c = { kind: q.type === 'focal' ? 'focal' : 'fold', pts: [], dirs: [], k0: ri, q0: q }; chains.set(q.chain, c); }
      c.pts.push(...q.p); c.dirs.push(...q.u); c.k1 = ri; c.q1 = q;
    }
  });
  // 焦线的头: 焦线在离眼最远 (最近) 处掉头时, 两条焦线停在同一层上, 中间那一小截夹在这一层与下一层之间 (层够不到).
  //   引擎的边界线里若有一条与焦线重合 (边界引擎把焦线当成接缝线 / 外缘), 沿它把两条焦线接成一条
  if (LN.length) {
    const fc = [...chains.entries()].filter(e => e[1].kind === 'focal');
    const tolL = r => Math.max(1e-6, 2e-3 * r), labOf = q => ({ lab: nearLines(q.p, k => k !== 'cusp', tolL(Math.hypot(q.p[0] - E[0], q.p[1] - E[1], q.p[2] - E[2]))) });
    const rev3 = a => { const o = []; for (let i = a.length - 3; i >= 0; i -= 3) o.push(a[i], a[i + 1], a[i + 2]); return o; };
    for (let x = 0; x < fc.length; x++) for (let y = x + 1; y < fc.length; y++) {
      const [ka, a] = fc[x], [kb, b] = fc[y];
      if (a.gone || b.gone) continue;
      for (const tail of [true, false]) {
        const ra = tail ? a.k1 : a.k0, rb = tail ? b.k1 : b.k0, nx = tail ? ra + 1 : ra - 1;
        if (ra !== rb || nx < 0 || nx >= rows.length) continue;
        const Q = linePath(labOf(tail ? a.q1 : a.q0), labOf(tail ? b.q1 : b.q0), rows[ra].r, rows[nx].r);
        if (Q === null || !Q.length) continue;
        const qp = [], qd = [];
        for (const p of Q) { qp.push(p[0], p[1], p[2]); qd.push(...unit([p[0] - E[0], p[1] - E[1], p[2] - E[2]])); }
        // a (头 → 尾) + 线 + b (尾 → 头); 接在头上的反过来
        const pts = tail ? [...a.pts, ...qp, ...rev3(b.pts)] : [...rev3(a.pts), ...qp, ...b.pts], dirs = tail ? [...a.dirs, ...qd, ...rev3(b.dirs)] : [...rev3(a.dirs), ...qd, ...b.dirs];
        a.pts = pts; a.dirs = dirs; b.gone = true;
        if (tail) { a.k1 = b.k0; a.q1 = b.q0; } else { a.k0 = b.k1; a.q0 = b.q1; const t0 = a.k0; a.k0 = a.k1; a.k1 = t0; const tq = a.q0; a.q0 = a.q1; a.q1 = tq; }
        stat.focalJoins = (stat.focalJoins || 0) + 1;
        break;
      }
    }
    for (const [k, c] of fc) if (c.gone) chains.delete(k);
  }
  // 折叠线只画延续得够长的: 径向跨度不到所在半径 15% 的 (且两头都不碰这一段的半径边界 —— 碰到的可能接着下一段) 不画.
  //   这类短链多是球壳与晕面相切 (鞍点) 处的假折点: 两条交叉的轮廓在切点重新连成 )(, 刚分开时双曲线顶点尖得像折点,
  //   半径稍一变就钝了 (3-4-5 Lowitz 入射面: r 3.158–3.535, 天空里却横扫 30°). 它们仍是网格的顶点, 只是不当成晕面的折痕画出来
  const rLo = rows.length ? rows[0].r : 0, rHi = rows.length ? rows[rows.length - 1].r : 0;
  // 链在天空里相邻两步的方向差 > 45° 处切开 (两条挨着的折叠线在灯前汇成一个折点时, 链会从一条跳到另一条), 取最长的一截
  const kinkCut = c => {
    const n = c.dirs.length / 3;
    if (n < 4) return c;
    const st = [];
    for (let i = 0; i + 1 < n; i++) st.push([c.dirs[3 * i + 3] - c.dirs[3 * i], c.dirs[3 * i + 4] - c.dirs[3 * i + 1], c.dirs[3 * i + 5] - c.dirs[3 * i + 2]]);
    const cuts = [0];
    let prev = -1;
    for (let i = 0; i < st.length; i++) {
      const l = Math.hypot(...st[i]);
      if (l < 1e-6) continue;
      if (prev >= 0 && dot(st[prev], st[i]) < 0.7 * Math.hypot(...st[prev]) * l) cuts.push(i + 1);
      prev = i;
    }
    if (cuts.length === 1) return c;
    cuts.push(n);
    let b = 0;
    for (let k = 1; k + 1 < cuts.length; k++) if (cuts[k + 1] - cuts[k] > cuts[b + 1] - cuts[b]) b = k;
    return { kind: c.kind, pts: c.pts.slice(3 * cuts[b], 3 * cuts[b + 1]), dirs: c.dirs.slice(3 * cuts[b], 3 * cuts[b + 1]) };
  };
  const special = [];
  for (const c0 of chains.values()) {
    const c = c0.kind === 'fold' ? kinkCut(c0) : c0;
    if (c.pts.length < 6) continue;
    if (c.kind === 'fold') {
      const n = c.pts.length / 3, ra = Math.hypot(c.pts[0] - E[0], c.pts[1] - E[1], c.pts[2] - E[2]), rb = Math.hypot(c.pts[3 * n - 3] - E[0], c.pts[3 * n - 2] - E[1], c.pts[3 * n - 1] - E[2]);
      const r0 = Math.min(ra, rb), r1 = Math.max(ra, rb), edge = r0 <= rLo * (1 + 1e-6) || r1 >= rHi * (1 - 1e-6);
      if (!edge && r1 - r0 < 0.15 * 0.5 * (r0 + r1)) { stat.foldsHidden = (stat.foldsHidden || 0) + 1; continue; }
    }
    special.push({ kind: c.kind, p: Float32Array.from(c.pts), d: Float32Array.from(c.dirs) });
    if (c.kind === 'fold') stat.folds++; else stat.focal++;
  }

  // ---- 自检与补洞 ("贴边检验"). 把所有三角形的边数一遍: 只属于一个三角形的边 (开口边) 应该都是晕面的真边界 —— 两端都是
  //   末端点 / 截断点 / 边界线上的点, 或都在这一段的最里 (最外) 一层上 (与相邻段相接). 其余的开口边是缺口: 拓扑变化处就近缝
  //   剩下的小洞 (鞍点: 两条弧在两层之间换了搭档, 中间剩一个四边形 …). 围成小闭环的 (或两头在边界上、隔得很近的一小串)
  //   补上 (最小面积三角化); 补不了的记进 stat.openDeg (天空里的总长, 度), 页面与测试拿它判这一段网格是否完整
  stat.openN = 0; stat.openDeg = 0; stat.holesFilled = 0;
  {
    // 位置相同的顶点并成一个 (沿边界线插的点, 相邻几条带各加了一次)
    const canon = new Int32Array(nv), km = new Map();
    for (let v = 0; v < nv; v++) {
      const k = VP[3 * v] + ',' + VP[3 * v + 1] + ',' + VP[3 * v + 2], c = km.get(k);
      if (c === undefined) { km.set(k, v); canon[v] = v; } else { canon[v] = c; VF[c] |= VF[v]; }
      // 截断球面上 (离灯 ≥ 截断距离) 与灯旁的顶点也算在边界上: 截断处晕面的边是截断球面的切口, 灯是锥顶
      const dl = Math.hypot(VP[3 * v] - L[0], VP[3 * v + 1] - L[1], VP[3 * v + 2] - L[2]);
      if (dl >= maxD * (1 - 2e-3) || dl < 3e-3 * D) VF[canon[v]] |= 1;
    }
    const K = 4194304, cnt = new Map();
    const add = (u, v) => { const k = u < v ? u * K + v : v * K + u; cnt.set(k, (cnt.get(k) || 0) + 1); };
    for (let i = 0; i < tris.length; i += 3) {
      const a = canon[tris[i]], b = canon[tris[i + 1]], c = canon[tris[i + 2]];
      if (a === b || b === c || c === a) continue;
      add(a, b); add(b, c); add(c, a);
    }
    const adj = new Map();
    const nb = (u, v) => { const l = adj.get(u); if (l) l.push(v); else adj.set(u, [v]); };
    for (const [k, n] of cnt) {
      if (n !== 1) continue;
      const a = Math.floor(k / K), b = k % K;
      if ((VF[a] & 1 && VF[b] & 1) || (VF[a] & VF[b] & 6)) continue;
      nb(a, b); nb(b, a);
    }
    const fillLoop = seq => {                              // 最小面积三角化 (动态规划; 顶点数 ≤ 28)
      const n = seq.length, X = seq.map(v => [VP[3 * v], VP[3 * v + 1], VP[3 * v + 2]]);
      const area = (i, j, k) => {
        const ux = X[j][0] - X[i][0], uy = X[j][1] - X[i][1], uz = X[j][2] - X[i][2], vx = X[k][0] - X[i][0], vy = X[k][1] - X[i][1], vz = X[k][2] - X[i][2];
        return Math.hypot(uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx);
      };
      const C = new Float64Array(n * n), Kx = new Int32Array(n * n).fill(-1);
      for (let len = 2; len < n; len++) for (let i = 0; i + len < n; i++) {
        const j = i + len;
        let best = Infinity, bk = -1;
        for (let k = i + 1; k < j; k++) { const c = C[i * n + k] + C[k * n + j] + area(i, k, j); if (c < best) { best = c; bk = k; } }
        C[i * n + j] = best; Kx[i * n + j] = bk;
      }
      const rec = (i, j) => { if (j - i < 2) return; const k = Kx[i * n + j]; tris.push(seq[i], seq[k], seq[j]); rec(i, k); rec(k, j); };
      rec(0, n - 1);
    };
    const seen = new Set();
    for (const v0 of adj.keys()) {
      if (seen.has(v0)) continue;
      const comp = [], st = [v0];
      seen.add(v0);
      while (st.length) { const v = st.pop(); comp.push(v); for (const w of adj.get(v)) if (!seen.has(w)) { seen.add(w); st.push(w); } }
      const deg1 = [];
      let simple = true, ne = 0;
      for (const v of comp) { const d = adj.get(v).length; ne += d; if (d === 1) deg1.push(v); else if (d !== 2) simple = false; }
      let loop = null;
      if (simple && comp.length >= 3 && comp.length <= 28 && (deg1.length === 0 || deg1.length === 2)) {
        const start = deg1.length ? deg1[0] : comp[0], seq = [start];
        for (let prev = -1, cur = start; seq.length <= comp.length;) {
          const l = adj.get(cur), nx = l[0] !== prev ? l[0] : l[1];
          if (nx === undefined || nx === start) break;
          seq.push(nx); prev = cur; cur = nx;
        }
        if (seq.length === comp.length) {
          let diam = 0;
          for (let i = 0; i < seq.length; i++) for (let j = i + 1; j < seq.length; j++) diam = Math.max(diam, vang(seq[i], seq[j]));
          // 一小串 (不闭合): 两头在边界上、隔得很近才补 (补上的那条弦成了新的边界边)
          if (diam <= MAXEDGE && (deg1.length === 0 || (VF[seq[0]] & 1 && VF[seq[seq.length - 1]] & 1 && vang(seq[0], seq[seq.length - 1]) <= 2.5 * GAPSKY))) loop = seq;
        }
      }
      if (loop) { fillLoop(loop); stat.holesFilled++; continue; }
      if (log) { let l = 0; for (const v of comp) for (const w of adj.get(v)) if (v < w) l += vang(v, w) / DEG; const v = comp[0], rr = Math.hypot(VP[3 * v] - E[0], VP[3 * v + 1] - E[1], VP[3 * v + 2] - E[2]); log(`open: ${ne / 2} edges ${l.toFixed(2)}° at r ${rr.toFixed(4)} sky ${(Math.atan2(VU[3 * v + 1], VU[3 * v]) / DEG).toFixed(2)}/${(Math.asin(VU[3 * v + 2]) / DEG).toFixed(2)}, ${comp.length} verts, deg1 ${deg1.length}, simple ${simple}`); }
      stat.openN += ne / 2;
      for (const v of comp) for (const w of adj.get(v)) if (v < w) stat.openDeg += vang(v, w) / DEG;
    }
  }

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
  if (log) log('stat ' + JSON.stringify(stat));
  return { pos, dir, h: hh, nrm, nTri: nT, lines, special, stat, rowR: rows.filter(r => r.arcs.length).map(r => r.r) };
}
