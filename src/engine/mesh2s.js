// 二自由度晕面的三角网格 —— 球壳法: 层 = 以眼为心的球壳 (发散度球壳, ρ = r / D) 与晕面的交线, 相邻两层在视方向里对应着缝.
//
// 为什么不用光锥层 (灯处仰角 h 的等值线, stitch2.js): 那一族切面以灯下的竖直线为轴, 眼睛是它的奇点 —— 光锥过眼睛的那一层
//   (h = Σ) 前后, 层线贴着眼睛的一截整段缩成一点 (离眼的距离 ∝ |h − Σ|), 相邻两层在这里谈不上对应; 而观测者天空正是
//   从这一点看出去的, 三维里几厘米的错连在天空里横跨几十度. 球壳以眼睛为心: 晕面在眼睛处是一个锥 (平行光晕弧上的锥),
//   每层球壳截出来的都是同一形状的线, 只是按半径缩放 —— 一直到 r → 0 都是正则的, 最里一层在天空里就是平行光弧.
//   同一层上 "位置" 与 "视方向" 一一对应 (都在同一个球面上), 所以层间的对应、距离、容差全部按角度量, 没有随尺度变的长度阈值.
//
// 输入: shellAt(r) → { arcs: [{ p, d, closed }], cusps: [[x, y, z] …] }
//   p / d = 这一层上一条弧各点的位置 / 单位视方向 (长度 3n 的数组, 或 [x, y, z] 的数组); 弧的两端是引擎给的精确末端点
//   (落在入射掠射 / 出射临界 … 边界上), closed = 闭合环; cusps = 这一层上的尖点. 引擎: divshell.js 的 computeShell2Dof.
// 做法:
//   1) 半径按对数等距取一组初始层; 每层的弧按天空张角抽稀 (≥ skyStep 留一个点, 转角大的、末端点、尖点总留).
//   2) 相邻两层配对: 先配末端点 (视方向互为最近, 且明显近于第二近的 —— 不含糊才配), 两端都配上的弧一段对一段;
//      闭合环按中心方向与长度配. 配上的两段按 "较短对角线" 拉链 (距离 = 视方向夹角), 两层配上的尖点是必经的横档
//      (尖点脊线成网格的边).
//   3) 配不全 (两层之间弧的条数 / 连法变了: 一条弧缩没了、两条在边界上接成一条、鞍点处换了连法 …): 在两层之间二分半径,
//      把变化夹到 r 的相对宽度 < eps 之内 (每步多算一层球壳, 不进网格), 只把紧挨变化前后的两层加进来. 这两层几乎重合,
//      之间配得上的照缝, 配不上的: 短的 (正在缩没的小弧 / 小环) 向自己的中心收口, 其余按逐点最近的对应分段缝.
//   4) 配全了但两层在天空里隔得远 (> skyGap) 就在中间加层, 直到够密或到行数上限.
//   灯在 r = D 的球壳上是晕面的另一个锥顶: 那里层线在天空里缩向灯的方向, 由 4) 自动加密.
//   最大距离截断 (离灯 maxD): 弧在截断球面处截断, 截口当末端点.
// 输出: { pos, dir, h, nrm (每个三角形三个顶点依次排), nTri, lines: [{ r, gap, curves: [{ p, d }] }] (各层的线, 画层线用), stat, rowR }
//   h = 各顶点的灯处仰角 (按位置算, 着色用).

const DEG = Math.PI / 180;

export function meshShells(shellAt, L, E, opt = {}) {
  const t0 = Date.now();
  const Dx = L[0] - E[0], Dy = L[1] - E[1], Dz = L[2] - E[2], D = Math.hypot(Dx, Dy, Dz);
  const maxD = opt.maxD || 80, rMin = opt.rMin || 1e-3 * D, rMax = opt.rMax || maxD + D;
  const N0 = opt.rows || 72, EPS = opt.eps || 1e-5, MAXROWS = opt.maxRows || 900, MAXCALLS = opt.maxCalls || 1600;
  // 预算 (顶点数): 天空里铺得大的曲面每层的线有几百度长, 按固定的 0.35° 取点、1.2° 加层会到几十万个三角形.
  //   先算 7 层看线有多长, 步长放宽到 "平均每层的点数 × 约 200 层 ≤ 预算", 层间距跟着步长走 (≥ 3 倍)
  const VBUDGET = opt.budget || 90000;
  let SKY = (opt.skyStep || 0.35) * DEG;
  {
    let tot = 0, cnt = 0;
    for (let k = 0; k < 7; k++) {
      const s0 = shellAt(rMin * Math.pow(rMax / rMin, (k + 0.5) / 7)) || { arcs: [] };
      let len = 0;
      for (const arc of s0.arcs || []) {
        const flat = !Array.isArray(arc.d[0]), n = flat ? arc.d.length / 3 : arc.d.length;
        for (let i = 1; i < n; i++) {
          const a = flat ? [arc.d[3 * i - 3], arc.d[3 * i - 2], arc.d[3 * i - 1]] : arc.d[i - 1], b = flat ? [arc.d[3 * i], arc.d[3 * i + 1], arc.d[3 * i + 2]] : arc.d[i];
          len += Math.acos(Math.max(-1, Math.min(1, a[0] * b[0] + a[1] * b[1] + a[2] * b[2])));
        }
      }
      if (len > 0) { tot += len; cnt++; }
    }
    if (cnt) SKY = Math.min(2 * DEG, Math.max(SKY, (tot / cnt) * 200 / VBUDGET));
  }
  const GAPSKY = Math.max((opt.skyGap || 1.2) * DEG, 3 * SKY), TOLN = 10 * DEG, CAPMAX = 1 * DEG;
  // 半径方向也要够密: 天空里张角很小的曲面 (几乎顺着视线铺开), 天空间距判据不触发, 三维里相邻两层却隔得远, 曲率被长条三角形抹平
  const RSTEP = opt.rStep || Infinity;          // 默认不用 (试过 5%: 不解决问题, 三角形数翻几倍)
  const MAXEDGE = Math.max(4 * GAPSKY, 5 * DEG), UNCOK = 1.5 * SKY;
  const PINL = opt.pinLamp || 1e-3;           // 灯所在的球壳 r = D 两侧各留这么宽 (相对) 不取层: 晕面在灯上是锥顶, 层线在这里线性地缩成一点
  const cancelled = opt.cancelled || null;
  // ---- 顶点表
  let VP = new Float64Array(3 * 65536), VU = new Float64Array(3 * 65536), nv = 0;
  const addV = (x, y, z, ux, uy, uz) => {
    if (3 * nv + 3 > VP.length) { const a = new Float64Array(2 * VP.length); a.set(VP); VP = a; const b = new Float64Array(2 * VU.length); b.set(VU); VU = b; }
    VP[3 * nv] = x; VP[3 * nv + 1] = y; VP[3 * nv + 2] = z; VU[3 * nv] = ux; VU[3 * nv + 1] = uy; VU[3 * nv + 2] = uz;
    return nv++;
  };
  const ang = (a, b) => { const c = VU[3 * a] * VU[3 * b] + VU[3 * a + 1] * VU[3 * b + 1] + VU[3 * a + 2] * VU[3 * b + 2]; return c >= 1 ? 0 : c <= -1 ? Math.PI : Math.acos(c); };
  const angU = (u, v) => { const c = u[0] * v[0] + u[1] * v[1] + u[2] * v[2]; return c >= 1 ? 0 : c <= -1 ? Math.PI : Math.acos(c); };
  const tris = [];
  // 天空里边长超过 MAXEDGE 的三角形不要: 层间对应都是近邻, 正常的三角形远小于它; 宁可留洞 (覆盖检查看得见), 不连远
  let nSkip = 0;
  const emit = (a, b, c) => {
    if (a === b || b === c || c === a) return;
    if (ang(a, b) > MAXEDGE || ang(b, c) > MAXEDGE || ang(c, a) > MAXEDGE) { nSkip++; return; }
    tris.push(a, b, c);
  };

  // ---- 一层: 引擎 → 抽稀 → 段 (两端是末端点的一条弧 / 一个闭合环) 与末端点
  let calls = 0;
  const cache = new Map();
  const prep = r => {
    const key = Math.round(Math.log(r) * 1e9);
    let row = cache.get(key);
    if (row) return row;
    if (cancelled && cancelled()) throw new Error('cancelled');
    calls++;
    const s = shellAt(r) || { arcs: [], cusps: [] };
    row = { r, segs: [], nodes: [], raw: [], lmin: Infinity, lmax: -Infinity };
    const cusp = (s.cusps || []).map(q => { const x = q[0] - E[0], y = q[1] - E[1], z = q[2] - E[2], l = Math.hypot(x, y, z) || 1; return [x / l, y / l, z / l]; });
    for (const arc of s.arcs || []) {
      const flat = !Array.isArray(arc.p[0]);
      let n = flat ? arc.p.length / 3 : arc.p.length, wrapped = false;
      if (n < 2) continue;
      // 绕了好几圈的弧: 刚出生的闭环在天空里比追踪步长还小时, 引擎认不出已经回到起点, 一直绕到步数上限
      //   (一自由度引擎报 stop = max_steps, 一个 0.3° 的小环绕了 4600°). 在第一次回到起点处截断, 当闭环
      {
        const dAt = i => flat ? [arc.d[3 * i], arc.d[3 * i + 1], arc.d[3 * i + 2]] : arc.d[i];
        const d0 = dAt(0), far0 = i => { const q = dAt(i); return Math.acos(Math.max(-1, Math.min(1, (q[0] * d0[0] + q[1] * d0[1] + q[2] * d0[2]) / (Math.hypot(q[0], q[1], q[2]) * Math.hypot(d0[0], d0[1], d0[2]) || 1)))); };
        let len = 0, mxStep = 0, mxFar = 0;
        for (let i = 1; i < n; i++) { const a0 = dAt(i - 1), a1 = dAt(i); const st = Math.acos(Math.max(-1, Math.min(1, (a0[0] * a1[0] + a0[1] * a1[1] + a0[2] * a1[2]) / (Math.hypot(...a0) * Math.hypot(...a1) || 1)))); len += st; if (st > mxStep) mxStep = st; }
        for (let i = 0; i < n; i += Math.max(1, n >> 9)) mxFar = Math.max(mxFar, far0(i));
        if (!arc.closed && mxFar > 0 && len > 4 * Math.PI * mxFar) {        // 周长超过直径 × 2π 的两倍: 至少绕了两圈
          let away = false;
          for (let i = 1; i < n; i++) {
            const f = far0(i);
            if (!away) { if (f > 0.5 * mxFar) away = true; continue; }
            if (f <= Math.max(1.5 * mxStep, 1e-9)) { n = i + 1; wrapped = true; break; }
          }
        }
      }
      const P = new Float64Array(3 * n), U = new Float64Array(3 * n), ok = new Uint8Array(n);
      for (let i = 0; i < n; i++) {
        const x = flat ? arc.p[3 * i] : arc.p[i][0], y = flat ? arc.p[3 * i + 1] : arc.p[i][1], z = flat ? arc.p[3 * i + 2] : arc.p[i][2];
        let ux = flat ? arc.d[3 * i] : arc.d[i][0], uy = flat ? arc.d[3 * i + 1] : arc.d[i][1], uz = flat ? arc.d[3 * i + 2] : arc.d[i][2];
        const ul = Math.hypot(ux, uy, uz) || 1; ux /= ul; uy /= ul; uz /= ul;
        P[3 * i] = x; P[3 * i + 1] = y; P[3 * i + 2] = z; U[3 * i] = ux; U[3 * i + 1] = uy; U[3 * i + 2] = uz;
        const dl = Math.hypot(x - L[0], y - L[1], z - L[2]);
        ok[i] = dl <= maxD ? 1 : 0;
        if (dl < row.lmin) row.lmin = dl; if (dl > row.lmax) row.lmax = dl;
      }
      // 截断球面内的连续段; 截口处线性插一个点
      let closed = !!arc.closed || wrapped;
      const pieces = [];
      for (let i = 0; i < n;) {
        if (!ok[i]) { i++; continue; }
        let j = i; while (j + 1 < n && ok[j + 1]) j++;
        pieces.push([i, j]); i = j + 1;
      }
      if (pieces.length !== 1 || pieces[0][0] !== 0 || pieces[0][1] !== n - 1) closed = false;
      for (const [i0, i1] of pieces) {
        const idx = [];
        const cut = (a, b) => {                        // a 在内 b 在外: 截口
          let lo = 0, hi = 1;
          const at = t => [P[3 * a] + t * (P[3 * b] - P[3 * a]), P[3 * a + 1] + t * (P[3 * b + 1] - P[3 * a + 1]), P[3 * a + 2] + t * (P[3 * b + 2] - P[3 * a + 2])];
          for (let it = 0; it < 30; it++) { const m = 0.5 * (lo + hi), q = at(m); if (Math.hypot(q[0] - L[0], q[1] - L[1], q[2] - L[2]) <= maxD) lo = m; else hi = m; }
          const q = at(lo), x = q[0] - E[0], y = q[1] - E[1], z = q[2] - E[2], l = Math.hypot(x, y, z) || 1;
          return addV(q[0], q[1], q[2], x / l, y / l, z / l);
        };
        const head = i0 > 0 ? cut(i0, i0 - 1) : -1, tail = i1 < n - 1 ? cut(i1, i1 + 1) : -1;
        // 抽稀: 与上一个保留点的夹角 ≥ SKY, 或这里转角 ≥ 12°, 或是尖点 (离引擎给的尖点最近的那个点)
        const isC = new Uint8Array(n);
        for (const c of cusp) {
          let bi = -1, bd = 2 * SKY;
          for (let i = i0; i <= i1; i++) { const d = angU(c, [U[3 * i], U[3 * i + 1], U[3 * i + 2]]); if (d < bd) { bd = d; bi = i; } }
          if (bi > i0 && bi < i1) isC[bi] = 1;
        }
        const vs = [], cs = [];
        if (head >= 0) vs.push(head);
        let last = -1;
        const cosT = Math.cos(12 * DEG);
        for (let i = i0; i <= i1; i++) {
          let keep = i === i0 || i === i1 || isC[i];
          if (!keep) {
            const d = angU([U[3 * last], U[3 * last + 1], U[3 * last + 2]], [U[3 * i], U[3 * i + 1], U[3 * i + 2]]);
            if (d >= SKY) keep = true;
            else {
              const ax = U[3 * i] - U[3 * i - 3], ay = U[3 * i + 1] - U[3 * i - 2], az = U[3 * i + 2] - U[3 * i - 1], bx = U[3 * i + 3] - U[3 * i], by = U[3 * i + 4] - U[3 * i + 1], bz = U[3 * i + 5] - U[3 * i + 2];
              const la = Math.hypot(ax, ay, az), lb = Math.hypot(bx, by, bz);
              if (la > 0 && lb > 0 && (ax * bx + ay * by + az * bz) / (la * lb) < cosT && d >= 0.1 * SKY) keep = true;
            }
          }
          if (!keep) continue;
          if (isC[i]) cs.push(vs.length);
          vs.push(addV(P[3 * i], P[3 * i + 1], P[3 * i + 2], U[3 * i], U[3 * i + 1], U[3 * i + 2]));
          last = i;
        }
        if (tail >= 0) vs.push(tail);
        // 闭合环: 首尾同一点只留一个; 引擎标成闭合、首尾却隔着一段的 (环在这里被边界切开了一个口) 当开弧
        let shut = closed;
        if (shut && vs.length > 3) { const g = ang(vs[0], vs[vs.length - 1]); if (g < 1e-9) vs.pop(); else if (g > SKY) shut = false; }
        if (vs.length < 2) continue;
        let len = 0;
        for (let i = 1; i < vs.length; i++) len += ang(vs[i - 1], vs[i]);
        if (shut) len += ang(vs[vs.length - 1], vs[0]);
        // 缩成一点的弧照留 (曲面的尖端: 弧的两个末端点在边界线交汇处并成一个, 前一层扇形连到这一点上)
        const cx = [0, 0, 0];
        for (const v of vs) { cx[0] += VU[3 * v]; cx[1] += VU[3 * v + 1]; cx[2] += VU[3 * v + 2]; }
        const cl = Math.hypot(cx[0], cx[1], cx[2]) || 1;
        const mid = vs[vs.length >> 1];
        row.segs.push({ v: vs, closed: shut, len, cusps: cs, ctr: [cx[0] / cl, cx[1] / cl, cx[2] / cl], mid: [VU[3 * mid], VU[3 * mid + 1], VU[3 * mid + 2]], n0: -1, n1: -1 });
      }
      row.raw.push({ p: Float32Array.from(P), d: Float32Array.from(U) });
    }
    // 末端点: 同一位置的 (两条弧在这里相接) 并成一个
    const nodeAt = v => {
      for (let i = 0; i < row.nodes.length; i++) if (ang(row.nodes[i].v, v) < 1e-5) return i;
      row.nodes.push({ v, u: [VU[3 * v], VU[3 * v + 1], VU[3 * v + 2]] });
      return row.nodes.length - 1;
    };
    for (const sg of row.segs) if (!sg.closed) { sg.n0 = nodeAt(sg.v[0]); sg.n1 = nodeAt(sg.v[sg.v.length - 1]); }
    cache.set(key, row);
    return row;
  };

  // ---- 两层配对
  const pair = (a, b) => {
    // 末端点: 互为最近, 且比各自第二近的明显近
    const near2 = (u, nodes) => { let i1 = -1, d1 = Infinity, d2 = Infinity; nodes.forEach((q, i) => { const d = angU(u, q.u); if (d < d1) { d2 = d1; d1 = d; i1 = i; } else if (d < d2) d2 = d; }); return [i1, d1, d2]; };
    const map = new Int32Array(a.nodes.length).fill(-1), back = new Int32Array(b.nodes.length).fill(-1);
    a.nodes.forEach((q, i) => {
      const [j, d1, d2] = near2(q.u, b.nodes);
      if (j < 0 || d1 > TOLN || d1 > 0.6 * d2) return;
      const [i2, e1, e2] = near2(b.nodes[j].u, a.nodes);
      if (i2 !== i || e1 > 0.6 * e2) return;
      map[i] = j; back[j] = i;
    });
    const usedB = new Uint8Array(b.segs.length), pairs = [], unA = [];
    a.segs.forEach((sa, ia) => {
      let best = -1, bd = Infinity, rev = false;
      if (sa.closed) {
        b.segs.forEach((sb, ib) => {
          if (usedB[ib] || !sb.closed) return;
          const d = angU(sa.ctr, sb.ctr), q = sb.len / sa.len;
          if (q < 0.4 || q > 2.5 || d > 0.5 * Math.max(sa.len, sb.len)) return;
          if (d < bd) { bd = d; best = ib; }
        });
      } else {
        const m0 = map[sa.n0], m1 = map[sa.n1];
        if (m0 >= 0 && m1 >= 0) b.segs.forEach((sb, ib) => {
          if (usedB[ib] || sb.closed) return;
          const fwd = sb.n0 === m0 && sb.n1 === m1, bwd = sb.n0 === m1 && sb.n1 === m0;
          if (!fwd && !bwd) return;
          const q = sb.len / sa.len;
          if (q < 0.25 || q > 4) return;
          const d = angU(sa.mid, sb.mid);
          if (d < bd) { bd = d; best = ib; rev = !fwd; }
        });
      }
      // 走向再按形状核一遍 (四分点对得近的那个方向): 两端落在同一个末端点上的弧 (绕了一圈回到原处, 首尾几乎相接)
      //   从末端点分不出走向, 引擎在相邻两层上给的方向又可能相反; 末端点说的走向与形状明显相悖的不配 (留给逐点对应)
      if (best >= 0 && !sa.closed) {
        const sb = b.segs[best], cost = r => { let c = 0; for (const f of [0.25, 0.5, 0.75]) { const i = Math.floor(f * (sa.v.length - 1)), j = Math.floor((r ? 1 - f : f) * (sb.v.length - 1)); c += ang(sa.v[i], sb.v[j]); } return c; };
        const cF = cost(false), cR = cost(true);
        if (sa.n0 === sa.n1) rev = cR < cF;
        else if ((rev ? cR : cF) > 3 * (rev ? cF : cR) + 3 * SKY) best = -1;
      }
      if (best >= 0) { usedB[best] = 1; pairs.push([ia, best, rev]); } else unA.push(ia);
    });
    // 段级配对: 末端点按点配不上的开弧 (短弧的两个末端点挨得近, 单看一个点分不清哪头对哪头). 两条段整体比: 正向 / 反向各算
    //   "两端夹角 + 四分点夹角" 之和, 取小的那个方向; 比次好的候选段明显小 (< 0.6 倍) 才配; 长度相称
    if (unA.length) {
      const at = (sg, f) => sg.v[Math.round(f * (sg.v.length - 1))];
      const cost = (sa, sb, rev) => { let c = 0; for (const f of [0, 0.25, 0.5, 0.75, 1]) c += ang(at(sa, f), at(sb, rev ? 1 - f : f)); return c; };
      const still = [];
      for (const ia of unA) {
        const sa = a.segs[ia];
        if (sa.closed) { still.push(ia); continue; }
        let best = -1, bc = Infinity, c2 = Infinity, brev = false;
        b.segs.forEach((sb, ib) => {
          if (usedB[ib] || sb.closed) return;
          const f = cost(sa, sb, false), r = cost(sa, sb, true), c = Math.min(f, r);
          if (c < bc) { c2 = bc; bc = c; best = ib; brev = r < f; } else if (c < c2) c2 = c;
        });
        const q = best >= 0 ? b.segs[best].len / Math.max(sa.len, 1e-12) : 0;
        if (best >= 0 && bc <= 5 * TOLN && bc < 0.6 * c2 && q > 0.25 && q < 4) { usedB[best] = 1; pairs.push([ia, best, brev]); } else still.push(ia);
      }
      unA.length = 0; unA.push(...still);
    }
    const unB = [];
    b.segs.forEach((_, ib) => { if (!usedB[ib]) unB.push(ib); });
    return { pairs, unA, unB, full: !unA.length && !unB.length };
  };

  // ---- 拉链: P、Q 两串顶点 (同向), 每步走视方向上较短的对角线; rungs = 必经的 [P 序号, Q 序号]
  const zip = (P, Q, rungs = null) => {
    let i = 0, j = 0, r = 0;
    while (i < P.length - 1 || j < Q.length - 1) {
      while (rungs && r < rungs.length && (rungs[r][0] < i || rungs[r][1] < j || (rungs[r][0] === i && rungs[r][1] === j))) r++;
      const ti = rungs && r < rungs.length ? rungs[r][0] : P.length - 1, tj = rungs && r < rungs.length ? rungs[r][1] : Q.length - 1;
      if (j === tj || (i < ti && ang(P[i + 1], Q[j]) <= ang(P[i], Q[j + 1]))) { emit(P[i], P[i + 1], Q[j]); i++; }
      else { emit(P[i], Q[j + 1], Q[j]); j++; }
    }
  };
  // 一对配上的段 → 两串同向的顶点 (闭合环: 转到最近的起点, 首点重复一次收口) 与尖点横档
  const strings = (sa, sb, rev) => {
    let P = sa.v, Q = rev ? sb.v.slice().reverse() : sb.v, cq = rev ? sb.cusps.map(k => sb.v.length - 1 - k).reverse() : sb.cusps;
    if (sa.closed) {
      let j0 = 0, bd = Infinity;
      for (let j = 0; j < Q.length; j++) { const d = ang(P[0], Q[j]); if (d < bd) { bd = d; j0 = j; } }
      // 走向: 两个方向各试, 取四分点对得近的
      const rot = (A, k, s) => A.map((_, i) => A[((k + s * i) % A.length + A.length) % A.length]);
      const cost = Qr => { let c = 0; for (const f of [0.25, 0.5, 0.75]) c += ang(P[Math.floor(f * P.length)], Qr[Math.floor(f * Qr.length)]); return c; };
      const q1 = rot(Q, j0, 1), q2 = rot(Q, j0, -1);
      Q = cost(q1) <= cost(q2) ? q1 : q2;
      P = P.concat([P[0]]); Q = Q.concat([Q[0]]);
      return { P, Q, rungs: null };
    }
    // 尖点横档: 两层的尖点按沿线先后一一对上 (个数相同、位置相近) 才用
    let rungs = null;
    if (sa.cusps.length && sa.cusps.length === cq.length) {
      rungs = [];
      for (let k = 0; k < cq.length; k++) { if (ang(P[sa.cusps[k]], Q[cq[k]]) > TOLN) { rungs = null; break; } rungs.push([sa.cusps[k], cq[k]]); }
    }
    return { P, Q, rungs };
  };
  // 配上的一对段在天空里隔多远: 照拉链的走法走一遍 (不出三角形), 一路上两边当前点之间的最大夹角.
  //   只比几个等分点会漏掉局部的错位 (一条线上多出一个小折回, 另一条没有)
  const zipGap = (P, Q, rungs) => {
    let i = 0, j = 0, r = 0, m = ang(P[0], Q[0]);
    while (i < P.length - 1 || j < Q.length - 1) {
      while (rungs && r < rungs.length && (rungs[r][0] < i || rungs[r][1] < j || (rungs[r][0] === i && rungs[r][1] === j))) r++;
      const ti = rungs && r < rungs.length ? rungs[r][0] : P.length - 1, tj = rungs && r < rungs.length ? rungs[r][1] : Q.length - 1;
      if (j === tj || (i < ti && ang(P[i + 1], Q[j]) <= ang(P[i], Q[j + 1]))) i++; else j++;
      const d = ang(P[i], Q[j]);
      if (d > m) m = d;
    }
    return m;
  };
  const pairGap = (a, b, q) => { const st = strings(a.segs[q[0]], b.segs[q[1]], q[2]); return zipGap(st.P, st.Q, st.rungs); };
  const sep = (a, b, pr) => { let m = 0; for (const c of corrOf(a, b, pr)) m = Math.max(m, zipGap(c.P, c.Q, c.rungs)); return m; };
  // 灯旁: 晕面在灯上是锥顶, 离灯近处沿半径方向弯得厉害 (天空里却都挤在灯的方向上, 隔不开), 两层的半径差
  //   不能超过离灯距离的 0.35 倍
  const lampDist = row => {
    if (row.rl === undefined) {
      let rl = Infinity;
      for (const sg of row.segs) for (const v of sg.v) { const d = Math.hypot(VP[3 * v] - L[0], VP[3 * v + 1] - L[1], VP[3 * v + 2] - L[2]); if (d < rl) rl = d; }
      row.rl = rl;
    }
    return row.rl;
  };
  const nearLamp = (a, b) => b.r - a.r > 0.35 * Math.max(Math.min(lampDist(a), lampDist(b)), 2e-3 * D);
  // 三维弦高: 两层之间配上的段, 沿拉链取 9 对对应点, 连线中点到中间那层球壳 (真的算出来) 上的线的距离.
  //   曲面几乎贴着球壳走时 (沿径向相邻两层只隔几十厘米, 沿曲面却隔得远), 天空间距和半径间距都看不出来, 三角形弯着.
  //   超过 1.5 cm + 0.15% r (r 取到灯、眼较近的距离) 就把中间那层加进来. 返回中间层 (要加) 或 null
  const corrOf = (a, b, pr) => {
    const out = pr.pairs.map(q => strings(a.segs[q[0]], b.segs[q[1]], q[2]));
    if (!pr.full) for (const [P, Q] of loosePlan(a, b, pr, LIM, false).strips) out.push({ P, Q, rungs: null });
    return out;
  };
  const sagCheck = (a, b, pr) => {
    const cs = corrOf(a, b, pr);
    // 半径只差 0.2% 以内的两层 (变化前后夹紧的那些) 弦高无从说起, 不另算中间层 —— 否则二分链上每一节都多算一层, 引擎调用数撑到上限
    if (!cs.length || b.r / a.r - 1 < 2e-3) return null;
    const m = prep(Math.sqrt(a.r * b.r));
    if (!m.segs.length) return null;
    const near = (x, y, z) => {
      let bd = Infinity;
      for (const sg of m.segs) { const w = sg.v; for (let k = 0; k + 1 < w.length; k++) {
        const i = w[k], j = w[k + 1], ax = VP[3 * i], ay = VP[3 * i + 1], az = VP[3 * i + 2], bx = VP[3 * j] - ax, by = VP[3 * j + 1] - ay, bz = VP[3 * j + 2] - az;
        const px = x - ax, py = y - ay, pz = z - az, l2 = bx * bx + by * by + bz * bz;
        let u = l2 > 0 ? (px * bx + py * by + pz * bz) / l2 : 0; u = u < 0 ? 0 : u > 1 ? 1 : u;
        const d = Math.hypot(px - u * bx, py - u * by, pz - u * bz); if (d < bd) bd = d; } }
      return bd;
    };
    for (const { P, Q } of cs) {
      for (let k = 1; k < 8; k++) {
        const i = P[Math.round(k * (P.length - 1) / 8)], j = Q[Math.round(k * (Q.length - 1) / 8)];
        const x = 0.5 * (VP[3 * i] + VP[3 * j]), y = 0.5 * (VP[3 * i + 1] + VP[3 * j + 1]), z = 0.5 * (VP[3 * i + 2] + VP[3 * j + 2]);
        const rr = Math.min(Math.hypot(x - E[0], y - E[1], z - E[2]), Math.hypot(x - L[0], y - L[1], z - L[2]));
        if (near(x, y, z) > 0.015 + 0.0015 * rr) return m;
      }
    }
    return null;
  };
  const zipPairs = (a, b, pr) => { for (const [ia, ib, rev] of pr.pairs) { const s = strings(a.segs[ia], b.segs[ib], rev); zip(s.P, s.Q, s.rungs); } };
  // 配不上的段 (引擎在相邻两层上把同一条线分成了不同的几段、多给了一遍, 或线的条数 / 连法真的变了): 不看分段, 逐点找对面
  //   视方向最近的点 (lim 之内), 同一段上脚点单调的一串一串缝; 再反过来把对面没盖到的补一遍. 返回方案:
  //   strips = 要缝的 [P, Q], unc = 两边没缝上的总长 (弧度), caps = 该向自己中心收口的小段
  const stat = { rows: 0, events: 0, caps: 0, loose: 0, looseDeg: 0, looseR: [], calls: 0, ms: 0, skipped: 0 };
  const loosePlan = (a, b, pr, lim, withCaps) => {
    const strips = [], caps = [];
    let unc = 0;
    const tiny = withCaps ? Math.min(CAPMAX, 30 * Math.sqrt(2 * Math.max(EPS, b.r / a.r - 1))) : -1;
    const split = (row, list) => list.filter(i => { if (row.segs[i].len > tiny) return true; caps.push([row, i]); return false; });
    const LA = split(a, pr.unA), LB = split(b, pr.unB);
    const seq = sg => sg.closed ? sg.v.concat([sg.v[0]]) : sg.v;
    const A = LA.map(i => seq(a.segs[i])), Bs = LB.map(i => seq(b.segs[i]));
    const covA = A.map(w => new Uint8Array(w.length)), covB = Bs.map(w => new Uint8Array(w.length));     // 边 k..k+1 缝上了没有
    const pass = (S, T, covS, covT, flip) => {
      S.forEach((w, si) => {
        const ft = w.map(v => { let bq = -1, bj = -1, bd = lim; T.forEach((x, q) => { for (let j = 0; j < x.length; j++) { const d = ang(v, x[j]); if (d < bd) { bd = d; bq = q; bj = j; } } }); return bq < 0 ? null : [bq, bj]; });
        for (let s0 = 0; s0 < ft.length;) {
          if (!ft[s0] || (s0 + 1 < w.length && covS[si][s0])) { s0++; continue; }
          let e = s0, dir = 0;
          while (e + 1 < ft.length && ft[e + 1] && ft[e + 1][0] === ft[s0][0] && !covS[si][e]) {
            const sg = Math.sign(ft[e + 1][1] - ft[e][1]);
            if ((sg && dir && sg !== dir) || Math.abs(ft[e + 1][1] - ft[e][1]) > 12) break;
            if (sg) dir = sg; e++;
          }
          if (e > s0) {
            const q = ft[s0][0], x = T[q], j0 = ft[s0][1], j1 = ft[e][1], st = j1 >= j0 ? 1 : -1, Q = [];
            for (let j = j0; j !== j1 + st; j += st) Q.push(x[j]);
            for (let j = Math.min(j0, j1); j < Math.max(j0, j1); j++) covT[q][j] = 1;
            for (let k = s0; k < e; k++) covS[si][k] = 1;
            strips.push(flip ? [Q, w.slice(s0, e + 1)] : [w.slice(s0, e + 1), Q]);
          }
          s0 = Math.max(e, s0 + 1);
        }
      });
    };
    if (A.length && Bs.length) { pass(A, Bs, covA, covB, false); pass(Bs, A, covB, covA, true); }
    // 没缝上的长度; 一条边也没缝上的段至少按 UNCOK 算 —— 缩成一点的段 (退化成一条曲线的晕面每层只有一个点) 长度是 0,
    //   只按长度算的话两层隔得再远也算 "缝全了", 中间一个三角形都没有
    const left = (W, cov) => { let u = 0; W.forEach((w, i) => { let ui = 0, any = false; for (let k = 0; k + 1 < w.length; k++) { if (!cov[i][k]) ui += ang(w[k], w[k + 1]); else any = true; } u += any ? ui : Math.max(ui, 1.01 * UNCOK); }); return u; };
    unc = left(A, covA) + left(Bs, covB);
    return { strips, caps, unc, r: a.r };
  };
  // 收口: 一般向这一小段自己的中心收; 灯所在球壳两侧的那一对层向灯收 (晕面在灯上是锥顶, 层线上的点与灯的连线就在面上)
  let lampV = -1;
  const applyLoose = (pl, toLamp = false) => {
    for (const [P, Q] of pl.strips) zip(P, Q);
    for (const [row, i] of pl.caps) {
      const sg = row.segs[i];
      let x = 0, y = 0, z = 0;
      for (const v of sg.v) { x += VP[3 * v]; y += VP[3 * v + 1]; z += VP[3 * v + 2]; }
      if (toLamp && lampV < 0) lampV = addV(L[0], L[1], L[2], Dx / D, Dy / D, Dz / D);
      const c = toLamp ? lampV : addV(x / sg.v.length, y / sg.v.length, z / sg.v.length, sg.ctr[0], sg.ctr[1], sg.ctr[2]);
      for (let k = 0; k + 1 < sg.v.length; k++) emit(sg.v[k], sg.v[k + 1], c);
      if (sg.closed) emit(sg.v[sg.v.length - 1], sg.v[0], c);
      stat.caps++;
    }
    if (pl.unc > UNCOK) { stat.loose++; stat.looseDeg += pl.unc / DEG; if (pl.r) stat.looseR.push(pl.r); }
  };
  // 两层之间能不能缝全: 配得全, 或配不上的段逐点对应后没缝上的不到 UNCOK
  const LIM = 2.5 * GAPSKY;
  const whole = (a, b, pr) => pr.full || loosePlan(a, b, pr, LIM, false).unc <= UNCOK;

  // ---- 层的序列: 初始对数等距, 再按配对情况加层
  // 缝一对层: 配上的段一段对一段; 但配上的两段若对应点在天空里隔得比 2 × 层间距还远 (两端的末端点对上了, 中间却不是
  //   同一条线: 一条多绕了一截, 或加层加到头了还隔得远), 不整段硬缝 —— 归到 "配不上的" 里逐点找近的对应, 够不着的留着
  const sewPair = (a, b, pr, withCaps, toLamp) => {
    const far = pr.pairs.filter(q => pairGap(a, b, q) > 2 * GAPSKY);
    if (far.length) pr = { pairs: pr.pairs.filter(q => !far.includes(q)), unA: pr.unA.concat(far.map(q => q[0])), unB: pr.unB.concat(far.map(q => q[1])), full: false };
    zipPairs(a, b, pr);
    if (!pr.full) applyLoose(loosePlan(a, b, pr, LIM, withCaps), toLamp);
  };
  // 预算按进了网格的层的顶点数算 (二分时试算的层不进网格, 不算)
  const rows = [], keep = row => { rows.push(row); usedV += row.segs.reduce((x, sg) => x + sg.v.length, 0); };
  let usedV = 0;
  const stitch = (a, b, depth) => {
    const pr = pair(a, b), ratio = b.r / a.r - 1;
    const room = rows.length < MAXROWS && calls < MAXCALLS && usedV < 2.5 * VBUDGET;
    const lampPair = a.r < D && b.r > D;
    const done = withCaps => sewPair(a, b, pr, withCaps, lampPair);
    // 灯所在的球壳两侧的那一对层: 晕面在灯上收成锥顶, 两侧的层线各自缩向灯的方向 (不是一一对应), 直接收口
    if (ratio <= EPS || !room || lampPair || (a.critLo && ratio < 3 * CRIT)) { done(true); return; }
    if (whole(a, b, pr)) {
      let m = null;
      if (depth < 12 && ratio > 4 * EPS && ((sep(a, b, pr) > GAPSKY || nearLamp(a, b) || ratio > RSTEP) || (m = sagCheck(a, b, pr)))) {
        if (!m) m = prep(Math.sqrt(a.r * b.r));
        keep(m); stitch(a, m, depth + 1); stitch(m, b, depth + 1);
        return;
      }
      done(false);
      return;
    }
    // 缝不全: 两层之间线的条数 / 连法变了. 二分半径把变化夹紧: lo 与 a 一路缝得全, hi 缝不全.
    //   途中 lo 与 hi 自己缝全了 → 其实没有变化, 只是隔得远时对应含糊, 加了这几层就够了
    const chain = [a];
    let lo = a, hi = b, real = true;
    while (hi.r / lo.r - 1 > EPS && calls < MAXCALLS) {
      let m = prep(Math.sqrt(lo.r * hi.r));
      if (suspect(m, lo) && !whole(lo, m, pair(lo, m))) m = redo(m);
      if (whole(lo, m, pair(lo, m))) { lo = m; chain.push(m); if (whole(lo, hi, pair(lo, hi))) { real = false; break; } } else hi = m;
    }
    // 夹紧之后外端还可疑 (末端点少了、中间没有临界半径): 用紧挨着的里端当种子重算外端; 这下缝得全了就不是变化
    if (real && suspect(hi, lo)) { redo(hi); if (whole(lo, hi, pair(lo, hi))) real = false; }
    if (real) { stat.events++; if (opt.log) opt.log(`event r ${lo.r.toFixed(5)} .. ${hi.r.toFixed(5)}: segs ${lo.segs.length} -> ${hi.segs.length}, nodes ${lo.nodes.length} -> ${hi.nodes.length}`); }
    // a → lo: 链上能直接缝全的尽量跳着连
    for (let i = 0; i < chain.length - 1;) {
      let j = chain.length - 1;
      while (j > i + 1 && !whole(chain[i], chain[j], pair(chain[i], chain[j]))) j--;
      keep(chain[j]);
      stitch(chain[i], chain[j], depth + 1);
      i = j;
    }
    if (hi !== b) keep(hi);
    if (real) sewPair(lo, hi, pair(lo, hi), true, false); else stitch(lo, hi, depth + 1);
    if (hi !== b) stitch(hi, b, depth + 1);
  };
  // 临界半径 (opt.critical): 层线的条数 / 连法只在这些半径上变 —— 弧的末端点落在三维边界线上, 边界线离眼的距离取极值
  //   (边界线与球壳相切) 或边界线起止的地方, 末端点成对地出现 / 消失, 一条弧分成两条或两条并成一条. 调用方从已经解出来的
  //   边界线直接算出这些半径 (criticalRadii), 这里在每个半径两侧各放一层 (相对 ± CRIT), 这一对层之间不再二分, 直接按
  //   "配得上的缝、缩没的小段收口、其余逐点对应" 处理; 没有预先给出的变化 (面内部的极值、鞍点) 仍由二分去夹
  const CRIT = 2e-4, crit = [];
  for (const rc of (opt.critical || []).slice().sort((p, q) => p - q)) {
    if (!(rc > rMin * 1.01 && rc < rMax * 0.99) || Math.abs(rc / D - 1) < 2 * PINL) continue;
    if (crit.length && rc / crit[crit.length - 1] - 1 < 4 * CRIT) continue;
    crit.push(rc);
  }
  const init = [], rs = [], allCrit = crit.concat([D]);
  for (let k = 0; k < N0; k++) {
    const r = rMin * Math.pow(rMax / rMin, k / (N0 - 1));
    if (Math.abs(r / D - 1) > PINL && !crit.some(rc => Math.abs(r / rc - 1) < 2 * CRIT)) rs.push(r);
  }
  if (rMin < D * (1 - PINL) && rMax > D * (1 + PINL)) rs.push(D * (1 - PINL), D * (1 + PINL));
  const critLo = new Set();
  for (const rc of crit) { const lo = rc * (1 - CRIT); rs.push(lo, rc * (1 + CRIT)); critLo.add(Math.round(Math.log(lo) * 1e9)); }
  rs.sort((p, q) => p - q);
  for (const r of rs) init.push(prep(r));
  // 截断 (离灯 maxD) 也是一条边界: 截断球面开始切到曲面 (某层弧上离灯最远的点越过 maxD) 和把它整个切掉 (最近的点也越过)
  //   的半径, 同样是层线拓扑变的地方 (截口成对出现 / 一段弧被切没). 用各层离灯距离的最大 / 最小值二分出来 (只算球壳, 便宜),
  //   和边界线给的临界半径一样处理
  {
    const extra = [];
    const has = row => row.raw.length > 0;
    for (let k = 0; k + 1 < init.length; k++) {
      const a = init[k], b = init[k + 1];
      if (!has(a) || !has(b)) continue;
      for (const key of ['lmax', 'lmin']) {
        if ((a[key] - maxD) * (b[key] - maxD) >= 0) continue;
        let lo = a.r, hi = b.r, flo = a[key] - maxD;
        for (let it = 0; it < 40 && hi / lo - 1 > CRIT && calls < MAXCALLS; it++) {
          const m = Math.sqrt(lo * hi), row = prep(m);
          if (!has(row)) break;
          const f = row[key] - maxD;
          if (f * flo > 0) { lo = m; flo = f; } else hi = m;
        }
        extra.push(Math.sqrt(lo * hi));
      }
    }
    if (extra.length) {
      for (const rc of extra) { const lo = rc * (1 - CRIT); rs.push(lo, rc * (1 + CRIT)); critLo.add(Math.round(Math.log(lo) * 1e9)); allCrit.push(rc); }
      rs.sort((p, q) => p - q);
      init.length = 0;
      for (const r of rs) init.push(prep(r));
      if (opt.log) opt.log(`truncation critical radii: ${extra.map(x => x.toFixed(4)).join(' ')}`);
    }
  }
  for (const row of init) if (critLo.has(Math.round(Math.log(row.r) * 1e9))) row.critLo = true;
  // 引擎漏弧的补救: 开弧的两个末端点在边界线上, 它只能在临界半径上消失 (边界线交汇 / 离眼距离取极值处) ——
  //   相邻两层之间末端点变少了、中间却没有已知的临界半径, 就是引擎的网格种子没找到 (剩下的长弧占了大头, 只按总长看不出来).
  //   这时用邻层的状态点当种子重算少的那一层 (shellAt.again), 多出来了才换. 线短了一成以上的也重算.
  //   两遍: 从外往里 (外层带回里层漏的), 从里往外 (反过来)
  const critBetween = (r1, r2) => { const lo = Math.min(r1, r2) * (1 - 2 * CRIT), hi = Math.max(r1, r2) * (1 + 2 * CRIT); return allCrit.some(rc => rc > lo && rc < hi); };
  const lenOfRow = row => row.segs.reduce((x, sg) => x + sg.len, 0);
  const suspect = (row, nb) => (row.nodes.length < nb.nodes.length && !critBetween(row.r, nb.r)) || lenOfRow(row) < 0.9 * lenOfRow(nb);
  // 原地更新 (init / rows / 二分链里引用的是同一个对象, 都看得到新结果)
  const redo = row => {
    if (!shellAt.again || calls >= MAXCALLS || !shellAt.again(row.r)) return row;
    const key = Math.round(Math.log(row.r) * 1e9);
    cache.delete(key);
    const nw = prep(row.r), keepLo = row.critLo;
    for (const k of Object.keys(row)) delete row[k];
    Object.assign(row, nw);
    if (keepLo) row.critLo = keepLo;
    cache.set(key, row);
    return row;
  };
  if (shellAt.again) {
    for (let k = init.length - 2; k >= 0; k--) if (suspect(init[k], init[k + 1])) init[k] = redo(init[k]);
    for (let k = 1; k < init.length; k++) if (suspect(init[k], init[k - 1])) init[k] = redo(init[k]);
    if (opt.log) opt.log(`after re-seeding: ${calls} engine calls`);
  }
  for (const r of init) keep(r);
  for (let k = 0; k + 1 < init.length; k++) {
    if (!init[k].segs.length && !init[k + 1].segs.length) continue;
    stitch(init[k], init[k + 1], 0);
  }
  rows.sort((p, q) => p.r - q.r);

  // ---- 输出: 每个三角形三个顶点依次排; 法向 = 共顶点的三角形里与本三角形朝向相近 (夹角 < 60°) 的面法向之和
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
        const i = cnt[v] + Math.floor(q * stp), j = inc[i], c = TN[3 * j] * TN[3 * t] + TN[3 * j + 1] * TN[3 * t + 1] + TN[3 * j + 2] * TN[3 * t + 2];
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
  stat.rows = rows.filter(r => r.segs.length).length; stat.calls = calls; stat.ms = Date.now() - t0; stat.skipped = nSkip;
  // 各层的线 (画层线用). gap = 这一层与里面一层在天空里隔多远 (弧度): 贴着眼睛的那些层在天空里几乎重合,
  //   画的时候按累计的 gap 挑着画, 免得平行光弧旁挤成一条粗带
  let prevRow = null;
  for (const row of rows) {
    if (!row.segs.length) continue;
    if (prevRow) { const pr = pair(prevRow, row); row.gap = pr.pairs.length ? sep(prevRow, row, pr) : GAPSKY; } else row.gap = Infinity;
    prevRow = row;
  }
  const lines = rows.filter(r => r.raw.length).map(r => ({ r: r.r, gap: r.gap === undefined ? GAPSKY : r.gap, curves: r.raw }));
  return { pos, dir, h: hh, nrm, nTri: nT, lines, stat, rowR: rows.filter(r => r.segs.length).map(r => r.r) };
}

// 验收: 球壳引擎偶尔在某些半径上漏弧 (大半径处整条焦散找不到, 个别区间少一段), 那里网格就是洞. 拿另一套引擎算的
//   光锥层线来查: 层线上抽 nSample 个点 (灯 / 眼 0.02 D 之内、截断处的不算), 量到网格 (点到三角形) 的距离,
//   超过 5 cm + 1% r (r = 到灯、眼较近者的距离) 的算漏. 返回 { n, miss, frac, missR (漏点离眼的距离), allR }
export function coverCheck(mesh, layers, L, E, maxD = 80, nSample = 800) {
  const P = mesh.pos, nT = P.length / 9, D = Math.hypot(L[0] - E[0], L[1] - E[1], L[2] - E[2]);
  let total = 0;
  for (const l of layers) for (const c of l.curves) total += (c.p || c).length / 3;
  if (!nT || !total) return { n: 0, miss: 0, frac: nT ? 0 : 1 };
  const CELL = Math.max(0.25, D / 40), grid = new Map(), big = [], key = (i, j, k) => (i + 2048) * 16777216 + (j + 2048) * 4096 + (k + 2048);
  for (let t = 0; t < nT; t++) {
    const o = 9 * t, lo = [0, 1, 2].map(c => Math.floor(Math.min(P[o + c], P[o + 3 + c], P[o + 6 + c]) / CELL)), hi = [0, 1, 2].map(c => Math.floor(Math.max(P[o + c], P[o + 3 + c], P[o + 6 + c]) / CELL));
    if ((hi[0] - lo[0] + 1) * (hi[1] - lo[1] + 1) * (hi[2] - lo[2] + 1) > 512) { big.push(t); continue; }     // 大三角形单列, 每次都查
    for (let i = lo[0]; i <= hi[0]; i++) for (let j = lo[1]; j <= hi[1]; j++) for (let k = lo[2]; k <= hi[2]; k++) { const q = key(i, j, k); let a = grid.get(q); if (!a) grid.set(q, a = []); a.push(t); }
  }
  const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]], dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  const ptTri = (p, a, b, c) => {
    const ab = sub(b, a), ac = sub(c, a), ap = sub(p, a), d1 = dot(ab, ap), d2 = dot(ac, ap);
    if (d1 <= 0 && d2 <= 0) return Math.hypot(...ap);
    const bp = sub(p, b), d3 = dot(ab, bp), d4 = dot(ac, bp);
    if (d3 >= 0 && d4 <= d3) return Math.hypot(...bp);
    const vc = d1 * d4 - d3 * d2;
    if (vc <= 0 && d1 >= 0 && d3 <= 0) { const v = d1 / (d1 - d3); return Math.hypot(ap[0] - v * ab[0], ap[1] - v * ab[1], ap[2] - v * ab[2]); }
    const cp = sub(p, c), d5 = dot(ab, cp), d6 = dot(ac, cp);
    if (d6 >= 0 && d5 <= d6) return Math.hypot(...cp);
    const vb = d5 * d2 - d1 * d6;
    if (vb <= 0 && d2 >= 0 && d6 <= 0) { const w = d2 / (d2 - d6); return Math.hypot(ap[0] - w * ac[0], ap[1] - w * ac[1], ap[2] - w * ac[2]); }
    const va = d3 * d6 - d5 * d4;
    if (va <= 0 && d4 - d3 >= 0 && d5 - d6 >= 0) { const w = (d4 - d3) / ((d4 - d3) + (d5 - d6)), bc = sub(c, b); return Math.hypot(bp[0] - w * bc[0], bp[1] - w * bc[1], bp[2] - w * bc[2]); }
    const den = 1 / (va + vb + vc), v = vb * den, w = vc * den;
    return Math.hypot(ap[0] - v * ab[0] - w * ac[0], ap[1] - v * ab[1] - w * ac[1], ap[2] - v * ab[2] - w * ac[2]);
  };
  const step = Math.max(1, Math.floor(total / nSample));
  let n = 0, miss = 0, k = 0;
  const missR = [], allR = [];
  for (const l of layers) for (const c of l.curves) {
    const q = c.p || c;
    for (let i = 0; i < q.length / 3; i++, k++) {
      if (k % step) continue;
      const X = [q[3 * i], q[3 * i + 1], q[3 * i + 2]], rl = Math.hypot(X[0] - L[0], X[1] - L[1], X[2] - L[2]), re = Math.hypot(X[0] - E[0], X[1] - E[1], X[2] - E[2]), rr = Math.min(rl, re);
      if (rr < 0.02 * D || rl > maxD - 0.5) continue;
      const tol = 0.05 + 0.01 * rr, cc = X.map(x => Math.floor(x / CELL));
      let dm = Infinity;
      for (let a = cc[0] - 1; a <= cc[0] + 1 && dm > tol; a++) for (let b = cc[1] - 1; b <= cc[1] + 1 && dm > tol; b++) for (let d = cc[2] - 1; d <= cc[2] + 1 && dm > tol; d++) {
        const g = grid.get(key(a, b, d)); if (!g) continue;
        for (const t of g) { const o = 9 * t, e = ptTri(X, [P[o], P[o + 1], P[o + 2]], [P[o + 3], P[o + 4], P[o + 5]], [P[o + 6], P[o + 7], P[o + 8]]); if (e < dm) { dm = e; if (dm <= tol) break; } }
      }
      // 离得比一格远的 (容差比格子大) 和大三角形: 逐个查
      for (let q = 0; q < big.length && dm > tol; q++) { const o = 9 * big[q], e = ptTri(X, [P[o], P[o + 1], P[o + 2]], [P[o + 3], P[o + 4], P[o + 5]], [P[o + 6], P[o + 7], P[o + 8]]); if (e < dm) dm = e; }
      if (dm > tol && tol > CELL) for (let a = cc[0] - 3; a <= cc[0] + 3 && dm > tol; a++) for (let b = cc[1] - 3; b <= cc[1] + 3 && dm > tol; b++) for (let d = cc[2] - 3; d <= cc[2] + 3 && dm > tol; d++) {
        if (Math.abs(a - cc[0]) <= 1 && Math.abs(b - cc[1]) <= 1 && Math.abs(d - cc[2]) <= 1) continue;
        const g = grid.get(key(a, b, d)); if (!g) continue;
        for (const t of g) { const o = 9 * t, e = ptTri(X, [P[o], P[o + 1], P[o + 2]], [P[o + 3], P[o + 4], P[o + 5]], [P[o + 6], P[o + 7], P[o + 8]]); if (e < dm) { dm = e; if (dm <= tol) break; } }
      }
      n++; allR.push(re); if (dm > tol) { miss++; missR.push(re); }
    }
  }
  return { n, miss, frac: n ? miss / n : 0, missR, allR };
}


// 球壳网格用到哪个半径为止: 漏点 (验收抽样里没被盖到的光锥层线点) 与没缝上的层都在它之外.
//   球壳引擎的缺口多在大半径处 (ρ ≳ 1.7 之外焦散找不到), 眼睛附近是完整的 —— 而眼睛附近正是光锥层缝合最弱的地方.
//   所以两者拼起来用: 半径 rc 以内用球壳网格, 以外用光锥层缝的网格 (沿球面 r = rc 裁开).
//   容许零星的漏点 (抽样容差): 从里往外数, 第 i 个漏点之前容许 1 + 0.2% × (它以内的抽样点数) 个.
//   返回 rc: Infinity = 全用球壳; 0 = 不用; 否则取不超过它的最大一层的半径 (网格正好止于一条层线)
export function shellCutoff(mesh, cv, D) {
  const bad = cv.missR.slice().sort((a, b) => a - b), all = cv.allR.slice().sort((a, b) => a - b);
  let rc = Infinity, j = 0;
  for (let i = 0; i < bad.length; i++) {
    while (j < all.length && all[j] < bad[i]) j++;
    if (i + 1 > 1 + 0.002 * j) { rc = 0.97 * bad[i]; break; }
  }
  for (const r of mesh.stat.looseR) rc = Math.min(rc, 0.97 * r);
  if (rc === Infinity) return Infinity;
  let best = 0;
  for (const r of mesh.rowR) if (r <= rc && r > best) best = r;
  return best >= 0.05 * D ? best : 0;
}

// 网格里半径 (到 E 的距离) 不超过 rc 的那些三角形
export function cutInside(mesh, E, rc) {
  const P = mesh.pos, nT = P.length / 9, keep = [], lim = rc * (1 + 1e-6);
  for (let t = 0; t < nT; t++) {
    let ok = true;
    for (let k = 0; k < 3 && ok; k++) { const o = 9 * t + 3 * k; if (Math.hypot(P[o] - E[0], P[o + 1] - E[1], P[o + 2] - E[2]) > lim) ok = false; }
    if (ok) keep.push(t);
  }
  const out = { pos: new Float32Array(9 * keep.length), dir: new Float32Array(9 * keep.length), nrm: new Float32Array(9 * keep.length), h: new Float32Array(3 * keep.length), nTri: keep.length };
  keep.forEach((t, i) => { out.pos.set(P.subarray(9 * t, 9 * t + 9), 9 * i); out.dir.set(mesh.dir.subarray(9 * t, 9 * t + 9), 9 * i); out.nrm.set(mesh.nrm.subarray(9 * t, 9 * t + 9), 9 * i); out.h.set(mesh.h.subarray(3 * t, 3 * t + 3), 3 * i); });
  return out;
}

// 网格在球面 |P − E| = rc 之外的部分: 整个在外的三角形照留, 跨球面的沿球面裁开 (交点按半径线性插值, 各属性一起插)
export function clipOutside(mesh, E, rc) {
  const P = mesh.pos, Dr = mesh.dir, Nr = mesh.nrm, H = mesh.h, nT = P.length / 9;
  const pos = [], dir = [], nrm = [], hh = [];
  const vert = (t, k) => { const o = 9 * t + 3 * k; return { p: [P[o], P[o + 1], P[o + 2]], d: [Dr[o], Dr[o + 1], Dr[o + 2]], n: [Nr[o], Nr[o + 1], Nr[o + 2]], h: H[3 * t + k], r: Math.hypot(P[o] - E[0], P[o + 1] - E[1], P[o + 2] - E[2]) }; };
  const mix = (a, b) => {
    const t = (rc - a.r) / (b.r - a.r), lp = (x, y) => [x[0] + t * (y[0] - x[0]), x[1] + t * (y[1] - x[1]), x[2] + t * (y[2] - x[2])], un = v => { const l = Math.hypot(v[0], v[1], v[2]) || 1; return [v[0] / l, v[1] / l, v[2] / l]; };
    return { p: lp(a.p, b.p), d: un(lp(a.d, b.d)), n: un(lp(a.n, b.n)), h: a.h + t * (b.h - a.h), r: rc };
  };
  const put = (a, b, c) => { for (const v of [a, b, c]) { pos.push(v.p[0], v.p[1], v.p[2]); dir.push(v.d[0], v.d[1], v.d[2]); nrm.push(v.n[0], v.n[1], v.n[2]); hh.push(v.h); } };
  for (let t = 0; t < nT; t++) {
    const v = [vert(t, 0), vert(t, 1), vert(t, 2)], out = v.map(q => q.r >= rc), n = out[0] + out[1] + out[2];
    if (n === 3) { put(v[0], v[1], v[2]); continue; }
    if (n === 0) continue;
    for (let k = 0; k < 3; k++) {
      const a = v[k], b = v[(k + 1) % 3], c = v[(k + 2) % 3];
      if (n === 1 && out[k]) { put(a, mix(a, b), mix(a, c)); break; }                         // 一个顶点在外: 留一个小三角形
      if (n === 2 && !out[k]) { const ab = mix(a, b), ac = mix(a, c); put(ab, b, c); put(ab, c, ac); break; }   // 一个顶点在内: 留一个四边形
    }
  }
  return { pos: Float32Array.from(pos), dir: Float32Array.from(dir), nrm: Float32Array.from(nrm), h: Float32Array.from(hh), nTri: hh.length / 3 };
}

export function joinMeshes(a, b) {
  const cat = (x, y) => { const o = new Float32Array(x.length + y.length); o.set(x); o.set(y, x.length); return o; };
  return { pos: cat(a.pos, b.pos), dir: cat(a.dir, b.dir), nrm: cat(a.nrm, b.nrm), h: cat(a.h, b.h), nTri: a.nTri + b.nTri };
}


// 球壳来源: 把 "算一层球壳" 包成 meshShells 要的 shellAt(r), 并做沿半径的延拓.
//   每层先不带种子算一遍 (引擎自己的网格种子, 快); 这一层的线比里外两个已算邻层里长的那个短了 8% 以上 (疑似漏了弧) 时,
//   拿里外两个邻层上各条弧的状态点 (每条取十几个, 总数不超过 60) 当额外种子再算一遍 (见 divshell.js traceCaustic), 线变长了才换.
//   compute(r, seeds) → { arcs: [{ p, d, closed }], cusps, states: [[状态点 …] …] }  (边界像没有 states, 也就不延拓)
//   返回的函数另带 again(r): 用现在的邻层带种子重算这一层, 线变长了才换, 返回换没换
export function shellSource(compute) {
  const done = [];                                    // { r, res, len } 按 r 排
  const lenOf = res => { let L = 0; for (const a of res.arcs) { const flat = !Array.isArray(a.d[0]), n = flat ? a.d.length / 3 : a.d.length; for (let i = 1; i < n; i++) { const p = flat ? [a.d[3 * i - 3], a.d[3 * i - 2], a.d[3 * i - 1]] : a.d[i - 1], q = flat ? [a.d[3 * i], a.d[3 * i + 1], a.d[3 * i + 2]] : a.d[i]; L += Math.acos(Math.max(-1, Math.min(1, p[0] * q[0] + p[1] * q[1] + p[2] * q[2]))); } } return L; };
  const near = r => {
    let lo = null, hi = null;
    for (const e of done) { if (e.r < r) lo = e; else if (e.r > r && !hi) hi = e; }
    return [lo, hi];
  };
  const seedsOf = nb => {
    const runs = [];
    for (const e of nb) if (e && e.res.states) for (const run of e.res.states) if (run.length) runs.push(run);
    if (!runs.length) return null;
    const per = Math.max(2, Math.min(16, Math.floor(60 / runs.length))), out = [];
    for (const run of runs) { const st = Math.max(1, Math.floor(run.length / per)); for (let i = st >> 1; i < run.length; i += st) out.push(run[i]); }
    return out;
  };
  const put = (r, res) => { const e = { r, res, len: lenOf(res) }; const i = done.findIndex(x => x.r >= r); if (i >= 0 && done[i].r === r) done[i] = e; else if (i < 0) done.push(e); else done.splice(i, 0, e); return e; };
  const shellAt = r => {
    const old = done.find(x => x.r === r);
    if (old) return old.res;
    let res = compute(r, null);
    // 和里外两个邻层里长的那个比 (二分时离得近的那个邻层常常就是外面空着的一层, 和它比就永远不补)
    const nb = near(r), len = lenOf(res), nl = Math.max(nb[0] ? nb[0].len : 0, nb[1] ? nb[1].len : 0);
    if (len < 0.92 * nl) { const sd = seedsOf(nb); if (sd) { const r2 = compute(r, sd); if (lenOf(r2) > len * 1.001) res = r2; } }
    return put(r, res).res;
  };
  shellAt.again = r => {
    const old = done.find(x => x.r === r), sd = seedsOf(near(r));
    if (!sd) return false;
    const res = compute(r, sd);
    if (old && lenOf(res) <= old.len * 1.001) return false;
    put(r, res);
    return true;
  };
  return shellAt;
}

// 临界半径: 三维边界线 / 脊线 (各是一条折线, 坐标数组长 3n) 上离眼距离 r 的局部极值, 以及每条线的两端.
//   球壳法里层线的拓扑只在这些半径上变 (见 meshShells 的 opt.critical). 灯 / 眼旁 (0.5% D 之内) 的不算.
//   极值带滞后: 从候选极值点往回退了 0.1% r (至少 2e-4 D) 才算真的转向 —— 边界线取样有细小的抖动.
//   显著程度 = 与前后相邻两个极值 (或线的端点) 之差的较小者; 半径相差不到 0.1% 的并成一个 (留显著的), 至多留 maxN 个
export function criticalRadii(feats, L, E, maxN = 64) {
  const D = Math.hypot(L[0] - E[0], L[1] - E[1], L[2] - E[2]), out = [];
  for (const f of feats) {
    const p = f.p || f, n = p.length / 3;
    if (n < 2) continue;
    const r = new Float64Array(n);
    for (let i = 0; i < n; i++) r[i] = Math.hypot(p[3 * i] - E[0], p[3 * i + 1] - E[1], p[3 * i + 2] - E[2]);
    const okAt = i => r[i] > 5e-3 * D && Math.hypot(p[3 * i] - L[0], p[3 * i + 1] - L[1], p[3 * i + 2] - L[2]) > 5e-3 * D;
    // 沿线的极值序列 (含两端)
    const ex = [0];
    let dir = 0, ext = 0;
    for (let i = 1; i < n; i++) {
      const hys = Math.max(1e-3 * r[ext], 2e-4 * D);
      if (dir >= 0 && r[i] > r[ext]) { ext = i; if (!dir && r[i] - r[0] > hys) dir = 1; continue; }
      if (dir <= 0 && r[i] < r[ext]) { ext = i; if (!dir && r[0] - r[i] > hys) dir = -1; continue; }
      if (Math.abs(r[i] - r[ext]) > hys) { if (dir && ext !== ex[ex.length - 1]) ex.push(ext); dir = r[i] > r[ext] ? 1 : -1; ext = i; }
    }
    if (ex[ex.length - 1] !== n - 1) ex.push(n - 1);
    ex.forEach((i, k) => {
      if (!okAt(i)) return;
      const end = k === 0 || k === ex.length - 1;
      const prom = end ? Infinity : Math.min(Math.abs(r[i] - r[ex[k - 1]]), Math.abs(r[i] - r[ex[k + 1]]));
      out.push([r[i], prom]);
    });
  }
  // 半径相差不到 0.1% 的并成一个
  out.sort((a, b) => a[0] - b[0]);
  const merged = [];
  for (const q of out) { const last = merged[merged.length - 1]; if (last && q[0] - last[0] < 1e-3 * q[0]) { if (q[1] > last[1]) { last[0] = q[0]; last[1] = q[1]; } } else merged.push(q.slice()); }
  merged.sort((a, b) => b[1] - a[1]);
  return merged.slice(0, maxN).map(q => q[0]);
}
