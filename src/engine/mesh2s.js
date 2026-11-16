// 球壳法网格的辅助函数 (网格本身见 shellmesh.js):
//   criticalRadii: 三维边界线 / 脊线离眼距离的极值与端点 —— 层线拓扑只在这些半径上变;
//   skyStepOf: 网格点距 (按平行光弧的天空长度与顶点预算);
//   coverCheck: 测试用, 一组曲线上的点有多少离网格太远.
//   (之前的球壳网格器 meshShells / shellSource 按 "天空最近点" 配对拉链, 已由 shellmesh.js 的特殊点骨架取代, 删去)

const DEG = Math.PI / 180;

// 临界半径: 三维边界线 / 脊线 (各是一条折线, 坐标数组长 3n) 上离眼距离 r 的局部极值, 以及每条线的两端.
//   球壳法里层线的拓扑只在这些半径上变 (shellmesh.js 在每个临界半径两侧各放一层). 灯 / 眼旁 (0.5% D 之内) 的不算.
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

// 网格点距 (弧度): 平行光弧 (r → 0 时的层线; Float32Array 折线, 点在以眼为心的大球面上) 的天空长度 × 约 200 层 ≤ 预算 (顶点数),
//   至少 0.35°、至多 maxDeg (默认 0.8°). 用平行光弧而不用某一段的层线: 分段并行时各段算出来的一样, 密度不在段间跳
export function skyStepOf(arcs, E, budget, maxDeg = 0.8) {
  let len = 0;
  for (const c of arcs || []) {
    const p = c.p || c;
    for (let i = 3; i + 2 < p.length; i += 3) {
      const ax = p[i - 3] - E[0], ay = p[i - 2] - E[1], az = p[i - 1] - E[2], bx = p[i] - E[0], by = p[i + 1] - E[1], bz = p[i + 2] - E[2];
      const la = Math.hypot(ax, ay, az), lb = Math.hypot(bx, by, bz);
      if (la > 0 && lb > 0) len += Math.acos(Math.max(-1, Math.min(1, (ax * bx + ay * by + az * bz) / (la * lb))));
    }
  }
  const DEG = Math.PI / 180;
  return Math.min(maxDeg * DEG, Math.max(0.35 * DEG, len * 200 / Math.max(1, budget)));
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
