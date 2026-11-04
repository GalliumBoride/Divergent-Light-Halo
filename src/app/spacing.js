// 三维折线的间距工具 (主线程与 Worker 共用): 到一层折线的距离查询, 补线 (夹在两层之间) 的逐点淡入与修剪. 用法见 main.js 的"补线".
export const FILL_LO = 1.4, FILL_HI = 2.4;   // 补线在两侧间距 1.4d* 处开始淡入, 2.4d* 处全显
const _segIx = new WeakMap();
// 到一层折线的距离函数, 只精确到 cap (更远的一律给 cap): 线段按 cap 宽的网格分桶 (整数哈希, 撞桶只多算几段), 查 3×3×3 邻格
export function distTo(curves, cap) {
  const e = _segIx.get(curves);
  if (e && e.cap === cap) return e.f;
  let ns = 0;
  for (const c of curves) ns += Math.max(0, (c.p || c).length / 3 - 1);
  const S6 = new Float64Array(6 * ns), map = new Map(), inv = 1 / cap;
  const key = (i, j, k) => (Math.imul(i, 73856093) ^ Math.imul(j, 19349663) ^ Math.imul(k, 83492791)) | 0;
  let s = 0;
  for (const c of curves) {
    const p = c.p || c;
    for (let i = 0; i + 5 < p.length; i += 3, s++) {
      for (let k = 0; k < 6; k++) S6[6 * s + k] = p[i + k];
      const n = Math.max(1, Math.ceil(Math.hypot(p[i + 3] - p[i], p[i + 4] - p[i + 1], p[i + 5] - p[i + 2]) * 2 * inv));
      let prev = null;
      for (let k = 0; k <= n; k++) {
        const t = k / n, kk = key(Math.floor((p[i] + (p[i + 3] - p[i]) * t) * inv), Math.floor((p[i + 1] + (p[i + 4] - p[i + 1]) * t) * inv),
          Math.floor((p[i + 2] + (p[i + 5] - p[i + 2]) * t) * inv));
        if (kk === prev) continue;
        prev = kk;
        const a = map.get(kk);
        if (!a) map.set(kk, [s]); else if (a[a.length - 1] !== s) a.push(s);
      }
    }
  }
  const f = (x, y, z) => {
    const cx = Math.floor(x * inv), cy = Math.floor(y * inv), cz = Math.floor(z * inv);
    let m2 = cap * cap;
    for (let i = -1; i <= 1; i++) for (let j = -1; j <= 1; j++) for (let k = -1; k <= 1; k++) {
      const a = map.get(key(cx + i, cy + j, cz + k));
      if (a) for (const q of a) {
        const o = 6 * q, ux = S6[o + 3] - S6[o], uy = S6[o + 4] - S6[o + 1], uz = S6[o + 5] - S6[o + 2];
        const wx = x - S6[o], wy = y - S6[o + 1], wz = z - S6[o + 2], l2 = ux * ux + uy * uy + uz * uz;
        const t = l2 > 0 ? Math.max(0, Math.min(1, (ux * wx + uy * wy + uz * wz) / l2)) : 0;
        const dx = wx - t * ux, dy = wy - t * uy, dz = wz - t * uz, d2 = dx * dx + dy * dy + dz * dz;
        if (d2 < m2) m2 = d2;
      }
    }
    return Math.sqrt(m2);
  };
  _segIx.set(curves, { cap, f });
  return f;
}
// 补的层 F (夹在 pa, pb 之间): 逐点不透明度 a = smoothstep(g; 1.4d*, 2.4d*), 只留 a > 0 的段 (连同两端各一个 a = 0 的点, 淡出到底)
export function trimFill(F, d) {
  const fa = distTo(F.pa.curves, FILL_HI * d), fb = distTo(F.pb.curves, FILL_HI * d), out = [];
  for (const c of F.curves) {
    const p = c.p, n = p.length / 3, al = new Float32Array(n);
    for (let i = 0; i < n; i++) {
      const g = fa(p[3 * i], p[3 * i + 1], p[3 * i + 2]) + fb(p[3 * i], p[3 * i + 1], p[3 * i + 2]);
      const t = Math.max(0, Math.min(1, (g / d - FILL_LO) / (FILL_HI - FILL_LO)));
      al[i] = t * t * (3 - 2 * t);
    }
    let i = 0;
    while (i < n) {
      if (!(al[i] > 0)) { i++; continue; }
      let a = i, b = i, top = 0;
      while (b + 1 < n && al[b + 1] > 0) b++;
      i = b + 1;
      for (let k = a; k <= b; k++) top = Math.max(top, al[k]);
      if (a > 0) a--;
      if (b + 1 < n) b++;
      if (b > a && top > 0.05) out.push({ p: p.slice(3 * a, 3 * b + 3), d: c.d.slice(3 * a, 3 * b + 3), a: al.slice(a, b + 1) });
    }
  }
  return out;
}
