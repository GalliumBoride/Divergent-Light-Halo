// 球壳法分段并行: 以眼为心的半径范围 [1e-3 D, maxDist + D] 按对数切成 k 段, 每段一个 Worker (各段的网格在交界那层上接上).
//   rHi > 0 (曲面所在的最远半径, 取自边界线) 时, 它以外 (多半没有曲面) 只占一段. 段的交界避开临界半径 (相对 ±0.3%)
//   和灯所在的球壳 (±0.5%) —— 那里两侧的一对层要落在同一段里. → [[r0, r1, 初始层数, 预算份额] …]
export function shellChunks(L, E, crit, k, rHi, maxDist) {
  const D = Math.hypot(L[0] - E[0], L[1] - E[1], L[2] - E[2]), rMin = 1e-3 * D, rMax = maxDist + D;
  const top = rHi > rMin * 10 && rHi * 1.2 < rMax * 0.8 ? rHi * 1.2 : rMax;
  const nIn = top < rMax ? Math.max(1, k - 1) : Math.max(1, k), cuts = [rMin];
  for (let i = 1; i < nIn; i++) cuts.push(rMin * Math.pow(top / rMin, i / nIn));
  if (top < rMax) cuts.push(top);
  cuts.push(rMax);
  for (let i = 1; i + 1 < cuts.length; i++) {
    let b = cuts[i];
    for (let it = 0; it < 4; it++) {
      const c = crit.find(rc => Math.abs(b / rc - 1) < 3e-3);
      if (c !== undefined) b = c * (b < c ? 1 - 3e-3 : 1 + 3e-3);
      else if (Math.abs(b / D - 1) < 5e-3) b = D * (b < D ? 1 - 5e-3 : 1 + 5e-3);
      else break;
    }
    cuts[i] = b;
  }
  const out = [], all = Math.log(rMax / rMin);
  for (let i = 0; i + 1 < cuts.length; i++) {
    if (!(cuts[i + 1] > cuts[i] * 1.001)) continue;
    const f = Math.log(cuts[i + 1] / cuts[i]) / all;
    out.push([cuts[i], cuts[i + 1], Math.max(6, Math.round(72 * f)) + 1, f]);
  }
  return out;
}
// 各段的网格 / 层线拼在一起 (还没回来的段跳过: 边算边画). 交界层两段各给一次, 留一条. 一段也没成 → null
export function mergeShell(parts) {
  const ok = parts.filter(d => d && d.ok && (d.nTri > 0 || (d.special && d.special.length)));
  if (!ok.length) return null;
  let n = 0, nh = 0;
  for (const d of ok) { n += d.pos.length; nh += d.h.length; }
  const pos = new Float32Array(n), dir = new Float32Array(n), nrm = new Float32Array(n), h = new Float32Array(nh);
  let o = 0, oh = 0, rows = 0, loose = 0, ms = 0, open = 0, filled = 0;
  for (const d of ok) {
    pos.set(d.pos, o); dir.set(d.dir, o); nrm.set(d.nrm, o); h.set(d.h, oh);
    o += d.pos.length; oh += d.h.length; rows += d.rows || 0; loose += d.loose || 0; ms = Math.max(ms, d.ms || 0);
    open += d.open || 0; filled += d.filled || 0;      // 网格自检 (shellmesh.js): 不在真边界上的开口边 (度) / 补上的小洞数
  }
  const all = [];
  for (const d of ok) for (const l of d.lines) all.push(l);
  all.sort((a, b) => a.r - b.r);
  const lines = all.filter((l, i) => i === 0 || Math.abs(l.r / all[i - 1].r - 1) > 1e-9);
  const special = [];
  for (const d of ok) for (const f of d.special || []) special.push(f);
  return { ok: true, pos, dir, nrm, h, nTri: n / 9, lines, special, rows, loose, ms, open, filled, parts: ok.length };
}
