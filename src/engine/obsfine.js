// 观测者焦散二自由度折叠的修正版串线 (obscaustic.js 保持与 C# 逐句一致, 这里只改"怎么连"):
//   原版沿 χ_内=0 曲线抽样到 ≤400 点, 并按 solveInner 返回的【序号】把内层解串成线 —
//   解的个数一变 (出现/消失), 序号就错位, 线会接到别的分支上 (长直弦、乱连).
//   这里: 抽样更密, 每个开着的段按"同一 α 支 + 内层角最接近 + 天空方向最接近"去接新点 (连续性匹配),
//   接不上的段收尾, 没人接的新点开新段. 输出与原版同格式 (Seg: p 为相对眼的位置, d 为视方向).

import * as OC from './obscaustic.js';

const wrapPi = a => a - 2 * Math.PI * Math.round(a / (2 * Math.PI));
const ang = (a, b) => { const c = a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; return Math.acos(c > 1 ? 1 : c < -1 ? -1 : c); };

export function fold2DofFine(c2, traceH = 0.01, maxPts = 1800, cancelled = null) {
  const probe = c2.inner(0.0);
  const f = OC.chiFieldOf(probe, c2.omegaInner());
  const curves = OC.traceZeros(f, traceH, 6000, 40, 96, 1440, cancelled);
  const segs = [];
  const close = o => { if (o.seg.p.length >= 2) segs.push(o.seg); };
  for (const cur of curves) {
    if (cancelled && cancelled()) break;
    const step = Math.max(1, Math.trunc(cur.length / Math.max(1, maxPts)));
    const dXi = traceH * step;                       // 相邻样本在 ξ 上的间距 (rad)
    const tolB = Math.max(0.12, 6 * dXi), tolD = Math.max(0.05, 8 * dXi);
    let open = [];
    for (let i = 0; i < cur.length; i += step) {
      if (cancelled && cancelled()) break;
      const xi = cur[i];
      const cand = [];
      for (const b of OC.solveInner(xi, c2)) {
        const cc = c2.inner(b);
        const brs = OC.pairAlpha(xi, cc);
        for (let k = 0; k < 2; k++)
          if (brs[k].vis && brs[k].r <= cc.maxDist) cand.push({ b, k, br: brs[k], cc });
      }
      const used = new Uint8Array(cand.length), next = [];
      for (const o of open) {
        let best = -1, bd = Infinity;
        for (let ci = 0; ci < cand.length; ci++) {
          const c = cand[ci];
          if (used[ci] || c.k !== o.k) continue;
          const db = Math.abs(wrapPi(c.b - o.b)), da = ang(c.br.d, o.d);
          if (db > tolB || da > tolD) continue;
          const sc = db / tolB + da / tolD;
          if (sc < bd) { bd = sc; best = ci; }
        }
        if (best < 0) { close(o); continue; }
        used[best] = 1;
        add(o, cand[best], xi);
        next.push(o);
      }
      for (let ci = 0; ci < cand.length; ci++) {
        if (used[ci]) continue;
        const c = cand[ci], o = { k: c.k, b: c.b, d: c.br.d, seg: new OC.Seg('K1-2dof', c.k) };
        add(o, c, xi);
        next.push(o);
      }
      open = next;
    }
    open.forEach(close);
  }
  return segs;
}

function add(o, c, xi) {
  const br = c.br, sg = o.seg;
  sg.xi.push(xi); sg.d.push(br.d); sg.alpha.push(br.alpha);
  sg.r.push(br.r); sg.t.push(br.t); sg.inner.push(c.b);
  sg.p.push([br.d[0] * br.r, br.d[1] * br.r, br.d[2] * br.r]);
  if (br.r > 0.98 * c.cc.maxDist) sg.truncated = true;
  o.b = c.b; o.d = br.d;
}
