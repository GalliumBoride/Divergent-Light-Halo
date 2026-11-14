// 页面的层安排搬到 node: 第一轮 nH = 72 层 (h 等分) + runDense 按三维层距等分布补 150 层 (main.js 的 layerGap / runDense 原样).
export function layerGap(A, B, how = 'mean', q = 0.75) {
  const pts = (L, cap) => {
    let n = 0; for (const c of L.curves) n += (c.p || c).length / 3;
    const st = Math.max(1, Math.ceil(n / cap)), out = [];
    let k = 0;
    for (const c of L.curves) { const p = c.p || c; for (let i = 0; i < p.length; i += 3, k++) if (k % st === 0) out.push(p[i], p[i + 1], p[i + 2]); }
    return out;
  };
  const one = (a, b) => {
    if (!a.length || !b.length) return 0;
    const ds = [];
    for (let i = 0; i < a.length; i += 3) {
      let m = Infinity;
      for (let j = 0; j < b.length; j += 3) { const dx = a[i] - b[j], dy = a[i + 1] - b[j + 1], dz = a[i + 2] - b[j + 2], d = dx * dx + dy * dy + dz * dz; if (d < m) m = d; }
      ds.push(Math.sqrt(m));
    }
    ds.sort((x, y) => x - y);
    return ds[Math.floor(q * (ds.length - 1))];
  };
  const g1 = one(pts(A, 150), pts(B, 1500));
  return how === 'from' ? g1 : 0.5 * (g1 + one(pts(B, 150), pts(A, 1500)));
}
export function pageLayers(EP, lp, el, D, low, z, th, surf, nH = 72, total = 150, nG = 101, maxD = 80) {
  const first = EP.lightCone2DofLayers(lp, el, D, 0, nH, nG, maxD, low, z, th, surf, 0, nH, null, true);
  const hOf = ih => -90 + 180 * ih / Math.max(1, nH - 1), clamp = (x, a, b) => Math.max(a, Math.min(b, x));
  let i0 = -1, i1 = -1;
  first.forEach((l, ih) => { if (l && l.curves.length) { if (i0 < 0) i0 = ih; i1 = ih; } });
  const all = [];
  first.forEach((l, ih) => { if (l && l.curves.length) all.push({ h: hOf(ih), curves: l.curves }); });
  if (i0 < 0) return all;
  const step = 180 / Math.max(1, nH - 1), kA = Math.max(0, i0 - 1), kB = Math.min(nH - 1, i1 + 1);
  const has = k => !!(first[k] && first[k].curves.length), g = [];
  for (let k = kA; k < kB; k++) g.push(has(k) && has(k + 1) ? layerGap(first[k], first[k + 1]) : -1);
  const known = g.filter(x => x > 0), gMed = known.length ? known.sort((x, y) => x - y)[Math.floor(known.length / 2)] : 1;
  for (let k = 0; k < g.length; k++) if (g[k] < 0) g[k] = gMed;
  const gSum = g.reduce((x, y) => x + y, 0) || 1, hs = [];
  g.forEach((gk, k) => { const n = Math.min(60, Math.round(total * gk / gSum)), h0 = hOf(kA + k); for (let j = 1; j <= n; j++) hs.push(clamp(h0 + step * j / (n + 1), -89.99, 89.99)); });
  for (const h of hs) { const c = EP.lightCone2DofLayers(lp, el, D, 0, 1, nG, maxD, low, z, th, surf, 0, 1, h, true)[0].curves; if (c.length) all.push({ h, curves: c }); }
  all.sort((a, b) => a.h - b.h);
  return all;
}
// 缝合 + 页面那一轮补层 (bad 层对的中点). feats: 特征线 (给缝合器); fixed: 不补层 (HF 给定的 h 列表已含补的层)
export function pageMesh(EP, ST, lp, el, D, low, z, th, surf, L, layers = null, feats = null, fixed = false) {
  layers = layers || pageLayers(EP, lp, el, D, low, z, th, surf);
  let m = ST.stitchLayers2(layers, L, D, feats);
  if (m.bad.length && !fixed) {
    for (const [a, b] of m.bad) { const h = 0.5 * (a + b), c = EP.lightCone2DofLayers(lp, el, D, 0, 1, 101, 80, low, z, th, surf, 0, 1, h, true)[0].curves; if (c.length) layers.push({ h, curves: c }); }
    layers.sort((a, b) => a.h - b.h);
    m = ST.stitchLayers2(layers, L, D, feats);
  }
  return { m, layers };
}
