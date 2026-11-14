// 二自由度曲面网格: 相邻两层层线缝成三角形 —— 一自由度网格"行 × 列"的做法搬过来.
//   行 = 每个 h 的层线 (端点已精确落在边界上, 折叠 / 尖点 / 接缝都已处理好); 列 = 沿层线的位置.
//   一自由度有共同的 α 对齐上下两行; 二自由度没有, 就用取向平面 (u, v) 里的最近点对应: 相邻两层 h 只差零点几度,
//   层曲线在 (u, v) 里只挪一点: 用 (u, v) 选对应的曲线与支 (两支在 (u, v) 里重合, 按三维分开, 见 foot), 沿线位置取三维最近
//   (见 descend). 对应好的两段用拉链式三角化缝起来 (每步走较短的那条对角线, Fuchs–Kedem–Uselton 的贪心版).
//   对不上的 (曲线在这两层之间出生 / 消亡 / 分叉) 留空; 下层没盖到的上层段再从上往下缝一遍.
// 特征线约束 (带特征线的轮廓缝合): 精确边界线、转折缝、尖点脊线、外缘必须是网格的边.
//   层线落在特征线上的顶点 (端点 / 折回点 / 折角) 是 "路标"; 相邻两层沿同一条特征线挨着的路标强制对应, 缝合在路标处断开
//   (两侧各缝各的, 两个路标的连线成了两段共用的边); 两个路标之间特征线弯得厉害时把线上的点插进网格 (见 arcOf).
//   不这样做: 层线在尖点脊线上折回处两臂在三维里重合, 贪心拉链把 A 的尖点连到 B 尖点前面的点上、再从回程臂扇回来,
//   脊线旁一排底边在脊线上的三角形缺口; 斜穿层线的折痕被三角形抄近路, 脊线一截埋在面下一截露出来.
// 输入 layers: [{ h, curves: [{ p, d, x, b }] }] (按 h 排好; p 位置, d 视方向, x 每点 (u, v), b 支); L 灯位置, D 灯眼距离;
//   feats: [{ p, d }] 特征线折线 (可省: 那就只按几何认折回点, 没有别的约束).
// 输出与一自由度网格同格式: { pos, dir, h, nrm } (三角形顶点不共享), 外加 nTri.

const TWO = 2 * Math.PI;
const wd = a => a - TWO * Math.round(a / TWO);
const CELL = 0.04, RUV = 0.12;          // (u, v) 分桶宽度与对应的最大距离 (弧度)

// 层线为画线加得很密 (灯旁、尖点旁), 网格用不着: 沿线保留与上一个保留点三维距离 ≥ smin(r)、或从眼看张角 ≥ 0.5°、或转角 ≥ 20° 的点,
//   端点总留. smin 随离灯距离 r 放宽 (层距本身 ∝ r, 远处按固定间距会切出细长三角形)
function thin(c, smin, must = null) {
  const p = c.p, d = c.d, x = c.x, n = p.length / 3, keep = [0];
  const cosA = Math.cos(0.5 * Math.PI / 180), cosT = Math.cos(20 * Math.PI / 180);
  // 折回点 (层线在尖点脊线上折回): 在原始的密折线上判 —— 尖点处 P(s) ≈ a s² + b s³, 两臂只在毫米尺度上才反向平行,
  //   抽稀到几厘米后夹角只剩七八十度. 一个顶点处转角 > 120°; 或折回落在两个取样点之间 (一小段横着的短边, 两头各转九十来度):
  //   短边前后两段反向, 取短边的前一头
  const seg = j => [p[3 * j + 3] - p[3 * j], p[3 * j + 4] - p[3 * j + 1], p[3 * j + 5] - p[3 * j + 2]];
  const cosv = (a, b) => (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) / (Math.hypot(a[0], a[1], a[2]) * Math.hypot(b[0], b[1], b[2]) || 1);
  const tip = new Uint8Array(n);
  for (let j = 1; j < n - 1; j++) {
    const a = seg(j - 1), b = seg(j);
    if (cosv(a, b) < -0.5) { if (!tip[j - 1]) tip[j] = 1; continue; }
    if (j + 2 < n) {
      const c2 = seg(j + 1), lb = Math.hypot(b[0], b[1], b[2]);
      if (cosv(a, c2) < -0.5 && lb < 0.5 * Math.min(Math.hypot(a[0], a[1], a[2]), Math.hypot(c2[0], c2[1], c2[2])) && !tip[j - 1]) tip[j] = 1;
    }
  }
  for (let j = 1; j < n - 1; j++) {
    const i = keep[keep.length - 1];
    const ds = Math.hypot(p[3 * j] - p[3 * i], p[3 * j + 1] - p[3 * i + 1], p[3 * j + 2] - p[3 * i + 2]);
    const ca = d[3 * j] * d[3 * i] + d[3 * j + 1] * d[3 * i + 1] + d[3 * j + 2] * d[3 * i + 2];
    const ax = p[3 * j] - p[3 * j - 3], ay = p[3 * j + 1] - p[3 * j - 2], az = p[3 * j + 2] - p[3 * j - 1];
    const bx = p[3 * j + 3] - p[3 * j], by = p[3 * j + 4] - p[3 * j + 1], bz = p[3 * j + 5] - p[3 * j + 2];
    const la = Math.hypot(ax, ay, az), lb = Math.hypot(bx, by, bz);
    const turn = la > 0 && lb > 0 && (ax * bx + ay * by + az * bz) / (la * lb) < cosT;
    if (ds >= smin(p[3 * i], p[3 * i + 1], p[3 * i + 2]) || ca < cosA || turn || tip[j] || (must && must[j])) keep.push(j);
  }
  if (n > 1) keep.push(n - 1);
  const tips = [];
  keep.forEach((j, k) => { if (tip[j]) tips.push(k); });
  if (keep.length === n) return { p, d, x, b: c.b, tips, keep };
  const P = new Float32Array(3 * keep.length), Dd = new Float32Array(3 * keep.length), X = new Float32Array(2 * keep.length);
  keep.forEach((j, k) => { P.set(p.subarray(3 * j, 3 * j + 3), 3 * k); Dd.set(d.subarray(3 * j, 3 * j + 3), 3 * k); X.set(x.subarray(2 * j, 2 * j + 2), 2 * k); });
  return { p: P, d: Dd, x: X, b: c.b, tips, keep };
}

// 同一层、同一支里首尾相接的曲线接成一条: 层线在转折缝 (可达带边, ε± 两半汇合) 处是两条, 端点重合;
//   分开缝的话 U 形尖端两侧的楔形谁也盖不到 (一排锯齿状的缺口)
function joinCurves(cs, tol) {
  const out = cs.map(c => ({ p: Array.from(c.p), d: Array.from(c.d), x: Array.from(c.x), b: c.b || 0 }));
  const end = (c, last) => { const k = last ? c.p.length - 3 : 0; return [c.p[k], c.p[k + 1], c.p[k + 2]]; };
  const near = (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]) < tol;
  const rev = c => {
    const n = c.p.length / 3, p = [], d = [], x = [];
    for (let j = n - 1; j >= 0; j--) { p.push(c.p[3 * j], c.p[3 * j + 1], c.p[3 * j + 2]); d.push(c.d[3 * j], c.d[3 * j + 1], c.d[3 * j + 2]); x.push(c.x[2 * j], c.x[2 * j + 1]); }
    c.p = p; c.d = d; c.x = x;
  };
  for (let again = true; again;) {
    again = false;
    for (let i = 0; i < out.length && !again; i++) for (let j = 0; j < out.length && !again; j++) {
      if (i === j || out[i].b !== out[j].b) continue;
      const a = out[i], b = out[j];
      if (near(end(a, true), end(b, true))) rev(b);
      else if (near(end(a, false), end(b, false))) rev(a);
      else if (near(end(a, false), end(b, true))) { out[i] = b; out[j] = a; continue; }
      if (!near(end(out[i], true), end(out[j], false))) continue;
      const A = out[i], Bc = out[j];
      A.p.push(...Bc.p.slice(3)); A.d.push(...Bc.d.slice(3)); A.x.push(...Bc.x.slice(2));
      out.splice(j, 1); again = true;
    }
  }
  return out.map(c => ({ p: Float32Array.from(c.p), d: Float32Array.from(c.d), x: Float32Array.from(c.x), b: c.b }));
}

// 特征线 (精确边界线 / 尖点脊线 / 转折缝 / 外缘) 的线段分桶: near(x, y, z, tol) → 容差内的各条特征线 [{ f, t (沿线弧长), d }]
function featIndex(feats, L, D) {
  const cell = Math.max(0.05 * D, 0.03), grid = new Map();
  const key = (i, j, k) => ((i + 4096) * 8192 + (j + 4096)) * 8192 + (k + 4096);
  const tol = r => 1.5e-4 * D + 3e-4 * r;
  const F = feats.map(f => {
    const p = f.p || f, n = p.length / 3, t = new Float64Array(n);
    for (let i = 1; i < n; i++) t[i] = t[i - 1] + Math.hypot(p[3 * i] - p[3 * i - 3], p[3 * i + 1] - p[3 * i - 2], p[3 * i + 2] - p[3 * i - 1]);
    return { p, t, d: f.d || null };
  });
  F.forEach((f, fi) => {
    const p = f.p;
    for (let i = 0; i + 1 < p.length / 3; i++) {
      const o = 3 * i, pad = 1.2 * tol(Math.max(Math.hypot(p[o] - L[0], p[o + 1] - L[1], p[o + 2] - L[2]), Math.hypot(p[o + 3] - L[0], p[o + 4] - L[1], p[o + 5] - L[2])));
      const lo = [0, 1, 2].map(a => Math.floor((Math.min(p[o + a], p[o + 3 + a]) - pad) / cell)), hi = [0, 1, 2].map(a => Math.floor((Math.max(p[o + a], p[o + 3 + a]) + pad) / cell));
      for (let a = lo[0]; a <= hi[0]; a++) for (let b = lo[1]; b <= hi[1]; b++) for (let c = lo[2]; c <= hi[2]; c++) {
        const k = key(a, b, c); let q = grid.get(k); if (!q) grid.set(k, q = []);
        q.push(fi, i);
      }
    }
  });
  const near = (x, y, z, tl) => {
    const q = grid.get(key(Math.floor(x / cell), Math.floor(y / cell), Math.floor(z / cell)));
    if (!q) return null;
    let out = null;
    for (let s = 0; s < q.length; s += 2) {
      const fi = q[s], i = q[s + 1], p = F[fi].p, o = 3 * i;
      const bx = p[o + 3] - p[o], by = p[o + 4] - p[o + 1], bz = p[o + 5] - p[o + 2], px = x - p[o], py = y - p[o + 1], pz = z - p[o + 2];
      const l2 = bx * bx + by * by + bz * bz;
      let u = l2 > 0 ? (px * bx + py * by + pz * bz) / l2 : 0; u = u < 0 ? 0 : u > 1 ? 1 : u;
      const d = Math.hypot(px - u * bx, py - u * by, pz - u * bz);
      if (d > tl) continue;
      if (!out) out = [];
      const e = out.find(w => w.f === fi), t = F[fi].t[i] + u * Math.sqrt(l2);
      if (!e) out.push({ f: fi, t, d }); else if (d < e.d) { e.d = d; e.t = t; }
    }
    return out;
  };
  return { near, tol, F };
}

export function stitchLayers2(layers, L, D, feats = null) {
  // 路标: 层线落在特征线上的顶点 (端点 / 折回点 / 折角). 同一条特征线上连着几个顶点都在容差内时取最近的一个.
  //   路标在抽稀时必须保留; 相邻两层沿同一条特征线挨着的路标强制连成网格的边 (见 pairMarks / sew)
  const FT = feats && feats.length ? featIndex(feats, L, D) : null;
  const landmarksOf = c => {
    const p = c.p, n = p.length / 3, hit = new Array(n);
    let any = false;
    for (let j = 0; j < n; j++) {
      const x = p[3 * j], y = p[3 * j + 1], z = p[3 * j + 2];
      const q = FT.near(x, y, z, FT.tol(Math.hypot(x - L[0], y - L[1], z - L[2])));
      if (q) { hit[j] = q; any = true; }
    }
    if (!any) return null;
    const lm = new Map();
    for (let j = 0; j < n; j++) if (hit[j]) for (const q of hit[j]) {
      if (j > 0 && hit[j - 1] && hit[j - 1].some(e => e.f === q.f)) continue;      // 不是这一串的开头
      let bj = j, bd = q.d, bt = q.t;
      for (let k = j + 1; k < n && hit[k]; k++) { const e = hit[k].find(w => w.f === q.f); if (!e) break; if (e.d < bd) { bd = e.d; bj = k; bt = e.t; } }
      let a = lm.get(bj); if (!a) lm.set(bj, a = []);
      a.push({ f: q.f, t: bt });
    }
    // 两条特征线在这里并在一起走 (外缘贴着接缝 …): 各自最近的顶点可能差一两个, 相距不到两倍容差的并成一个路标
    const js = [...lm.keys()].sort((a, b) => a - b);
    for (let i = 1, k = 0; i < js.length; i++) {
      const a = js[k], b = js[i], x = p[3 * a], y = p[3 * a + 1], z = p[3 * a + 2];
      if (b - a <= 4 && Math.hypot(p[3 * b] - x, p[3 * b + 1] - y, p[3 * b + 2] - z) < 2 * FT.tol(Math.hypot(x - L[0], y - L[1], z - L[2]))) { lm.get(a).push(...lm.get(b)); lm.delete(b); }
      else k = i;
    }
    return lm;
  };
  // 全部顶点编号; 每条曲线记首顶点号
  let nv = 0;
  const smin = (x, y, z) => Math.max(D / 200, 0.015 * Math.hypot(x - L[0], y - L[1], z - L[2]));
  const rows = layers.map(l => ({ h: l.h, cs: joinCurves(l.curves.filter(c => c.x && c.p.length >= 6), 1e-4 * D)
    .map(c0 => {
      const lm0 = FT ? landmarksOf(c0) : null;
      let must = null;
      if (lm0) { must = new Uint8Array(c0.p.length / 3); for (const j of lm0.keys()) must[j] = 1; }
      const c = thin(c0, smin, must), o = { n: c.p.length / 3, base: nv, c, lm: [] };
      if (lm0) c.keep.forEach((j, k) => { const a = lm0.get(j); if (a) o.lm.push({ j: k, fs: a }); });
      nv += o.n; return o;
    }) }));
  // 顶点表留出余量: 特征线上的点可以插进网格 (见 arcOf), 最多每条特征线每个点一次
  let nvT = nv;
  const nvCap = nv + (FT ? FT.F.reduce((a, f) => a + f.p.length / 3, 0) : 0);
  const VP = new Float64Array(3 * nvCap), VD = new Float32Array(3 * nvCap), VH = new Float32Array(nvCap);
  for (const r of rows) for (const o of r.cs) { VP.set(o.c.p, 3 * o.base); VD.set(o.c.d, 3 * o.base); VH.fill(r.h, o.base, o.base + o.n); }
  const d3 = (i, j) => Math.hypot(VP[3 * i] - VP[3 * j], VP[3 * i + 1] - VP[3 * j + 1], VP[3 * i + 2] - VP[3 * j + 2]);
  // 每个顶点沿所在曲线的累计弧长
  const ARC = new Float64Array(nv);
  for (const r of rows) for (const o of r.cs) for (let j = 1; j < o.n; j++) ARC[o.base + j] = ARC[o.base + j - 1] + d3(o.base + j - 1, o.base + j);
  // 层间允许的三维距离 (只挡荒唐的远距对应): 层距大致 ∝ 离灯距离, 远处一瓣锥映射沿层线方向拉得很长, 斜着的对应也合法
  //   (拉链按三维较短的对角线走, 自己把斜切补正), 所以放得很宽: 0.5 r + 0.05 D
  const gap = i => 0.5 * Math.hypot(VP[3 * i] - L[0], VP[3 * i + 1] - L[1], VP[3 * i + 2] - L[2]) + 0.05 * D;
  const tris = [];

  // 折回点 (层线在尖点脊线上折回: 相邻两段夹角 > 120°). 折回处两臂在三维里几乎重合, 拉链按 "较短对角线" 走会把 A 的尖点
  //   连到 B 尖点前面的点上、再从 A 的回程臂扇回来 —— 两层尖点之间的弦不成网格的边, 脊线旁留一排底边在脊线上的三角形缺口.
  //   做法同带特征线的轮廓缝合: 相邻两层的折回点配对 (取向平面里最近; 尖点在取向平面里随 h 连续移动), 配上的强制对应,
  //   缝合在这里断成两段 (每臂各缝各的), 两层尖点的连线就成了两段共用的边.
  for (const r of rows) for (const o of r.cs) o.tips = o.c.tips;
  // 配上对的两个路标之间那段特征线: 弯得厉害时 (远处边界线几乎顺着层线走, 相邻两层的路标沿线隔十几米, 弦离真边界几十厘米)
  //   把这一段上的点插进网格, 接在 A 侧那一串的后面 (或前面) 一起拉链 —— 网格的边就沿着真边界走, 折痕两侧共用这些点.
  //   arcOf(va, vb) → 从 va 到 vb (都不含) 的新顶点号; 弦高不到两倍容差的不插
  const linkF = new Map(), arcs = new Map();
  const arcOf = (va, vb) => {
    const lo = Math.min(va, vb), hi = Math.max(va, vb), key = lo * 67108864 + hi;
    let a = arcs.get(key);
    if (a === undefined) {
      a = [];
      const info = linkF.get(key);
      if (info && FT) {
        const f = FT.F[info.f], n = f.p.length / 3, t0 = Math.min(info.tLo, info.tHi), t1 = Math.max(info.tLo, info.tHi), eps = 1e-3 * (t1 - t0);
        const ax = VP[3 * lo], ay = VP[3 * lo + 1], az = VP[3 * lo + 2], bx = VP[3 * hi] - ax, by = VP[3 * hi + 1] - ay, bz = VP[3 * hi + 2] - az, l2 = bx * bx + by * by + bz * bz;
        let i0 = 0, i1 = n - 1;
        { let a0 = 0, b0 = n - 1; while (b0 - a0 > 1) { const m = (a0 + b0) >> 1; if (f.t[m] <= t0 + eps) a0 = m; else b0 = m; } i0 = f.t[a0] > t0 + eps ? a0 : b0; }
        const idx = [];
        let sag = 0, lx = ax, ly = ay, lz = az;
        if (info.tLo > info.tHi) { lx = VP[3 * hi]; ly = VP[3 * hi + 1]; lz = VP[3 * hi + 2]; }
        for (let i = i0; i < n && f.t[i] < t1 - eps; i++) {
          const x = f.p[3 * i], y = f.p[3 * i + 1], z = f.p[3 * i + 2];
          const u = l2 > 0 ? Math.max(0, Math.min(1, ((x - ax) * bx + (y - ay) * by + (z - az) * bz) / l2)) : 0;
          sag = Math.max(sag, Math.hypot(x - ax - u * bx, y - ay - u * by, z - az - u * bz));
          if (Math.hypot(x - lx, y - ly, z - lz) >= Math.min(smin(x, y, z), 0.25 * (t1 - t0))) { idx.push(i); lx = x; ly = y; lz = z; }
        }
        // 最后一个离终点太近的去掉
        if (idx.length) { const i = idx[idx.length - 1], e = info.tLo > info.tHi ? lo : hi; if (Math.hypot(f.p[3 * i] - VP[3 * e], f.p[3 * i + 1] - VP[3 * e + 1], f.p[3 * i + 2] - VP[3 * e + 2]) < 0.5 * Math.min(smin(f.p[3 * i], f.p[3 * i + 1], f.p[3 * i + 2]), 0.25 * (t1 - t0))) idx.pop(); }
        if (sag > 2 * FT.tol(Math.hypot(ax - L[0], ay - L[1], az - L[2])) && idx.length && nvT + idx.length <= nvCap) {
          if (info.tLo > info.tHi) idx.reverse();           // 从 lo 到 hi
          for (const i of idx) {
            const v = nvT++, w = (f.t[i] - info.tLo) / (info.tHi - info.tLo);
            VP[3 * v] = f.p[3 * i]; VP[3 * v + 1] = f.p[3 * i + 1]; VP[3 * v + 2] = f.p[3 * i + 2];
            if (f.d) { VD[3 * v] = f.d[3 * i]; VD[3 * v + 1] = f.d[3 * i + 1]; VD[3 * v + 2] = f.d[3 * i + 2]; }
            else { let x = 0, y = 0, z = 0; x = VD[3 * lo] + w * (VD[3 * hi] - VD[3 * lo]); y = VD[3 * lo + 1] + w * (VD[3 * hi + 1] - VD[3 * lo + 1]); z = VD[3 * lo + 2] + w * (VD[3 * hi + 2] - VD[3 * lo + 2]); const l = Math.hypot(x, y, z) || 1; VD[3 * v] = x / l; VD[3 * v + 1] = y / l; VD[3 * v + 2] = z / l; }
            VH[v] = VH[lo] + w * (VH[hi] - VH[lo]);
            a.push(v);
          }
        }
      }
      arcs.set(key, a);
    }
    return va === lo ? a : a.slice().reverse();
  };

  // 路标配对: 每条特征线上, 两层的路标按沿线弧长排, 挨着的一对 (一层一个) 就是同一条边界线在相邻两层上的两个点.
  //   同一层同一位置的几个 (两条曲线在这里收头) 并成一组. 一个路标前后各可配一个 (特征线在这一层上折返时).
  //   返回 [A 顶点 → B 成员们, B 顶点 → A 成员们], 成员 { ci, j, v }
  const pairMarks = (A, B) => {
    const byF = new Map();
    const add = (row, side) => row.cs.forEach((o, ci) => {
      for (const m of o.lm) for (const e of m.fs) { let a = byF.get(e.f); if (!a) byF.set(e.f, a = []); a.push({ t: e.t, side, ci, j: m.j, v: o.base + m.j }); }
    });
    add(A, 0); add(B, 1);
    const fw = new Map(), bw = new Map();
    const link = (map, ga, gb, fid) => {
      for (const e of ga.m) {
        let a = map.get(e.v); if (!a) map.set(e.v, a = []);
        for (const w of gb.m) {
          if (!a.some(q => q.v === w.v)) a.push(w);
          const lo = Math.min(e.v, w.v), k = lo * 67108864 + Math.max(e.v, w.v);
          if (!linkF.has(k)) linkF.set(k, { f: fid, tLo: e.v === lo ? e.t : w.t, tHi: e.v === lo ? w.t : e.t });
        }
      }
    };
    for (const [fid, a] of byF) {
      a.sort((x, y) => x.t - y.t);
      const g = [];
      // 同一层上沿特征线挨得很近的几个路标 (几条曲线在这里相交 / 收头, 或两条不同的层线先后在近旁穿过同一条脊线) 并成一组:
      //   组与组配对, 组内谁配谁留给 sew 按 "邻点脚点在哪条曲线上" 定 (同一张面上的才相连)
      for (const e of a) {
        const q = g[g.length - 1];
        if (q && q.side === e.side && e.t - q.t1 < 8 * FT.tol(2 * gap(e.v) - 0.1 * D)) { q.m.push(e); q.t1 = e.t; } else g.push({ side: e.side, t: e.t, t1: e.t, m: [e] });
      }
      const tie = (x, y) => {
        if (x.side === y.side) return;
        const ga = x.side === 0 ? x : y, gb = x.side === 0 ? y : x;
        if (Math.abs(y.t - x.t1) > Math.max(gap(ga.m[0].v), gap(gb.m[0].v))) return;
        link(fw, ga, gb, fid); link(bw, gb, ga, fid);
      };
      for (let i = 0; i + 1 < g.length; i++) {
        tie(g[i], g[i + 1]);
        // 两层在这里是同一个点 (捏合点: 各层都穿过它): 沿线排在谁前谁后是随机的, 它的两个身份都要与前后的路标配
        if (g[i].side !== g[i + 1].side && d3(g[i].m[0].v, g[i + 1].m[0].v) < 2e-4 * D) {
          if (i + 2 < g.length) tie(g[i], g[i + 2]);
          if (i > 0) tie(g[i - 1], g[i + 1]);
        }
      }
    }
    return [fw, bw];
  };
  const pairTips = (A, B, fw, bw) => {
    const cands = [];
    A.cs.forEach((oa, ca) => {
      for (const k of oa.tips) {
        const u = oa.c.x[2 * k], v = oa.c.x[2 * k + 1], q = oa.base + k;
        let bj = -1, bcb = -1, bd = RUV;
        B.cs.forEach((ob, cb) => {
          for (const j of ob.tips) {
            const dd = Math.hypot(wd(ob.c.x[2 * j] - u), wd(ob.c.x[2 * j + 1] - v));
            if (dd < bd) { bd = dd; bj = j; bcb = cb; }
          }
        });
        if (bj < 0) continue;
        // 三维也得挨着: 不超过 A 的尖点到 B 那条曲线距离的 8 倍 (脊线斜穿层线时尖点间距比层距大)
        const ob = B.cs[bcb];
        let dB = Infinity;
        for (let j = 0; j < ob.n; j++) { const d = d3(q, ob.base + j); if (d < dB) dB = d; }
        if (d3(q, ob.base + bj) <= Math.min(gap(q), 8 * dB + 0.02 * D)) cands.push([bd, ca, k, bcb, bj]);
      }
    });
    cands.sort((a, b) => a[0] - b[0]);
    for (const [, ca, k, cb, j] of cands) {
      const qa = A.cs[ca].base + k, qb = B.cs[cb].base + j;
      if (fw.has(qa) || bw.has(qb)) continue;
      fw.set(qa, [{ ci: cb, j, v: qb }]); bw.set(qb, [{ ci: ca, j: k, v: qa }]);
    }
  };

  // 一层的顶点按 (u, v) 分桶 (两支一起: 同一取向的两支在 (u, v) 里重合, 由三维距离分开)
  const index = row => {
    const m = new Map();
    row.cs.forEach((o, ci) => {
      const x = o.c.x;
      for (let j = 0; j < o.n; j++) {
        const iu = Math.floor((x[2 * j] - TWO * Math.floor(x[2 * j] / TWO)) / CELL), iv = Math.floor((x[2 * j + 1] - TWO * Math.floor(x[2 * j + 1] / TWO)) / CELL);
        const k = iu * 4096 + iv;
        let a = m.get(k); if (!a) m.set(k, a = []);
        a.push(ci, j);
      }
    });
    return m;
  };
  const NC = Math.ceil(TWO / CELL), R = Math.ceil(RUV / CELL);
  // (u, v) 点 (三维点 q) 在另一层上的对应顶点 { ci, j } (RUV 以内没有则 null). 主要按 (u, v) 找, 不按三维:
  //   U 形层线 (转折缝处 ε± 两半) 在三维里两臂贴得很近, 按三维最近会串到另一臂, 三角形横跨 U 的内侧、本该盖的带反而空着;
  //   在 (u, v) 里两臂是分开的 (ζ 一正一负). 但锥映射的两支在 (u, v) 里重合 (同一取向两个落点): 每支各取 (u, v) 最近的
  //   (与最近者差 TIE 以内的再按三维挑; TIE 不能大: U 尖处两臂在 (u, v) 里也只差零点几度), 两支之间再按三维挑 ——
  //   层线整段落在折叠线 (滚角 60°) 上时两支的 (u, v) 是交叠的同一条线段, 只按 (u, v) 会一半对到另一支;
  //   光锥穿过眼睛的那层前后层线整个从一支换到另一支, 也靠两支之间按三维挑接上.
  //   由近及远一圈圈查格子, 最近者比下一圈还近 BTOL 以上就停
  const TIE = 1e-4, BTOL = 0.02;
  const foot = (row, idx, u, v, q) => {
    const iu = Math.floor((u - TWO * Math.floor(u / TWO)) / CELL), iv = Math.floor((v - TWO * Math.floor(v / TWO)) / CELL);
    const cand = [];
    let bd = RUV;
    const scan = (a, c) => {
      const e = idx.get(((iu + a + NC) % NC) * 4096 + ((iv + c + NC) % NC));
      if (!e) return;
      for (let t = 0; t < e.length; t += 2) {
        const x = row.cs[e[t]].c.x, j = e[t + 1], dd = Math.hypot(wd(x[2 * j] - u), wd(x[2 * j + 1] - v));
        if (dd <= bd + BTOL) { cand.push(e[t], j, dd); if (dd < bd) bd = dd; }
      }
    };
    for (let rr = 0; rr <= R; rr++) {
      if (rr === 0) scan(0, 0);
      else for (let a = -rr; a <= rr; a++) { scan(a, -rr); scan(a, rr); if (Math.abs(a) < rr) { scan(-rr, a); scan(rr, a); } }
      if (cand.length && bd + BTOL <= rr * CELL) break;
    }
    const bb = [Infinity, Infinity], br = ci => row.cs[ci].c.b ? 1 : 0;
    for (let t = 0; t < cand.length; t += 3) if (cand[t + 2] <= bd + BTOL && cand[t + 2] < bb[br(cand[t])]) bb[br(cand[t])] = cand[t + 2];
    let best = null, b3 = Infinity;
    for (let t = 0; t < cand.length; t += 3) {
      if (cand[t + 2] > bb[br(cand[t])] + TIE) continue;
      const d = d3(q, row.cs[cand[t]].base + cand[t + 1]);
      if (d < b3) { b3 = d; best = { ci: cand[t], j: cand[t + 1] }; }
    }
    return best;
  };
  // (u, v) 选定曲线之后, 沿这条曲线走到三维最近点 (局部下降). 取向 → 位置的映射可以剪切得很厉害: 层线落在折叠线上时
  //   u 只差 2° 就从灯旁走到 80 m 外, 相邻两层 u 又整体挪了 0.7°, "同一取向" 在两层上隔几十米. (u, v) 管选哪条线、哪一支,
  //   沿线位置交给三维
  const descend = (row, f, q) => {
    const o = row.cs[f.ci];
    let j = f.j, dj = d3(q, o.base + j);
    for (;;) {
      const a = j > 0 ? d3(q, o.base + j - 1) : Infinity, b = j < o.n - 1 ? d3(q, o.base + j + 1) : Infinity;
      if (a < dj && a <= b) { j--; dj = a; } else if (b < dj) { j++; dj = b; } else break;
    }
    return j === f.j ? f : { ci: f.ci, j };
  };
  // (u, v) 里找不到时的后备: 三维最近. 取向平面里的层线可能整段换到对称的另一组取向 (如 Lowitz 太阳过水平面,
  //   ζ ↔ 180° − ζ, 同一张三维曲面), (u, v) 对不上而三维里挨着; 那里曲面又常贴着水平的锥面, 相邻两层在三维里隔得远,
  //   距离只能按宽的 gap 限; 另要求从眼看的方向差 < 2° —— 真正相邻两层的对应点在天空里总挨着, 连错的两处在天空里隔得很远
  const COS2 = Math.cos(2 * Math.PI / 180);
  const foot3 = (row, q) => {
    let best = null, bd = gap(q);
    row.cs.forEach((o, ci) => {
      for (let j = 0; j < o.n; j++) {
        const v = o.base + j;
        if (VD[3 * q] * VD[3 * v] + VD[3 * q + 1] * VD[3 * v + 1] + VD[3 * q + 2] * VD[3 * v + 2] < COS2) continue;
        const d = d3(q, v);
        if (d < bd) { bd = d; best = { ci, j }; }
      }
    });
    return best;
  };
  // 拉链: P、Q 两串顶点号 (同向), 每步走较短的对角线; 三角形边长超出层间允许距离的不要
  const zip = (P, Q) => {
    let i = 0, j = 0;
    const add = (a, b, c) => {
      const g = Math.max(gap(a), gap(b), gap(c));
      if (d3(a, b) > g || d3(b, c) > g || d3(c, a) > g) return;
      tris.push(a, b, c);
    };
    while (i < P.length - 1 || j < Q.length - 1) {
      if (j === Q.length - 1 || (i < P.length - 1 && d3(P[i + 1], Q[j]) <= d3(P[i], Q[j + 1]))) { add(P[i], P[i + 1], Q[j]); i++; }
      else { add(P[i], Q[j + 1], Q[j]); j++; }
    }
  };
  // 把 A 层一条曲线 (或其中 [k0, k1] 一段) 对到 B 层上, 一段段缝起来; covered: B 层被盖到的顶点
  let miss = 0;                       // 这一对层里找不到对应的点数 (多了就在中间补一层, 见 bad)
  const sew = (A, oa, k0, k1, B, idxB, covered, orph = null, lim = gap, pair = null) => {
    const x = oa.c.x, ft = [];
    for (let k = k0; k <= k1; k++) {
      const f0 = foot(B, idxB, x[2 * k], x[2 * k + 1], oa.base + k), f = f0 ? descend(B, f0, oa.base + k) : foot3(B, oa.base + k);
      ft.push(f && d3(oa.base + k, B.cs[f.ci].base + f.j) < lim(oa.base + k) ? f : null);
      if (!f && orph) orph.push(k);
      if (!ft[ft.length - 1]) miss++;
    }
    // 配对的路标: 强制对应, 缝合在这里断开. 一段以它收尾时用 fL、从它起头时用 fR (B 层在这里是两条曲线收头, 或前后各配了一个时两者不同):
    //   取邻点脚点所在曲线上、离那个脚点最近的成员. 再把紧挨着它、脚点落到它另一侧去的点夹回来 (两臂重合处三维最近点会串臂;
    //   特征线斜穿层线时两层的路标沿线错开很远, 中间这些点本该扇形连到路标上)
    let cut = null, fL = null, fR = null;
    if (pair) {
      const run = new Map();            // B 的曲线 → 这一段在它上面的走向 (首尾脚点序号之差的符号)
      for (let i = 0; i < ft.length; i++) if (ft[i]) { const q = run.get(ft[i].ci); if (!q) run.set(ft[i].ci, [ft[i].j, ft[i].j]); else q[1] = ft[i].j; }
      const dirOf = ci => { const q = run.get(ci); return q && q[1] < q[0] ? -1 : 1; };
      for (let c = 0; c < ft.length; c++) {
        const M = pair.get(oa.base + k0 + c);
        if (!M) continue;
        if (!cut) { cut = new Uint8Array(ft.length); fL = new Array(ft.length); fR = new Array(ft.length); }
        // 成员里只看邻点脚点所在曲线上的, 按位置分组 (同一位置的几个序号: B 的层线在这里走了两遍 / 闭合线的首尾) ——
        //   每组取序号离邻点脚点最近的. 只有一组就是它; 有几组 (位置不同: 特征线在 A 上折返、前后各配了一个, 或两条特征线
        //   在 A 上并在一点、在 B 上分开) 时, 收尾的取沿走向最靠前的一组、起头的取最靠后的一组, 中间那段 B 扇形连到这个路标上
        const va = oa.base + k0 + c, tl2 = FT ? 2 * FT.tol(2 * gap(va) - 0.1 * D) : 1e-4 * D;
        const pick = (nb, side) => {
          let best = null, bk = Infinity;
          if (nb) {
            const reps = [];
            for (const m of M) if (m.ci === nb.ci) {
              const q = reps.find(r => d3(r.v, m.v) < tl2);
              if (!q) reps.push(m); else if (Math.abs(m.j - nb.j) < Math.abs(q.j - nb.j)) reps[reps.indexOf(q)] = m;
            }
            if (reps.length) {
              let r0 = reps[0];
              for (const r of reps) if (Math.abs(r.j - nb.j) < Math.abs(r0.j - nb.j)) r0 = r;
              best = r0;
              const sg = side * dirOf(nb.ci);
              for (const r of reps) if (sg * r.j < sg * best.j && Math.abs(ARC[r.v] - ARC[r0.v]) <= 3 * gap(va)) best = r;
            }
          }
          if (!best) { for (const m of M) { const d = d3(va, m.v); if (d < bk) { bk = d; best = m; } } }
          return best;
        };
        const nbL = c > 0 ? (cut[c - 1] ? fR[c - 1] : ft[c - 1]) : null, nbR = c + 1 < ft.length ? ft[c + 1] : null;
        cut[c] = 1;
        fL[c] = pick(nbL || nbR, 1);
        fR[c] = pick(nbR || nbL, -1);
        if (!ft[c]) miss--;
        ft[c] = fL[c];
        if (fL[c].ci === fR[c].ci && fL[c].j !== fR[c].j && d3(fL[c].v, fR[c].v) >= 1e-4 * D) {
          const ob = B.cs[fL[c].ci];
          if (Math.abs(ARC[ob.base + fR[c].j] - ARC[ob.base + fL[c].j]) <= 3 * gap(va)) {
            const Q = [], st = fR[c].j > fL[c].j ? 1 : -1;
            for (let j = fL[c].j; j !== fR[c].j + st; j += st) { Q.push(ob.base + j); if (covered) covered[ob.base + j] = 1; }
            zip([va], Q);
          }
        }
      }
      // 走向按路标近旁十来个点的脚点趋势定 (首尾之差): B 折回成两臂、A 只有一臂时, 整条上的走向一半朝前一半朝后, 不能用
      const trend = (i0, i1, ci) => {
        let a = -1, b = -1;
        for (let i = i0; i <= i1; i++) if (i >= 0 && i < ft.length && !cut[i] && ft[i] && ft[i].ci === ci) { if (a < 0) a = i; b = i; }
        return a < 0 || ft[a].j === ft[b].j ? dirOf(ci) : ft[b].j > ft[a].j ? 1 : -1;
      };
      // 夹回来的范围有上限: 沿 A 离路标不超过 3 × (路标到它对应点的距离 + 一个点距) —— 两层的路标错开多远, 中间要扇过去的就多长;
      //   路标与对应点重合 (灯上的捏合点) 时一个也不夹 (走向一旦判反, 整条曲线都会被夹到那一点上)
      if (cut) for (let c = 0; c < ft.length; c++) if (cut[c]) {
        const vc = oa.base + k0 + c, far = (i, m) => Math.abs(ARC[oa.base + k0 + i] - ARC[vc]) > 3 * (d3(vc, m.v) + d3(vc, oa.base + k0 + (i < c ? c - 1 : c + 1)));
        let m = fL[c], sg = trend(c - 12, c - 1, m.ci);
        for (let i = c - 1; i >= 0 && !cut[i] && ft[i] && ft[i].ci === m.ci && sg * (ft[i].j - m.j) > 0 && !far(i, m); i--) ft[i] = m;
        m = fR[c]; sg = trend(c + 1, c + 12, m.ci);
        for (let i = c + 1; i < ft.length && !cut[i] && ft[i] && ft[i].ci === m.ci && sg * (ft[i].j - m.j) < 0 && !far(i, m); i++) ft[i] = m;
      }
    }
    const fAt = (i, asEnd) => cut && cut[i] ? (asEnd ? fL[i] : fR[i]) : ft[i];
    // 曲线末尾一段对不上, 而旁边已对上的点正好对到 B 那条曲线的端点: 两层层线都斜着停在同一条边界上, 相邻两层的端点
    //   沿边隔得远, 中间留一个楔形 (一排锯齿). 把这一段扇形连到 B 的端点补上 (同一自由度网格在边界上用精确区间端点收边)
    const fan = (ks, ke, f) => {
      const ob = B.cs[f.ci];
      if (f.j !== 0 && f.j !== ob.n - 1) return;
      const apex = ob.base + f.j, rl = v => Math.hypot(VP[3 * v] - L[0], VP[3 * v + 1] - L[1], VP[3 * v + 2] - L[2]);
      for (let k = ks; k < ke; k++) {
        const a = oa.base + k0 + k, b = a + 1, g = 0.5 * Math.max(rl(a), rl(b), rl(apex)) + 0.05 * D;
        if (d3(a, b) <= g && d3(b, apex) <= g && d3(apex, a) <= g) tris.push(a, b, apex);
      }
    };
    const first = ft.findIndex(f => f);
    if (first > 0 && k0 === 0) fan(0, first, ft[first]);
    if (first >= 0 && k1 === oa.n - 1) {
      let last = ft.length - 1;
      while (!ft[last]) last--;
      if (last < ft.length - 1) fan(last, ft.length - 1, ft[last]);
    }
    // 中间一段对不上, 两侧对上的又都是 B 层曲线的端点、而且这两个端点在三维里是同一点: B 层在这里穿过一个捏合点
    //   (层线经过灯, 在灯点断成两段, 中间一截映到灯后看不见). 晕面在这里收成一点, 这一段扇形连到它
    const endOf = f => { const ob = B.cs[f.ci]; return f.j === 0 || f.j === ob.n - 1 ? ob.base + f.j : -1; };
    for (let s1 = 1; s1 < ft.length - 1; s1++) {
      if (ft[s1] || !ft[s1 - 1]) continue;
      let s2 = s1;
      while (s2 + 1 < ft.length && !ft[s2 + 1]) s2++;
      if (s2 + 1 >= ft.length) break;
      const e1 = endOf(ft[s1 - 1]), e2 = endOf(ft[s2 + 1]);
      if (e1 >= 0 && e2 >= 0 && d3(e1, e2) < 1e-3 * D)
        for (let k = s1 - 1; k <= s2; k++) tris.push(oa.base + k0 + k, oa.base + k0 + k + 1, e1);
      s1 = s2;
    }
    const doneL = cut ? new Uint8Array(ft.length) : null, doneR = cut ? new Uint8Array(ft.length) : null;
    let s = 0;
    while (s < ft.length) {
      const f0 = fAt(s, false);
      if (!f0) { s++; continue; }
      // 连续一段: 对到同一条曲线, 脚点单调 (同向); 相邻两脚点之间跳过的 B 弧长不超过 3 × (A 上的点距 + 层间距)
      //   (U 形层线的尖端: B 的尖端比 A 的靠外时, A 绕过尖端那几点的最近点从 B 的一臂跳到另一臂, 跳过的正是 B 的尖端 ——
      //   拉链会把 A 的尖点扇形连到它; 按三维点距判断就断开, 尖端一带一排缺口). 路标处一定断开: 这一段以它收尾, 下一段从它起头
      let e = s, dir = 0, fe = f0;
      while (e + 1 < ft.length) {
        if (cut && cut[e] && e > s) break;
        const fn = fAt(e + 1, true);
        if (!fn || fn.ci !== f0.ci) break;
        const sg = Math.sign(fn.j - fe.j);
        if (sg && dir && sg !== dir) break;
        const ob = B.cs[f0.ci].base, a = oa.base + k0 + e;
        // 跳过的 B 弧长按当地实际层距 (点到脚点的距离) 限: 按宽的 gap 限时, B 在边界处折回成两臂、A 只有一臂的地方,
        //   A 的脚点在两臂之间来回跳, 每跳一次就从一个点扇出十几米长的三角形横跨 U 形内侧. 路标的强制对应不受此限
        if (!(cut && (cut[e] || cut[e + 1])) && Math.abs(ARC[ob + fn.j] - ARC[ob + fe.j]) > 3 * (d3(a, a + 1) + Math.max(d3(a, ob + fe.j), d3(a + 1, ob + fn.j)))) break;
        if (sg) dir = sg;
        e++; fe = fn;
      }
      if (e > s) {
        const ob = B.cs[f0.ci], P = [], Q = [];
        for (let k = s; k <= e; k++) P.push(oa.base + k0 + k);
        const j0 = f0.j, j1 = fe.j, st = j1 >= j0 ? 1 : -1;
        for (let j = j0; j !== j1 + st; j += st) { Q.push(ob.base + j); if (covered) covered[ob.base + j] = 1; }
        if (cut && (cut[s] || cut[e])) {
          // 起头 / 收尾在路标上: 两个路标之间那段特征线上插的点接到哪一串 —— 它顺着谁的走向就接在对面那串上
          //   (特征线斜穿层线: A 的路标落后, 它后面一段 A 对着的是这段特征线而不是 B, 特征线就是 B 那串的前缀; 反之亦然).
          //   接错了一串会成发卡形, 拉链从一个点扇出去盖住整条曲线
          const ta0 = P[0], tb0 = Q[0], an0 = P[1], bn0 = Q.length > 1 ? Q[1] : -1;
          const ta1 = P[P.length - 1], tb1 = Q[Q.length - 1], an1 = P[P.length - 2], bn1 = Q.length > 1 ? Q[Q.length - 2] : -1;
          const cosAt = (o, a, b) => {
            const ax = VP[3 * a] - VP[3 * o], ay = VP[3 * a + 1] - VP[3 * o + 1], az = VP[3 * a + 2] - VP[3 * o + 2], bx = VP[3 * b] - VP[3 * o], by = VP[3 * b + 1] - VP[3 * o + 1], bz = VP[3 * b + 2] - VP[3 * o + 2];
            return (ax * bx + ay * by + az * bz) / (Math.hypot(ax, ay, az) * Math.hypot(bx, by, bz) || 1);
          };
          const onQ = (ta, tb, an, bn) => cosAt(ta, tb, an) > (bn < 0 ? -2 : cosAt(tb, ta, bn));
          if (cut[e]) { const a = arcOf(ta1, tb1); if (a.length) { if (onQ(ta1, tb1, an1, bn1)) for (let i = a.length - 1; i >= 0; i--) Q.push(a[i]); else P.push(...a); } }
          if (cut[s]) { const a = arcOf(ta0, tb0); if (a.length) { if (onQ(ta0, tb0, an0, bn0)) Q.unshift(...a); else P.unshift(...a.slice().reverse()); } }
        }
        if (cut) { if (cut[s]) doneR[s] = 1; if (cut[e]) doneL[e] = 1; }
        zip(P, Q);
      }
      s = cut && cut[e] && e > s ? e : e + 1;
    }
    // 路标的强制对应没连上的 (邻点的脚点越过了它 —— 脚点趋势被两臂来回跳搅乱时夹不回来, 上面的一段就在它前面断了):
    //   单独把 (邻点, 路标) 与 B 上 [邻点的脚点 .. 路标的对应点] 这一小段缝上, 不问走向
    if (cut) for (let c = 0; c < ft.length; c++) if (cut[c]) {
      for (const sd of [-1, 1]) {
        const i = c + sd, m = sd < 0 ? fL[c] : fR[c];
        if ((sd < 0 ? doneL[c] : doneR[c]) || i < 0 || i >= ft.length) continue;
        const f = cut[i] ? (sd < 0 ? fR[i] : fL[i]) : ft[i];
        if (!f || f.ci !== m.ci) continue;
        const ob = B.cs[m.ci], va = oa.base + k0 + c, vn = oa.base + k0 + i;
        if (Math.abs(ARC[ob.base + f.j] - ARC[ob.base + m.j]) > 3 * (d3(va, vn) + gap(va))) continue;
        const Q = [], st = m.j >= f.j ? 1 : -1;
        for (let j = f.j; j !== m.j + st; j += st) { Q.push(ob.base + j); if (covered) covered[ob.base + j] = 1; }
        zip([vn, va], Q);
      }
    }
  };

  // 两层对不上的 (曲线在这两层之间变化太大: 换了一组对称取向、离开灯的方向转了个角 …) 报给调用方, 在中间补算一层再缝
  //   (同一自由度网格: 相邻两行区间对不上就插一行)
  const bad = [], idxC = new Map();
  const idxOf = t => { let m = idxC.get(t); if (!m) { idxC.set(t, m = index(rows[t])); idxC.delete(t - 4); } return m; };
  for (let r = 0; r + 1 < rows.length; r++) {
    const A = rows[r], B = rows[r + 1];
    miss = 0;
    const nA = A.cs.reduce((s, o) => s + o.n, 0), nB = B.cs.reduce((s, o) => s + o.n, 0);
    {
      const idxB = idxOf(r + 1), idxA = idxOf(r), covered = new Uint8Array(nv), [tipF, tipB] = pairMarks(A, B);
      pairTips(A, B, tipF, tipB);
      for (const oa of A.cs) {
        const orph = [];
        sew(A, oa, 0, oa.n - 1, B, idxB, covered, orph, gap, tipF);
        // B 层在这里整段没有曲线 (连 (u, v) 候选都没有): 隔一两层再找. 光锥穿过眼睛的那几层 (h ≈ −Σ) 前后,
        //   沿灯–眼连线延长线的那段层线会连着一两层缺失 (层线计算在这里本来就不全), 中间留一道细缝.
        //   隔层的对应按跨过的 h 收紧三维距离 (0.05 D + 10 r Δh), 免得跨过真的空洞
        const save = miss;
        for (let a = 0; a < orph.length;) {
          let b = a;
          while (b + 1 < orph.length && orph[b + 1] === orph[b] + 1) b++;
          if (b > a) {
            const k0 = Math.max(0, orph[a] - 1), k1 = Math.min(oa.n - 1, orph[b] + 1);
            for (let t = r + 2; t <= Math.min(r + 3, rows.length - 1); t++) {
              const dh = (rows[t].h - A.h) * Math.PI / 180, o2 = [];
              sew(A, oa, k0, k1, rows[t], idxOf(t), null, o2, i => 0.05 * D + 10 * dh * Math.hypot(VP[3 * i] - L[0], VP[3 * i + 1] - L[1], VP[3 * i + 2] - L[2]));
              if (o2.length < orph[b] - orph[a] + 1) break;
            }
          }
          a = b + 1;
        }
        miss = save;
      }
      // 上层没盖到的段 (各带一个已盖到的邻点接上): 从上往下再缝
      for (const ob of B.cs) {
        let j = 0;
        while (j < ob.n) {
          if (covered[ob.base + j]) { j++; continue; }
          let e = j;
          while (e + 1 < ob.n && !covered[ob.base + e + 1]) e++;
          sew(B, ob, Math.max(0, j - 1), Math.min(ob.n - 1, e + 1), A, idxA, null, null, gap, tipB);
          j = e + 1;
        }
      }
    }
    if (miss > 8 && miss > 0.05 * Math.min(nA, nB) && B.h - A.h > 0.02) bad.push([A.h, B.h]);
  }

  // 法向: 相邻三角形面积加权; 三角形的绕向不统一 (上下两遍缝的方向相反), 按与已累加的同号再加 (NA: 每顶点一个, 只用来定正反).
  //   每个三角形的每个角再各算一个 (折痕角, 60°): 只与法向夹角 60° 以内的相邻三角形平均 —— 折叠处两张面共用一条层线上的顶点,
  //   尖点脊线 / 转折缝两侧的面也在这里相交成折痕; 每顶点只有一个法向时两张面的明暗互相带歪, 面上一块块暗斑
  const NA = new Float64Array(3 * nvT), nT = tris.length / 3;
  const TN = new Float64Array(3 * nT);
  for (let t = 0; t < nT; t++) {
    const a = tris[3 * t], b = tris[3 * t + 1], c = tris[3 * t + 2];
    const ux = VP[3 * b] - VP[3 * a], uy = VP[3 * b + 1] - VP[3 * a + 1], uz = VP[3 * b + 2] - VP[3 * a + 2];
    const vx = VP[3 * c] - VP[3 * a], vy = VP[3 * c + 1] - VP[3 * a + 1], vz = VP[3 * c + 2] - VP[3 * a + 2];
    TN[3 * t] = uy * vz - uz * vy; TN[3 * t + 1] = uz * vx - ux * vz; TN[3 * t + 2] = ux * vy - uy * vx;
    for (const v of [a, b, c]) {
      const sg = NA[3 * v] * TN[3 * t] + NA[3 * v + 1] * TN[3 * t + 1] + NA[3 * v + 2] * TN[3 * t + 2] < 0 ? -1 : 1;
      for (let k = 0; k < 3; k++) NA[3 * v + k] += sg * TN[3 * t + k];
    }
  }
  const pos = new Float32Array(9 * nT), dir = new Float32Array(9 * nT), hh = new Float32Array(3 * nT), nrm = new Float32Array(9 * nT);
  const off = new Int32Array(nvT + 1);
  for (let i = 0; i < tris.length; i++) off[tris[i] + 1]++;
  for (let v = 0; v < nvT; v++) off[v + 1] += off[v];
  const inc = new Int32Array(tris.length), fill = off.slice(0, nvT);
  for (let i = 0; i < tris.length; i++) inc[fill[tris[i]]++] = (i / 3) | 0;
  const TL = new Float64Array(nT);
  for (let t = 0; t < nT; t++) TL[t] = Math.hypot(TN[3 * t], TN[3 * t + 1], TN[3 * t + 2]);
  for (let t = 0; t < nT; t++) for (let k = 0; k < 3; k++) {
    const v = tris[3 * t + k], o = 9 * t + 3 * k;
    let x = 0, y = 0, z = 0;
    for (let i = off[v]; i < off[v + 1]; i++) {
      const j = inc[i], c = TN[3 * j] * TN[3 * t] + TN[3 * j + 1] * TN[3 * t + 1] + TN[3 * j + 2] * TN[3 * t + 2];
      if (j !== t && Math.abs(c) < 0.5 * TL[j] * TL[t]) continue;
      const sg = c < 0 ? -1 : 1;
      x += sg * TN[3 * j]; y += sg * TN[3 * j + 1]; z += sg * TN[3 * j + 2];
    }
    if (x * NA[3 * v] + y * NA[3 * v + 1] + z * NA[3 * v + 2] < 0) { x = -x; y = -y; z = -z; }
    const nl = Math.hypot(x, y, z) || 1;
    pos[o] = VP[3 * v]; pos[o + 1] = VP[3 * v + 1]; pos[o + 2] = VP[3 * v + 2];
    dir[o] = VD[3 * v]; dir[o + 1] = VD[3 * v + 1]; dir[o + 2] = VD[3 * v + 2];
    nrm[o] = x / nl; nrm[o + 1] = y / nl; nrm[o + 2] = z / nl;
    hh[3 * t + k] = VH[v];
  }
  const used = new Uint8Array(nv);
  for (const v of tris) used[v] = 1;
  let nu = 0; for (let i = 0; i < nv; i++) nu += used[i];
  return { pos, dir, h: hh, nrm, nTri: nT, cover: nv ? nu / nv : 0, bad };
}
