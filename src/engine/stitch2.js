// 二自由度曲面网格: 相邻两层层线缝成三角形 —— 一自由度网格"行 × 列"的做法搬过来.
//   行 = 每个 h 的层线 (端点已精确落在边界上, 折叠 / 尖点 / 接缝都已处理好); 列 = 沿层线的位置.
//   一自由度有共同的 α 对齐上下两行; 二自由度没有, 就用取向平面 (u, v) 里的最近点对应: 相邻两层 h 只差零点几度,
//   层曲线在 (u, v) 里只挪一点: 用 (u, v) 选对应的曲线与支 (两支在 (u, v) 里重合, 按三维分开, 见 foot), 沿线位置取三维最近
//   (见 descend). 对应好的两段用拉链式三角化缝起来
//   (每步走较短的那条对角线, Fuchs–Kedem–Uselton 的贪心版). 对不上的 (曲线在这两层之间出生 / 消亡 / 分叉) 留空;
//   下层没盖到的上层段再从上往下缝一遍; 曲线末尾斜着停在边界上的楔形用扇形补上.
// 输入 layers: [{ h, curves: [{ p, d, x, b }] }] (按 h 排好; p 位置, d 视方向, x 每点 (u, v), b 支); L 灯位置, D 灯眼距离.
// 输出与一自由度网格同格式: { pos, dir, h, nrm } (三角形顶点不共享), 外加 nTri.

const TWO = 2 * Math.PI;
const wd = a => a - TWO * Math.round(a / TWO);
const CELL = 0.04, RUV = 0.12;          // (u, v) 分桶宽度与对应的最大距离 (弧度)

// 层线为画线加得很密 (灯旁、尖点旁), 网格用不着: 沿线保留与上一个保留点三维距离 ≥ smin(r)、或从眼看张角 ≥ 0.5°、或转角 ≥ 20° 的点,
//   端点总留. smin 随离灯距离 r 放宽 (层距本身 ∝ r, 远处按固定间距会切出细长三角形)
function thin(c, smin) {
  const p = c.p, d = c.d, x = c.x, n = p.length / 3, keep = [0];
  const cosA = Math.cos(0.5 * Math.PI / 180), cosT = Math.cos(20 * Math.PI / 180);
  for (let j = 1; j < n - 1; j++) {
    const i = keep[keep.length - 1];
    const ds = Math.hypot(p[3 * j] - p[3 * i], p[3 * j + 1] - p[3 * i + 1], p[3 * j + 2] - p[3 * i + 2]);
    const ca = d[3 * j] * d[3 * i] + d[3 * j + 1] * d[3 * i + 1] + d[3 * j + 2] * d[3 * i + 2];
    const ax = p[3 * j] - p[3 * j - 3], ay = p[3 * j + 1] - p[3 * j - 2], az = p[3 * j + 2] - p[3 * j - 1];
    const bx = p[3 * j + 3] - p[3 * j], by = p[3 * j + 4] - p[3 * j + 1], bz = p[3 * j + 5] - p[3 * j + 2];
    const la = Math.hypot(ax, ay, az), lb = Math.hypot(bx, by, bz);
    const turn = la > 0 && lb > 0 && (ax * bx + ay * by + az * bz) / (la * lb) < cosT;
    if (ds >= smin(p[3 * i], p[3 * i + 1], p[3 * i + 2]) || ca < cosA || turn) keep.push(j);
  }
  if (n > 1) keep.push(n - 1);
  if (keep.length === n) return c;
  const P = new Float32Array(3 * keep.length), Dd = new Float32Array(3 * keep.length), X = new Float32Array(2 * keep.length);
  keep.forEach((j, k) => { P.set(p.subarray(3 * j, 3 * j + 3), 3 * k); Dd.set(d.subarray(3 * j, 3 * j + 3), 3 * k); X.set(x.subarray(2 * j, 2 * j + 2), 2 * k); });
  return { p: P, d: Dd, x: X, b: c.b };
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

export function stitchLayers2(layers, L, D) {
  // 全部顶点编号; 每条曲线记首顶点号
  let nv = 0;
  const smin = (x, y, z) => Math.max(D / 200, 0.015 * Math.hypot(x - L[0], y - L[1], z - L[2]));
  const rows = layers.map(l => ({ h: l.h, cs: joinCurves(l.curves.filter(c => c.x && c.p.length >= 6), 1e-4 * D)
    .map(c0 => { const c = thin(c0, smin), o = { n: c.p.length / 3, base: nv, c }; nv += o.n; return o; }) }));
  const VP = new Float64Array(3 * nv), VD = new Float32Array(3 * nv), VH = new Float32Array(nv);
  for (const r of rows) for (const o of r.cs) { VP.set(o.c.p, 3 * o.base); VD.set(o.c.d, 3 * o.base); VH.fill(r.h, o.base, o.base + o.n); }
  const d3 = (i, j) => Math.hypot(VP[3 * i] - VP[3 * j], VP[3 * i + 1] - VP[3 * j + 1], VP[3 * i + 2] - VP[3 * j + 2]);
  // 每个顶点沿所在曲线的累计弧长
  const ARC = new Float64Array(nv);
  for (const r of rows) for (const o of r.cs) for (let j = 1; j < o.n; j++) ARC[o.base + j] = ARC[o.base + j - 1] + d3(o.base + j - 1, o.base + j);
  // 层间允许的三维距离 (只挡荒唐的远距对应): 层距大致 ∝ 离灯距离, 远处一瓣锥映射沿层线方向拉得很长, 斜着的对应也合法
  //   (拉链按三维较短的对角线走, 自己把斜切补正), 所以放得很宽: 0.5 r + 0.05 D
  const gap = i => 0.5 * Math.hypot(VP[3 * i] - L[0], VP[3 * i + 1] - L[1], VP[3 * i + 2] - L[2]) + 0.05 * D;
  const tris = [];

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
  const sew = (A, oa, k0, k1, B, idxB, covered, orph = null, lim = gap) => {
    const x = oa.c.x, ft = [];
    for (let k = k0; k <= k1; k++) {
      const f0 = foot(B, idxB, x[2 * k], x[2 * k + 1], oa.base + k), f = f0 ? descend(B, f0, oa.base + k) : foot3(B, oa.base + k);
      ft.push(f && d3(oa.base + k, B.cs[f.ci].base + f.j) < lim(oa.base + k) ? f : null);
      if (!f && orph) orph.push(k);
      if (!ft[ft.length - 1]) miss++;
    }
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
    let s = 0;
    while (s < ft.length) {
      if (!ft[s]) { s++; continue; }
      // 连续一段: 对到同一条曲线, 脚点单调 (同向); 相邻两脚点之间跳过的 B 弧长不超过 3 × (A 上的点距 + 层间距)
      //   (U 形层线的尖端: B 的尖端比 A 的靠外时, A 绕过尖端那几点的最近点从 B 的一臂跳到另一臂, 跳过的正是 B 的尖端 ——
      //   拉链会把 A 的尖点扇形连到它; 按三维点距判断就断开, 尖端一带一排缺口)
      let e = s, dir = 0;
      while (e + 1 < ft.length && ft[e + 1] && ft[e + 1].ci === ft[s].ci) {
        const dj = ft[e + 1].j - ft[e].j, sg = Math.sign(dj);
        if (sg && dir && sg !== dir) break;
        const ob = B.cs[ft[s].ci].base, a = oa.base + k0 + e;
        if (Math.abs(ARC[ob + ft[e + 1].j] - ARC[ob + ft[e].j]) > 3 * (d3(a, a + 1) + gap(a))) break;
        if (sg) dir = sg;
        e++;
      }
      if (e > s) {
        const ob = B.cs[ft[s].ci], P = [], Q = [];
        for (let k = s; k <= e; k++) P.push(oa.base + k0 + k);
        const j0 = ft[s].j, j1 = ft[e].j, st = j1 >= j0 ? 1 : -1;
        for (let j = j0; j !== j1 + st; j += st) { Q.push(ob.base + j); if (covered) covered[ob.base + j] = 1; }
        zip(P, Q);
      }
      s = e + 1;
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
      const idxB = idxOf(r + 1), idxA = idxOf(r), covered = new Uint8Array(nv);
      for (const oa of A.cs) {
        const orph = [];
        sew(A, oa, 0, oa.n - 1, B, idxB, covered, orph);
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
          sew(B, ob, Math.max(0, j - 1), Math.min(ob.n - 1, e + 1), A, idxA, null);
          j = e + 1;
        }
      }
    }
    if (miss > 8 && miss > 0.05 * Math.min(nA, nB) && B.h - A.h > 0.02) bad.push([A.h, B.h]);
  }

  // 每顶点法向: 相邻三角形面积加权; 三角形的绕向不统一 (上下两遍缝的方向相反), 按与已累加的同号再加
  const NA = new Float64Array(3 * nv), nT = tris.length / 3;
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
  for (let t = 0; t < nT; t++) for (let k = 0; k < 3; k++) {
    const v = tris[3 * t + k], o = 9 * t + 3 * k;
    const nl = Math.hypot(NA[3 * v], NA[3 * v + 1], NA[3 * v + 2]) || 1;
    for (let q = 0; q < 3; q++) { pos[o + q] = VP[3 * v + q]; dir[o + q] = VD[3 * v + q]; nrm[o + q] = NA[3 * v + q] / nl; }
    hh[3 * t + k] = VH[v];
  }
  const used = new Uint8Array(nv);
  for (const v of tris) used[v] = 1;
  let nu = 0; for (let i = 0; i < nv; i++) nu += used[i];
  return { pos, dir, h: hh, nrm, nTri: nT, cover: nv ? nu / nv : 0, bad };
}
