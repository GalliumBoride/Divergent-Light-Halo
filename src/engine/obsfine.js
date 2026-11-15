// 观测者焦散二自由度折叠的修正版串线 (obscaustic.js 保持与 C# 逐句一致, 这里只改"怎么连"):
//   原版沿 χ_内=0 曲线抽样到 ≤400 点, 并按 solveInner 返回的【序号】把内层解串成线 —
//   解的个数一变 (出现/消失), 序号就错位, 线会接到别的分支上 (长直弦、乱连).
//   这里: 抽样更密, 每个开着的段按"同一 α 支 + 内层角最接近 + 天空方向最接近"去接新点 (连续性匹配),
//   接不上的段收尾, 没人接的新点开新段. 输出与原版同格式 (Seg: p 为相对眼的位置, d 为视方向).

import * as OC from './obscaustic.js';

const wrapPi = a => a - 2 * Math.PI * Math.round(a / (2 * Math.PI));
const ang = (a, b) => { const c = a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; return Math.acos(c > 1 ? 1 : c < -1 ? -1 : c); };

//   段端: 一段接不上下一个样本 (那一支在两个样本之间消失: 两支 α 在重根处汇合, 或走出可见 / 有效区, 或超出最大距离) 时,
//   在两个样本之间沿 χ=0 曲线二分 30 次, 把最后一个接得上的点补进段里; 新开的段同样朝前一个样本二分补上起点.
//   重根附近位置随参数呈 √ 变化, 相邻两个样本可以隔好几米 —— 不补的话两支在汇合点两侧各差一截, 焦散线中间空一段.
//   段内相邻两点隔得远 (> 基线 / 60) 时在零线上取中点加密 (最多 5 层).
export function fold2DofFine(c2, traceH = 0.01, maxPts = 1800, cancelled = null) {
  const probe = c2.inner(0.0);
  const f = OC.chiFieldOf(probe, c2.omegaInner());
  const curves = OC.traceZeros(f, traceH, 6000, 40, 96, 1440, cancelled);
  const segs = [];
  const close = o => { if (o.seg.p.length >= 2) segs.push(o.seg); };
  const candAt = xi => {
    const cand = [];
    for (const b of OC.solveInner(xi, c2)) {
      const cc = c2.inner(b);
      const brs = OC.pairAlpha(xi, cc);
      for (let k = 0; k < 2; k++)
        if (brs[k].vis && brs[k].r <= cc.maxDist) cand.push({ b, k, br: brs[k], cc });
    }
    return cand;
  };
  const mid = (a, b) => {
    const nu = cr3(a, b), m = Math.hypot(nu[0], nu[1], nu[2]);
    if (!(m > 1e-14)) return null;
    return OC.projectTransverse(unit([a[0] + b[0], a[1] + b[1], a[2] + b[2]]), [nu[0] / m, nu[1] / m, nu[2] / m], f, 3e-3);
  };
  const smax = Math.hypot(c2.dvec[0], c2.dvec[1], c2.dvec[2]) / 60;
  for (const cur of curves) {
    if (cancelled && cancelled()) break;
    const step = Math.max(1, Math.trunc(cur.length / Math.max(1, maxPts)));
    const dXi = traceH * step;                       // 相邻样本在 ξ 上的间距 (rad)
    const tolB = Math.max(0.12, 6 * dXi), tolD = Math.max(0.05, 8 * dXi);
    // cand 里与状态 (k, b, d) 接得上的最好的一个 (used: 已被别的段占了的)
    //   wide: 二分 / 加密时放宽 (内层角 4 倍, 天空方向到 0.5 rad): 重根旁、眼旁天空方向一步能转几十度, 按样本间距定的
    //   连续性门限会把 "还在、只是走得快" 判成 "没了"
    const match = (cand, k, b, d, used, wide = false) => {
      let best = -1, bd = Infinity;
      const tB = wide ? 4 * tolB : tolB, tD = wide ? Math.max(0.5, tolD) : tolD;
      for (let ci = 0; ci < cand.length; ci++) {
        const c = cand[ci];
        if ((used && used[ci]) || c.k !== k) continue;
        const db = Math.abs(wrapPi(c.b - b)), da = ang(c.br.d, d);
        if (db > tB || da > tD) continue;
        const sc = db / tolB + da / tolD;
        if (sc < bd) { bd = sc; best = ci; }
      }
      return best;
    };
    // 零线上 lo (接得上状态 k, b, d) 与 hi (接不上) 之间最后一个接得上的点
    const edge = (lo, hi, k, b, d) => {
      let best = null;
      for (let it = 0; it < 30; it++) {
        const m = mid(lo, hi);
        if (!m) break;
        const cand = candAt(m), ci = match(cand, k, b, d, null, true);
        if (ci >= 0) { lo = m; best = cand[ci]; b = best.b; d = best.br.d; } else hi = m;
      }
      return best ? { xi: lo, c: best } : null;
    };
    // 段 o 的末点与 (xiB, cB) 之间隔得远时在零线上补点 (不含两端)
    const fill = (o, xiB, cB, depth) => {
      const q = o.seg.p[o.seg.p.length - 1], rB = cB.br.r, dB = cB.br.d;
      if (depth <= 0 || Math.hypot(q[0] - dB[0] * rB, q[1] - dB[1] * rB, q[2] - dB[2] * rB) <= smax) return;
      const m = mid(o.xiL, xiB);
      if (!m) return;
      const cand = candAt(m);
      let ci = match(cand, o.k, o.b, o.d, null, true);
      if (ci < 0) ci = match(cand, o.k, cB.b, dB, null, true);
      if (ci < 0) return;
      const cm = cand[ci];
      fill(o, m, cm, depth - 1); add(o, cm, m); fill(o, xiB, cB, depth - 1);
    };
    let open = [], prevXi = null;
    for (let i = 0; i < cur.length; i += step) {
      if (cancelled && cancelled()) break;
      const xi = cur[i];
      const cand = candAt(xi);
      const used = new Uint8Array(cand.length), next = [];
      for (const o of open) {
        const best = match(cand, o.k, o.b, o.d, used);
        if (best < 0) {
          // 按门限接不上: 二分看这一支到底走到哪. 一直走到了这个样本跟前、并且与这里的某个解接得上 → 没断, 只是走得快
          const e = edge(o.xiL, xi, o.k, o.b, o.d);
          if (e) { fill(o, e.xi, e.c, 14); add(o, e.c, e.xi); }
          const b2 = e && ang(e.xi, xi) < 1e-3 * dXi ? match(cand, o.k, o.b, o.d, used, true) : -1;
          if (b2 < 0) { close(o); continue; }
          used[b2] = 1; add(o, cand[b2], xi); next.push(o);
          continue;
        }
        used[best] = 1;
        fill(o, xi, cand[best], 14);
        add(o, cand[best], xi);
        next.push(o);
      }
      for (let ci = 0; ci < cand.length; ci++) {
        if (used[ci]) continue;
        const c = cand[ci], o = { k: c.k, b: c.b, d: c.br.d, seg: new OC.Seg('K1-2dof', c.k) };
        const e = prevXi ? edge(xi, prevXi, c.k, c.b, c.br.d) : null;      // 朝前一个样本二分: 这一段真正的起点
        if (e) { add(o, e.c, e.xi); fill(o, xi, c, 14); }
        add(o, c, xi);
        next.push(o);
      }
      open = next; prevXi = xi;
    }
    open.forEach(close);
  }
  return segs;
}

// ---------------------------------------------------------------- 一自由度: 段端二分到精确终点
// 原版 (curvesFromField) 在可见窗口里取 400 点, 段停在最后一个可见样本上. 段端多是两支 α 汇合处 (重根, 位置随参数呈 √ 变化)
// 或截断处, 最后一个样本离真正终点可以差一米 (晕面上两段焦散中间留一截空). 这里同样追零线、配对,
// 但每个 可见 ↔ 不可见 的相邻样本之间沿 χ=0 曲线二分 30 次, 把精确终点补进段里 → 两支在汇合点接上.
const cr3 = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const unit = a => { const m = Math.hypot(a[0], a[1], a[2]); return m > 1e-300 ? [a[0] / m, a[1] / m, a[2] / m] : a; };

function curvesFromFieldFine(c, f, kind, extraValid, traceH, coarse, windowPts, cancelled, maxSegs = 4000) {
  const segs = [];
  const visOf = (xi, bi) => {
    const br = OC.pairAlpha(xi, c)[bi];
    return br.visGeom && (!extraValid || extraValid(xi)) && br.r <= c.maxDist ? br : null;
  };
  // 零线上 a (可见) 与 b (不可见) 之间的最后一个可见点
  const edge = (a, b, bi) => {
    let lo = a, hi = b, blo = visOf(a, bi);
    for (let it = 0; it < 30; it++) {
      const nu = unit(cr3(lo, hi));
      const m = OC.projectTransverse(unit([lo[0] + hi[0], lo[1] + hi[1], lo[2] + hi[2]]), nu, f, 3e-3);
      const bm = visOf(m, bi);
      if (bm) { lo = m; blo = bm; } else hi = m;
    }
    return { xi: lo, br: blo };
  };
  const curves = OC.traceZeros(f, traceH, 6000, 40, 96, 1440, cancelled);
  for (const cur of curves) {
    if (cancelled && cancelled()) break;
    if (segs.length >= maxSegs) break;
    const cd = OC.densify(cur, f, coarse);
    const vis = cd.map(xi => { const b = OC.pairAlpha(xi, c), ex = !extraValid || extraValid(xi); return [b[0].visGeom && ex, b[1].visGeom && ex]; });
    for (let bi = 0; bi < 2; bi++) {
      let i0 = -1;
      for (let i = 0; i <= cd.length; i++) {
        const on = i < cd.length && vis[i][bi];
        if (on && i0 < 0) i0 = i;
        if (!on && i0 >= 0) {
          const a = Math.max(0, i0 - 2), b = Math.min(cd.length - 1, i - 1 + 2);
          const sub = cd.slice(a, b + 1);
          let arc = 0;
          for (let q = 0; q + 1 < sub.length; q++) arc += Math.acos(Math.max(-1, Math.min(1, sub[q][0] * sub[q + 1][0] + sub[q][1] * sub[q + 1][1] + sub[q][2] * sub[q + 1][2])));
          const fine = OC.densify(sub, f, Math.max(arc / windowPts, 1e-7));
          let sg = null, prev = null, prevBr = null;
          const push = (xi, br) => {
            if (!sg) sg = new OC.Seg(kind, bi);
            sg.xi.push(xi); sg.d.push(br.d); sg.alpha.push(br.alpha); sg.r.push(br.r); sg.t.push(br.t);
            sg.p.push([br.d[0] * br.r, br.d[1] * br.r, br.d[2] * br.r]);
            if (br.r > 0.98 * c.maxDist) sg.truncated = true;
          };
          const close = () => { if (sg && sg.p.length >= 2) segs.push(sg); sg = null; };
          for (const xi of fine) {
            const br = visOf(xi, bi);
            if (br && !prevBr && prev) { const e = edge(xi, prev, bi); push(e.xi, e.br); }
            if (!br && prevBr) { const e = edge(prev, xi, bi); push(e.xi, e.br); close(); }
            if (br) push(xi, br);
            prev = xi; prevBr = br;
          }
          close();
          i0 = -1;
        }
      }
    }
  }
  return segs;
}

export function fold1DofFine(c, traceH = 0.008, cancelled = null) {
  const f = OC.chiFieldOf(c, c.what);
  return curvesFromFieldFine(c, f, 'K1', xi => OC.chain(xi, c, false).valid, traceH, 3e-4, 400, cancelled);
}

export function seam1DofFine(c, which, traceH = 0.008, cancelled = null) {
  const f = OC.seamFieldOf(c, which);
  return curvesFromFieldFine(c, f, 'K5:' + which, xi => OC.seamValid(xi, c, which, 1e-6), traceH, 3e-4, 400, cancelled, 600);
}

function add(o, c, xi) {
  const br = c.br, sg = o.seg;
  sg.xi.push(xi); sg.d.push(br.d); sg.alpha.push(br.alpha);
  sg.r.push(br.r); sg.t.push(br.t); sg.inner.push(c.b);
  sg.p.push([br.d[0] * br.r, br.d[1] * br.r, br.d[2] * br.r]);
  if (br.r > 0.98 * c.cc.maxDist) sg.truncated = true;
  o.b = c.b; o.d = br.d; o.xiL = xi;
}
