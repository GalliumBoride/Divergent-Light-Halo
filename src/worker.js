// 计算 Worker: 收 { kind, args }, 回 { ok, data } (类型化数组走 transfer). 任务种类:
//   coneH / coneS   一自由度光锥 h 扫 / ψ 扫 的一段切片 (原窗体的采样点)
//   mesh1           一自由度晕面三角网格 ((h, α) 网格 + 边界二分裁剪)
//   bdry1           一自由度物理边界 (边界引擎 v4) + 平行光弧
//   obs1 / obs2     观测者视角发散焦散 (一 / 二自由度)
//   layers2         二自由度光锥的一段 h 层 (焦散或选定的平行光曲线 → 锥映射)
//   layers2h        二自由度光锥的指定 h 层 (只要显示曲线; 有曲线的 h 区间里加密用)
//   edges2          二自由度曲面边界线 (末端点身份 h-延拓) + 平行光弧
//   rims2           二自由度焦散面的外缘 (焦散∩有效区边, 沿 h 精确延拓); 选的不是焦散面时返回 null
//   cone1 / cone2   锥切片 (一自由度闭式 / 二自由度单层)
//   shell1 / shell2 发散度球壳 (DivShell 闭式)

import { setApex, baseNormalFromCompact, baseNormalTable, lpFaceFromCompact } from './engine/crystal.js';
import { LP, setIndex } from './engine/lp.js';
import { lightConeH, lightConeSwap, parallelArc, divergentOnCone } from './engine/cone.js';
import { lightConeBoundaries } from './engine/bdry4.js';
import { surfaceMesh1 } from './engine/mesh1.js';
import { lightCone2DofLayers, parallelArc2Dof, edges2Dof, resolveSurface, causticRims2Dof } from './engine/epdiv.js';
import * as OC from './engine/obscaustic.js';
import { fold2DofFine, fold1DofFine, seam1DofFine } from './engine/obsfine.js';
import { computeShell, computeShell2Dof } from './engine/divshell.js';

function transferOf(obj, list = []) {
  if (obj && typeof obj === 'object') {
    if (ArrayBuffer.isView(obj)) { if (!list.includes(obj.buffer)) list.push(obj.buffer); return list; }
    for (const k in obj) transferOf(obj[k], list);
  }
  return list;
}

// 公共配置: 锥晶顶角 / 折射率 / 部分反射 (窗体二自由度恒用 ModePartial = true) / 文字语言
function setup(a) {
  setApex(a.apexUp, a.apexLo);
  LP.upperApexDeg = a.apexUp; LP.lowerApexDeg = a.apexLo;
  LP.partial = true;
  setIndex(a.n);
  OC.setLang(a.lang === 'en' ? 'en' : 'zh');
}

const f32of = pts => {          // [[x,y,z]…] → Float32Array
  const a = new Float32Array(pts.length * 3);
  for (let i = 0; i < pts.length; i++) { a[3 * i] = pts[i][0]; a[3 * i + 1] = pts[i][1]; a[3 * i + 2] = pts[i][2]; }
  return a;
};

// 焦散段 → { p: 位置, d: 精确视方向 }
const segPD = sg => ({ p: sg.pts3(), d: f32of(sg.d) });
function obsPack(fold, seams, asym, note) {
  const ridge = [], ridgeInf = [], seam = [];
  for (const sg of fold || []) {
    if (sg.p.length < 2) continue;
    (sg.truncated ? ridgeInf : ridge).push(segPD(sg));
  }
  for (const sg of seams || []) if (sg.p.length >= 2) seam.push(segPD(sg));
  return { ridge, ridgeInf, seam, asym: asym || [], note: note || '' };
}

const shift = (segs, eyeH) => { for (const sg of segs) for (const p of sg.p) p[2] += eyeH; };

// χ 在整个有效域上只剩舍入噪声 (Fibonacci 球面 4000 点, 有效点 |χ| 全 < 1e-12) ⟹ 退化
function chiDegenerate(c) {
  const f = OC.chiFieldOf(c, c.what), o = { ok: false };
  let n = 0, mx = 0;
  const N = 4000;
  for (let i = 0; i < N; i++) {
    const ph = Math.acos(1 - 2 * (i + 0.5) / N), th = Math.PI * (1 + Math.sqrt(5)) * (i + 0.5);
    const v = f([Math.sin(ph) * Math.cos(th), Math.sin(ph) * Math.sin(th), Math.cos(ph)], o);
    if (!o.ok || !Number.isFinite(v)) continue;
    n++; mx = Math.max(mx, Math.abs(v));
  }
  return n > 20 && mx < 1e-12;
}

const RUN = {
  coneH: a => lightConeH(a.fp, a.el, a.dist, a.eye, a.n, a.nH, a.nA, a.maxD, a.zeta, a.theta, a.from, a.to),
  coneS: a => lightConeSwap(a.fp, a.el, a.dist, a.eye, a.n, a.nA, a.nPsi, a.maxD, a.zeta, a.theta, a.from, a.to),
  mesh1: a => surfaceMesh1(a.fp, a.el, a.dist, a.eye, a.n, a.maxD, a.zeta, a.theta, { budget: a.budget }),

  bdry1: a => {
    const b = lightConeBoundaries(a.fp, a.el, a.dist, a.eye, a.n, a.maxD, a.zeta, a.theta);
    b.parArc = parallelArc(a.fp, a.el, a.zeta, a.theta, a.n, [0, 0, a.eye]);
    return b;
  },

  obs1: a => {
    const er = a.el * Math.PI / 180;
    const dvec = [a.dist * Math.cos(er), 0.0, a.dist * Math.sin(er) - a.eye];
    const c = OC.make1Dof(a.fp, baseNormalTable(), a.n, 'alpha', 0, a.zeta * Math.PI / 180, a.theta * Math.PI / 180,
      dvec, [0, 0, a.eye], a.maxD);
    // 退化光路 (如板晶 1-4-6-2: 天球映射与扫描轴对易, χ ≡ 0): 两个引擎都只追得出舍入噪声 → 不画, 给说明
    if (chiDegenerate(c)) return obsPack([], [], null, a.lang === 'en'
      ? 'χ ≡ 0 for this path and orientation (the sky map commutes with the sweep): degenerate, no fold caustic.'
      : '这条光路在此取向下 χ ≡ 0 (天球映射与扫描转轴对易): 退化, 没有折叠焦散.');
    const fold = fold1DofFine(c, 0.008);                  // 段端二分到精确终点 (两支在汇合处接上), 见 engine/obsfine.js
    let seams = null;
    if (a.seams) { seams = []; for (const w of OC.seamNames(c)) seams.push(...seam1DofFine(c, w, 0.008)); }
    const asym = a.asym ? OC.retroCircles(c, Math.max(50.0, a.maxD * 2.5)) : null;
    shift(fold, a.eye);
    if (seams) shift(seams, a.eye);
    const note = fold.length === 0 ? OC.diagnose(c, 60000, 0.008).verdict : '';
    return obsPack(fold, seams, asym, note);
  },

  obs2: a => {
    const er = a.el * Math.PI / 180;
    const dvec = [a.dist * Math.cos(er), 0.0, a.dist * Math.sin(er) - a.eye];
    const c2 = new OC.Ctx2({
      faces: a.fp, baseN: baseNormalTable(), n: a.n,
      mode: a.lowitz ? OC.Mode2.Lowitz : OC.Mode2.Column,
      fixedRad: (a.lowitz ? a.theta : a.zeta) * Math.PI / 180,
      dvec, eyePos: [0, 0, a.eye], maxDist: a.maxD,
    });
    const fold = fold2DofFine(c2, 0.0025, 4000);      // 修正版串线 (连续性匹配 + 更密), 见 engine/obsfine.js
    shift(fold, a.eye);
    const note = fold.length === 0
      ? (a.lang === 'en'
        ? '2-DOF: χ_inner=0 and χ_outer=0 have no common solution (or all solutions are invisible) — no 2-DOF fold caustic for this path/orientation.'
        : '二自由度: χ_内=0 与 χ_外=0 未同时可解 (或解全部不可见) — 该光路/取向下没有二自由度折叠焦散。')
      : '';
    return obsPack(fold, null, null, note);
  },

  layers2: a => {
    const lp = a.fp.map(lpFaceFromCompact);
    const surf = resolveSurface(lp, a.surf);
    return lightCone2DofLayers(lp, a.el, a.dist, a.eye, a.nH, a.nG, a.maxD, a.lowitz, a.zeta, a.theta, surf, a.from, a.to, null, true);
  },

  layers2h: a => {
    const lp = a.fp.map(lpFaceFromCompact);
    const surf = resolveSurface(lp, a.surf);
    return a.hs.map(h => ({ h, curves: lightCone2DofLayers(lp, a.el, a.dist, a.eye, 1, a.nG, a.maxD, a.lowitz, a.zeta, a.theta, surf, 0, 1, h, true)[0].curves }));
  },

  rims2: a => {
    const lp = a.fp.map(lpFaceFromCompact);
    if (resolveSurface(lp, a.surf) !== 'caustic') return null;
    return causticRims2Dof(lp, a.el, a.dist, a.eye, a.maxD, a.lowitz, a.zeta, a.theta, a.nG);
  },

  edges2: a => {
    const lp = a.fp.map(lpFaceFromCompact);
    const surf = resolveSurface(lp, a.surf);
    const e = edges2Dof(lp, a.el, a.dist, a.eye, a.maxD, a.lowitz, a.zeta, a.theta, surf, null, true);
    e.parArc = parallelArc2Dof(lp, a.el, a.eye, a.nG, a.lowitz, a.zeta, a.theta, surf, true);
    e.surf = surf;
    return e;
  },

  cone1: a => divergentOnCone(a.fp, a.h, a.el, a.dist, a.eye, a.zeta, a.theta, a.n),

  cone2: a => {
    // 二自由度的"锥切片" = 光锥计算里仰角 h 的那一层 (焦散或选定曲线经锥映射)
    const lp = a.fp.map(lpFaceFromCompact);
    const surf = resolveSurface(lp, a.surf);
    const h = Math.max(-89.99, Math.min(89.99, a.h));
    return lightCone2DofLayers(lp, a.el, a.dist, a.eye, 1, a.nG, a.maxD, a.lowitz, a.zeta, a.theta, surf, 0, 1, h, true)[0];
  },

  shell1: a => {
    const nr = a.fp.map(baseNormalFromCompact);
    return shellPack1(computeShell(nr, a.zeta, a.theta, a.n, a.L, a.E, a.r, true, 0, 360));
  },

  shell2: a => {
    const nr = a.fp.map(baseNormalFromCompact);
    // 这条光路上没有选中的那条边界 (如无反弹光路选了 "反弹1·接缝") → 与晕面一致, 回退到焦散
    let sf = a.surfFence;
    const nBounce = nr.length - 2;
    if (sf && sf.length > 4 && /^\d+$/.test(sf.substring(4)) && parseInt(sf.substring(4), 10) > nBounce) sf = null;
    const wantC = a.selCaustic || (sf === null && !a.selBdr);
    const wantF = nm => nm === sf || (nm === 'entry' ? a.bdrEntry : nm.startsWith('seam') ? a.bdrSeam : a.bdrTir);
    const r = computeShell2Dof(nr, a.lowitz ? 'lowitz' : 'column', a.lowitz ? a.theta : a.zeta, a.n, a.L, a.E, a.r, true, wantC, wantF);
    return shellPack2(r, wantC);
  },
};

function shellPack1(res) {
  const curves = [];
  for (let k = 0; k < res.polylines.length; k++) curves.push({ p: f32of(res.polylines[k]), d: f32of(res.polylineDirs[k]) });
  const ends = [], endDirs = [], folds = [], foldDirs = [];
  for (const e of res.ends) if (e.p && e.dir && e.inAlphaRange) { ends.push(e.p); endDirs.push(e.dir); }
  for (const a of res.arcs) for (let k = 0; k < a.foldP.length; k++) { folds.push(a.foldP[k]); foldDirs.push(a.foldDir[k]); }
  const cusps = res.cuspP.slice();
  return {
    curves, kinds: curves.map(() => 0),
    ends: f32of(ends), endDirs: f32of(endDirs), folds: f32of(folds), foldDirs: f32of(foldDirs),
    cusps: f32of(cusps), ok: res.ok, note: res.ok ? '' : res.note,
    nArcs: curves.length, nEnds: ends.length, nFolds: folds.length,
  };
}

function shellPack2(res, wantC) {
  const curves = [], kinds = [];
  let nCau = 0, nFen = 0;
  for (const cu of res.curves) {
    curves.push({ p: f32of(cu.p), d: f32of(cu.dir) });
    kinds.push(cu.kind === 0 ? 5 : cu.kind);      // 5 = 二自由度焦散; 1..4 = 边界像 (掠入 / 出射临界 / 接缝 / 蓝边)
    if (cu.kind === 0) nCau++; else nFen++;
  }
  return {
    curves, kinds,
    ends: f32of(res.endP), endDirs: f32of(res.endDir), folds: new Float32Array(0), foldDirs: new Float32Array(0),
    cusps: f32of(res.cuspP), ok: !res.note, note: res.note || '', wantC,
    nArcs: curves.length, nCau, nFen, nEnds: res.endP.length,
  };
}

onmessage = e => {
  const { kind, args } = e.data;
  try {
    setup(args);
    const data = RUN[kind](args);
    postMessage({ ok: true, data }, transferOf(data));
  } catch (err) {
    postMessage({ ok: false, error: String((err && err.stack) || err) });
  }
};
