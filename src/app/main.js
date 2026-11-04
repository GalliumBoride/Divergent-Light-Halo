// 发散光晕 3D — 界面主控: 状态 → 计算任务 (Worker 池, 按通道取消) → 显示列表 → 3D 场景 + 观测者天空.
// 计算全在 Worker 里, 与原窗体逐点一致 (见 tools/): 一自由度 = 光锥 h 扫 + ψ 扫 + 边界引擎 v4;
// 二自由度 = 光点法焦散 × 光锥映射 + 末端点 h-延拓; 另有观测者焦散 (ObsCaustic) 与发散度球壳 (DivShell).

import { Pool } from './pool.js';
import { SkyView, LENS_COUNT, fovMax } from './sky.js';
import { Scene3D } from './scene3d.js';
import { CrystalPreview } from './crystal.js';
import { t, term, termKeys, setLang, getLang, detectLang, applyI18n } from './i18n.js';
import { PRESETS } from './presets.js';
import { distTo, FILL_LO, FILL_HI } from './spacing.js';
import { parsePath, parseMultiPaths, faceNumFromCompact } from '../engine/crystal.js';
import { apexDegFromMiller, millerFromApexDeg } from '../engine/lp.js';

/* global WORKER_SRC */
const $ = id => document.getElementById(id);
const clamp = (x, a, b) => Math.max(a, Math.min(b, x));
const DEG = Math.PI / 180;

// ============================================================ 配色 (白底, matplotlib 风格)
const C = {
  pts: [70, 130, 180], entry: [214, 39, 40], mirror: [148, 103, 189], tir: [31, 119, 180],
  seam: [44, 160, 44], turn: [140, 86, 75], cusp: [227, 0, 227], caustic: [255, 127, 14],
  cSeam: [107, 142, 35], asym: [110, 110, 110], parArc: [184, 134, 11], ref: [120, 120, 120],
  cone: [214, 39, 40], shell: [0, 128, 96], shellC: [230, 120, 0], blueEdge: [23, 160, 190],
  end: [20, 20, 20], fold: [255, 127, 14], cuspM: [227, 0, 227], rim: [70, 70, 70],
};
const PATHCOL = [[214, 39, 40], [31, 119, 180], [44, 160, 44], [148, 103, 189], [255, 127, 14], [140, 86, 75], [23, 190, 207], [188, 189, 34]];
const VIRIDIS = [[68, 1, 84], [72, 40, 120], [62, 74, 137], [49, 104, 142], [38, 130, 142], [31, 158, 137], [53, 183, 121], [109, 205, 89], [180, 222, 44], [253, 231, 37]];
function cmap(u) {
  u = clamp(u, 0, 1) * (VIRIDIS.length - 1);
  const i = Math.min(VIRIDIS.length - 2, Math.floor(u)), f = u - i, a = VIRIDIS[i], b = VIRIDIS[i + 1];
  return [a[0] + (b[0] - a[0]) * f, a[1] + (b[1] - a[1]) * f, a[2] + (b[2] - a[2]) * f];
}

// ============================================================ 状态
const DEFAULTS = {
  elev: 10, dist: 15, eye: 0,
  paths: '3-5', mode: 'parry', zeta: -90, theta: 0, colZeta: -90, lowTheta: 0, surf2: 'caustic',
  n: 1.31, apexUp: 56, apexLo: 56, q: 0, maxDist: 80, pyr: false,
  lay: { surf: true, pts: false, colorBy: 'h', bEntry: true, bMirror: true, bTir: true, bSeam: true, bTurn: true, bCusp: true, bRim: true,
         caustic: false, cSeam: false, cAsym: true, parArc: true, ref22: false, foci: true, ground: true, mirror: true },
  cone: { on: false, h: 22, alpha: 30 },
  shell: { on: false, kappa: -0.301, r: 7.5, master: 'kappa', sphere: true, alpha: 12 },
  view: 'split', single: 'sky', sky: { az: 0, el: 10, fov: 180, lens: 2 }, lang: null,
  panel: true, rt: 'layers', crystalOpen: false,
};
const STORE = 'divergent-light-halo.v2';
function loadState() {
  let s = null;
  try { s = JSON.parse(localStorage.getItem(STORE) || 'null'); } catch (e) { s = null; }
  const d = structuredClone(DEFAULTS);
  if (s && typeof s === 'object') {
    for (const k of Object.keys(d)) {
      if (s[k] === undefined) continue;
      if (d[k] && typeof d[k] === 'object' && !Array.isArray(d[k])) Object.assign(d[k], s[k]);
      else d[k] = s[k];
    }
  }
  return d;
}
const S = loadState();
let saveTimer = 0;
function save() {
  clearTimeout(saveTimer);
  saveTimer = setTimeout(() => { try { localStorage.setItem(STORE, JSON.stringify(S)); } catch (e) { /* 存不了就算了 */ } }, 300);
}

const MODES = ['plate', 'parry', 'custom', 'column', 'lowitz'];
const is2 = () => S.mode === 'column' || S.mode === 'lowitz';
function orient() {
  switch (S.mode) {
    case 'plate': return { dof: 1, zeta: 0, theta: 0 };
    case 'parry': return { dof: 1, zeta: -90, theta: 0 };
    case 'custom': return { dof: 1, zeta: S.zeta, theta: S.theta };
    case 'column': return { dof: 2, lowitz: false, zeta: S.colZeta, theta: 0 };
    default: return { dof: 2, lowitz: true, zeta: 0, theta: S.lowTheta };
  }
}
const scene = () => {
  const e = S.elev * DEG;
  return { L: [S.dist * Math.cos(e), 0, S.dist * Math.sin(e)], E: [0, 0, S.eye] };
};
const baseline = () => { const { L, E } = scene(); return Math.hypot(L[0] - E[0], L[1] - E[1], L[2] - E[2]); };
const codes = () => parseMultiPaths(S.paths);
const SURF_IDS = ['caustic', 'entry', 'exit', 'deg1', 'blue1', 'deg2', 'blue2'];

// ============================================================ 计算
const NW = clamp((navigator.hardwareConcurrency || 4) - 1, 2, 12);
const pool = new Pool(WORKER_SRC, NW);
const common = () => ({ n: S.n, apexUp: S.apexUp, apexLo: S.apexLo, lang: getLang() });
function ranges(n, k) {
  const out = [], step = Math.max(1, Math.ceil(n / k));
  for (let a = 0; a < n; a += step) out.push([a, Math.min(n, a + step)]);
  return out;
}

let R = null;            // 主结果 { key, paths: [ … ], t0, ms, done, total }
let RO = null;           // 焦散 { key, paths: [ … ] }
let RC = null;           // 锥切片 { key, data, dof }
let RS = null;           // 球壳 { key, paths: [ … ], ms }
const busy = { main: false, dense: false, fill: false, obs: false, cone: false, shell: false };
const errs = [];

function mainKey() {
  const o = orient();
  return JSON.stringify([S.elev, S.dist, S.eye, codes(), o, o.dof === 2 ? S.surf2 : '', S.n, S.apexUp, S.apexLo, S.q, S.maxDist]);
}
function obsKey() { return mainKey() + '|' + S.lay.cSeam + S.lay.cAsym; }
function coneKey() { return mainKey() + '|' + S.cone.h; }
function shellKey() {
  const r = shellRadius();
  const sel = is2() ? [S.lay.caustic, S.lay.bEntry, S.lay.bTir, S.lay.bSeam] : [];
  return mainKey() + '|' + r + '|' + sel.join(',');
}

function runMain() {
  const key = mainKey();
  if (R && R.key === key) return;
  const cs = codes(), o = orient(), q = S.q;
  const tasks = [], route = [];
  R = { key, paths: [], t0: performance.now(), ms: 0, done: 0, total: 0, dof: o.dof, dDone: 0, dTotal: 0 };
  pool.cancel('dense'); pool.cancel('fill'); busy.dense = busy.fill = false;
  errs.length = 0;
  const K = Math.max(2, NW * 2);
  cs.forEach((code, ip) => {
    const fp = parsePath(code);
    const base = { ...common(), fp, el: S.elev, dist: S.dist, eye: S.eye, maxD: S.maxDist, zeta: o.zeta, theta: o.theta, lowitz: !!o.lowitz };
    const P = { code, fp, dof: o.dof, h: [], s: [], layers: [], dense: [], nH: 0, mesh: null, bdry: null, edges: null, rims2: undefined, cloud: null, t0: performance.now(), ms: 0, left: 0 };
    R.paths.push(P);
    if (o.dof === 1) {
      const nH = [320, 640, 1280][q], nA = [360, 720, 1440][q], nAS = [360, 720, 1440][q], nPs = [720, 1440, 2880][q];
      tasks.push({ kind: 'mesh1', args: { ...base, budget: [60000, 160000, 400000][q] } }); route.push([ip, 'mesh']);
      tasks.push({ kind: 'bdry1', args: base }); route.push([ip, 'bdry']);
      for (const [a, b] of ranges(nAS, K)) { tasks.push({ kind: 'coneS', args: { ...base, nA: nAS, nPsi: nPs, from: a, to: b } }); route.push([ip, 's', a]); }
      for (const [a, b] of ranges(nH, K)) { tasks.push({ kind: 'coneH', args: { ...base, nH, nA, from: a, to: b } }); route.push([ip, 'h', a]); }
    } else {
      const nH = [72, 120, 180][q], nG = [101, 141, 181][q];
      P.nH = nH;
      tasks.push({ kind: 'edges2', args: { ...base, surf: S.surf2, nG } }); route.push([ip, 'edges']);
      tasks.push({ kind: 'rims2', args: { ...base, surf: S.surf2, nG } }); route.push([ip, 'rims']);
      for (const [a, b] of ranges(nH, K)) { tasks.push({ kind: 'layers2', args: { ...base, surf: S.surf2, nH, nG, from: a, to: b } }); route.push([ip, 'layers', a]); }
    }
  });
  R.total = tasks.length;
  for (const r of route) R.paths[r[0]].left++;
  busy.main = tasks.length > 0;
  const myR = R;
  setStatus();
  pool.run('main', tasks, (i, d) => {
    if (R !== myR) return;
    const [ip, part, at] = route[i], P = R.paths[ip];
    if (part === 'bdry') P.bdry = d;
    else if (part === 'edges') P.edges = d;
    else if (part === 'rims') P.rims2 = d || null;
    else if (part === 'mesh') { P.mesh = d; P.meshCol = null; }
    else if (part === 'h') d.forEach((a, k) => { P.h[at + k] = a; });
    else if (part === 's') d.forEach((a, k) => { P.s[at + k] = a; });
    else if (part === 'layers') { d.forEach((a, k) => { P.layers[at + k] = a; }); P.lver = (P.lver || 0) + 1; }
    if (part !== 'bdry' && part !== 'edges' && part !== 'mesh' && part !== 'rims') P.cloud = null;
    if (--P.left === 0) P.ms = performance.now() - P.t0;
    R.done++;
    scheduleRefresh();
    setStatus();
  }).then(res => {
    if (res === null || R !== myR) return;
    R.ms = performance.now() - R.t0;
    busy.main = false;
    refresh(); setStatus();
    frame3d();
    runDense(myR);
  }).catch(err => {
    if (R !== myR) return;
    busy.main = false; errs.push(String(err.message || err)); console.error(err);
    refresh(); setStatus();
  });
}

// 两层曲线 (各若干条折线) 在三维里相隔多远: A 上取样点到 B 上最近点的距离, 取 q 分位 (默认 75%: 层在灯旁汇聚、在远处散开,
//   取偏高的分位才反映稀疏的那一侧). how = 'mean': 两个方向取平均 (分配补层用); 'from': 只量 A → B (A 是否整条都贴着 B)
function layerGap(A, B, how = 'mean', q = 0.75) {
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

// 二自由度: 原窗体的 h 层在 −90°..90° 上等分, 真正有曲线的往往只占一小段, 而且等 h 间隔的层在三维里疏密悬殊
//   (h → ±90° 处光锥退化成竖线, 相邻层几乎重合; 近水平处层在远处散得很开). 第一轮算完后按【三维间距等分布】补层:
//   每个原间隔的间距 g_k 由相邻两层的曲线量出 (layerGap), 补的层数 ∝ g_k (de Boor 等分布), 稀的地方多补, 挤的地方不补.
//   曲线出现 / 消失的间隔 (只有一侧有曲线) 按邻近间隔的间距补, 把收口处夹紧.
function runDense(myR) {
  if (myR.dof !== 2) return;
  const o = orient(), tasks = [], route = [];
  const K = Math.max(2, NW * 2);
  myR.paths.forEach((P, ip) => {
    let i0 = -1, i1 = -1;
    P.layers.forEach((l, ih) => { if (l && l.curves.length) { if (i0 < 0) i0 = ih; i1 = ih; } });
    if (i0 < 0) return;
    const step = 180 / Math.max(1, P.nH - 1);
    const kA = Math.max(0, i0 - 1), kB = Math.min(P.nH - 1, i1 + 1);
    const has = k => !!(P.layers[k] && P.layers[k].curves.length);
    const g = [];
    for (let k = kA; k < kB; k++) g.push(has(k) && has(k + 1) ? layerGap(P.layers[k], P.layers[k + 1]) : -1);
    const known = g.filter(x => x > 0), gMed = known.length ? known.sort((x, y) => x - y)[Math.floor(known.length / 2)] : 1;
    for (let k = 0; k < g.length; k++) if (g[k] < 0) g[k] = gMed;            // 出现 / 消失的间隔
    const total = [150, 240, 360][S.q], gSum = g.reduce((x, y) => x + y, 0) || 1;
    const hs = [];
    g.forEach((gk, k) => {
      const n = Math.min(60, Math.round(total * gk / gSum));
      const h0 = layerH(P, kA + k);
      for (let j = 1; j <= n; j++) hs.push(clamp(h0 + step * j / (n + 1), -89.99, 89.99));
    });
    const base = { ...common(), fp: P.fp, el: S.elev, dist: S.dist, eye: S.eye, maxD: S.maxDist, zeta: o.zeta, theta: o.theta,
      lowitz: !!o.lowitz, surf: S.surf2, nG: [101, 141, 181][S.q] };
    const n = Math.max(1, Math.ceil(hs.length / K));
    for (let a = 0; a < hs.length; a += n) { tasks.push({ kind: 'layers2h', args: { ...base, hs: hs.slice(a, a + n) } }); route.push(ip); }
  });
  if (!tasks.length) return;
  myR.dDone = 0; myR.dTotal = tasks.length;
  busy.dense = true; setStatus();
  pool.run('dense', tasks, (i, d) => {
    if (R !== myR) return;
    const P = myR.paths[route[i]];
    P.dense.push(...d); P.lver = (P.lver || 0) + 1;
    myR.dDone++;
    scheduleRefresh(); setStatus();
  }).then(res => {
    if (res === null || R !== myR) return;
    busy.dense = false; myR.ms = performance.now() - myR.t0;
    refresh(); setStatus();
    runFill(myR);
  }).catch(err => {
    if (R !== myR) return;
    busy.dense = false; errs.push(String(err.message || err)); console.error(err);
    refresh(); setStatus();
  });
}

// 整条画的层 (按 h 排好的 [{ h, curves }]). 全部算完后去掉多余的层: 整条 (95% 的点) 都贴着上一条留下的层
//   (不到中位间距的 0.3 倍). 只量"这一层 → 上一层"一个方向: 上一层可能很短 (只剩灯旁一小段, 层都在那里汇聚), 反方向量出来很近,
//   却不说明这一层多余 —— 它在远处也许铺开一大片. 以前取两个方向的小者, 就把这样的层连同后面一串都丢了.
function baseLayers(P) {
  const ver = (P.lver || 0) + (busy.dense ? 0.5 : 0);
  if (P.base && P.base.ver === ver) return P.base.list;
  const all = [];
  P.layers.forEach((l, ih) => { if (l) all.push({ h: layerH(P, ih), curves: l.curves }); });
  all.push(...P.dense);
  all.sort((a, b) => a.h - b.h);
  let show = all;
  if (!busy.dense && all.length > 8) {
    const nz = all.filter(l => l.curves.length), gaps = [];
    for (let k = 0; k + 1 < nz.length; k++) gaps.push(layerGap(nz[k], nz[k + 1]));
    const med = gaps.slice().sort((x, y) => x - y)[Math.floor(gaps.length / 2)] || 0;
    show = [];
    let last = null;
    nz.forEach((l, k) => {
      if (last && k + 1 < nz.length && layerGap(l, last, 'from', 0.95) < 0.3 * med) return;
      show.push(l); last = l;
    });
  }
  P.base = { ver, list: show };
  return show;
}

// ------------------------------------------------------------ 补线 (只画在 3D 里)
// 层线 = 仰角 h 的光锥与晕面的交线, 每条都收向灯, 相邻两层的三维间距大致 ∝ 离灯距离 (灯旁几厘米, 远处几十厘米),
//   整条的层线怎么排都压不平. 做法同地形图的补充等高线 / 等间距流线 (Jobard–Lefer): 整条层线照旧画; 相邻两层 A, B 张开处
//   在中间仰角补算一层 (那一段上它就是该 h 的真层线), 只画张开的那一段, 而且逐点淡入: 不透明度 = smoothstep(g; 1.4d*, 2.4d*),
//   g = 到 A 与到 B 的距离之和. 线的"墨量" ≈ Σ不透明度 / 间距 在过渡带里保持在 0.71–0.83 / d* (硬切则在 2d* 处跳一倍,
//   补线起止处一道道白条). 新层再与两侧比, 逐级对分 (≤ FILL_ROUNDS 级). d* = 整条层线逐点间距中位数的 0.75 倍
//   (墨量对上原来的中位疏密). 天空里的疏密按张角算, 是另一回事, 补线不画进天空.
const FILL_ROUNDS = 3;
// d*: 相邻整层逐点间距 (A 上的点到 B 的距离) 中位数的 0.75 倍
function fillScale(base) {
  const nz = base.filter(l => l.curves.length);
  if (nz.length < 3) return 0;
  const g = [];
  for (let k = 0; k + 1 < nz.length; k++) g.push(layerGap(nz[k], nz[k + 1]));
  g.sort((x, y) => x - y);
  const cap = 3 * (g[Math.floor(g.length / 2)] || 0), ds = [];
  if (!(cap > 0)) return 0;
  for (let k = 0; k + 1 < nz.length; k++) {
    const f = distTo(nz[k + 1].curves, cap);
    let n = 0; for (const c of nz[k].curves) n += (c.p || c).length / 3;
    const st = Math.max(1, Math.floor(n / 100));
    let i = 0;
    for (const c of nz[k].curves) { const p = c.p || c; for (let j = 0; j < p.length; j += 3, i++) if (i % st === 0) ds.push(f(p[j], p[j + 1], p[j + 2])); }
  }
  ds.sort((x, y) => x - y);
  return 0.75 * ds[Math.floor(ds.length / 2)];
}
// 折线们 pieces 上离层 X 超过 1.4d* (补线开始淡入) 的那部分有多长 (每隔一点量, 只用来排先后)
function gapLength(pieces, X, d) {
  const f = distTo(X.curves, FILL_HI * d);
  let L = 0;
  for (const c of pieces) {
    const p = c.p || c;
    for (let i = 0; i + 8 < p.length; i += 6)
      if (f(p[i], p[i + 1], p[i + 2]) > FILL_LO * d) L += Math.hypot(p[i + 6] - p[i], p[i + 7] - p[i + 1], p[i + 8] - p[i + 2]);
  }
  return L;
}
function runFill(myR, round = 1) {
  if (myR.dof !== 2 || R !== myR) return;
  const o = orient(), tasks = [], route = [];
  const K = Math.max(2, NW * 2);
  myR.paths.forEach((P, ip) => {
    if (round === 1) { P.fill = []; P.fillD = fillScale(baseLayers(P)); P.fver = (P.fver || 0) + 1; }
    const d = P.fillD;
    if (!(d > 0)) return;
    // 这一轮要比的相邻对: 第一轮 = 相邻整层 (两个方向都量); 之后 = 上一轮补的层与它两侧 (只量它画出来的那部分)
    const pairs = [];
    if (round === 1) {
      const nz = baseLayers(P).filter(l => l.curves.length);
      for (let k = 0; k + 1 < nz.length; k++) {
        const A = nz[k], B = nz[k + 1];
        pairs.push({ A, B, L: gapLength(A.curves, B, d) + gapLength(B.curves, A, d) });
      }
    } else for (const F of P.fill) if (F.lv === round - 1 && F.kept.length) {
      pairs.push({ A: F.pa, B: F, L: gapLength(F.kept, F.pa, d) });
      pairs.push({ A: F, B: F.pb, L: gapLength(F.kept, F.pb, d) });
    }
    // 第一轮不设上限 (随机漏掉几个间隔就是一道道白条); 往下几轮总数封顶 = 整层间隔数的 1.5 倍, 按张开的长度排先后
    if (round === 1) P.fillCap = Math.round(1.5 * pairs.length);
    const cand = pairs.filter(q => q.L > 3 * d && Math.abs(q.B.h - q.A.h) > 1e-4)
      .sort((a, b) => b.L - a.L).slice(0, Math.max(round === 1 ? pairs.length : 0, P.fillCap - P.fill.length));
    if (!cand.length) return;
    const base = { ...common(), fp: P.fp, el: S.elev, dist: S.dist, eye: S.eye, maxD: S.maxDist, zeta: o.zeta, theta: o.theta,
      lowitz: !!o.lowitz, surf: S.surf2, nG: [61, 85, 109][S.q], d };   // 种子网格粗一些: 补线只画在层张开的大片区域, 比一格还小的闭环漏了无妨
    // 修剪 (逐点到两侧的距离) 也在 Worker 里做, 主线程只收结果: 两侧层的折线随任务送过去
    const pos = X => X.curves.map(c => c.p || c);
    const n = Math.max(1, Math.ceil(cand.length / K));
    for (let a = 0; a < cand.length; a += n) {
      const part = cand.slice(a, a + n);
      tasks.push({ kind: 'fill2', args: { ...base, items: part.map(q => ({ h: 0.5 * (q.A.h + q.B.h), pa: pos(q.A), pb: pos(q.B) })) } });
      route.push([ip, part]);
    }
  });
  if (!tasks.length) { if (busy.fill) myR.ms = performance.now() - myR.t0; busy.fill = false; refresh(); setStatus(); return; }
  if (round === 1) { myR.dDone = 0; myR.dTotal = 0; }
  myR.dTotal += tasks.length;
  busy.fill = true; setStatus();
  pool.run('fill', tasks, (i, data) => {
    if (R !== myR) return;
    const [ip, part] = route[i], P = myR.paths[ip];
    data.forEach((l, k) => { if (l.curves.length) P.fill.push({ h: l.h, curves: l.curves, kept: l.kept, lv: round, pa: part[k].A, pb: part[k].B }); });
    P.fver++;
    myR.dDone++;
    scheduleRefresh(); setStatus();
  }).then(res => {
    if (res === null || R !== myR) return;
    if (round < FILL_ROUNDS) runFill(myR, round + 1);
    else { busy.fill = false; myR.ms = performance.now() - myR.t0; refresh(); setStatus(); }
  }).catch(err => {
    if (R !== myR) return;
    busy.fill = false; errs.push(String(err.message || err)); console.error(err);
    refresh(); setStatus();
  });
}

function runObs() {
  if (!S.lay.caustic) { pool.cancel('obs'); busy.obs = false; return; }
  const key = obsKey();
  if (RO && RO.key === key) return;
  const cs = codes(), o = orient();
  RO = { key, paths: cs.map(() => null), t0: performance.now() };
  const myR = RO;
  const tasks = cs.map(code => ({
    kind: o.dof === 1 ? 'obs1' : 'obs2',
    args: { ...common(), fp: parsePath(code), el: S.elev, dist: S.dist, eye: S.eye, maxD: S.maxDist,
            zeta: o.zeta, theta: o.theta, lowitz: !!o.lowitz, seams: S.lay.cSeam, asym: S.lay.cAsym },
  }));
  busy.obs = tasks.length > 0; setStatus();
  pool.run('obs', tasks, (i, d) => { if (RO === myR) { RO.paths[i] = d; scheduleRefresh(); } })
    .then(r => { if (r !== null && RO === myR) { busy.obs = false; RO.ms = performance.now() - RO.t0; refresh(); setStatus(); } })
    .catch(err => { if (RO === myR) { busy.obs = false; errs.push(String(err.message || err)); refresh(); setStatus(); } });
}

function runCone() {
  if (!S.cone.on) { pool.cancel('cone'); busy.cone = false; RC = null; return; }
  const key = coneKey();
  if (RC && RC.key === key) return;
  const cs = codes(), o = orient();
  if (!cs.length) { RC = null; return; }
  RC = { key, data: null, dof: o.dof };
  const myR = RC;
  const args = { ...common(), fp: parsePath(cs[0]), h: S.cone.h, el: S.elev, dist: S.dist, eye: S.eye, maxD: S.maxDist,
                 zeta: o.zeta, theta: o.theta, lowitz: !!o.lowitz, surf: S.surf2, nG: [101, 141, 181][S.q] };
  busy.cone = true; setStatus();
  pool.run('cone', [{ kind: o.dof === 1 ? 'cone1' : 'cone2', args }]).then(r => {
    if (r === null || RC !== myR) return;
    RC.data = r[0]; busy.cone = false; refresh(); setStatus();
  }).catch(err => { if (RC === myR) { busy.cone = false; errs.push(String(err.message || err)); setStatus(); } });
}

function shellRadius() {
  const D = baseline();
  return S.shell.master === 'r' ? S.shell.r : D * Math.pow(10, S.shell.kappa);
}
function runShell() {
  if (!S.shell.on) { pool.cancel('shell'); busy.shell = false; RS = null; return; }
  const key = shellKey();
  if (RS && RS.key === key) return;
  const cs = codes(), o = orient(), { L, E } = scene(), r = shellRadius();
  RS = { key, paths: cs.map(() => null), t0: performance.now(), r };
  const myR = RS;
  const sf = S.surf2 === 'caustic' ? null : S.surf2.replace('deg', 'seam');
  const tasks = cs.map(code => ({
    kind: o.dof === 1 ? 'shell1' : 'shell2',
    args: { ...common(), fp: parsePath(code), zeta: o.zeta, theta: o.theta, lowitz: !!o.lowitz, L, E, r,
            surfFence: sf, selCaustic: S.lay.caustic, selBdr: S.lay.bEntry || S.lay.bTir || S.lay.bSeam,
            bdrEntry: S.lay.bEntry, bdrTir: S.lay.bTir, bdrSeam: S.lay.bSeam },
  }));
  busy.shell = tasks.length > 0; setStatus();
  pool.run('shell', tasks, (i, d) => { if (RS === myR) RS.paths[i] = d; })
    .then(res => { if (res !== null && RS === myR) { busy.shell = false; RS.ms = performance.now() - RS.t0; refresh(); setStatus(); } })
    .catch(err => { if (RS === myR) { busy.shell = false; errs.push(String(err.message || err)); refresh(); setStatus(); } });
}

let schedTimer = 0;
function schedule(immediate) {
  clearTimeout(schedTimer);
  const go = () => {
    if (!codes().length) { R = null; RO = null; RC = null; RS = null; for (const ch of ['main', 'dense', 'fill', 'obs', 'cone', 'shell']) pool.cancel(ch); busy.main = busy.dense = busy.fill = busy.obs = busy.cone = busy.shell = false; refresh(); setStatus(); return; }
    runMain(); runObs(); runCone(); runShell(); refresh(); setStatus();
  };
  if (immediate) go(); else schedTimer = setTimeout(go, 140);
}

// ============================================================ 显示列表
function pathCloud(P) {
  if (P.cloud) return P.cloud;
  const parts = [];
  for (const a of P.h) if (a) parts.push(a);
  for (const a of P.s) if (a) parts.push(a);
  for (const l of P.layers) if (l && l.verts && l.verts.length) parts.push(l.verts);
  let n = 0;
  for (const a of parts) n += a.length;
  const out = new Float32Array(n);
  let k = 0;
  for (const a of parts) { out.set(a, k); k += a.length; }
  P.cloud = { pts: out, cols: null, colKey: '' };
  return P.cloud;
}

// 采样点的视方向 (天空 GPU 层用). 采样点只有位置, 只能由 P−E 求; 与眼重合的点给零向量 (着色器剔除)
function cloudDirs(cl) {
  if (cl.dir) return cl.dir;
  const { E } = scene(), P = cl.pts, n = P.length / 3, D = new Float32Array(3 * n);
  for (let i = 0; i < n; i++) {
    const dx = P[3 * i] - E[0], dy = P[3 * i + 1] - E[1], dz = P[3 * i + 2] - E[2], r = Math.hypot(dx, dy, dz);
    if (r > 1e-9) { D[3 * i] = dx / r; D[3 * i + 1] = dy / r; D[3 * i + 2] = dz / r; }
  }
  cl.dir = D;
  return D;
}

// 着色: h = 灯处入射仰角 (光锥层), ρ = r/D 发散度 (lg)
function colorize(cl, mode, ip, rng) {
  const key = mode + '|' + ip + '|' + (rng ? rng.join(',') : '');
  if (cl.colKey === key) return cl.cols;
  const P = cl.pts, n = P.length / 3;
  if (mode === 'mono' || mode === 'path') { cl.cols = null; cl.colKey = key; return null; }
  const cols = new Uint8Array(3 * n);
  const { L, E } = scene(), D = baseline();
  for (let i = 0; i < n; i++) {
    const v = mode === 'h' ? hOf(P, i, L) : rhoOf(P, i, E, D);
    const c = cmap((v - rng[0]) / Math.max(1e-9, rng[1] - rng[0]));
    cols[3 * i] = c[0]; cols[3 * i + 1] = c[1]; cols[3 * i + 2] = c[2];
  }
  cl.cols = cols; cl.colKey = key;
  return cols;
}
function hOf(P, i, L) {
  const dx = P[3 * i] - L[0], dy = P[3 * i + 1] - L[1], dz = P[3 * i + 2] - L[2];
  return -Math.asin(clamp(dz / Math.max(1e-12, Math.hypot(dx, dy, dz)), -1, 1)) / DEG;
}
function rhoOf(P, i, E, D) {
  const r = Math.hypot(P[3 * i] - E[0], P[3 * i + 1] - E[1], P[3 * i + 2] - E[2]);
  return Math.log10(Math.max(1e-6, r / D));
}
// 一条光路的晕面位置 (网格顶点优先, 没有网格就用采样点) — 取景 / 着色范围 / 锥面大小用
function surfPts(P) {
  if (P.mesh && P.mesh.pos.length) return P.mesh.pos;
  return pathCloud(P).pts;
}
function colorRange(mode) {
  if (!R || (mode !== 'h' && mode !== 'rho')) return null;
  const { L, E } = scene(), D = baseline();
  let lo = Infinity, hi = -Infinity;
  const take = v => { if (v < lo) lo = v; if (v > hi) hi = v; };
  for (const P of R.paths) {
    if (mode === 'h' && P.mesh && P.mesh.h.length) {
      const a = P.mesh.h, st = Math.max(1, Math.floor(a.length / 40000));
      for (let i = 0; i < a.length; i += st) take(a[i]);
    } else {
      const pts = surfPts(P), n = pts.length / 3, st = Math.max(1, Math.floor(n / 40000));
      for (let i = 0; i < n; i += st) take(mode === 'h' ? hOf(pts, i, L) : rhoOf(pts, i, E, D));
    }
    if (mode === 'h') {
      P.layers.forEach((l, ih) => { if (l && l.curves.length) take(layerH(P, ih)); });
      for (const l of P.dense) if (l.curves.length) take(l.h);
    }
  }
  if (!(hi > lo)) { lo = (lo === Infinity ? 0 : lo) - 1; hi = lo + 2; }
  return [lo, hi];
}
const layerH = (P, ih) => -90 + 180 * ih / Math.max(1, P.nH - 1);

// 网格顶点颜色 (按 h 用网格自带的精确 h; 按 ρ 由位置算)
function meshColors(P, mode, rng) {
  if (mode === 'mono' || mode === 'path' || !rng) return null;
  const key = mode + '|' + rng.join(',');
  if (P.meshCol && P.meshCol.key === key) return P.meshCol.cols;
  const pos = P.mesh.pos, n = pos.length / 3, cols = new Uint8Array(3 * n);
  const { E } = scene(), D = baseline(), span = Math.max(1e-9, rng[1] - rng[0]);
  for (let i = 0; i < n; i++) {
    const v = mode === 'h' ? P.mesh.h[i] : rhoOf(pos, i, E, D);
    const c = cmap((v - rng[0]) / span);
    cols[3 * i] = c[0]; cols[3 * i + 1] = c[1]; cols[3 * i + 2] = c[2];
  }
  P.meshCol = { key, cols };
  return cols;
}

// 锥面画多大: 晕面离灯最远处 (不超过截断距离), 没有晕面时取灯距
function coneReach() {
  const { L } = scene();
  let m = 0;
  if (R) for (const P of R.paths) {
    const a = surfPts(P), n = a.length / 3, st = Math.max(1, Math.floor(n / 5000));
    for (let i = 0; i < n; i += st) m = Math.max(m, Math.hypot(a[3 * i] - L[0], a[3 * i + 1] - L[1], a[3 * i + 2] - L[2]));
  }
  return Math.min(S.maxDist, Math.max(S.dist, m * 1.1));
}
// 眼看灯的仰角 (眼高不为 0 时与 Σ 略有不同)
function lampElevFromEye() {
  const { L, E } = scene();
  return Math.atan2(L[2] - E[2], Math.hypot(L[0] - E[0], L[1] - E[1])) / DEG;
}

// 二自由度晕面的外缘: 把相邻 h 层的层曲线端点连起来. 精确边界引擎只认得平行光曲线的末端点身份 (曲线∩边界、
//   转折、尖点), 而晕面还会在别处收口 —— 那些末端点映到灯后 / 眼后看不见, 或者曲线挨着链失效区结束 (如 Lowitz 4-1).
//   层曲线本身已二分到真正的端点 (见 epdiv 修正版), 端点随 h 连成的线就是实际画出来的晕面的边.
//   不算端点的: 灯 / 眼上的捏合点; 同层两段首尾相接处 (锥面折叠或取向坐标接缝, 曲线在这里光滑穿过, 不是边).
// 焦散面: 用引擎沿 h 精确延拓的外缘 (rims2); 还没算完时先不画 (避免先出一版层端点连线再跳变).
// 选的是某条边界面: 用层端点连线 (这类面的层端点已二分到精确端点, 连起来就在边界上).
function rimsOf(P) {
  if (P.rims2) return P.rims2;
  if (P.rims2 === undefined && S.surf2 === 'caustic') return [];
  return layerRims(P);
}
function layerRims(P) {
  const ver = (P.lver || 0) + (P.edges ? 0.5 : 0);
  if (P.rims && P.rims.v === ver) return P.rims.segs;
  const { L, E } = scene(), D = baseline();
  const all = [];
  P.layers.forEach((l, ih) => { if (l && l.curves.length) all.push({ h: layerH(P, ih), curves: l.curves }); });
  for (const l of P.dense) if (l.curves.length) all.push(l);
  all.sort((a, b) => a.h - b.h);
  const near = (a, b) => Math.hypot(a.x - b.x, a.y - b.y, a.z - b.z);
  // 精确边界线的点 (加密过, 间距 ~0.1 m) 放进网格哈希: 落在精确边界上的端点不再画外缘 (避免两条线贴着走成"双层边")
  const cell = Math.max(0.25, 0.02 * D), hash = new Map();
  const hk = (x, y, z) => Math.floor(x / cell) + ',' + Math.floor(y / cell) + ',' + Math.floor(z / cell);
  if (P.edges) for (const k of ['entry', 'tir', 'seam', 'turn', 'cusp']) for (const sg of P.edges[k] || []) {
    const p = sg.p || sg;
    for (let i = 0; i < p.length; i += 3) { const key = hk(p[i], p[i + 1], p[i + 2]); let a = hash.get(key); if (!a) hash.set(key, a = []); a.push(p[i], p[i + 1], p[i + 2]); }
  }
  const onExact = q => {
    if (!hash.size) return false;
    const tolE = Math.max(0.1, 0.02 * Math.hypot(q.x - L[0], q.y - L[1], q.z - L[2]));
    const cx = Math.floor(q.x / cell), cy = Math.floor(q.y / cell), cz = Math.floor(q.z / cell), R = Math.ceil(tolE / cell);
    for (let i = -R; i <= R; i++) for (let j = -R; j <= R; j++) for (let k = -R; k <= R; k++) {
      const a = hash.get((cx + i) + ',' + (cy + j) + ',' + (cz + k));
      if (a) for (let t = 0; t < a.length; t += 3) if (Math.hypot(a[t] - q.x, a[t + 1] - q.y, a[t + 2] - q.z) < tolE) return true;
    }
    return false;
  };
  const ends = all.map(l => {
    const raw = [];
    for (const c of l.curves) {
      const p = c.p, d = c.d, n = p.length / 3;
      if (n < 2 || !d) continue;
      for (const i of [0, n - 1]) raw.push({ x: p[3 * i], y: p[3 * i + 1], z: p[3 * i + 2], dx: d[3 * i], dy: d[3 * i + 1], dz: d[3 * i + 2] });
    }
    const tol = 2e-4 * D, Lp = { x: L[0], y: L[1], z: L[2] }, Ep = { x: E[0], y: E[1], z: E[2] };
    return raw.filter((q, i) => near(q, Lp) > 1e-3 * D && near(q, Ep) > 1e-3 * D && !raw.some((r, j) => j !== i && near(q, r) < tol)
      && !onExact(q));
  });
  // 逐层贪心配对: 还开着的链尾 (最多空过 2 层) 与新一层的端点, 按距离从近到远配, 两边都没配过且在上限内才连
  const chains = [];
  let open = [];
  ends.forEach((layer, k) => {
    const pairs = [];
    for (let ci = 0; ci < open.length; ci++) {
      const c = open[ci], q = c.pts[c.pts.length - 1];
      const dh = Math.abs(all[k].h - all[c.k].h);
      const r = Math.hypot(q.x - L[0], q.y - L[1], q.z - L[2]);
      const cap = Math.max(0.1 * D, 0.3 * r) * Math.min(1.5, dh / 0.5);
      layer.forEach((e, j) => { const dd = near(q, e); if (dd < cap) pairs.push([dd, ci, j]); });
    }
    pairs.sort((x, y) => x[0] - y[0]);
    const usedC = new Uint8Array(open.length), usedE = new Uint8Array(layer.length);
    for (const [, ci, j] of pairs) {
      if (usedC[ci] || usedE[j]) continue;
      usedC[ci] = usedE[j] = 1;
      open[ci].pts.push(layer[j]); open[ci].k = k;
    }
    const next = [];
    open.forEach((c, ci) => { if (usedC[ci] || k - c.k < 3) next.push(c); else chains.push(c); });
    layer.forEach((e, j) => { if (!usedE[j]) next.push({ k, pts: [e] }); });
    open = next;
  });
  chains.push(...open);
  const segs = [];
  for (const c of chains) {
    if (c.pts.length < 2) continue;
    const p3 = [], d3 = [];
    for (const q of c.pts) { p3.push(q.x, q.y, q.z); d3.push(q.dx, q.dy, q.dz); }
    segs.push({ p: Float32Array.from(p3), d: Float32Array.from(d3) });
  }
  P.rims = { v: ver, segs };
  return segs;
}

function buildDL() {
  const { L, E } = scene();
  const dl = { eye: E, lamp: L, mirror: S.lay.mirror, meshes: [], clouds: [], lines: [], marks: [], cone: null, shell: null };
  const lay = S.lay;
  const line = (segs, col, width, extra = {}) => { if (segs && segs.length) dl.lines.push({ segs, col, width, ...extra }); };
  const rng = colorRange(lay.colorBy);
  dl.colorRange = rng;
  if (R) R.paths.forEach((P, ip) => {
    const pcol = lay.colorBy === 'path' ? PATHCOL[ip % PATHCOL.length] : C.pts;
    if (lay.surf && P.mesh && P.mesh.pos.length)
      dl.meshes.push({ pos: P.mesh.pos, dir: P.mesh.dir, nrm: P.mesh.nrm, cols: meshColors(P, lay.colorBy, rng), col: pcol });
    // 平面退化 (灯眼等高的板晶柱面光路): 晕面整块躺在灯眼所在的水平面里, 天空里全在地平线上 — 画成地平线上的弧
    if (lay.surf && P.mesh && P.mesh.sky) line(P.mesh.sky, pcol, 3, { skyOnly: true });
    const ob = P.mesh && P.mesh.outline;            // 平面退化时的边界线 (边界引擎在这里也退化, 用网格的精确区间端点)
    if (ob) { if (lay.bEntry) line(ob.entry, C.entry, 2); if (lay.bMirror) line(ob.mirror, C.mirror, 2); if (lay.bTir) line(ob.tir, C.tir, 2); }
    if (lay.surf && P.dof === 2) {
      // 二自由度: 每个 h 层是一条 (或几条) 平行光曲线经锥映射 → 按层着色的线 (数据与配色不变就复用同一组对象)
      const lk = [lay.colorBy, rng ? rng.join(',') : '', P.lver || 0, P.fver || 0, pcol, busy.dense ? 1 : 0].join('|');
      if (!P.layerLines || P.layerLines.key !== lk) {
        const D = baseline(), span = rng ? Math.max(1e-9, rng[1] - rng[0]) : 1;
        const colOf = (h, cp) => rng && lay.colorBy === 'h' ? cmap((h - rng[0]) / span)
          : rng && lay.colorBy === 'rho' ? cmap((rhoOf(cp, Math.floor(cp.length / 6), E, D) - rng[0]) / span) : pcol;
        const segs = [], segCols = [], fSegs = [], fCols = [], fAl = [];
        for (const l of baseLayers(P)) for (const c of l.curves) { segs.push(c); segCols.push(colOf(l.h, c.p || c)); }
        for (const F of P.fill || []) for (const c of F.kept) { fSegs.push(c); fCols.push(colOf(F.h, c.p)); fAl.push(c.a); }
        P.layerLines = { key: lk, segs, segCols, fSegs, fCols, fAl };
      }
      const LL = P.layerLines, wl = P.dense.length ? 1 : 1.3;
      if (LL.segs.length) dl.lines.push({ segs: LL.segs, segCols: LL.segCols, col: pcol, width: wl, alpha: 0.85 });
      if (LL.fSegs.length) dl.lines.push({ segs: LL.fSegs, segCols: LL.fCols, vAlpha: LL.fAl, col: pcol, width: wl, alpha: 0.85, only3d: true });
    }
    if (lay.pts) {
      const cl = pathCloud(P);
      if (cl.pts.length) {
        const cols = colorize(cl, lay.colorBy, ip, rng);
        dl.clouds.push({ pts: cl.pts, dir: cloudDirs(cl), cols, col: lay.colorBy === 'path' ? PATHCOL[ip % PATHCOL.length] : C.pts, size: 1.6, alpha: cols ? 0.6 : 0.45 });
      }
    }
    if (P.bdry) {
      if (lay.bEntry) line(P.bdry.entry, C.entry, 2);
      if (lay.bMirror) line(P.bdry.mirror, C.mirror, 2);
      if (lay.bTir) line(P.bdry.tir, C.tir, 2);
      if (lay.parArc) line(P.bdry.parArc, C.parArc, 1.6, { dash: [6, 4], skyOnly: true });
    }
    if (P.dof === 2 && lay.bRim) line(rimsOf(P), C.rim, 1.6);
    if (P.edges) {
      if (lay.bEntry) line(P.edges.entry, C.entry, 2);
      if (lay.bTir) line(P.edges.tir, C.tir, 2);
      if (lay.bSeam) line(P.edges.seam, C.seam, 2);
      if (lay.bTurn) line(P.edges.turn, C.turn, 2);
      if (lay.bCusp) line(P.edges.cusp, C.cusp, 2);
      if (lay.parArc) line(P.edges.parArc, C.parArc, 1.6, { dash: [6, 4], skyOnly: true });
    }
  });
  if (RO && lay.caustic) for (const o of RO.paths) {
    if (!o) continue;
    line(o.ridge, C.caustic, 2.8);
    line(o.ridgeInf, C.caustic, 2.4, { dash: [7, 4] });
    if (lay.cSeam) line(o.seam, C.cSeam, 1.8);
    if (lay.cAsym) line(o.asym, C.asym, 1.2, { dash: [4, 4] });
  }
  if (RC && RC.data && S.cone.on) {
    const d = RC.data;
    if (RC.dof === 1) { if (d.crystal) dl.marks.push({ pts: d.crystal, col: C.cone, kind: 'pt', size: 1.6 }); }
    else line(d.curves, C.cone, 2.4);
  }
  if (S.cone.on) dl.cone = { h: S.cone.h, R: coneReach(), alpha: S.cone.alpha / 100 * 0.35 };
  if (RS && S.shell.on) {
    for (const p of RS.paths) {
      if (!p) continue;
      p.curves.forEach((c, k) => {
        const kd = p.kinds[k];
        const col = kd === 0 ? C.shell : kd === 5 ? C.shellC : kd === 1 ? C.entry : kd === 2 ? C.tir : kd === 3 ? C.seam : C.blueEdge;
        line([c], col, kd >= 1 && kd <= 4 ? 1.8 : 2.6);
      });
      if (p.ends.length) dl.marks.push({ pts: p.ends, dirs: p.endDirs, col: [255, 255, 255], kind: 'dot', size: 3.6 });
      if (p.folds.length) dl.marks.push({ pts: p.folds, dirs: p.foldDirs, col: C.fold, kind: 'tri', size: 4 });
      if (p.cusps.length) dl.marks.push({ pts: p.cusps, col: C.cuspM, kind: 'star', size: 4.4 });
    }
  }
  if (S.shell.on) dl.shell = { r: RS ? RS.r : shellRadius(), sphere: S.shell.sphere, alpha: S.shell.alpha / 100 };
  return dl;
}

// ============================================================ 视图
const sky = new SkyView($('cvSky'), $('glSky'));
const v3 = new Scene3D($('host3d'));
let baseKey = '';
function refresh() {
  rafRefresh = 0;
  const dl = buildDL();
  // 3D 底座 (灯/眼/地面) 只在场景尺度或开关变化时重建
  const { L, E } = scene();
  const extent = Math.max(S.dist, S.maxDist * 0.6) + 10;
  const bk = JSON.stringify([L, E, extent, S.lay.foci, S.lay.ground, getLang()]);
  if (bk !== baseKey) {
    baseKey = bk;
    v3.setBase({ lamp: L, eye: E, extent, foci: S.lay.foci, ground: S.lay.ground, lampText: t('lab.lamp'), eyeText: t('lab.eye') });
  }
  v3.setScene(dl);
  sky.lampLabel = t('lab.lampDir');
  sky.refCircle = S.lay.ref22 ? 22 : 0;
  sky.setScene(dl);
  renderColorbar(dl.colorRange);
  renderLayers(); renderResults(); renderProbes(); renderCrumbs();
  const none = !codes().length;
  $('note3d').hidden = !none; $('note3d').textContent = none ? t('note.noPath') : '';
}
// 3D 取景: 灯、眼与晕面点云 (抽样, 去掉 2% 离群) 的包围盒
function frame3d() {
  const { L, E } = scene();
  const xs = [L[0], E[0]], ys = [L[1], E[1]], zs = [L[2], E[2]];
  if (R) for (const P of R.paths) {
    const parts = [surfPts(P)];
    for (const l of P.layers) if (l) for (const c of l.curves) parts.push(c.p || c);
    for (const a of parts) {
      const n = a.length / 3, st = Math.max(1, Math.floor(n / 4000));
      for (let i = 0; i < n; i += st) { xs.push(a[3 * i]); ys.push(a[3 * i + 1], -a[3 * i + 1]); zs.push(a[3 * i + 2]); }
    }
  }
  const q = (v, f) => { const s = v.slice().sort((x, y) => x - y); return s[Math.min(s.length - 1, Math.max(0, Math.round(f * (s.length - 1))))]; };
  const lo = [q(xs, 0.005), q(ys, 0.005), q(zs, 0.005)], hi = [q(xs, 0.995), q(ys, 0.995), q(zs, 0.995)];
  for (let k = 0; k < 3; k++) { lo[k] = Math.min(lo[k], L[k], E[k]); hi[k] = Math.max(hi[k], L[k], E[k]); }
  v3.frame([(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, (lo[2] + hi[2]) / 2], 0.5 * Math.hypot(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]));
}
let rafRefresh = 0, lastRefresh = 0;
function scheduleRefresh() {
  if (rafRefresh) return;
  const wait = Math.max(0, 180 - (performance.now() - lastRefresh));
  rafRefresh = setTimeout(() => { lastRefresh = performance.now(); refresh(); }, wait);
}

function renderColorbar(rng) {
  const cb = $('colorbar');
  const mode = S.lay.colorBy;
  if (!rng || !(S.lay.pts || S.lay.surf) || (mode !== 'h' && mode !== 'rho')) { cb.hidden = true; return; }
  cb.hidden = false;
  const cv = $('cvBar'), g = cv.getContext('2d');
  for (let x = 0; x < cv.width; x++) { const c = cmap(x / (cv.width - 1)); g.fillStyle = `rgb(${c[0] | 0},${c[1] | 0},${c[2] | 0})`; g.fillRect(x, 0, 1, cv.height); }
  $('cbLab').textContent = mode === 'h' ? t('cb.h') : t('cb.rho');
  const f = v => (mode === 'h' ? v.toFixed(1) + '°' : v.toFixed(2));
  $('cbMin').textContent = f(rng[0]); $('cbMax').textContent = f(rng[1]);
}

// ============================================================ 状态栏 / 面包屑 / 结果
function setStatus() {
  const st = $('status'), tx = $('statusText');
  const any = busy.main || busy.dense || busy.fill || busy.obs || busy.cone || busy.shell;
  st.classList.toggle('busy', any);
  st.classList.toggle('err', errs.length > 0);
  $('prog3d').classList.toggle('on', any); $('progSky').classList.toggle('on', any);
  if (errs.length) { tx.textContent = t('st.err', { e: errs[0].split('\n')[0].slice(0, 120) }); return; }
  if (!codes().length) { tx.textContent = t('st.noPath'); return; }
  if (busy.main && R) { tx.textContent = t('st.main', { d: R.done, n: R.total }); return; }
  if ((busy.dense || busy.fill) && R) { tx.textContent = t(busy.dense ? 'st.dense' : 'st.fill', { d: R.dDone, n: R.dTotal }); return; }
  if (any) { tx.textContent = t(busy.obs ? 'st.obs' : busy.shell ? 'st.shell' : 'st.cone'); return; }
  if (R) tx.textContent = t('st.done', { ms: Math.round(R.ms) });
}

function renderCrumbs() {
  const cs = codes();
  const mode = S.mode === 'parry' ? 'Parry' : S.mode === 'lowitz' ? 'Lowitz' : t('mode.' + S.mode);
  $('crumbs').innerHTML = cs.length
    ? `<b>${cs.map(esc).join(' ; ')}</b><span class="sep">·</span>${esc(mode)}<span class="sep">·</span>Σ ${S.elev}°<span class="sep">·</span>D ${S.dist} m`
    : '';
}
const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const fmtN = n => n >= 1e6 ? (n / 1e6).toFixed(2) + 'M' : n >= 1e4 ? (n / 1e3).toFixed(1) + 'k' : String(n);
const nPts = a => (a || []).reduce((s, x) => s + (x.p || x).length / 3, 0);

function renderResults() {
  const box = $('results');
  if (!R || !R.paths.length) { box.innerHTML = `<p class="rintro">${t('res.empty')}</p>`; return; }
  let html = '';
  R.paths.forEach((P, ip) => {
    const col = PATHCOL[ip % PATHCOL.length];
    const np = pathCloud(P).pts.length / 3;
    const rows = [];
    if (P.mesh) rows.push([t('res.mesh'), t('res.meshLine', { n: fmtN(P.mesh.nTri), ms: P.mesh.ms })]);
    if (P.dof === 2) rows.push([t('res.layers'), t('res.layersLine', {
      n: P.layers.filter(l => l && l.curves.length).length + P.dense.filter(l => l.curves.length).length,
      c: P.layers.reduce((s, l) => s + (l ? l.curves.length : 0), 0) + P.dense.reduce((s, l) => s + l.curves.length, 0) })]);
    rows.push([t('res.points'), fmtN(np)]);
    if (P.bdry) rows.push([t('res.bdry'), `${P.bdry.entry.length} / ${P.bdry.mirror.length} / ${P.bdry.tir.length}` + `<span class="lcount"> ${t('res.bdryKinds1')}</span>`]);
    if (P.edges) rows.push([t('res.bdry'), `${P.edges.entry.length} / ${P.edges.tir.length} / ${P.edges.seam.length} / ${P.edges.turn.length} / ${P.edges.cusp.length}` + `<span class="lcount"> ${t('res.bdryKinds2')}</span>`]);
    if (P.edges && P.edges.surf && P.edges.surf !== S.surf2) rows.push([t('res.surf'), t('res.surfFallback')]);
    const o = RO && S.lay.caustic ? RO.paths[ip] : null;
    if (o) rows.push([t('res.caustic'), `${o.ridge.length + o.ridgeInf.length}` + (S.lay.cSeam ? ` · ${t('res.seams')} ${o.seam.length}` : '')]);
    const sh = RS && S.shell.on ? RS.paths[ip] : null;
    if (sh) rows.push([t('res.shell'), t('res.shellLine', { a: sh.nArcs, e: sh.nEnds, f: sh.nFolds || 0 })]);
    const ms = P.left > 0 || busy.dense ? '…' : Math.round(P.ms) + ' ms';
    html += `<div class="rpath"><header><i class="pathdot" style="background:rgb(${col})"></i><code>${esc(faceStr(P.code))}</code><span>${P.dof === 1 ? t('res.dof1') : t('res.dof2')}</span><span class="rt">${ms}</span></header><dl>`;
    for (const [k, v] of rows) html += `<dt>${k}</dt><dd>${v}</dd>`;
    html += '</dl>';
    if (o && o.note) html += `<div class="rnote">${esc(engineNote(o.note))}</div>`;
    if (sh && sh.note) html += `<div class="rnote">${esc(sh.note)}</div>`;
    if (!P.left && np === 0 && !(P.mesh && P.mesh.nTri) && !(P.bdry && nPts(P.bdry.entry) + nPts(P.bdry.tir) > 0)) html += `<div class="rnote">${t('res.nothing')}</div>`;
    html += '</div>';
  });
  for (const e of errs) html += `<div class="rnote err">${esc(e.split('\n')[0])}</div>`;
  box.innerHTML = html;
  $('resBadge').hidden = !(html.includes('rnote') && S.rt !== 'res');
}
// 引擎的诊断文字原样来自 C# 窗体: 去掉 markdown 加粗, 窗体勾选名换成本页的图层名
const engineNote = s => s.replace(/\*\*/g, '').replace(/"含缝像"|"With seams"/g, `「${t('lay.cSeam')}」`);
const faceStr = code => { const f = parsePath(code); return f ? f.map(faceNumFromCompact).join('-') : code; };

// ============================================================ 图层面板
const LAYERS = [
  { g: 'lg.surface' },
  { id: 'surf', sw: 'fill', term: 'surface', colorMode: true },
  { id: 'pts', sw: 'dots', col: C.pts, term: 'pts' },
  { g: 'lg.bdry' },
  { id: 'bEntry', col: C.entry, term: 'bEntry' },
  { id: 'bMirror', col: C.mirror, term: 'bMirror', dof: 1 },
  { id: 'bTir', col: C.tir, term: 'bTir' },
  { id: 'bSeam', col: C.seam, term: 'bSeam', dof: 2 },
  { id: 'bTurn', col: C.turn, term: 'bTurn', dof: 2 },
  { id: 'bCusp', col: C.cusp, term: 'bCusp', dof: 2 },
  { id: 'bRim', col: C.rim, term: 'bRim', dof: 2 },
  { g: 'lg.caustic' },
  { id: 'caustic', col: C.caustic, w: 3.2, term: 'caustic' },
  { id: 'cSeam', col: C.cSeam, term: 'cSeam', dof: 1, sub: 'caustic' },
  { id: 'cAsym', col: C.asym, dash: true, term: 'cAsym', dof: 1, sub: 'caustic' },
  { g: 'lg.ref' },
  { id: 'parArc', col: C.parArc, dash: true, term: 'parArc' },
  { id: 'ref22', col: C.ref, dash: true, w: 1.5, term: 'ref22' },
  { id: 'foci', sw: 'none' },
  { id: 'ground', sw: 'none' },
  { id: 'mirror', sw: 'none', term: 'mirror' },
];
function layerCount(id) {
  if (!R) return '';
  let n = 0;
  const add = a => { if (a) n += a.length; };
  for (const P of R.paths) {
    if (P.bdry) { if (id === 'bEntry') add(P.bdry.entry); if (id === 'bMirror') add(P.bdry.mirror); if (id === 'bTir') add(P.bdry.tir); if (id === 'parArc') add(P.bdry.parArc); }
    const ob = P.mesh && P.mesh.outline;
    if (ob) { if (id === 'bEntry') add(ob.entry); if (id === 'bMirror') add(ob.mirror); if (id === 'bTir') add(ob.tir); }
    if (P.dof === 2 && id === 'bRim') add(rimsOf(P));
    if (P.edges) { if (id === 'bEntry') add(P.edges.entry); if (id === 'bTir') add(P.edges.tir); if (id === 'bSeam') add(P.edges.seam); if (id === 'bTurn') add(P.edges.turn); if (id === 'bCusp') add(P.edges.cusp); if (id === 'parArc') add(P.edges.parArc); }
  }
  if (RO && S.lay.caustic) for (const o of RO.paths) if (o) { if (id === 'caustic') { add(o.ridge); add(o.ridgeInf); } if (id === 'cSeam') add(o.seam); if (id === 'cAsym') add(o.asym); }
  if (id === 'pts') { let k = 0; for (const P of R.paths) k += pathCloud(P).pts.length / 3; return fmtN(k); }
  if (id === 'surf') {
    let k = 0, c = 0;
    for (const P of R.paths) { if (P.mesh) k += P.mesh.nTri; for (const l of P.layers) if (l) c += l.curves.length; for (const l of P.dense) c += l.curves.length; }
    return k ? fmtN(k) + ' △' : c ? String(c) : '';
  }
  if (['bEntry', 'bMirror', 'bTir', 'bSeam', 'bTurn', 'bCusp', 'bRim', 'caustic', 'cSeam', 'cAsym', 'parArc'].includes(id)) return String(n);
  return '';
}
let layersBuilt = '';
function renderLayers() {
  const box = $('layers'), dof = is2() ? 2 : 1;
  const sig = getLang() + dof;
  if (layersBuilt !== sig) {
    layersBuilt = sig;
    let html = '';
    for (const ly of LAYERS) {
      if (ly.g) { html += `<div class="lgroup">${t(ly.g)}</div>`; continue; }
      const swCls = ly.sw || (ly.dash ? 'dash' : '');
      const style = ly.col ? `--c:rgb(${ly.col});--w:${ly.w || 2.5}px` : '';
      html += `<label class="lrow${ly.sub ? ' lsub' : ''}" data-ly="${ly.id}"><input type="checkbox" data-lay="${ly.id}">`
        + `<i class="lsw ${swCls}" style="${style}"></i><span class="lname">${t('lay.' + ly.id)}</span><span class="lcount" data-cnt="${ly.id}"></span>`
        + (ly.term ? `<button class="qmark" data-term="${ly.term}" aria-label="?">?</button>` : '') + '</label>';
      if (ly.colorMode) html += `<div class="color-mode"><span class="lcount">${t('lay.colorBy')}</span><select id="selColor">`
        + ['h', 'rho', 'mono', 'path'].map(m => `<option value="${m}">${t('color.' + m)}</option>`).join('') + '</select></div>';
    }
    box.innerHTML = html;
    for (const cb of box.querySelectorAll('input[data-lay]')) cb.addEventListener('change', () => {
      S.lay[cb.dataset.lay] = cb.checked; save();
      const id = cb.dataset.lay;
      if (id === 'caustic' || id === 'cSeam' || id === 'cAsym') runObs();
      if (is2() && S.shell.on && ['caustic', 'bEntry', 'bTir', 'bSeam'].includes(id)) runShell();
      refresh(); setStatus();
    });
    $('selColor').addEventListener('change', e => { S.lay.colorBy = e.target.value; save(); refresh(); });
    bindTerms(box);
  }
  for (const ly of LAYERS) {
    if (!ly.id) continue;
    const row = box.querySelector(`[data-ly="${ly.id}"]`), cb = row.querySelector('input');
    cb.checked = !!S.lay[ly.id];
    const na = ly.dof && ly.dof !== dof;
    row.classList.toggle('na', !!na);
    row.classList.toggle('off', !cb.checked);
    row.hidden = !!(ly.sub && !S.lay[ly.sub]) || !!na;
    const c = box.querySelector(`[data-cnt="${ly.id}"]`); if (c) c.textContent = layerCount(ly.id);
  }
  $('selColor').value = S.lay.colorBy;
}

// ============================================================ 左栏控件
function bindPair(sl, inp, get, set, opt = {}) {
  const toSl = opt.toSl || (v => v), fromSl = opt.fromSl || (v => v);
  const sync = () => { const v = get(); sl.value = toSl(v); inp.value = opt.fmt ? opt.fmt(v) : v; };
  sl.addEventListener('input', () => { set(fromSl(+sl.value)); inp.value = opt.fmt ? opt.fmt(get()) : get(); changed(); });
  inp.addEventListener('change', () => {
    let v = parseFloat(inp.value);
    if (!Number.isFinite(v)) { sync(); return; }
    v = clamp(v, +inp.min, +inp.max);
    set(v); sync(); changed();
  });
  sync();
  return sync;
}
let afterChange = null;
function changed() { save(); drawSceneFig(); renderModeUI(); syncShellPair(); schedule(); if (afterChange) afterChange(); }

const distToSl = d => Math.round(1000 * Math.log(d) / Math.log(300));
const slToDist = v => { const d = Math.exp(v / 1000 * Math.log(300)); return d < 20 ? Math.round(d * 2) / 2 : Math.round(d); };
const syncers = [];
syncers.push(bindPair($('slElev'), $('inElev'), () => S.elev, v => { S.elev = v; }));
syncers.push(bindPair($('slDist'), $('inDist'), () => S.dist, v => { S.dist = clamp(v, 1, 300); }, { toSl: distToSl, fromSl: slToDist }));
syncers.push(bindPair($('slEye'), $('inEye'), () => S.eye, v => { S.eye = Math.round(v * 100) / 100; }));
syncers.push(bindPair($('slZeta'), $('inZeta'), () => (S.mode === 'column' ? S.colZeta : S.zeta), v => { if (S.mode === 'column') S.colZeta = v; else S.zeta = v; }));
syncers.push(bindPair($('slTheta'), $('inTheta'), () => (S.mode === 'lowitz' ? S.lowTheta : S.theta), v => { if (S.mode === 'lowitz') S.lowTheta = v; else S.theta = v; }));

function drawSceneFig() {
  const svg = $('sceneFig');
  const W = 300, H = 96, gy = 82, ex = 34;
  const e = clamp(S.elev, -80, 80) * DEG;
  const len = 150 + 60 * Math.min(1, Math.log10(S.dist) / Math.log10(300));
  const eyeY = gy - 6 - Math.min(14, S.eye * 6);
  let lx = ex + len * Math.cos(e), ly = eyeY - len * Math.sin(e);
  const sc = Math.min(1, (eyeY - 10) / Math.max(1, eyeY - ly), (gy + 8 - eyeY) / Math.max(1, ly - eyeY + 1e-9) * (ly > eyeY ? 1 : 99));
  lx = ex + (lx - ex) * Math.min(1, Math.abs(sc)); ly = eyeY + (ly - eyeY) * Math.min(1, Math.abs(sc));
  const ar = 26, a0 = 0, a1 = -e;
  const arc = `M ${ex + ar} ${eyeY} A ${ar} ${ar} 0 0 ${e > 0 ? 0 : 1} ${ex + ar * Math.cos(a1)} ${eyeY + ar * Math.sin(a1)}`;
  svg.innerHTML = `<line class="gnd" x1="0" y1="${gy}" x2="${W}" y2="${gy}"/>`
    + `<line class="base" x1="${ex}" y1="${eyeY}" x2="${W - 8}" y2="${eyeY}" stroke-dasharray="2 3"/>`
    + `<line class="ray" x1="${ex}" y1="${eyeY}" x2="${lx}" y2="${ly}"/>`
    + (Math.abs(S.elev) > 0.5 ? `<path class="arc" d="${arc}"/>` : '')
    + `<text x="${ex + ar + 4}" y="${eyeY + (e > 0 ? -4 : 12)}">Σ ${S.elev}°</text>`
    + `<text x="${(ex + lx) / 2 - 10}" y="${(eyeY + ly) / 2 + (e > 0 ? -6 : 14)}">D ${S.dist} m</text>`
    + `<circle class="eyec" cx="${ex}" cy="${eyeY}" r="4"/><text x="${ex - 26}" y="${eyeY + 4}">${t('lab.eye')}</text>`
    + (S.eye > 0 ? `<line class="base" x1="${ex}" y1="${eyeY + 4}" x2="${ex}" y2="${gy}"/>` : '')
    + `<circle class="lampc" cx="${lx}" cy="${ly}" r="5.5"/><text x="${lx + 8}" y="${ly + 4}">${t('lab.lamp')}</text>`;
  void a0;
}

// ---------- 光路 + 晶体 ----------
const cry = new CrystalPreview($('cvCrystal'), f => {
  const parts = S.paths.split(/[;,\n]/);
  const last = parts.length - 1;
  const cur = parts[last].trim();
  parts[last] = (cur ? cur + '-' : '') + f;
  S.paths = parts.map(s => s.trim()).filter((s, i) => s.length || i === last).join('; ');
  $('inPath').value = S.paths; onPathChanged();
});
function onPathChanged() {
  const cs = S.paths.split(/[;,\n]/).map(s => s.trim()).filter(s => s.length);
  const bad = cs.filter(c => parsePath(c) === null);
  $('pathMsg').textContent = bad.length ? t('path.bad', { p: bad.join(', ') }) : '';
  $('inPath').classList.toggle('invalid', bad.length > 0);
  renderChain();
  const f0 = parsePath(cs[0] || '');
  cry.set([1, 1, 1, 1, 1, 1], 1, f0 ? f0.map(faceNumFromCompact) : []);
  changed();
}
function renderChain() {
  const cs = codes();
  $('chain').innerHTML = cs.map((c, ip) => {
    const f = parsePath(c).map(faceNumFromCompact), col = PATHCOL[ip % PATHCOL.length];
    return `<div class="chain-row"><i class="pathdot" style="background:rgb(${col})"></i>` + f.map((x, i) => {
      const role = i === 0 ? 'in' : i === f.length - 1 ? 'out' : 'mid';
      return (i ? '<span class="arrow">→</span>' : '') + `<span class="face ${role}"><b>${x}</b><small>${t('role.' + role)}</small></span>`;
    }).join('') + '</div>';
  }).join('');
}
$('inPath').addEventListener('input', () => { S.paths = $('inPath').value; onPathChanged(); });
$('btnBack').addEventListener('click', () => {
  S.paths = S.paths.replace(/\s*[-;,]?\s*\d+\s*$/, '').replace(/[;,]\s*$/, '');
  $('inPath').value = S.paths; onPathChanged();
});
$('btnNewPath').addEventListener('click', () => { if (S.paths.trim()) { S.paths = S.paths.trim() + '; '; $('inPath').value = S.paths; $('inPath').focus(); } });
$('btnClear').addEventListener('click', () => { S.paths = ''; $('inPath').value = ''; onPathChanged(); });
for (const b of $('crySeg').querySelectorAll('button')) b.addEventListener('click', () => { S.pyr = b.dataset.cry === '1'; renderSeg('crySeg', S.pyr ? '1' : '0', 'cry'); cry.setShape(S.pyr, S.apexUp, S.apexLo); save(); });

// ---------- 取向 ----------
for (const b of $('modes').querySelectorAll('button')) b.addEventListener('click', () => { S.mode = b.dataset.mode; syncers.forEach(f => f()); changed(); });
function renderModeUI() {
  for (const b of $('modes').querySelectorAll('button')) { b.classList.toggle('on', b.dataset.mode === S.mode); b.setAttribute('aria-checked', b.dataset.mode === S.mode); }
  $('modeDesc').innerHTML = t('mdesc.' + S.mode);
  const showZ = S.mode === 'custom' || S.mode === 'column', showT = S.mode === 'custom' || S.mode === 'lowitz';
  $('rowZeta').hidden = !showZ; $('rowTheta').hidden = !showT;
  $('lbZeta').textContent = t(S.mode === 'column' ? 'o.zetaCol' : 'o.zeta');
  $('lbTheta').textContent = t(S.mode === 'lowitz' ? 'o.thetaLow' : 'o.theta');
  $('rowSurf').hidden = !is2();
}
function fillSurf() {
  $('selSurf').innerHTML = SURF_IDS.map(id => `<option value="${id}">${t('surf.' + id)}</option>`).join('');
  $('selSurf').value = S.surf2;
}
$('selSurf').addEventListener('change', e => { S.surf2 = e.target.value; changed(); });

// ---------- 更多选项 ----------
function renderSeg(id, val, attr) {
  for (const b of $(id).querySelectorAll('button')) { const on = b.dataset[attr] === val; b.classList.toggle('on', on); b.setAttribute('aria-checked', on); }
}
for (const b of $('qSeg').querySelectorAll('button')) b.addEventListener('click', () => { S.q = +b.dataset.q; renderSeg('qSeg', String(S.q), 'q'); renderMoreTags(); changed(); });
const numIn = (id, get, set) => {
  const el = $(id);
  el.value = get();
  el.addEventListener('change', () => {
    const v = parseFloat(el.value);
    if (!Number.isFinite(v) || v < +el.min || v > +el.max) { el.classList.add('invalid'); return; }
    el.classList.remove('invalid'); set(v); renderMoreTags(); changed();
  });
};
numIn('inMaxD', () => S.maxDist, v => { S.maxDist = v; });
numIn('inN', () => S.n, v => { S.n = v; });
numIn('inApexUp', () => S.apexUp, v => { S.apexUp = v; cry.setShape(S.pyr, S.apexUp, S.apexLo); renderMiller(); });
numIn('inApexLo', () => S.apexLo, v => { S.apexLo = v; cry.setShape(S.pyr, S.apexUp, S.apexLo); renderMiller(); });
$('btnNReset').addEventListener('click', () => { S.n = 1.31; $('inN').value = 1.31; renderMoreTags(); changed(); });
function renderMiller() {
  for (const [id, a] of [['millerUp', S.apexUp], ['millerLo', S.apexLo]]) {
    const [p, q] = millerFromApexDeg(a);
    $(id).textContent = `≈ (${p} 0 ${-p} ${q}) ${apexDegFromMiller(p, q).toFixed(1)}°`;
  }
}
function renderMoreTags() {
  const tags = [];
  if (S.q) tags.push(t(S.q === 1 ? 'more.qMid' : 'more.qHigh'));
  if (S.maxDist !== 80) tags.push(`${S.maxDist} m`);
  if (S.n !== 1.31) tags.push(`n ${S.n}`);
  $('moreTags').innerHTML = tags.map(x => `<span>${esc(x)}</span>`).join('');
}

// ---------- 常见晕弧 ----------
function renderPresets() {
  $('presetsPop').innerHTML = PRESETS.map(p => `<button class="preset" data-id="${p.id}"><b>${t('pre.' + p.id)}</b><code>${p.path}</code><small>${t('pre.' + p.id + '.d')}</small></button>`).join('');
  for (const b of $('presetsPop').querySelectorAll('.preset')) b.addEventListener('click', () => applyPreset(PRESETS.find(p => p.id === b.dataset.id)));
}
function applyPreset(p) {
  S.paths = p.path; S.mode = p.mode; S.elev = p.elev; S.dist = p.dist;
  if (p.mode === 'column') S.colZeta = -90;
  if (p.mode === 'lowitz') S.lowTheta = 0;
  S.sky.az = p.sky[0]; S.sky.el = p.sky[1]; S.sky.fov = p.sky[2];
  $('inPath').value = S.paths;
  syncers.forEach(f => f()); syncSkyUI(); sky.setView(S.sky);
  togglePresets(false);
  onPathChanged();
}
function togglePresets(on) {
  $('presetsPop').hidden = !on; $('btnPresets').setAttribute('aria-expanded', on);
}
$('btnPresets').addEventListener('click', e => { e.stopPropagation(); togglePresets($('presetsPop').hidden); });
document.addEventListener('click', e => { if (!$('presetsPop').hidden && !$('presetsPop').contains(e.target)) togglePresets(false); });

// ============================================================ 视图切换 / 天空工具条 / 3D 镜头
// 舞台太窄 (< 820 px) 或窄屏页签时不并排: 只显示一个视图 (S.single 记住最后单看的那个)
const SPLIT_MIN = 820;
let narrowStage = false;
const effView = () => (S.view === 'split' && narrowStage ? S.single : S.view);
function renderView() {
  const ev = effView();
  renderSeg('viewSeg', ev, 'view');
  const v = $('views');
  v.classList.toggle('v-3d', ev === '3d'); v.classList.toggle('v-sky', ev === 'sky');
  v3.setActive(ev !== 'sky');
  if (typeof renderTabs === 'function') renderTabs();
  requestAnimationFrame(() => { sky.draw(); v3.dirty = true; });
}
for (const b of $('viewSeg').querySelectorAll('button')) b.addEventListener('click', () => {
  S.view = b.dataset.view;
  if (S.view !== 'split') S.single = S.view;
  save(); renderView();
});
new ResizeObserver(() => {
  const w = $('stage').getBoundingClientRect().width, nw = w > 0 && w < SPLIT_MIN;
  if (nw !== narrowStage) { narrowStage = nw; $('app').classList.toggle('narrow-stage', nw); renderView(); }
}).observe($('stage'));

// 右栏: 收起 / 三页
function renderPanel() {
  $('app').classList.toggle('no-read', !S.panel);
  $('btnPanel').setAttribute('aria-pressed', S.panel);
  $('read').dataset.rt = S.rt;
  for (const b of $('readTabs').querySelectorAll('button')) b.classList.toggle('on', b.dataset.rt === S.rt);
  $('resBadge').hidden = !($('results').innerHTML.includes('rnote') && S.rt !== 'res');
  requestAnimationFrame(() => { sky.draw(); v3.dirty = true; });
}
$('btnPanel').addEventListener('click', () => { S.panel = !S.panel; save(); renderPanel(); });
for (const b of $('readTabs').querySelectorAll('button')) b.addEventListener('click', () => { S.rt = b.dataset.rt; save(); renderPanel(); });

// 导出菜单
function toggleMenu(on) { $('exportMenu').hidden = !on; $('btnExport').setAttribute('aria-expanded', on); }
$('btnExport').addEventListener('click', e => { e.stopPropagation(); toggleMenu($('exportMenu').hidden); });
document.addEventListener('click', e => { if (!$('exportMenu').hidden && !$('exportMenu').contains(e.target)) toggleMenu(false); });

// 晶体点选 (默认收起)
$('crystalDet').open = !!S.crystalOpen;
$('crystalDet').addEventListener('toggle', () => { S.crystalOpen = $('crystalDet').open; save(); if (S.crystalOpen) requestAnimationFrame(() => cry.draw()); });
for (const b of document.querySelectorAll('[data-cam]')) b.addEventListener('click', () => v3.preset(b.dataset.cam));

function fillLens() {
  $('selLens').innerHTML = Array.from({ length: LENS_COUNT }, (_, i) => `<option value="${i}">${t('lens.' + i)}</option>`).join('');
  $('selLens').value = S.sky.lens;
}
function syncSkyUI() { $('selLens').value = S.sky.lens; $('inFov').value = Math.round(S.sky.fov); $('inFov').max = fovMax(S.sky.lens); }
$('selLens').addEventListener('change', e => { S.sky.lens = +e.target.value; S.sky.fov = Math.min(S.sky.fov, fovMax(S.sky.lens)); if (S.sky.lens === 5) S.sky.el = 0; sky.setView(S.sky); syncSkyUI(); save(); });
$('inFov').addEventListener('change', () => { const v = parseFloat($('inFov').value); if (Number.isFinite(v)) { S.sky.fov = clamp(v, 5, fovMax(S.sky.lens)); sky.setView(S.sky); syncSkyUI(); save(); } });
for (const b of document.querySelectorAll('[data-look]')) b.addEventListener('click', () => {
  const k = b.dataset.look, lampEl = lampElevFromEye();
  if (k === 'lamp') { S.sky.az = 0; S.sky.el = clamp(lampEl, -90, 90); }
  else if (k === 'anti') { S.sky.az = 180; S.sky.el = clamp(-lampEl + 20, -90, 90); }
  else { S.sky.az = 0; S.sky.el = 90; }
  sky.setView(S.sky); save();
});
sky.onView = V => { Object.assign(S.sky, { az: V.az, el: V.el, fov: V.fov }); $('inFov').value = Math.round(V.fov); save(); };
sky.onHover = h => {
  $('readout').innerHTML = h ? `${t('ro.az')} <b>${h.az.toFixed(1)}°</b>  ${t('ro.el')} <b>${h.el.toFixed(1)}°</b><br>${t('ro.fromLamp')} <b>${h.ang.toFixed(2)}°</b>` : '';
};

// ============================================================ 切片工具
syncers.push(bindPair($('slConeH'), $('inConeH'), () => S.cone.h, v => { S.cone.h = v; }));
$('ckCone').addEventListener('change', e => { S.cone.on = e.target.checked; save(); runCone(); renderProbes(); refresh(); setStatus(); });
$('slConeA').addEventListener('input', e => { S.cone.alpha = +e.target.value; save(); refresh(); });
$('ckShell').addEventListener('change', e => { S.shell.on = e.target.checked; save(); syncShellPair(); runShell(); renderProbes(); refresh(); setStatus(); });
$('ckSphere').addEventListener('change', e => { S.shell.sphere = e.target.checked; save(); refresh(); });
$('slShellA').addEventListener('input', e => { S.shell.alpha = +e.target.value; save(); refresh(); });
$('slKappa').addEventListener('input', e => { S.shell.kappa = +e.target.value / 100; S.shell.master = 'kappa'; syncShellPair(); save(); schedule(); });
$('inKappa').addEventListener('change', e => {
  const v = parseFloat(e.target.value); if (!Number.isFinite(v)) return;
  S.shell.kappa = clamp(v, -3, 3); S.shell.master = 'kappa'; syncShellPair(); save(); schedule();
});
$('inRad').addEventListener('change', e => {
  const v = parseFloat(e.target.value); if (!(v > 0)) return;
  S.shell.r = v; S.shell.master = 'r'; syncShellPair(); save(); schedule();
});
// κ 与 r 互为主从: 用户最后改的那个是主量, 灯/眼变动时保持它, 回填另一个
function syncShellPair() {
  const D = baseline();
  if (!(D > 1e-9)) return;
  if (S.shell.master === 'r') S.shell.kappa = Math.round(clamp(Math.log10(S.shell.r / D), -3, 3) * 1000) / 1000;
  else S.shell.r = Math.round(D * Math.pow(10, S.shell.kappa) * 1000) / 1000;
  $('inKappa').value = S.shell.kappa.toFixed(3);
  $('slKappa').value = Math.round(S.shell.kappa * 100);
  $('inRad').value = S.shell.r;
}
function renderProbes() {
  $('ckCone').checked = S.cone.on; $('coneBody').classList.toggle('off', !S.cone.on);
  $('slConeA').value = S.cone.alpha;
  $('ckShell').checked = S.shell.on; $('shellBody').classList.toggle('off', !S.shell.on);
  $('ckSphere').checked = S.shell.sphere; $('slShellA').value = S.shell.alpha;
  // 锥切片状态
  let cs = '';
  if (S.cone.on) {
    if (!RC || !RC.data) cs = t('probe.computing');
    else if (RC.dof === 1) cs = t('probe.coneStat1', { n: RC.data.n });
    else cs = t('probe.coneStat2', { n: RC.data.curves.length, p: fmtN(nPts(RC.data.curves)) });
    if (RC && RC.data && !(RC.dof === 1 ? RC.data.n : RC.data.curves.length) && Math.abs(S.cone.h - lampElevFromEye()) < 0.05)
      cs += '\n' + t('probe.coneEye');
    if (codes().length > 1) cs += '\n' + t('probe.firstPath', { p: faceStr(codes()[0]) });
  }
  $('coneStatus').textContent = cs;
  let ss = '';
  if (S.shell.on) {
    const D = baseline(), r = shellRadius();
    ss = t('probe.shellHead', { rho: (r / D).toPrecision(4), r: r.toPrecision(5), D: D.toPrecision(5) });
    if (!RS || RS.paths.some(p => !p)) ss += '\n' + t('probe.computing');
    else {
      let a = 0, e = 0, f = 0, notes = [];
      RS.paths.forEach((p, i) => { a += p.nArcs; e += p.nEnds; f += p.nFolds || 0; if (p.note) notes.push(faceStr(codes()[i]) + ': ' + p.note); });
      ss += '\n' + t('probe.shellStat', { a, e, f, ms: Math.round(RS.ms || 0) });
      if (a === 0) ss += '\n' + t('probe.shellNone');
      if (is2()) ss += '\n' + t('probe.shell2Hint');
      for (const n of notes.slice(0, 2)) ss += '\n⚠ ' + n;
    }
  }
  $('shellStatus').textContent = ss;
}

// ============================================================ 术语提示 / 术语表
const tip = $('tip');
function showTip(btn) {
  const tm = term(btn.dataset.term);
  if (!tm) return;
  tip.innerHTML = `<b>${tm.title}</b>${tm.body}`;
  tip.hidden = false;
  const r = btn.getBoundingClientRect(), tw = tip.offsetWidth, th = tip.offsetHeight;
  let x = r.left + r.width / 2 - tw / 2, y = r.bottom + 8;
  x = clamp(x, 8, innerWidth - tw - 8);
  if (y + th > innerHeight - 8) y = r.top - th - 8;
  tip.style.left = x + 'px'; tip.style.top = y + 'px';
}
function bindTerms(root) {
  for (const b of root.querySelectorAll('.qmark[data-term]')) {
    b.addEventListener('mouseenter', () => showTip(b));
    b.addEventListener('mouseleave', () => { tip.hidden = true; });
    b.addEventListener('click', e => { e.preventDefault(); e.stopPropagation(); if (tip.hidden) showTip(b); else tip.hidden = true; });
  }
}
document.addEventListener('scroll', () => { tip.hidden = true; }, true);
$('btnGlossary').addEventListener('click', () => {
  $('glBody').innerHTML = termKeys().map(k => { const tm = term(k); return `<div class="gl-item"><b>${tm.title}</b>${tm.body}</div>`; }).join('');
  $('glossary').showModal();
});
$('btnGlClose').addEventListener('click', () => $('glossary').close());

// ============================================================ 导出
function stamp() {
  const d = new Date(), p = n => String(n).padStart(2, '0');
  return `${d.getFullYear()}${p(d.getMonth() + 1)}${p(d.getDate())}_${p(d.getHours())}${p(d.getMinutes())}`;
}
function download(name, blob) {
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob); a.download = name;
  document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(a.href), 4000);
}
const fileTag = () => codes().map(c => faceStr(c).replace(/-/g, '')).join('_') || 'none';
$('btnPng').addEventListener('click', () => {
  toggleMenu(false);
  const parts = [], ev = effView();
  if (ev !== 'sky') parts.push(v3.snapshot());
  if (ev !== '3d') parts.push(sky.snapshot());
  const H = Math.max(...parts.map(c => c.height)), pad = 16, head = 44;
  const W = parts.reduce((s, c) => s + c.width, 0) + pad * (parts.length + 1);
  const cv = document.createElement('canvas'); cv.width = W; cv.height = H + head + pad;
  const g = cv.getContext('2d');
  g.fillStyle = '#fff'; g.fillRect(0, 0, cv.width, cv.height);
  let x = pad;
  for (const c of parts) { g.drawImage(c, x, head); g.strokeStyle = '#d0d5dc'; g.strokeRect(x + 0.5, head + 0.5, c.width - 1, c.height - 1); x += c.width + pad; }
  g.fillStyle = '#1d2129'; g.font = '600 20px system-ui, "Microsoft YaHei", sans-serif'; g.textBaseline = 'middle';
  const mode = S.mode === 'parry' ? 'Parry' : S.mode === 'lowitz' ? 'Lowitz' : t('mode.' + S.mode);
  g.fillText(`${codes().map(faceStr).join(' ; ')}  ·  ${mode}  ·  Σ ${S.elev}°  ·  D ${S.dist} m  ·  ${t('s1.eye')} ${S.eye} m`, pad, head / 2);
  cv.toBlob(b => download(`divergent_${fileTag()}_${stamp()}.png`, b));
});
$('btnCsv').addEventListener('click', () => {
  toggleMenu(false);
  if (!R) return;
  const { L, E } = scene(), D = baseline(), o = orient();
  const lines = [];
  lines.push('# Divergent Light Halo — data export');
  lines.push(`# lamp elev=${S.elev} deg, dist=${S.dist} m, L=(${L.map(v => v.toFixed(6)).join(',')}); eye height=${S.eye} m`);
  lines.push(`# mode=${S.mode} zeta=${o.zeta} theta=${o.theta}${o.dof === 2 ? ' surf=' + S.surf2 : ''}; n=${S.n}; maxDist=${S.maxDist} m; quality=${S.q}`);
  lines.push('# section,path,part,index,x,y,z,az_deg,el_deg,dist_lamp,dist_eye,h_deg,lg_rho');
  // dv: 精确视方向 (引擎给的, 靠近眼也准); 没有就由 P−E 求
  const row = (sec, path, part, idx, a, i, dv) => {
    const x = a[3 * i], y = a[3 * i + 1], z = a[3 * i + 2];
    const r = Math.hypot(x - E[0], y - E[1], z - E[2]);
    const [ux, uy, uz] = dv ? [dv[3 * i], dv[3 * i + 1], dv[3 * i + 2]] : [(x - E[0]) / Math.max(r, 1e-12), (y - E[1]) / Math.max(r, 1e-12), (z - E[2]) / Math.max(r, 1e-12)];
    const az = -Math.atan2(uy, ux) / DEG, el = Math.asin(clamp(uz, -1, 1)) / DEG;
    const rl = Math.hypot(x - L[0], y - L[1], z - L[2]);
    lines.push([sec, path, part, idx, x.toFixed(5), y.toFixed(5), z.toFixed(5), az.toFixed(4), el.toFixed(4), rl.toFixed(4), r.toFixed(4),
      hOf(a, i, L).toFixed(4), Math.log10(Math.max(1e-9, r / D)).toFixed(5)].join(','));
  };
  const dumpSegs = (sec, path, segs) => (segs || []).forEach((s, k) => { const a = s.p || s; for (let i = 0; i < a.length / 3; i++) row(sec, path, k, i, a, i, s.d); });
  R.paths.forEach((P, ip) => {
    const pc = faceStr(P.code);
    P.h.forEach((a, k) => { if (a) for (let i = 0; i < a.length / 3; i++) row('surface', pc, 'h-scan', k, a, i); });
    P.s.forEach((a, k) => { if (a) for (let i = 0; i < a.length / 3; i++) row('surface', pc, 'psi-scan', k, a, i); });
    P.layers.forEach((l, k) => { if (l && l.verts) for (let i = 0; i < l.verts.length / 3; i++) row('surface', pc, 'layer', k, l.verts, i); });
    if (P.mesh) {                     // 网格顶点 (三角形不共享顶点, 这里去重)
      const seen = new Set(), a = P.mesh.pos;
      for (let i = 0; i < a.length / 3; i++) {
        const k = a[3 * i] + ',' + a[3 * i + 1] + ',' + a[3 * i + 2];
        if (seen.has(k)) continue;
        seen.add(k); row('mesh', pc, 'vertex', seen.size - 1, a, i, P.mesh.dir);
      }
    }
    if (P.bdry) { dumpSegs('entry', pc, P.bdry.entry); dumpSegs('mirror', pc, P.bdry.mirror); dumpSegs('exit', pc, P.bdry.tir); }
    if (P.edges) { dumpSegs('entry', pc, P.edges.entry); dumpSegs('exit', pc, P.edges.tir); dumpSegs('seam', pc, P.edges.seam); dumpSegs('turn', pc, P.edges.turn); dumpSegs('cusp', pc, P.edges.cusp); }
    const ob = RO && S.lay.caustic ? RO.paths[ip] : null;
    if (ob) { dumpSegs('caustic', pc, ob.ridge); dumpSegs('caustic_unbounded', pc, ob.ridgeInf); dumpSegs('seam_image', pc, ob.seam); }
    const sh = RS && S.shell.on ? RS.paths[ip] : null;
    if (sh) dumpSegs('shell', pc, sh.curves);
  });
  download(`divergent_${fileTag()}_${stamp()}.csv`, new Blob([lines.join('\n') + '\n'], { type: 'text/csv' }));
});

// ============================================================ 语言 / 窄屏页签 / 启动
function applyLang(code) {
  S.lang = setLang(code); save();
  applyI18n();
  renderSeg('langSeg', S.lang, 'lang');
  fillLens(); fillSurf(); renderPresets(); renderModeUI(); renderMiller(); renderMoreTags(); renderChain(); drawSceneFig();
  layersBuilt = ''; baseKey = '';
  RO = null;                                   // 焦散诊断文字随语言
  schedule(true);
}
for (const b of $('langSeg').querySelectorAll('button')) b.addEventListener('click', () => applyLang(b.dataset.lang));
// 窄屏页签: 设置 | 3D | 天空 | 显示 (视图切换并进页签)
function renderTabs() {
  const tab = $('app').dataset.tab, ev = effView();
  for (const x of $('tabs').querySelectorAll('button'))
    x.classList.toggle('on', x.dataset.tab === tab && (!x.dataset.v || x.dataset.v === ev));
}
for (const b of $('tabs').querySelectorAll('button')) b.addEventListener('click', () => {
  $('app').dataset.tab = b.dataset.tab;
  if (b.dataset.v) { S.view = b.dataset.v; S.single = b.dataset.v; save(); }
  renderView(); renderTabs();
});
$('app').dataset.tab = 'stage';
new ResizeObserver(() => sky.draw()).observe($('paneSky'));
// 中等宽度屏幕: 右栏默认收起 (浮在视图上)
if (localStorage.getItem(STORE) === null && innerWidth <= 1400) S.panel = false;
renderPanel();

$('inPath').value = S.paths;
renderSeg('qSeg', String(S.q), 'q'); renderSeg('crySeg', S.pyr ? '1' : '0', 'cry');
cry.setShape(S.pyr, S.apexUp, S.apexLo);
const f0 = parsePath(codes()[0] || '');
cry.set([1, 1, 1, 1, 1, 1], 1, f0 ? f0.map(faceNumFromCompact) : []);
bindTerms(document);
syncShellPair();
sky.setView(S.sky);
renderView();
applyLang(detectLang(S.lang));
syncSkyUI();

// 调试/自动化截图用 (控制台可用): 视图对象与状态
window.__dlh = { v3, sky, S, refresh, pool, get R() { return R; } };
