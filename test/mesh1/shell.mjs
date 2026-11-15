// 一自由度网格 (球壳法: divshell.js computeShell 的弧 + mesh2s.js 缝合) 的校验:
//   node test/mesh1/shell.mjs <光路> <Σ> <D> <ζ> <θ>
//   (1) 网格 → 晕面: 同 phys.mjs 的物理校验 (三角形重心处是否存在晶体方位把灯光折进眼睛);
//   (2) 晕面 → 网格: 拿 (h, α) 网格 (mesh1.js, 另一套算法, 已单独校验过) 的顶点当真值, 量到球壳网格的距离 ——
//       球壳引擎若在某些半径上漏了弧, 这里看得见. 灯 / 眼 0.02 D 之内、截断处的不算.
import path from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';
const ROOT = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../..');
const eng = f => pathToFileURL(path.join(ROOT, 'src/engine', f)).href;
const BD = await import(eng('bdry4.js')), CR = await import(eng('crystal.js')), LPm = await import(eng('lp.js')), DS = await import(eng('divshell.js')), MS = await import(eng('mesh2s.js')), M1 = await import(eng('mesh1.js'));
const PH = await import('./phys.mjs');

export function shellMesh1(code, el, D, zeta, theta, opt = {}) {
  CR.setApex(56, 56);
  const nr = CR.parsePath(code).map(CR.baseNormalFromCompact), er = el * Math.PI / 180, L = [D * Math.cos(er), 0, D * Math.sin(er)], E = [0, 0, 0];
  const compute = (r, seeds) => {
    const res = DS.computeShell(nr, zeta, theta, 1.31, L, E, r, true, 0, 360, { seeds: process.env.NOSEED ? null : seeds });
    const cusps = res.cuspP.slice();
    for (const a of res.arcs) for (const q of a.foldP) cusps.push(q);
    return { arcs: res.arcs.filter(a => a.p.length >= 2).map(a => ({ p: a.p, d: a.dir, closed: a.closed })), cusps, states: res.arcs.map(a => a.w) };
  };
  // 临界半径: 已经解出来的三维边界线 (入射掠射 / 镜面掠射 / 出射临界) 离眼距离的极值与端点 —— 弧的末端点在边界线上,
  //   层线的条数 / 连法只在这些半径上变 (曲面的尖端 = 两条边界线的交汇点, 弧在那里缩成一点)
  let critical = [];
  if (!process.env.NOCRIT) {
    LPm.LP.partial = true; LPm.setIndex(1.31);
    const b = BD.lightConeBoundaries(CR.parsePath(code), el, D, 0, 1.31, 80, zeta, theta);
    critical = MS.criticalRadii([...b.entry, ...b.mirror, ...b.tir].map(c => c.p || c), L, E);
  }
  return MS.meshShells(MS.shellSource(compute), L, E, { maxD: 80, critical, log: process.env.LOG ? console.log : null, ...opt });
}

export function check1(code, el, D, zeta, theta) {
  const m = shellMesh1(code, el, D, zeta, theta);
  const f = PH.physCheck(code, el, D, zeta, theta, 160000, 80, m);
  // 反向: mesh1 的顶点到球壳网格
  CR.setApex(56, 56); LPm.LP.partial = true; LPm.setIndex(1.31);
  const ref = M1.surfaceMesh1(CR.parsePath(code), el, D, 0, 1.31, 80, zeta, theta, { budget: 60000 });
  const er = el * Math.PI / 180, L = [D * Math.cos(er), 0, D * Math.sin(er)];
  const cv = ref.nTri ? MS.coverCheck(m, [{ curves: [{ p: ref.pos }] }], L, [0, 0, 0], 80, 3000) : { n: 0, miss: 0, frac: 0 };
  return { m, f, cv, refTri: ref.nTri };
}

if (process.argv[1] && path.resolve(process.argv[1]) === fileURLToPath(import.meta.url)) {
  const [code = '3-5', el = '10.5', D = '15', zeta = '-138', theta = '60'] = process.argv.slice(2);
  const r = check1(code, +el, +D, +zeta, +theta), m = r.m;
  console.log(`${code} Σ${el} D${D} ζ${zeta} θ${theta}: 球壳法 ${m.nTri} 三角形 ${r.f.area.toFixed(2)} m², ${m.stat.rows} 层 (引擎算了 ${m.stat.calls} 层), 变化 ${m.stat.events}, 没缝上 ${m.stat.loose} (${m.stat.looseDeg.toFixed(1)}°), ${m.stat.ms} ms`);
  console.log(`  网格→晕面: 偏差 > 3 cm + 0.3% r 的面积 ${r.f.bad2.toFixed(3)} m², > 10 cm + 1% r 的 ${r.f.bad10.toFixed(3)} m²`);
  console.log(`  晕面→网格: (h, α) 网格 ${r.refTri} 三角形的顶点抽 ${r.cv.n} 个, 离球壳网格 > 5 cm + 1% r 的 ${r.cv.miss} 个 (${(100 * r.cv.frac).toFixed(2)}%)`);
}
