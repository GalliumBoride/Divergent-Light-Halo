// 二自由度网格 (球壳法, engine/mesh2s.js —— 页面用的就是它) 的几何校验: 同 phys.mjs 的双向、三维 + 天空两遍.
//   真值取光锥法的层线和精确边界线: 两套引擎互相独立 (球壳轮廓来自 divshell.js, 层线来自 epdiv.js), 所以这也是两套引擎的对账.
//   node test/mesh2/shell.mjs code mode el D surf fix
//     env: NT, NL, TOP (同 phys.mjs), ROWS = 初始层数, EPS, SKYGAP, SKYSTEP, LOG=1 (打印层间的变化), NOSEED=1 (不做沿半径的延拓), NOCRIT=1 (不给临界半径),
//       GRID=1 (延拓也做网格找种子), NOMOVE=1 (不沿半径搬运, 只用种子追踪) —— 对比用; PAGE=1: 与页面一致 (不给临界半径, 点距按平行光弧)
//   读 test/mesh2/cache/ 里 edges.mjs 存下的场景 (层线与边界线).
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
const HERE = path.dirname(fileURLToPath(import.meta.url));
const B = process.env.ENG || new URL('../../src/engine/', import.meta.url).href;
const CR = await import(B + 'crystal.js'), DS = await import(B + 'divshell.js'), MS = await import(process.env.MESH2S || (B + 'mesh2s.js'));
const EP = await import(B + 'epdiv.js'), LPm = await import(B + 'lp.js');
const PH = await import('./phys.mjs');

// 场景 → 球壳法网格. surf: caustic (焦散) 或一条边界的名字 (entry / exit / seam1 …)
export function shellMesh(code, mode, el, Ds, surfName, fixAng, opt = {}) {  // eslint-disable-line
  CR.setApex(56, 56);
  const nr = CR.parsePath(code).map(CR.baseNormalFromCompact), D = +Ds, er = +el * Math.PI / 180, L = [D * Math.cos(er), 0, D * Math.sin(er)], E = [0, 0, 0];
  const wantC = surfName === 'caustic', fam = mode === 'lowitz' ? 'lowitz' : 'column';
  const compute = (r, o) => {
    const res = DS.computeShell2Dof(nr, fam, +fixAng, 1.31, L, E, r, true, wantC, nm => !wantC && nm === surfName,
      { seeds: process.env.NOSEED ? null : o.seeds, noGrid: !process.env.GRID && o.noGrid, from: process.env.NOMOVE ? null : o.from, keepStates: true, skyStep: +(process.env.SKYSTEPE || 2) });
    return { arcs: res.curves.filter(c => (wantC ? c.kind === 0 : c.name === surfName) && c.p.length >= 2).map(c => ({ p: c.p, d: c.dir, closed: c.closed })), cusps: res.cuspP, focal: wantC ? res.focalP : [], states: res.states || null, transported: !!res.transported };
  };
  // 临界半径取自缓存里的边界线 (页面里取自 edges2 / rims2 的结果)
  let critical = [];
  // PAGE=1: 与页面一致 —— 不给临界半径 (页面不等边界线), 点距按平行光弧与预算 60000 定
  if (process.env.PAGE) {
    LPm.LP.partial = true; LPm.setIndex(1.31);
    const lp = CR.parsePath(code).map(CR.lpFaceFromCompact), low = mode === 'lowitz';
    if (!process.env.PAGESKY0) opt = { ...opt, skyFixed: MS.skyStepOf(EP.parallelArc2Dof(lp, +el, 0, 101, low, low ? 0 : +fixAng, low ? +fixAng : 0, EP.resolveSurface(lp, surfName), true), E, 60000, +(process.env.SKYMAX || 0.8)) };
    if (!process.env.PAGECRIT) return MS.meshShells(MS.shellSource(compute), L, E, { maxD: 80, critical: [], ...opt });
  }
  const cname = `${process.env.CACHE || HERE + '/cache'}/${[code, mode, el, Ds, surfName, fixAng].join('_')}.json`;
  if (!process.env.NOCRIT && fs.existsSync(cname)) {
    const dec = s2 => { const b = Buffer.from(s2, 'base64'); return new Float32Array(b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength)); };
    critical = MS.criticalRadii(JSON.parse(fs.readFileSync(cname, 'utf8')).feats.map(f => dec(f.p)), L, E);
  }
  return MS.meshShells(MS.shellSource(compute), L, E, { maxD: 80, critical, ...opt });
}

if (process.argv[1] && path.resolve(process.argv[1]) === fileURLToPath(import.meta.url)) {
  const [code = '4-6', mode = 'lowitz', el = '10.5', Ds = '15', surfName = 'caustic', fixAng = '0'] = process.argv.slice(2);
  const opt = {};
  if (process.env.ROWS) opt.rows = +process.env.ROWS;
  if (process.env.EPS) opt.eps = +process.env.EPS;
  if (process.env.SKYGAP) opt.skyGap = +process.env.SKYGAP;
  if (process.env.SKYSTEP) opt.skyStep = +process.env.SKYSTEP;
  if (process.env.LOG) opt.log = console.log;
  const m = shellMesh(code, mode, el, Ds, surfName, fixAng, opt);
  console.log(`球壳法: ${m.nTri} 三角形, ${m.stat.rows} 层 (引擎算了 ${m.stat.calls} 层: 全搜 ${m.stat.full}, 延拓 ${m.stat.cont}, 其中搬运 ${m.stat.move}; 全搜里: 长度不对重搜 ${m.stat.re}, again ${m.stat.again}, 无邻层 ${m.stat.first}; 弦高查 ${m.stat.sagN} 次 (加层 ${m.stat.sagHit}), 天空间距加层 ${m.stat.gapSplit}, 二分 ${m.stat.bis}), 变化 ${m.stat.events} 处, 收口 ${m.stat.caps}, 没缝上的段 ${m.stat.loose} (${m.stat.looseDeg.toFixed(1)}°), ${m.stat.ms} ms`);
  const r = PH.physCheck2(code, mode, el, Ds, surfName, fixAng, +(process.env.NT || 500), +(process.env.NL || 40), m);
  if (!r) { console.log('no cache: run edges.mjs first'); process.exit(2); }
  PH.printCheck([code, mode, el, Ds, surfName, fixAng].join(' '), r, +(process.env.TOP || 5));
}
