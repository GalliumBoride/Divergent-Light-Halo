// 二自由度网格 (球壳法, engine/shellmesh.js —— 页面用的就是它) 的几何校验: 同 phys.mjs 的双向、三维 + 天空两遍.
//   真值取光锥法的层线和精确边界线: 两套引擎互相独立 (球壳轮廓来自 divshell.js, 层线来自 epdiv.js), 所以这也是两套引擎的对账.
//   网格的骨架用缓存里的三维边界线 (页面里取自 edges2 / rims2), 与页面一致; 点距按平行光弧与预算 60000 定 (同页面).
//   node test/mesh2/shell.mjs code mode el D surf fix
//     env: NT, NL, TOP (同 phys.mjs), ROWS = 初始层数, LOG=1 (打印加密与拓扑变化), NOLINES=1 (不给三维边界线, 对比用)
//   读 test/mesh2/cache/ 里 edges.mjs 存下的场景 (层线与边界线).
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
const HERE = path.dirname(fileURLToPath(import.meta.url));
const B = process.env.ENG || new URL('../../src/engine/', import.meta.url).href;
const CR = await import(B + 'crystal.js'), DS = await import(B + 'divshell.js'), MS = await import(B + 'mesh2s.js'), SM = await import(B + 'shellmesh.js');
const EP = await import(B + 'epdiv.js'), LPm = await import(B + 'lp.js');
const PH = await import('./phys.mjs');

const dec = s2 => { const b = Buffer.from(s2, 'base64'); return new Float32Array(b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength)); };

// 场景 → 球壳法网格. surf: caustic (焦散) 或一条边界的名字 (entry / exit / seam1 …)
export function shellMesh(code, mode, el, Ds, surfName, fixAng, opt = {}) {
  CR.setApex(56, 56); LPm.LP.partial = true; LPm.setIndex(1.31);
  const nr = CR.parsePath(code).map(CR.baseNormalFromCompact), D = +Ds, er = +el * Math.PI / 180, L = [D * Math.cos(er), 0, D * Math.sin(er)], E = [0, 0, 0];
  const wantC = surfName === 'caustic', fam = mode === 'lowitz' ? 'lowitz' : 'column';
  const compute = (r, o) => {
    const res = DS.computeShell2Dof(nr, fam, +fixAng, 1.31, L, E, r, true, wantC, nm => !wantC && nm === surfName, { seeds: o.seeds, keepStates: true, skyStep: 2 });
    return { arcs: res.curves.filter(c => (wantC ? c.kind === 0 : c.name === surfName) && c.p.length >= 1).map(c => ({ p: c.p, d: c.dir, closed: c.closed })), focal: wantC ? res.focalP : [], states: res.states || null };
  };
  const cname = `${process.env.CACHE || HERE + '/cache'}/${[code, mode, el, Ds, surfName, fixAng].join('_')}.json`;
  let lines = [];
  if (!process.env.NOLINES && fs.existsSync(cname)) lines = JSON.parse(fs.readFileSync(cname, 'utf8')).feats.map(f => ({ p: dec(f.p), kind: f.g === 'cusp' ? 'cusp' : f.g }));
  const lp = CR.parsePath(code).map(CR.lpFaceFromCompact), low = mode === 'lowitz';
  const sky = MS.skyStepOf(EP.parallelArc2Dof(lp, +el, 0, 101, low, low ? 0 : +fixAng, low ? +fixAng : 0, EP.resolveSurface(lp, surfName), true), E, 60000);
  return SM.meshShellSurface(compute, L, E, { maxD: 80, lines, skyStep: sky, ...opt });
}

if (process.argv[1] && path.resolve(process.argv[1]) === fileURLToPath(import.meta.url)) {
  const [code = '4-6', mode = 'lowitz', el = '10.5', Ds = '15', surfName = 'caustic', fixAng = '0'] = process.argv.slice(2);
  const opt = {};
  if (process.env.ROWS) opt.rows = +process.env.ROWS;
  if (process.env.LOG) opt.log = console.log;
  const m = shellMesh(code, mode, el, Ds, surfName, fixAng, opt);
  const st = m.stat;
  console.log(`球壳法: ${m.nTri} 三角形, ${st.rows} 层 (引擎算了 ${st.calls} 层), 拓扑变化 ${st.events} 处, 收口 ${st.caps}, 没缝上的段 ${st.loose} (${st.looseDeg.toFixed(1)}°), 折叠线 ${st.folds}, 焦线 ${st.focal}, 末端点有标签 ${st.labelled} / 无 ${st.unlabelled}, 边界线插点 ${st.lineVerts}, ${st.ms} ms`);
  const r = PH.physCheck2(code, mode, el, Ds, surfName, fixAng, +(process.env.NT || 500), +(process.env.NL || 40), m);
  if (!r) { console.log('no cache: run edges.mjs first'); process.exit(2); }
  PH.printCheck([code, mode, el, Ds, surfName, fixAng].join(' '), r, +(process.env.TOP || 5));
}
