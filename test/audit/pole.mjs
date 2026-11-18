// 极点线对账 (engine/divshell.js poleLines2Dof): node test/audit/pole.mjs
//   极点线是闭式直接算的三维边界线; 球壳引擎 (computeShell2Dof) 在每层球壳上另外用二分找出边界像的角点末端.
//   两者是互相独立的两条计算路线: 极点线与球壳 r 的交点, 应当就是那一层上 "本边界|另一条边界" 的角点末端.
//   逐场景、逐半径比: 引擎报的这类末端到最近交点的角距 (应 < 0.01°), 交点数 = 这类末端数; 线整条落在 "h = 定值" 上
//   (Lowitz: 过灯的水平面).
import assert from 'node:assert';
const B = process.env.ENG || new URL('../../src/engine/', import.meta.url).href;
const DS = await import(B + 'divshell.js'), CR = await import(B + 'crystal.js');
CR.setApex(56, 56);
const DEG = Math.PI / 180;
// [code, fam, 灯仰角, 灯距, 定角, 曲面, 另一条边界, [半径 …]]
const CASES = [
  ['3-1-5-8', 'lowitz', 10.5, 15, 0, 'deg1', 'entry', [12, 20, 40]],
  ['3-1-5-8', 'lowitz', 87.5, 4.5, 0, 'entry', 'seam1', [15.8, 15.805]],
  ['3-1-5-8', 'lowitz', 34.5, 12.5, 0, 'entry', 'seam1', [26, 30, 40]],
];
let bad = 0;
for (const [code, fam, el, D, fix, surf, other, rs] of CASES) {
  const nr = CR.parsePath(code).map(CR.baseNormalFromCompact), er = el * DEG, L = [D * Math.cos(er), 0, D * Math.sin(er)], E = [0, 0, 0], sf = surf.replace('deg', 'seam');
  const lines = DS.poleLines2Dof(nr, fam, fix, 1.31, L, E, true, 80, sf);
  let zdev = 0, n = 0;
  for (const l of lines) for (let i = 0; i < l.p.length; i += 3) { zdev = Math.max(zdev, Math.abs(l.p[i + 2] - L[2])); n++; }
  let worst = 0, cnt = '';
  for (const r of rs) {
    const res = DS.computeShell2Dof(nr, fam, fix, 1.31, L, E, r, true, false, nm => nm === sf);
    const X = [];
    for (const l of lines) {
      const m = l.p.length / 3, rr = i => Math.hypot(l.p[3 * i], l.p[3 * i + 1], l.p[3 * i + 2]);
      for (let i = 0; i + 1 < m; i++) {
        const a = rr(i) - r, b = rr(i + 1) - r;
        if (a === 0 || a * b < 0) { const t = a / (a - b), q = [0, 1, 2].map(c => l.p[3 * i + c] + t * (l.p[3 * i + 3 + c] - l.p[3 * i + c])), len = Math.hypot(...q); X.push(q.map(x => x / len)); }
      }
    }
    // 这一层上落在 "过灯的水平面" 里的这类角点末端 (同名的角点也可以不在极点上: 那些由光锥法的边界线管)
    const mine = res.fenceEndDir.filter((d, i) => { const k = res.fenceEndKind[i].split('|'); return k.includes(sf) && k.includes(other) && Math.abs(r * d[2] - L[2]) < 1e-5 * r; });
    for (const d of mine) { let bd = Infinity; for (const x of X) bd = Math.min(bd, Math.acos(Math.min(1, d[0] * x[0] + d[1] * x[1] + d[2] * x[2]))); worst = Math.max(worst, bd / DEG); }
    cnt += ` r=${r}: ${mine.length}/${X.length}`;
    if (mine.length !== X.length) worst = Infinity;
  }
  const ok = lines.length > 0 && zdev < 1e-4 && worst < 0.01;
  if (!ok) bad++;
  console.log(`${ok ? '  ' : '✗ '}${code} ${fam} Σ${el} D${D} ${surf}: ${lines.length} 条极点线 ${n} 点, 离 "过灯的水平面" 最远 ${zdev.toExponential(1)} m; 角点末端/交点${cnt}; 最大角距 ${worst.toExponential(2)}°`);
}
// 极点不是角点的光路 / 取向: 不该有极点线. 3-4-5 (3 → 5 隔着 4 反射 = 平行板加一面镜子): 入射掠射与出射临界是同一条曲线,
//   两条边界在极点都取零, 却不是角点 (边界像光滑地穿过那里, 不在那里结束)
for (const [code, fam, fix, surf] of [['3-5-6-8', 'lowitz', 0, 'blue2'], ['3-4-5', 'lowitz', 0, 'entry'], ['3-1-5', 'column', -90, 'entry']]) {
  const nr = CR.parsePath(code).map(CR.baseNormalFromCompact), er = 34.5 * DEG, L = [12.5 * Math.cos(er), 0, 12.5 * Math.sin(er)];
  const n = DS.poleLines2Dof(nr, fam, fix, 1.31, L, [0, 0, 0], true, 80, surf.replace('deg', 'seam')).length;
  if (n) bad++;
  console.log(`${n ? '✗ ' : '  '}${code} ${fam} ${surf}: ${n} 条极点线 (应为 0)`);
}
console.log(bad ? `${bad} 项不合格` : '全部合格');
assert.strictEqual(bad, 0);
