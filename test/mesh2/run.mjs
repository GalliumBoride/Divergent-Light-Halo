// 二自由度曲面网格贴边测试全套: node test/mesh2/run.mjs [并行数=6]
//   逐行跑 suite.txt 里的场景 (edges.mjs), 汇总 "特征线上相邻两层的路标是否连成了网格的边" —— 目标: missV / missE / teeth 全 0.
//   第一次要算每个场景的层和精确边界线 (几分钟到十几分钟), 之后读缓存, 一分钟左右.
import { spawn } from 'node:child_process';
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
const HERE = path.dirname(fileURLToPath(import.meta.url));
const scenes = fs.readFileSync(path.join(HERE, 'suite.txt'), 'utf8').split('\n').map(s => s.trim()).filter(Boolean);
const K = +(process.argv[2] || 6), outs = new Array(scenes.length);
let next = 0;
const one = () => new Promise(res => {
  const i = next++; if (i >= scenes.length) return res();
  const p = spawn(process.execPath, [path.join(HERE, 'edges.mjs'), ...scenes[i].split(' ')], { env: process.env });
  let o = ''; p.stdout.on('data', d => o += d); p.stderr.on('data', d => o += d);
  p.on('close', () => { outs[i] = o; res(one()); });
});
await Promise.all(Array.from({ length: K }, one));
const tot = { links: 0, missV: 0, missE: 0, teeth: 0, caps: 0, skips: 0 };
console.log('scene'.padEnd(46) + 'links  missV missE teeth  caps skips   where (missV/missE/teeth)');
outs.forEach((o, i) => {
  const x = { links: 0, missV: 0, missE: 0, teeth: 0, caps: 0, skips: 0 }, by = [];
  for (const l of o.split('\n')) {
    const m = l.match(/^\s+(\w+)\s+marks\s+(\d+) links\s+(\d+) missV\s+(\d+) missE\s+(\d+) teeth\s+(\d+)\/(\d+)\s+caps\s+(\d+) \(([\d.]+) m\) skips (\d+)/);
    if (m) { x.links += +m[3]; x.missV += +m[4]; x.missE += +m[5]; x.teeth += +m[6]; x.caps += +m[8]; x.skips += +m[10]; if (+m[4] + +m[5] + +m[6]) by.push(`${m[1]} ${m[4]}/${m[5]}/${m[6]}`); }
  }
  for (const k in tot) tot[k] += x[k];
  const ok = /tris \d+/.test(o);
  console.log(scenes[i].padEnd(46) + String(x.links).padStart(5) + String(x.missV).padStart(7) + String(x.missE).padStart(6) + String(x.teeth).padStart(6) + String(x.caps).padStart(6) + String(x.skips).padStart(6) + '   ' + (ok ? by.join('  ') : '!! ' + o.split('\n').slice(0, 2).join(' ').slice(0, 120)));
});
console.log(`TOTAL  links ${tot.links}  missV ${tot.missV}  missE ${tot.missE}  teeth ${tot.teeth}  (caps ${tot.caps}, skips ${tot.skips}: 没有纳入约束的两类, 见 README)`);
process.exit(tot.missV + tot.missE + tot.teeth > 0 ? 1 : 0);
