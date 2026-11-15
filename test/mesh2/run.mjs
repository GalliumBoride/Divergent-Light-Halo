// 二自由度曲面网格贴边测试全套: node test/mesh2/run.mjs [并行数=6]
//   逐行跑 suite.txt 里的场景, 两项测试:
//     edges.mjs —— 特征线上相邻两层的路标是否连成了网格的边 (missV / missE / teeth, 目标全 0);
//     open.mjs  —— 网格的开口边是否都落在特征线上 (不在的 = 洞 / 缺口, 按总长计: layer 沿层线的, cross 跨层的).
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
  p.on('close', () => {
    const q = spawn(process.execPath, [path.join(HERE, 'open.mjs'), ...scenes[i].split(' ')], { env: { ...process.env, N: '0' } });
    q.stdout.on('data', d => o += d); q.stderr.on('data', d => o += d);
    q.on('close', () => { outs[i] = o; res(one()); });
  });
});
await Promise.all(Array.from({ length: K }, one));
const tot = { links: 0, missV: 0, missE: 0, teeth: 0, caps: 0, skips: 0, layer: 0, cross: 0 };
console.log('scene'.padEnd(46) + 'links  missV missE teeth  caps skips  open: layer m  cross m   where (missV/missE/teeth)');
outs.forEach((o, i) => {
  const x = { links: 0, missV: 0, missE: 0, teeth: 0, caps: 0, skips: 0, layer: 0, cross: 0 }, by = [];
  for (const l of o.split('\n')) {
    const m = l.match(/^\s+(\w+)\s+marks\s+(\d+) links\s+(\d+) missV\s+(\d+) missE\s+(\d+) teeth\s+(\d+)\/(\d+)\s+caps\s+(\d+) \(([\d.]+) m\) skips (\d+)/);
    if (m) { x.links += +m[3]; x.missV += +m[4]; x.missE += +m[5]; x.teeth += +m[6]; x.caps += +m[8]; x.skips += +m[10]; if (+m[4] + +m[5] + +m[6]) by.push(`${m[1]} ${m[4]}/${m[5]}/${m[6]}`); }
  }
  const mb = o.match(/BAD\s+layer \d+ \(([\d.]+) m\)\s+cross \d+ \(([\d.]+) m\)/);
  if (mb) { x.layer = +mb[1]; x.cross = +mb[2]; }
  for (const k in tot) tot[k] += x[k];
  const ok = /tris \d+/.test(o);
  console.log(scenes[i].padEnd(46) + String(x.links).padStart(5) + String(x.missV).padStart(7) + String(x.missE).padStart(6) + String(x.teeth).padStart(6) + String(x.caps).padStart(6) + String(x.skips).padStart(6) + x.layer.toFixed(1).padStart(14) + x.cross.toFixed(1).padStart(9) + '   ' + (ok ? by.join('  ') : '!! ' + o.split('\n').slice(0, 2).join(' ').slice(0, 120)));
});
console.log(`TOTAL  links ${tot.links}  missV ${tot.missV}  missE ${tot.missE}  teeth ${tot.teeth}  (caps ${tot.caps}, skips ${tot.skips})   open edges off features: layer ${tot.layer.toFixed(0)} m, cross ${tot.cross.toFixed(0)} m`);
process.exit(tot.missV + tot.missE + tot.teeth > 0 ? 1 : 0);
