// 二自由度网格 (球壳法) 几何校验全套: node test/mesh2/shellrun.mjs [并行数=4]
//   逐行跑 suite.txt 里的场景 (test/mesh2/shell.mjs): 球壳法建网格 → 双向校验 (同 phys.mjs). 每个场景一行:
//     假面 (网格 → 晕面): 抽查的三角形里重心离真实层线 > 10 cm + 1% r 的面积 (m²), 天空里偏差 > 1° 的面积 (平方度);
//     漏 (晕面 → 网格): 另算的光锥层线、精确边界线 / 脊线里离网格 > 10 cm + 1% r 的长度 (m), 天空里离网格 > 1° 的长度 (度).
import { spawn } from 'node:child_process';
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
const HERE = path.dirname(fileURLToPath(import.meta.url));
const scenes = fs.readFileSync(path.join(HERE, 'suite.txt'), 'utf8').split('\n').map(s => s.trim()).filter(Boolean);
const K = +(process.argv[2] || 4), outs = new Array(scenes.length);
let next = 0;
const one = () => new Promise(res => {
  const i = next++; if (i >= scenes.length) return res();
  const p = spawn(process.execPath, [path.join(HERE, 'shell.mjs'), ...scenes[i].split(' ')], { env: { ...process.env, TOP: '1', NT: process.env.NT || '300', NL: process.env.NL || '24' } });
  let o = ''; p.stdout.on('data', d => o += d); p.stderr.on('data', d => o += d);
  p.on('close', () => { outs[i] = o; res(one()); });
});
await Promise.all(Array.from({ length: K }, one));
console.log('scene'.padEnd(38) + ' 三角形   层   秒   假面m² 假面°²  层线漏m   漏°  边界漏m   漏°');
const tot = [0, 0, 0, 0, 0, 0];
let clean = 0;
outs.forEach((o, i) => {
  const a = o.match(/球壳法: (\d+) 三角形, (\d+) 层.*?, (\d+) ms/);
  const f = o.match(/三维: 查.*?>10cm ([\d.]+) m²/), s = o.match(/天空: 查.*?>1° ([\d.]+) 平方度/);
  const L = [...o.matchAll(/晕面→网格 .*?: 三维 [\d.]+ m 里漏 >3cm [\d.]+ m, >10cm ([\d.]+) m; 天空 \d+° 里漏 >0.3° [\d.]+°, >1° ([\d.]+)°/g)];
  if (!a || !f || !s || L.length < 2) { console.log(scenes[i].padEnd(38) + ' (failed) ' + o.trim().split('\n').slice(-1)[0].slice(0, 90)); return; }
  const v = [+f[1], +s[1], +L[0][1], +L[0][2], +L[1][1], +L[1][2]];
  v.forEach((x, k) => { tot[k] += x; });
  if (v.every(x => x === 0)) clean++;
  console.log(scenes[i].padEnd(38) + a[1].padStart(7) + a[2].padStart(5) + (a[3] / 1000).toFixed(0).padStart(5) + v[0].toFixed(2).padStart(9) + v[1].toFixed(1).padStart(7) + v[2].toFixed(2).padStart(9) + v[3].toFixed(1).padStart(6) + v[4].toFixed(2).padStart(9) + v[5].toFixed(1).padStart(6));
});
console.log('TOTAL'.padEnd(55) + tot[0].toFixed(2).padStart(9) + tot[1].toFixed(1).padStart(7) + tot[2].toFixed(2).padStart(9) + tot[3].toFixed(1).padStart(6) + tot[4].toFixed(2).padStart(9) + tot[5].toFixed(1).padStart(6));
console.log(`${clean} / ${scenes.length} 个场景六项全为 0`);
