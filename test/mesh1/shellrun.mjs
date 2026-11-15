// 一自由度球壳法全套: node test/mesh1/shellrun.mjs [并行数=4]
//   13 条光路 × 14 组灯位 / 取向 (同 phys.mjs 用过的那一批), 每个场景跑 test/mesh1/shell.mjs:
//   耗时、网格→晕面 (物理校验的偏差面积), 晕面→网格 ((h, α) 网格顶点没盖到的比例). 环境变量 LIST="code el D ζ θ;…" 换场景.
import { spawn } from 'node:child_process';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
const HERE = path.dirname(fileURLToPath(import.meta.url));
let scenes = (process.env.LIST || '').split(';').map(s => s.trim()).filter(Boolean);
if (!scenes.length)
  for (const code of ['3-5', '1-3', '3-1-5', '3-7', '1-4', '3-2', '4-1-8', '3-8-5', '3-6', '1-3-2', '13-5', '3-25', '3-1-7'])
    for (const cfg of ['10.5 15 -138 60', '20 15 0 0', '5 15 90 0', '30 15 90 30', '45 10 60 20', '2 15 0 0', '60 15 -40 75', '15 15 0 30', '25 15 90 15', '8 15 35 10', '40 15 120 45', '12 30 -100 5', '70 15 20 50', '1 15 90 0'])
      scenes.push(code + ' ' + cfg);
const K = +(process.argv[2] || 4), outs = new Array(scenes.length), ms = new Array(scenes.length);
let next = 0;
const one = () => new Promise(res => {
  const i = next++; if (i >= scenes.length) return res();
  const t0 = Date.now(), p = spawn(process.execPath, [path.join(HERE, 'shell.mjs'), ...scenes[i].split(' ')], { env: process.env });
  let o = ''; p.stdout.on('data', d => o += d); p.stderr.on('data', d => o += d);
  p.on('close', () => { outs[i] = o; ms[i] = Date.now() - t0; res(one()); });
});
await Promise.all(Array.from({ length: K }, one));
let n = 0, clean = 0, bad10 = 0, miss = 0, pts = 0, tmax = 0, ttot = 0;
console.log('scene'.padEnd(26) + ' 三角形  面积m²  球壳ms  假面>10cm m²   漏盖');
outs.forEach((o, i) => {
  const a = o.match(/球壳法 (\d+) 三角形 ([\d.]+) m².*?, (\d+) ms/), f = o.match(/> 10 cm \+ 1% r 的 ([\d.]+) m²/), c = o.match(/顶点抽 (\d+) 个, 离球壳网格 > 5 cm \+ 1% r 的 (\d+) 个/);
  if (!a || !f || !c) { console.log(scenes[i].padEnd(26) + ' (failed) ' + o.trim().split('\n').slice(-1)[0].slice(0, 100)); return; }
  n++; bad10 += +f[1]; miss += +c[2]; pts += +c[1]; ttot += +a[3]; tmax = Math.max(tmax, +a[3]);
  const ok = +f[1] === 0 && +c[2] <= 0.002 * +c[1]; if (ok) clean++;
  if (!ok || process.env.ALL) console.log(scenes[i].padEnd(26) + a[1].padStart(7) + (+a[2]).toFixed(1).padStart(8) + a[3].padStart(8) + (+f[1]).toFixed(3).padStart(13) + `${c[2]}/${c[1]}`.padStart(11));
});
console.log(`${n} 个场景跑完: 全对 (假面 0、漏盖 ≤ 0.2%) ${clean} 个; 假面合计 ${bad10.toFixed(3)} m², 漏盖 ${miss}/${pts}; 球壳法耗时合计 ${(ttot / 1000).toFixed(0)} s, 最慢 ${(tmax / 1000).toFixed(1)} s`);
