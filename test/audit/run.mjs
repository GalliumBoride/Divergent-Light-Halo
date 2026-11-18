// 批量贴边审计: node test/audit/run.mjs [scenes.txt] [并行数=5] [--strict]
//   每行一个场景: code mode el D surf fix (同 meshaudit.mjs 的参数; # 开头是注释). 默认跑 test/audit/scenes.txt.
//   env: SRC (同 meshaudit.mjs), OUT = 存每个场景完整输出的目录.  --strict: 有不合格的场景时退出码 1
//   合格: 缺口 ≤ 0.5° (天空里的总长), 段间缝 ≤ 0.05°, 边界线没贴上 ≤ 0.5°, 没有失败的段;
//         降自由度的光路没有边界线可比, 看网格器的自检 (开口 ≤ 0.5°)
import { spawn } from 'node:child_process';
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
const HERE = path.dirname(fileURLToPath(import.meta.url));
const argv = process.argv.slice(2).filter(a => a !== '--strict'), strict = process.argv.includes('--strict');
const scenes = fs.readFileSync(argv[0] || path.join(HERE, 'scenes.txt'), 'utf8').split(/\r?\n/).map(s => s.trim()).filter(s => s && !s.startsWith('#'));
const K = +(argv[1] || 5), outs = new Array(scenes.length), times = new Array(scenes.length);
let next = 0;
const one = () => new Promise(res => {
  const i = next++; if (i >= scenes.length) return res();
  const t0 = Date.now();
  const p = spawn(process.execPath, [path.join(HERE, 'meshaudit.mjs'), ...scenes[i].split(/\s+/)], { env: { ...process.env, TOP: process.env.TOP || '3' } });
  let o = ''; p.stdout.on('data', d => o += d); p.stderr.on('data', d => o += d);
  p.on('close', () => { outs[i] = o; times[i] = Date.now() - t0; res(one()); });
});
await Promise.all(Array.from({ length: K }, one));
if (process.env.OUT) { fs.mkdirSync(process.env.OUT, { recursive: true }); scenes.forEach((s, i) => fs.writeFileSync(path.join(process.env.OUT, s.replace(/[^\w.-]+/g, '_') + '.txt'), outs[i])); }
console.log('场景'.padEnd(36) + '三角形'.padStart(8) + ' 网格s 最慢s' + '  缺口条/天空°' + '  段间缝°' + ' 边界线没贴上m/°' + '  丢边' + ' 自检开口°/补洞');
let bad = 0, sumDef = 0, sumSeam = 0, sumBl = 0, sumT = 0, maxT = 0;
scenes.forEach((s, i) => {
  const m = (outs[i] || '').match(/^AUDIT (.*)$/m);
  if (!m) { bad++; console.log(s.padEnd(38) + ' (失败) ' + (outs[i] || '').trim().split('\n').slice(-1)[0].slice(0, 100)); return; }
  const a = JSON.parse(m[1]);
  const ok = (a.rank1 ? (a.selfOpen || 0) <= 0.5 : a.defSky <= 0.5 && a.blMissSky <= 0.5) && a.seamSky <= 0.05 && !a.failed;
  if (!ok) bad++;
  sumDef += a.defSky; sumSeam += a.seamSky; sumBl += a.blMissSky; sumT += a.tSum; maxT = Math.max(maxT, a.tMax);
  console.log((ok ? '  ' : '✗ ') + s.padEnd(36) + String(a.nTri).padStart(8) + (a.tSum / 1000).toFixed(1).padStart(6) + (a.tMax / 1000).toFixed(1).padStart(6) + `${a.defN}/${a.defSky.toFixed(1)}`.padStart(14) + a.seamSky.toFixed(1).padStart(9) + `${a.blMiss.toFixed(2)}/${a.blMissSky.toFixed(1)}`.padStart(15) + String(a.skipped).padStart(7) + `${(a.selfOpen || 0).toFixed(1)}/${a.selfFilled || 0}`.padStart(12) + (a.failed ? '  失败段 ' + a.failed : ''));
});
console.log(`${scenes.length} 个场景: 不合格 ${bad} 个 (缺口 > 0.5° / 段间缝 > 0.05° / 边界线没贴上 > 0.5° / 有段失败); 缺口合计 ${sumDef.toFixed(1)}°, 段间缝 ${sumSeam.toFixed(1)}°, 边界线没贴上 ${sumBl.toFixed(1)}°; 网格用时合计 ${(sumT / 1000).toFixed(1)} s, 最慢一段 ${(maxT / 1000).toFixed(1)} s`);
if (strict && bad) process.exit(1);
