// 常见晕弧在路灯光下的样子: 光路 + 取向 + 灯仰角/灯距 + 天空视角 [方位, 仰角, 视场]. 名字与说明在 strings.js (pre.<id>).
// 面号是本工程 (发散光窗体) 的逆时针约定; 不对称光路 (3-1-5 等) 打开「左右镜像补全」即得另一半.

export const PRESETS = [
  { id: 'sundog', path: '3-5', mode: 'plate', elev: 10, dist: 15, sky: [0, 10, 160] },
  { id: 'cza', path: '1-3', mode: 'plate', elev: 25, dist: 15, sky: [0, 60, 180] },
  { id: 'cha', path: '3-2', mode: 'plate', elev: 60, dist: 15, sky: [0, 20, 180] },
  { id: 'parry', path: '3-5', mode: 'parry', elev: 15, dist: 15, sky: [0, 25, 140] },
  { id: 'tangent', path: '3-5', mode: 'column', elev: 15, dist: 20, sky: [0, 15, 160] },
  { id: 'circumscribed', path: '3-5', mode: 'column', elev: 40, dist: 20, sky: [0, 40, 160] },
  { id: 'supralateral', path: '3-1', mode: 'column', elev: 20, dist: 20, sky: [0, 45, 180] },
  { id: 'infralateral', path: '1-3', mode: 'column', elev: 20, dist: 20, sky: [0, 15, 180] },
  { id: 'wegener', path: '3-1-5', mode: 'column', elev: 10, dist: 20, sky: [180, 15, 180] },
  { id: 'lowitz', path: '3-5', mode: 'lowitz', elev: 20, dist: 15, sky: [0, 20, 140] },
  { id: 'tape', path: '4-1', mode: 'parry', elev: 15, dist: 15, sky: [0, 30, 180] },
  { id: 'hastings', path: '3-1-5', mode: 'parry', elev: 10, dist: 15, sky: [0, 20, 180] },
];
