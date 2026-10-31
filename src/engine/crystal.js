// 晶体面、光路解析、取向法向 — Divergentlightform.cs 常量区 / ParsePath / BuildCanonNormals 的 JS 移植.
// 紧凑面索引 (ParsePath 输出, 法向数组下标): 1-8 → 0-7, 13-18 → 8-13 (上锥), 23-28 → 14-19 (下锥).
// 取向 R(α,ζ,θ) = Rz(α)·Ry(ζ)·Rz(θ); 世界系 Z 朝上, 灯在 XZ 平面: L = (D cosΣ, 0, D sinΣ), 眼 E = (0,0,h).

const S32 = Math.sqrt(3) / 2.0;
export const CRYSTAL_NORMALS = [
  [0, 0, 1],          // 面1: 底面 顶
  [0, 0, -1],         // 面2: 底面 底
  [1, 0, 0],          // 面3: 柱面
  [0.5, S32, 0],      // 面4
  [-0.5, S32, 0],     // 面5
  [-1, 0, 0],         // 面6
  [-0.5, -S32, 0],    // 面7
  [0.5, -S32, 0],     // 面8
];

// 锥晶顶角 (锥面 13-18 上 / 23-28 下), 度 (= 窗体 s_upApexDeg / s_loApexDeg)
export const APEX = { up: 56.0, lo: 56.0 };
export function setApex(up, lo) { APEX.up = up; APEX.lo = lo; }

export function faceNumFromCompact(ci) {
  return ci < 8 ? ci + 1 : ci < 14 ? 13 + (ci - 8) : 23 + (ci - 14);
}

// 晶体系基准法向 (含锥面: 方位同柱面, 抬起 δ=顶角/2, 上锥 +z / 下锥 -z)
export function baseNormalFromCompact(ci) {
  if (ci < 8) return CRYSTAL_NORMALS[ci].slice();
  const p = (ci < 14 ? ci - 8 : ci - 14) + 2;                  // 柱面 3..8 的行号 2..7
  const delta = (ci < 14 ? APEX.up : APEX.lo) * Math.PI / 360.0;
  const cd = Math.cos(delta), sd = Math.sin(delta);
  return [cd * CRYSTAL_NORMALS[p][0], cd * CRYSTAL_NORMALS[p][1], (ci < 14 ? 1.0 : -1.0) * sd];
}

/** 20 个紧凑面的基准法向表 (观测者焦散引擎用). */
export function baseNormalTable() {
  const b = [];
  for (let f = 0; f < 20; f++) b.push(baseNormalFromCompact(f));
  return b;
}

/** 面号 → 紧凑索引; 不认识的面号 → -1. */
export function compactOfFace(f) {
  return (f >= 1 && f <= 8) ? f - 1
       : (f >= 13 && f <= 18) ? 8 + (f - 13)
       : (f >= 23 && f <= 28) ? 14 + (f - 23) : -1;
}

/** "3-1-5" → 紧凑索引数组; 无效 → null (同 ParsePath: 按 '-' 或空格分隔, 至少两面). */
export function parsePath(code) {
  if (code == null) return null;
  const parts = String(code).split(/[-\s]+/).filter(s => s.length > 0);
  const faces = [];
  for (const p of parts) {
    if (!/^[+]?\d+$/.test(p.trim())) return null;
    const ci = compactOfFace(parseInt(p, 10));
    if (ci < 0) return null;
    faces.push(ci);
  }
  return faces.length >= 2 ? faces : null;
}

/** 多条光路: 分号 / 逗号 / 换行分隔 (窗体只认分号); 返回有效的光路代号. */
export function parseMultiPaths(text) {
  if (text == null) return [];
  return String(text).split(/[;,\n]+/).map(s => s.trim()).filter(s => s.length > 0 && parsePath(s) !== null);
}

/** 光路代号规范化: "3 1 5" → "3-1-5". */
export function pathName(code) {
  const f = parsePath(code);
  return f ? f.map(faceNumFromCompact).join('-') : String(code);
}

/** 固定 (ζ,θ) 的规范法向 m_f = Ry(ζ)·Rz(θ)·n_f (20 个紧凑面). 与 h、α 无关. */
export function buildCanonNormals(zetaDeg, thetaDeg) {
  const z = zetaDeg * Math.PI / 180, t = thetaDeg * Math.PI / 180;
  const cz = Math.cos(z), sz = Math.sin(z), ct = Math.cos(t), st = Math.sin(t);
  const m = new Array(20);
  for (let f = 0; f < 20; f++) {
    const b = baseNormalFromCompact(f);
    const x1 = ct * b[0] - st * b[1], y1 = st * b[0] + ct * b[1], z1 = b[2];   // Rz(θ)
    const x2 = cz * x1 + sz * z1, y2 = y1, z2 = -sz * x1 + cz * z1;             // Ry(ζ)
    m[f] = [x2, y2, z2];
  }
  return m;
}

/** 由规范法向施加 Rz(α) 得世界法向 (ca=cosα, sa=sinα). */
export function worldNormalsFromCanon(canon, ca, sa) {
  const wn = new Array(canon.length);
  for (let f = 0; f < canon.length; f++) {
    const x = canon[f][0], y = canon[f][1], z = canon[f][2];
    wn[f] = [ca * x - sa * y, sa * x + ca * y, z];
  }
  return wn;
}

export function getAllWorldNormals(alphaDeg, zetaDeg, thetaDeg) {
  const a = alphaDeg * Math.PI / 180;
  return worldNormalsFromCanon(buildCanonNormals(zetaDeg, thetaDeg), Math.cos(a), Math.sin(a));
}

// 紧凑索引 → HaloLP 面号. 本窗体柱面方位 CCW (面4=+60°), HaloLP 是 CW → 4↔8, 5↔7 镜像;
// 取向映射: al=α, ze=ζ, ρ=θ−90°
const LP_M = [0, 1, 2, 3, 8, 7, 6, 5, 4];
export function lpFaceFromCompact(ci) {
  const f = faceNumFromCompact(ci);
  if (f <= 8) return LP_M[f];
  return (f < 20 ? 10 : 20) + LP_M[f % 10];
}

/** 灯/眼位置 (世界系). */
export function lampEye(lampElevDeg, lampDist, eyeHeight) {
  const e = lampElevDeg * Math.PI / 180.0;
  return { L: [lampDist * Math.cos(e), 0, lampDist * Math.sin(e)], E: [0, 0, eyeHeight] };
}
