// 光锥法 — Divergentlightform.cs 的 ComputeHaloSurfaceLightCone (h 扫) / ComputeHaloSurfaceLightConeSwap (ψ 扫) /
// ComputeParallelArc / ComputeDivergentOnCone 的 JS 移植, 运算次序与 C# 一致 (逐点对账见 tools/).
//
// 【问题】灯 L 与眼 E 为两焦点; 光自 L 出发在 P 进晶体, 按光路折射/反射后到 E. 固定 (ζ,θ), 扫 α, 求所有 P (三维"晕面").
// 【锥上消元】入射仰角 h、方位 ψ: û = Rz(ψ)·(−cos h, 0, −sin h); 偏折绕天顶对称 ⟹ 出射 d_div = Rz(ψ)·d(α, h).
// 【共面】E − L = t·û + s·d_div (t,s>0) ⟺ F(ψ) = Pc sinψ + Qc cosψ − Rc = 0, 闭式两支;
// 【定位】c = û·d_div, t = (c(d·𝐃) − û·𝐃)/(1−c²), s = (c(û·𝐃) − d·𝐃)/(1−c²), P = L + t û.
// h 扫在 ε₀→0 时 h 窄带会漏采, 由 ψ 扫 (每个 (α,ψ) 二分求 h) 补齐, 两者取并集.
//
// 切片可以分给多个 Worker: 每个函数只算 [from, to) 区间, 结果按切片序号拼接即与 C# 完全相同.

import { buildCanonNormals, worldNormalsFromCanon, getAllWorldNormals } from './crystal.js';

const f32 = Math.fround;

// 平行光 Snell 链: 入射传播方向 (dix,diy,diz) → 归一化出射方向, 无效 → null
function snell(wn, facePath, dix, diy, diz, eta, n2, refN) {
  const last = facePath.length - 1;
  const nIn = wn[facePath[0]], nOut = wn[facePath[last]];
  const cosI = -(dix * nIn[0] + diy * nIn[1] + diz * nIn[2]);
  if (cosI < 1e-10) return null;
  const s2 = eta * eta * (1 - cosI * cosI);
  if (s2 > 1) return null;
  const cosT = Math.sqrt(1 - s2);
  const c1 = eta * cosI - cosT;
  let ix = eta * dix + c1 * nIn[0];
  let iy = eta * diy + c1 * nIn[1];
  let iz = eta * diz + c1 * nIn[2];
  for (let k = 1; k < last; k++) {
    const nR = wn[facePath[k]];
    const dn = ix * nR[0] + iy * nR[1] + iz * nR[2];
    if (dn <= 0) return null;
    ix -= 2 * dn * nR[0]; iy -= 2 * dn * nR[1]; iz -= 2 * dn * nR[2];
  }
  const cosI2 = ix * nOut[0] + iy * nOut[1] + iz * nOut[2];
  if (cosI2 < 1e-10) return null;
  const s2_2 = n2 * (1 - cosI2 * cosI2);
  if (s2_2 > 1) return null;
  const cosT2 = Math.sqrt(1 - s2_2);
  const c2 = cosT2 - refN * cosI2;
  const dx = refN * ix + c2 * nOut[0];
  const dy = refN * iy + c2 * nOut[1];
  const dz = refN * iz + c2 * nOut[2];
  const dn2 = Math.sqrt(dx * dx + dy * dy + dz * dz);
  if (dn2 < 1e-10) return null;
  return [dx / dn2, dy / dn2, dz / dn2];
}

/** TryConePoint: 成功返回 [Px,Py,Pz], 否则 null. */
export function tryConePoint(psi, cosH, sinH, dx, dy, dz, Dx, Dy, Dz, L, maxDist) {
  const cp = Math.cos(psi), sp = Math.sin(psi);
  const ux = -cosH * cp, uy = -cosH * sp, uz = -sinH;      // û
  const dxDiv = dx * cp - dy * sp;                          // d_div = Rz(ψ)·d
  const dyDiv = dx * sp + dy * cp;
  const dzDiv = dz;
  const cVal = ux * dxDiv + uy * dyDiv + uz * dzDiv;
  const det = 1 - cVal * cVal;
  if (Math.abs(det) < 1e-10) return null;
  const uD = ux * Dx + uy * Dy + uz * Dz;
  const dD = dxDiv * Dx + dyDiv * Dy + dzDiv * Dz;
  const t = (cVal * dD - uD) / det;
  const s = (cVal * uD - dD) / det;
  if (t < 1e-6 || s < 1e-6) return null;
  const Px = L[0] + t * ux, Py = L[1] + t * uy, Pz = L[2] + t * uz;
  const rx = Px - L[0], ry = Py - L[1], rz = Pz - L[2];
  if (rx * rx + ry * ry + rz * rz > maxDist * maxDist) return null;
  return [Px, Py, Pz];
}

function addConePoint(psi, cosH, sinH, dx, dy, dz, Dx, Dy, Dz, L, maxDist, verts) {
  const P = tryConePoint(psi, cosH, sinH, dx, dy, dz, Dx, Dy, Dz, L, maxDist);
  if (P) verts.push(f32(P[0]), f32(P[1]), f32(P[2]));
}

/** 场景常量 (两种扫描共用). */
export function coneScene(lampElev, lampDist, eyeHeight, refN) {
  const elevRad = lampElev * Math.PI / 180.0;
  const Lx = lampDist * Math.cos(elevRad), Lz = lampDist * Math.sin(elevRad);
  return {
    L: [Lx, 0, Lz], E: [0, 0, eyeHeight],
    Dx: Lx, Dy: 0, Dz: Lz - eyeHeight,
    eta: 1.0 / refN, n2: refN * refN, refN,
  };
}

/**
 * h 扫: 切片 ih ∈ [from, to) (共 nH 个, h = −90..90), 每个切片扫 nAlpha 个 α, 闭式 ψ 两支.
 * 返回 slices[ih−from] = Float32Array (x,y,z 连排) 或 null.
 */
export function lightConeH(facePath, lampElev, lampDist, eyeHeight, refN, nH, nAlpha, maxDist,
  zetaDeg, thetaDeg, from = 0, to = nH) {
  const sc = coneScene(lampElev, lampDist, eyeHeight, refN);
  const { L, Dx, Dy, Dz, eta, n2 } = sc;
  const last = facePath.length - 1;
  const PSI_FULL_SCAN = 361;
  const canon = buildCanonNormals(zetaDeg, thetaDeg);
  const wnByAlpha = new Array(nAlpha);
  for (let iA = 0; iA < nAlpha; iA++) {
    const aRad = 2 * Math.PI * iA / nAlpha;
    wnByAlpha[iA] = worldNormalsFromCanon(canon, Math.cos(aRad), Math.sin(aRad));
  }
  const out = [];
  for (let ih = from; ih < to; ih++) {
    const hDeg = -90.0 + 180.0 * ih / Math.max(1, nH - 1);
    const hRad = hDeg * Math.PI / 180.0;
    const cosH = Math.cos(hRad), sinH = Math.sin(hRad);
    const dixP = -cosH, diyP = 0, dizP = -sinH;
    const sliceVerts = [];
    for (let iA = 0; iA < nAlpha; iA++) {
      const wn = wnByAlpha[iA];
      const nIn = wn[facePath[0]], nOut = wn[facePath[last]];
      const cosI = -(dixP * nIn[0] + diyP * nIn[1] + dizP * nIn[2]);
      if (cosI < 1e-10) continue;
      const s2 = eta * eta * (1 - cosI * cosI);
      if (s2 > 1) continue;
      const cosT = Math.sqrt(1 - s2);
      const c1 = eta * cosI - cosT;
      let ix = eta * dixP + c1 * nIn[0];
      let iy = eta * diyP + c1 * nIn[1];
      let iz = eta * dizP + c1 * nIn[2];
      let okSnell = true;
      for (let k = 1; k < last; k++) {
        const nR = wn[facePath[k]];
        const dn = ix * nR[0] + iy * nR[1] + iz * nR[2];
        if (dn <= 0) { okSnell = false; break; }
        ix -= 2 * dn * nR[0]; iy -= 2 * dn * nR[1]; iz -= 2 * dn * nR[2];
      }
      if (!okSnell) continue;
      const cosI2 = ix * nOut[0] + iy * nOut[1] + iz * nOut[2];
      if (cosI2 < 1e-10) continue;
      const s2_2 = n2 * (1 - cosI2 * cosI2);
      if (s2_2 > 1) continue;
      const cosT2 = Math.sqrt(1 - s2_2);
      const c2 = cosT2 - refN * cosI2;
      let dx = refN * ix + c2 * nOut[0];
      let dy = refN * iy + c2 * nOut[1];
      let dz = refN * iz + c2 * nOut[2];
      const dn2 = Math.sqrt(dx * dx + dy * dy + dz * dz);
      if (dn2 < 1e-10) continue;
      dx /= dn2; dy /= dn2; dz /= dn2;

      // 共面 F(ψ)=Pc·sinψ+Qc·cosψ−Rc=0 → 闭式两支 (cos h 缩放形式, 全 h 无奇点)
      const Pc = sinH * (Dx * dx + Dy * dy) - cosH * Dx * dz;
      const Qc = sinH * (Dx * dy - Dy * dx) + cosH * Dy * dz;
      const Rc = cosH * Dz * dy;
      const mag = Math.sqrt(Pc * Pc + Qc * Qc);
      if (mag < 1e-12) {
        if (Math.abs(Rc) >= 1e-12) continue;                 // 无解
        for (let i = 0; i < PSI_FULL_SCAN; i++) {             // 2D 塌缩: ψ 连续族, 全扫
          const psiF = -Math.PI + 2 * Math.PI * i / (PSI_FULL_SCAN - 1);
          addConePoint(psiF, cosH, sinH, dx, dy, dz, Dx, Dy, Dz, L, maxDist, sliceVerts);
        }
        continue;
      }
      let rhs = Rc / mag;
      if (rhs > 1.0 + 1e-9 || rhs < -1.0 - 1e-9) continue;
      if (rhs > 1.0) rhs = 1.0; else if (rhs < -1.0) rhs = -1.0;
      const phi = Math.atan2(Qc, Pc);
      const asr = Math.asin(rhs);
      addConePoint(asr - phi, cosH, sinH, dx, dy, dz, Dx, Dy, Dz, L, maxDist, sliceVerts);
      if (rhs < 1.0 - 1e-9 && rhs > -1.0 + 1e-9)
        addConePoint(Math.PI - asr - phi, cosH, sinH, dx, dy, dz, Dx, Dy, Dz, L, maxDist, sliceVerts);
    }
    out.push(sliceVerts.length > 0 ? Float32Array.from(sliceVerts) : null);
  }
  return out;
}

/**
 * ψ 扫 (swap): 切片 iA ∈ [from, to) (共 nAlpha 个), 每个 α 扫 nPsi 个 ψ, F(h) 粗扫 (1°) 零穿越 → 二分求 h.
 * ε₀ = 0 (灯与眼同高且灯仰角 0) 时 F ≡ 0, 直接对 h = 0 解 (t,s).
 */
export function lightConeSwap(facePath, lampElev, lampDist, eyeHeight, refN, nAlpha, nPsi, maxDist,
  zetaDeg, thetaDeg, from = 0, to = nAlpha) {
  const sc = coneScene(lampElev, lampDist, eyeHeight, refN);
  const { L, Dx, Dy, Dz, eta, n2 } = sc;
  const Lz = L[2];
  const isEps0 = Math.abs(lampElev) < 0.01 && Math.abs(eyeHeight - Lz) < 0.001;
  const H_COARSE_N = 181, H_BISECT = 25;
  const hCoarse = new Float64Array(H_COARSE_N);
  for (let i = 0; i < H_COARSE_N; i++) hCoarse[i] = -90.0 + 180.0 * i / (H_COARSE_N - 1);
  const canon = buildCanonNormals(zetaDeg, thetaDeg);
  const out = [];
  for (let iA = from; iA < to; iA++) {
    const aRad = 2 * Math.PI * iA / nAlpha;
    const wn = worldNormalsFromCanon(canon, Math.cos(aRad), Math.sin(aRad));
    const snellAt = hDeg => {
      const hRad = hDeg * Math.PI / 180.0;
      const cH = Math.cos(hRad), sH = Math.sin(hRad);
      const d = snell(wn, facePath, -cH, 0, -sH, eta, n2, refN);
      return d === null ? null : [d[0], d[1], d[2], cH, sH];
    };
    const computeF = (snl, psi) => {
      if (snl === null) return NaN;
      const dx = snl[0], dy = snl[1], dz = snl[2], cH = snl[3], sH = snl[4];
      const Pc = sH * (Dx * dx + Dy * dy) - cH * Dx * dz;
      const Qc = sH * (Dx * dy - Dy * dx) + cH * Dy * dz;
      const Rc = cH * Dz * dy;
      return Pc * Math.sin(psi) + Qc * Math.cos(psi) - Rc;
    };
    const alphaVerts = [];
    const snellCache = new Array(H_COARSE_N);
    for (let ih = 0; ih < H_COARSE_N; ih++) snellCache[ih] = snellAt(hCoarse[ih]);
    for (let iP = 0; iP < nPsi; iP++) {
      const psi = -Math.PI + 2 * Math.PI * iP / nPsi;
      if (isEps0) {
        const snl0 = snellAt(0.0);
        if (snl0 !== null)
          addConePoint(psi, snl0[3], snl0[4], snl0[0], snl0[1], snl0[2], Dx, Dy, Dz, L, maxDist, alphaVerts);
        continue;
      }
      let prevF = computeF(snellCache[0], psi);
      for (let ih = 1; ih < H_COARSE_N; ih++) {
        const curF = computeF(snellCache[ih], psi);
        if (!Number.isNaN(prevF) && !Number.isNaN(curF) && prevF * curF < 0) {
          let lo = hCoarse[ih - 1], hi = hCoarse[ih], flo = prevF;
          for (let bi = 0; bi < H_BISECT; bi++) {
            const mid = 0.5 * (lo + hi);
            const Fm = computeF(snellAt(mid), psi);
            if (Number.isNaN(Fm)) break;
            if (flo * Fm < 0) hi = mid;
            else { lo = mid; flo = Fm; }
          }
          const snl = snellAt(0.5 * (lo + hi));
          if (snl !== null)
            addConePoint(psi, snl[3], snl[4], snl[0], snl[1], snl[2], Dx, Dy, Dz, L, maxDist, alphaVerts);
        }
        prevF = curF;
      }
    }
    out.push(alphaVerts.length > 0 ? Float32Array.from(alphaVerts) : null);
  }
  return out;
}

/** 平行光弧: 把灯当太阳 (仰角 Σ) 的平行光晕, 投到远天球 (眼位 − 200 m·d_out). α 步长 0.1°. */
export function parallelArc(faces, lampElevDeg, zetaDeg, thetaDeg, refN, E) {
  const result = [];
  const elevRad = lampElevDeg * Math.PI / 180.0;
  const dix = -Math.cos(elevRad), diy = 0, diz = -Math.sin(elevRad);
  const eta = 1.0 / refN, n2 = refN * refN, skyR = 200.0;
  let seg = [];
  const flush = () => { if (seg.length >= 6) result.push(Float32Array.from(seg)); seg = []; };
  const canon = buildCanonNormals(zetaDeg, thetaDeg);
  for (let ia = 0; ia <= 3600; ia++) {
    const a = ia * 0.1 * Math.PI / 180;
    const wn = worldNormalsFromCanon(canon, Math.cos(a), Math.sin(a));
    const d = snell(wn, faces, dix, diy, diz, eta, n2, refN);
    if (d === null) { flush(); continue; }
    seg.push(f32(E[0] - skyR * d[0]), f32(E[1] - skyR * d[1]), f32(E[2] - skyR * d[2]));
  }
  flush();
  return result;
}

/**
 * 锥切片: 灯处仰角 h 的入射锥上, 对每个 α' 解析求 ψ, 得晶体位置与观测者视方向 (投到 200 m 天球).
 * 返回 { crystal: Float32Array | null, view: Float32Array | null } (各一条点列, 同 C# 单段输出).
 */
export function divergentOnCone(faces, coneElevDeg, lampElevDeg, lampDist, eyeHeight, zetaDeg, thetaDeg, refN) {
  const lampElevRad = lampElevDeg * Math.PI / 180.0;
  const Lx = lampDist * Math.cos(lampElevRad), Lz = lampDist * Math.sin(lampElevRad);
  const L = [Lx, 0, Lz], O = [0, 0, eyeHeight];
  const Dx = L[0] - O[0], Dy = L[1] - O[1], Dz = L[2] - O[2];
  const hRad = coneElevDeg * Math.PI / 180.0;
  const cosH = Math.cos(hRad), sinH = Math.sin(hRad);
  const res = { crystal: null, view: null, n: 0 };
  if (Math.abs(sinH) < 1e-10) return res;                   // h≈0 退化
  const k = cosH / sinH;                                   // cot h
  const dix = -Math.cos(hRad), diy = 0, diz = -Math.sin(hRad);
  const eta = 1.0 / refN, n2 = refN * refN;
  const crystalPts = [], viewPts = [];
  const canon = buildCanonNormals(zetaDeg, thetaDeg);
  for (let ia = 0; ia < 3601; ia++) {
    const a = ia * 0.1 * Math.PI / 180;
    const wn = worldNormalsFromCanon(canon, Math.cos(a), Math.sin(a));
    const d = snell(wn, faces, dix, diy, diz, eta, n2, refN);
    if (d === null) continue;
    const dx = d[0], dy = d[1], dz = d[2];
    const A = dx - k * dz, B = dy;
    const P = Dx * A + Dy * B, Q = Dx * B - Dy * A, R = k * Dz * B;
    const mag = Math.sqrt(P * P + Q * Q);
    if (mag < 1e-15) continue;
    const ratio = R / mag;
    if (Math.abs(ratio) > 1.0) continue;
    const baseAngle = Math.atan2(P, Q);
    const delta = Math.acos(Math.max(-1, Math.min(1, ratio)));
    const sols = [baseAngle + delta, baseAngle - delta];
    for (let iSol = 0; iSol < 2; iSol++) {
      const psi = sols[iSol];
      const phi = psi + Math.PI;
      const cosPsi = Math.cos(psi), sinPsi = Math.sin(psi);
      const lhs_x = A * cosPsi - B * sinPsi, lhs_y = A * sinPsi + B * cosPsi;
      const rhs_x = k * Dz * cosPsi - Dx, rhs_y = k * Dz * sinPsi - Dy;
      let s;
      if (Math.abs(lhs_x) > Math.abs(lhs_y)) s = rhs_x / lhs_x;
      else if (Math.abs(lhs_y) > 1e-15) s = rhs_y / lhs_y;
      else continue;
      const t = (s * dz + Dz) / sinH;
      if (t < 0.01 || s < 0.01) continue;
      const cosPhi = Math.cos(phi), sinPhi = Math.sin(phi);
      const Px = L[0] + t * cosH * cosPhi, Py = L[1] + t * cosH * sinPhi, Pz = L[2] + t * (-sinH);
      const dist2 = (Px - O[0]) * (Px - O[0]) + (Py - O[1]) * (Py - O[1]) + (Pz - O[2]) * (Pz - O[2]);
      if (dist2 > 40000) continue;
      crystalPts.push(f32(Px), f32(Py), f32(Pz));
      const vx = Px - O[0], vy = Py - O[1], vz = Pz - O[2];
      const vn = Math.sqrt(vx * vx + vy * vy + vz * vz);
      if (vn > 1e-10) {
        const skyR = 200.0;
        viewPts.push(f32(O[0] + skyR * vx / vn), f32(O[1] + skyR * vy / vn), f32(O[2] + skyR * vz / vn));
      }
    }
  }
  if (crystalPts.length >= 3) res.crystal = Float32Array.from(crystalPts);
  if (viewPts.length >= 3) res.view = Float32Array.from(viewPts);
  res.n = crystalPts.length / 3;
  return res;
}

export { getAllWorldNormals };
