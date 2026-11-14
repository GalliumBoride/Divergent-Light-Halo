// 观测者天空视图: 站在眼睛处看出去. 晕面上每个晶体位置 P 映成方向 −d (出射光传播方向的反向 — 精确, 不靠 P−E 相减),
// 再按镜头投影到画面. 两层画布: 底层 WebGL 画晕面网格与采样点 (GPU 逐顶点投影, 拖动不卡);
// 上层 2D 画布画网格、线、标记、文字 (透明底). 镜头 (SkyProj) 与 LightPointTrace / 原窗体一致:
//   方位 az 为配置方位角 (俯视顺时针增), az=0 朝灯; 画面右手 = 面朝视线的右手.

const DEG = Math.PI / 180;
const FONT = '"PingFang SC", "Microsoft YaHei UI", "Microsoft YaHei", "Noto Sans SC", system-ui, sans-serif';

// ============================================================ 天球镜头
// lens: 0 直线 (透视) / 1 鱼眼等积 / 2 鱼眼等距 / 3 鱼眼体视 / 4 鱼眼正交 / 5 等距圆柱 (全景)
export const LENS_COUNT = 6;
export function fovMax(lens) {
  switch (lens) { case 0: return 170; case 3: return 355; case 4: return 180; default: return 360; }
}

export class SkyProj {
  // frame = { w, h }: 投影按这块居中的画框取景 (载入底图时 = 照片按长宽比等比放进视图), 鱼眼像圈内切画框、不留边;
  //   不给就按整个视图 (鱼眼留 18 px 边)
  set(az, el, fov, lens, W, H, frame = null) {
    const FW = frame ? frame.w : W, FH = frame ? frame.h : H;
    this.lens = lens;
    const a = -az * DEG, e = el * DEG;
    const ce = Math.cos(e), se = Math.sin(e), ca = Math.cos(a), sa = Math.sin(a);
    const ax = [ce * ca, ce * sa, se];
    const up = [-se * ca, -se * sa, ce];
    this.ax = ax; this.up = up;
    this.rt = [ax[1] * up[2] - ax[2] * up[1], ax[2] * up[0] - ax[0] * up[2], ax[0] * up[1] - ax[1] * up[0]];
    this.hf = Math.max(0.5, fov) * Math.PI / 360.0;
    this.cosHf = Math.cos(Math.min(this.hf, Math.PI));
    this.azRad = -az * DEG;
    const hfLin = Math.min(this.hf, 85.0 * DEG);
    this.linFocal = 0.5 * Math.sqrt(FW * FW + FH * FH) / Math.tan(hfLin);
    this.eqScale = Math.min(FW / 2.0, FH) / Math.PI;
    this.cx = W / 2; this.cy = H / 2; this.W = W; this.H = H;
    this.rpx = Math.min(FW, FH) / 2 - (frame ? 0 : 18);
  }

  /** 方向 → 屏幕 (this.x, this.y); 不可见返回 false */
  project(Dx, Dy, Dz) {
    const cx = this.cx, cy = this.cy, rpx = this.rpx, lens = this.lens;
    if (lens === 5) {
      let lon = Math.atan2(Dy, Dx) - this.azRad;
      lon -= 2 * Math.PI * Math.floor((lon + Math.PI) / (2 * Math.PI));
      const lat = Math.asin(Dz > 1 ? 1 : (Dz < -1 ? -1 : Dz));
      this.x = cx - lon * this.eqScale; this.y = cy - lat * this.eqScale;
      return true;
    }
    const ax = this.ax;
    let cc = Dx * ax[0] + Dy * ax[1] + Dz * ax[2];
    if (cc > 1) cc = 1; else if (cc < -1) cc = -1;
    if (lens === 0) { if (cc < 0.0087) return false; }
    else if (cc < this.cosHf) return false;
    const px = Dx - cc * ax[0], py = Dy - cc * ax[1], pz = Dz - cc * ax[2];
    let pn = Math.sqrt(px * px + py * py + pz * pz);
    const th = Math.atan2(pn, cc), hf = this.hf;        // atan2: 小角度也精确 (acos 在 1 附近丢精度)
    if (pn < 1e-15) pn = 1;
    const rt = this.rt, up = this.up;
    const tr = (px * rt[0] + py * rt[1] + pz * rt[2]) / pn;
    const tu = (px * up[0] + py * up[1] + pz * up[2]) / pn;
    let r;
    switch (lens) {
      case 0: r = this.linFocal * Math.tan(th); break;
      case 1: r = rpx * Math.sin(th / 2) / Math.max(1e-9, Math.sin(hf / 2)); break;
      case 3: r = rpx * Math.tan(th / 2) / Math.max(1e-9, Math.tan(hf / 2)); break;
      case 4:
        if (th > Math.PI / 2) return false;
        r = rpx * Math.sin(th) / Math.max(1e-9, Math.sin(Math.min(hf, Math.PI / 2)));
        break;
      default: r = rpx * th / hf;
    }
    this.x = cx + r * tr; this.y = cy - r * tu;
    return true;
  }

  /** 屏幕 → 方向 (单位向量); 画面外返回 null */
  unproject(x, y) {
    const lens = this.lens;
    if (lens === 5) {
      const lon = (this.cx - x) / this.eqScale, lat = (this.cy - y) / this.eqScale;
      if (Math.abs(lon) > Math.PI || Math.abs(lat) > Math.PI / 2) return null;
      const wl = lon + this.azRad;
      return [Math.cos(lat) * Math.cos(wl), Math.cos(lat) * Math.sin(wl), Math.sin(lat)];
    }
    const dx = x - this.cx, dy = this.cy - y, r = Math.hypot(dx, dy), hf = this.hf, rpx = this.rpx;
    let th;
    switch (lens) {
      case 0: th = Math.atan(r / this.linFocal); break;
      case 1: { const s = r * Math.sin(hf / 2) / rpx; if (s > 1) return null; th = 2 * Math.asin(s); break; }
      case 3: th = 2 * Math.atan(r * Math.tan(hf / 2) / rpx); break;
      case 4: { const s = r * Math.sin(Math.min(hf, Math.PI / 2)) / rpx; if (s > 1) return null; th = Math.asin(s); break; }
      default: th = r * hf / rpx;
    }
    if (lens !== 0 && th > hf + 1e-9) return null;
    if (th > Math.PI) return null;
    const tr = r > 0 ? dx / r : 0, tu = r > 0 ? dy / r : 0;
    const s = Math.sin(th), c = Math.cos(th), ax = this.ax, rt = this.rt, up = this.up;
    return [c * ax[0] + s * (tr * rt[0] + tu * up[0]), c * ax[1] + s * (tr * rt[1] + tu * up[1]), c * ax[2] + s * (tr * rt[2] + tu * up[2])];
  }

  /** 相邻两点是否连线: 夹角不大且等距圆柱下不跨接缝 */
  connect(ax_, ay, az, bx, by, bz, cosMax) {
    if (ax_ * bx + ay * by + az * bz < cosMax) return false;
    if (this.lens === 5) {
      let la = Math.atan2(ay, ax_) - this.azRad, lb = Math.atan2(by, bx) - this.azRad;
      la -= 2 * Math.PI * Math.floor((la + Math.PI) / (2 * Math.PI));
      lb -= 2 * Math.PI * Math.floor((lb + Math.PI) / (2 * Math.PI));
      if (Math.abs(la - lb) > Math.PI) return false;
    }
    return true;
  }
}

/** 世界方向 → (配置方位°, 仰角°) */
export function azElOf(d) {
  let az = -Math.atan2(d[1], d[0]) / DEG;
  if (az <= -180) az += 360;
  return [az, Math.asin(Math.max(-1, Math.min(1, d[2]))) / DEG];
}

export const css = (c, a = 1) => a >= 1 ? `rgb(${c[0] | 0},${c[1] | 0},${c[2] | 0})` : `rgba(${c[0] | 0},${c[1] | 0},${c[2] | 0},${a})`;

// ============================================================ GPU 层: 晕面网格 + 采样点
const VS_COMMON = `
precision highp float;
uniform vec3 uAx, uUp, uRt;
uniform float uHf, uLin, uRpx, uEq, uCx, uCy, uAz, uLens, uW, uH, uMirror;
float wrapPi(float a) { return a - 6.283185307179586 * floor((a + 3.141592653589793) / 6.283185307179586); }
vec2 proj(vec3 D, vec3 C) {
  if (uLens > 4.5) {                       // 等距圆柱: 经度按三角形重心展开, 跨接缝的三角形不撕裂
    float lonC = wrapPi(atan(C.y, C.x) - uAz);
    float lon = lonC + wrapPi(atan(D.y, D.x) - uAz - lonC);
    float lat = asin(clamp(D.z, -1.0, 1.0));
    return vec2(uCx - lon * uEq, uCy - lat * uEq);
  }
  float c = dot(D, uAx);
  vec3 pv = D - c * uAx;
  float pn = length(pv);
  float th = atan(pn, c);
  if (pn < 1e-12) pn = 1.0;
  float tr = dot(pv, uRt) / pn, tu = dot(pv, uUp) / pn, r;
  if (uLens < 0.5) r = uLin * tan(min(th, 1.55));
  else if (uLens < 1.5) r = uRpx * sin(th * 0.5) / max(1e-9, sin(uHf * 0.5));
  else if (uLens < 2.5) r = uRpx * th / uHf;
  else if (uLens < 3.5) r = uRpx * tan(min(th, 3.1) * 0.5) / max(1e-9, tan(uHf * 0.5));
  else r = uRpx * sin(min(th, 1.5707963)) / max(1e-9, sin(min(uHf, 1.5707963)));
  return vec2(uCx + r * tr, uCy - r * tu);
}
bool culled(vec3 C) {
  if (dot(C, C) < 0.25) return true;       // 重心置零 = 该三角形在天空上太大 (见 setMeshes) 或该点无方向
  if (uLens > 4.5) return false;
  float c = dot(C, uAx);
  if (uLens < 0.5) return c < 0.02;
  if (uLens > 3.5) return c < 0.0;
  return atan(length(C - c * uAx), c) > min(uHf + 0.35, 3.05);
}
vec4 toClip(vec2 px) { return vec4(2.0 * px.x / uW - 1.0, 1.0 - 2.0 * px.y / uH, 0.5, 1.0); }
`;
const VS_MESH = VS_COMMON + `
attribute vec3 aDir, aCen, aCol;
varying vec3 vCol;
void main() {
  vec3 D = vec3(aDir.x, aDir.y * uMirror, aDir.z), C = vec3(aCen.x, aCen.y * uMirror, aCen.z);
  vCol = aCol;
  if (culled(C)) { gl_Position = vec4(2.0, 2.0, 2.0, 1.0); return; }
  gl_Position = toClip(proj(D, C));
}`;
const VS_PTS = VS_COMMON + `
attribute vec3 aDir, aCol;
uniform float uPs;
varying vec3 vCol;
void main() {
  vec3 D = vec3(aDir.x, aDir.y * uMirror, aDir.z);
  vCol = aCol;
  gl_PointSize = uPs;
  if (culled(D)) { gl_Position = vec4(2.0, 2.0, 2.0, 1.0); return; }
  gl_Position = toClip(proj(D, D));
}`;
const FS = `
precision highp float;
varying vec3 vCol;
uniform float uAlpha, uRound, uDpr, uCx, uCy, uRpx, uH, uPt;
void main() {
  vec2 p = vec2(gl_FragCoord.x / uDpr, uH - gl_FragCoord.y / uDpr);
  if (uRound > 0.5 && dot(p - vec2(uCx, uCy), p - vec2(uCx, uCy)) > uRpx * uRpx) discard;
  if (uPt > 0.5 && dot(gl_PointCoord - 0.5, gl_PointCoord - 0.5) > 0.25) discard;
  gl_FragColor = vec4(vCol, uAlpha);
}`;
// 底图: 一个带纹理的四边形 (屏幕像素坐标), 不透明, 画在晕面下面
const VS_BG = `
precision highp float;
attribute vec2 aPos, aUv;
uniform float uW, uH;
varying vec2 vUv;
void main() { vUv = aUv; gl_Position = vec4(2.0 * aPos.x / uW - 1.0, 1.0 - 2.0 * aPos.y / uH, 0.9, 1.0); }`;
const FS_BG = `
precision highp float;
varying vec2 vUv;
uniform sampler2D uTex;
void main() { gl_FragColor = vec4(texture2D(uTex, vUv).rgb, 1.0); }`;

class SkyGL {
  constructor(canvas) {
    this.cv = canvas;
    const gl = canvas.getContext('webgl', { antialias: true, preserveDrawingBuffer: true, alpha: false });
    this.gl = gl;
    if (!gl) return;
    const mk = (src, type) => { const s = gl.createShader(type); gl.shaderSource(s, src); gl.compileShader(s); if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) console.error(gl.getShaderInfoLog(s)); return s; };
    const prog = (vs, fs = FS) => {
      const p = gl.createProgram();
      gl.attachShader(p, mk(vs, gl.VERTEX_SHADER)); gl.attachShader(p, mk(fs, gl.FRAGMENT_SHADER));
      gl.linkProgram(p);
      if (!gl.getProgramParameter(p, gl.LINK_STATUS)) console.error(gl.getProgramInfoLog(p));
      const u = {};
      const n = gl.getProgramParameter(p, gl.ACTIVE_UNIFORMS);
      for (let i = 0; i < n; i++) { const a = gl.getActiveUniform(p, i); u[a.name] = gl.getUniformLocation(p, a.name); }
      return { p, u, aDir: gl.getAttribLocation(p, 'aDir'), aCen: gl.getAttribLocation(p, 'aCen'), aCol: gl.getAttribLocation(p, 'aCol'),
        aPos: gl.getAttribLocation(p, 'aPos'), aUv: gl.getAttribLocation(p, 'aUv') };
    };
    this.pm = prog(VS_MESH); this.pp = prog(VS_PTS); this.pb = prog(VS_BG, FS_BG);
    this.bgTex = null; this.bgQuad = gl.createBuffer(); this.bgUv = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, this.bgUv); gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([0, 0, 1, 0, 0, 1, 1, 1]), gl.STATIC_DRAW);
    this.meshes = new Map();       // 数据数组 → GPU 缓冲
    this.clouds = new Map();
  }

  _buf(data) { const gl = this.gl, b = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, b); gl.bufferData(gl.ARRAY_BUFFER, data, gl.STATIC_DRAW); return b; }

  /** 底图纹理 (画布 / 位图, 已缩到纹理上限以内); null = 清除 */
  setBgImage(src) {
    const gl = this.gl;
    if (!gl) return;
    if (this.bgTex) { gl.deleteTexture(this.bgTex); this.bgTex = null; }
    if (!src) return;
    const tex = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, tex);
    gl.pixelStorei(gl.UNPACK_PREMULTIPLY_ALPHA_WEBGL, false);
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, src);
    // 任意尺寸 (WebGL1 非 2 的幂): 不做 mipmap, 边缘夹紧
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    this.bgTex = tex;
  }
  get maxTex() { return this.gl ? this.gl.getParameter(this.gl.MAX_TEXTURE_SIZE) : 4096; }
  _free(e) { for (const k of ['dir', 'cen', 'col']) if (e[k]) this.gl.deleteBuffer(e[k]); }

  /** 网格: [{ dir: Float32Array (每顶点), cols: Uint8Array | null, col: [r,g,b], subdiv? }] (三角形不共享顶点) */
  setMeshes(list) {
    if (!this.gl) return;
    const keep = new Set();
    for (const m0 of list) {
      keep.add(m0.dir);
      const colKey = m0.cols || String(m0.col);
      let e = this.meshes.get(m0.dir);
      if (e && e.colKey === colKey) continue;
      // subdiv: 天球上张角大的三角形细分 (方向取球面中点, 一分为四, 直到都小于 25°) 而不剔除 ——
      //   二自由度曲面在灯点收拢, 扇形三角形一个顶点在灯上、另两个在二十几度外, 是真的面, 剔掉就是一道从灯拉出去的白带
      const m = m0.subdiv ? { ...m0, ...subdivSky(m0.dir, m0.cols, Math.cos(25 * Math.PI / 180)) } : m0;
      if (e && m0.subdiv) { this._free(e); e = null; }
      if (!e) {
        // 每顶点再带一份三角形重心方向: 着色器据此统一剔除 / 展开全景经度. 顶点间张角 > 25° 的三角形
        // (逐顶点投影后直边会严重偏离真实像) 重心置零 = 剔除; 网格加密后正常不会出现
        const D = m.dir, n = D.length / 3, C = new Float32Array(D.length), cosBig = Math.cos(25 * Math.PI / 180);
        for (let t = 0; t < n / 3; t++) {
          const o = 9 * t;
          let x = D[o] + D[o + 3] + D[o + 6], y = D[o + 1] + D[o + 4] + D[o + 7], z = D[o + 2] + D[o + 5] + D[o + 8];
          const l = Math.hypot(x, y, z) || 1; x /= l; y /= l; z /= l;
          let ok = l > 1e-6;
          for (let k = 0; k < 3 && ok; k++) ok = D[o + 3 * k] * x + D[o + 3 * k + 1] * y + D[o + 3 * k + 2] * z > cosBig;
          if (!ok) x = y = z = 0;
          for (let k = 0; k < 3; k++) { C[o + 3 * k] = x; C[o + 3 * k + 1] = y; C[o + 3 * k + 2] = z; }
        }
        e = { n, dir: this._buf(D), cen: this._buf(C), col: null };
        this.meshes.set(m0.dir, e);
      } else if (e.col) this.gl.deleteBuffer(e.col);
      const n = e.n, col = new Float32Array(3 * n);
      if (m.cols) for (let i = 0; i < 3 * n; i++) col[i] = m.cols[i] / 255;
      else for (let i = 0; i < n; i++) { col[3 * i] = m.col[0] / 255; col[3 * i + 1] = m.col[1] / 255; col[3 * i + 2] = m.col[2] / 255; }
      e.col = this._buf(col); e.colKey = colKey;
    }
    for (const [k, e] of this.meshes) if (!keep.has(k)) { this._free(e); this.meshes.delete(k); }
    this.meshList = list.map(m => this.meshes.get(m.dir));
  }

  /** 采样点: [{ dir: Float32Array, cols: Uint8Array | null, col }] */
  setClouds(list) {
    if (!this.gl) return;
    const keep = new Set();
    for (const c of list) {
      keep.add(c.dir);
      const colKey = c.cols || String(c.col);
      let e = this.clouds.get(c.dir);
      if (e && e.colKey === colKey) continue;
      if (!e) { e = { n: c.dir.length / 3, dir: this._buf(c.dir), col: null }; this.clouds.set(c.dir, e); }
      else if (e.col) this.gl.deleteBuffer(e.col);
      const n = e.n, col = new Float32Array(3 * n);
      if (c.cols) for (let i = 0; i < 3 * n; i++) col[i] = c.cols[i] / 255;
      else for (let i = 0; i < n; i++) { col[3 * i] = c.col[0] / 255; col[3 * i + 1] = c.col[1] / 255; col[3 * i + 2] = c.col[2] / 255; }
      e.col = this._buf(col); e.colKey = colKey;
    }
    for (const [k, e] of this.clouds) if (!keep.has(k)) { this._free(e); this.clouds.delete(k); }
    this.cloudList = list.map(c => this.clouds.get(c.dir));
  }

  // bgQuad: 底图四角的屏幕坐标 [左上, 右上, 左下, 右下] (不画底图时 null); haloA: 晕本身的不透明度倍数 (「晕透明」)
  render(p, V, dpr, mirror, bgQuad = null, haloA = 1) {
    const gl = this.gl;
    if (!gl) return;
    gl.viewport(0, 0, this.cv.width, this.cv.height);
    gl.clearColor(1, 1, 1, 1); gl.clearDepth(1);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    if (bgQuad && this.bgTex) {
      const pr = this.pb;
      gl.useProgram(pr.p);
      gl.uniform1f(pr.u.uW, p.W); gl.uniform1f(pr.u.uH, p.H);
      gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, this.bgTex); gl.uniform1i(pr.u.uTex, 0);
      for (let i = 0; i < 4; i++) gl.disableVertexAttribArray(i);       // 上一帧网格 / 点留下的属性数组 (缓冲可能已释放)
      gl.bindBuffer(gl.ARRAY_BUFFER, this.bgQuad); gl.bufferData(gl.ARRAY_BUFFER, bgQuad, gl.DYNAMIC_DRAW);
      gl.enableVertexAttribArray(pr.aPos); gl.vertexAttribPointer(pr.aPos, 2, gl.FLOAT, false, 0, 0);
      gl.bindBuffer(gl.ARRAY_BUFFER, this.bgUv);
      gl.enableVertexAttribArray(pr.aUv); gl.vertexAttribPointer(pr.aUv, 2, gl.FLOAT, false, 0, 0);
      gl.disable(gl.DEPTH_TEST); gl.disable(gl.BLEND);
      gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
      gl.disableVertexAttribArray(pr.aPos); gl.disableVertexAttribArray(pr.aUv);
    }
    if (!(haloA > 0.002)) return;
    const round = V.lens >= 1 && V.lens <= 4 && V.fov < 359 ? 1 : 0;
    const setU = (pr) => {
      const u = pr.u;
      gl.uniform3fv(u.uAx, p.ax); gl.uniform3fv(u.uUp, p.up); gl.uniform3fv(u.uRt, p.rt);
      gl.uniform1f(u.uHf, p.hf); gl.uniform1f(u.uLin, p.linFocal); gl.uniform1f(u.uRpx, p.rpx); gl.uniform1f(u.uEq, p.eqScale);
      gl.uniform1f(u.uCx, p.cx); gl.uniform1f(u.uCy, p.cy); gl.uniform1f(u.uAz, p.azRad); gl.uniform1f(u.uLens, V.lens);
      gl.uniform1f(u.uW, p.W); gl.uniform1f(u.uH, p.H); gl.uniform1f(u.uRound, round); gl.uniform1f(u.uDpr, dpr);
    };
    const sides = mirror ? [1, -1] : [1];
    // 网格: 同深度 + LESS ⟹ 每个像素只上色一次 (两片重叠处不叠加变暗), 再按固定透明度混合到白底
    if (this.meshList && this.meshList.length) {
      const pr = this.pm;
      gl.useProgram(pr.p); setU(pr);
      gl.uniform1f(pr.u.uAlpha, 0.62 * haloA); gl.uniform1f(pr.u.uPt, 0);
      gl.enable(gl.DEPTH_TEST); gl.depthFunc(gl.LESS);
      gl.enable(gl.BLEND); gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
      for (const ys of sides) {
        gl.uniform1f(pr.u.uMirror, ys);
        for (const e of this.meshList) {
          if (!e) continue;
          bind(gl, pr.aDir, e.dir); bind(gl, pr.aCen, e.cen); bind(gl, pr.aCol, e.col);
          gl.drawArrays(gl.TRIANGLES, 0, e.n);
        }
      }
      gl.disable(gl.DEPTH_TEST);
    }
    if (this.cloudList && this.cloudList.length) {
      const pr = this.pp;
      gl.useProgram(pr.p); setU(pr);
      gl.uniform1f(pr.u.uAlpha, 0.55 * haloA); gl.uniform1f(pr.u.uPt, 1); gl.uniform1f(pr.u.uPs, 2.2 * dpr);
      gl.enable(gl.BLEND); gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
      for (const ys of sides) {
        gl.uniform1f(pr.u.uMirror, ys);
        for (const e of this.cloudList) {
          if (!e) continue;
          bind(gl, pr.aDir, e.dir); bind(gl, pr.aCol, e.col);
          gl.drawArrays(gl.POINTS, 0, e.n);
        }
      }
    }
    gl.disable(gl.BLEND);
  }
}
// 天球上张角大的三角形细分: 任一顶点离重心方向超过 acos(cosBig) 就按三边的球面中点一分为四 (颜色取平均), 最多 4 层.
//   张角超过 60° 的不细分 (交给剔除): 那么大的多半是连错的, 不是灯旁收拢的扇形
function subdivSky(D, cols, cosBig) {
  const cosHuge = Math.cos(60 * Math.PI / 180);
  const dir = [], col = cols ? [] : null;
  const big = (a, b, c) => {
    let x = a[0] + b[0] + c[0], y = a[1] + b[1] + c[1], z = a[2] + b[2] + c[2];
    const l = Math.hypot(x, y, z);
    if (!(l > 1e-6)) return false;                 // 退化的不管 (交给剔除)
    x /= l; y /= l; z /= l;
    return [a, b, c].some(v => v[0] * x + v[1] * y + v[2] * z < cosBig);
  };
  const mid = (a, b) => { const x = a[0] + b[0], y = a[1] + b[1], z = a[2] + b[2], l = Math.hypot(x, y, z) || 1; return [x / l, y / l, z / l]; };
  const cm = (a, b) => [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2];
  const emit = (a, b, c, ca, cb, cc, dep) => {
    if (dep < 4 && big(a, b, c)) {
      const ab = mid(a, b), bc = mid(b, c), ca_ = mid(c, a);
      const cab = col && cm(ca, cb), cbc = col && cm(cb, cc), cca = col && cm(cc, ca);
      emit(a, ab, ca_, ca, cab, cca, dep + 1); emit(ab, b, bc, cab, cb, cbc, dep + 1);
      emit(ca_, bc, c, cca, cbc, cc, dep + 1); emit(ab, bc, ca_, cab, cbc, cca, dep + 1);
      return;
    }
    dir.push(...a, ...b, ...c);
    if (col) col.push(...ca, ...cb, ...cc);
  };
  for (let o = 0; o + 8 < D.length; o += 9) {
    const v = k => [D[o + 3 * k], D[o + 3 * k + 1], D[o + 3 * k + 2]], c = k => cols ? [cols[o + 3 * k], cols[o + 3 * k + 1], cols[o + 3 * k + 2]] : null;
    const a = v(0), b = v(1), cc = v(2), dmin = Math.min(a[0] * b[0] + a[1] * b[1] + a[2] * b[2], b[0] * cc[0] + b[1] * cc[1] + b[2] * cc[2], cc[0] * a[0] + cc[1] * a[1] + cc[2] * a[2]);
    emit(a, b, cc, c(0), c(1), c(2), dmin < cosHuge ? 99 : 0);
  }
  return { dir: Float32Array.from(dir), cols: col ? Uint8Array.from(col, x => Math.round(x)) : null };
}

function bind(gl, loc, buf) {
  if (loc < 0) return;
  gl.bindBuffer(gl.ARRAY_BUFFER, buf);
  gl.enableVertexAttribArray(loc);
  gl.vertexAttribPointer(loc, 3, gl.FLOAT, false, 0, 0);
}

// ============================================================ 视图
/**
 * 数据 (setScene) 由 main 组装: { eye, lamp, mirror,
 *   meshes: [{ dir: Float32Array (三角形顶点视方向), cols, col }], clouds: [{ dir, cols, col }],
 *   lines: [{ segs: [{ p, d } | Float32Array(位置)…], col, segCols?, width, dash, alpha, skyOnly }],
 *   marks: [{ pts: Float32Array (位置), dirs?: Float32Array, col, kind, size }] }
 */
export class SkyView {
  constructor(canvas, glCanvas) {
    this.cv = canvas; this.ctx = canvas.getContext('2d');
    this.glc = glCanvas;
    this.gpu = new SkyGL(glCanvas);
    this.proj = new SkyProj();
    this.V = { az: 0, el: 10, fov: 180, lens: 2 };
    this.S = null;
    this.hover = null;
    this.onView = null;
    this.onHover = null;
    this.onBg = null;              // 拖图 / 滚轮改了底图摆放 → 回调同步面板
    this.dragging = false;
    // 底图 (同原窗体「发散光 BSDF 渲染」的底图拟合): 载入后整个投影按照片长宽比的画框取景 (显示与否都不变,
    //   开关底图时模拟不跳); 照片在画框里: scale 1 = 铺满画框, ox/oy 单位 = 半个画框 (oy 向上为正), rot 度 (逆时针为正)
    this.bgImg = null;             // { iw, ih } 原图像素 (纹理可能缩过, 比例不变)
    this.bg = { show: true, edit: false, scale: 1, ox: 0, oy: 0, rot: 0 };
    this.haloA = 1;                // 晕本身 (网格 / 采样点 / 二自由度层线) 的不透明度倍数
    this._bind();
  }

  setView(v) { Object.assign(this.V, v); this.V.fov = Math.min(this.V.fov, fovMax(this.V.lens)); this.draw(); }
  setBgImage(src, iw, ih) { this.gpu.setBgImage(src); this.bgImg = src ? { iw, ih } : null; this.draw(); }
  setBg(o) { Object.assign(this.bg, o); this.cv.style.cursor = this.bg.edit && this.bgImg ? 'move' : ''; this.draw(); }
  setHaloAlpha(a) { this.haloA = a; this.draw(); }
  get maxTex() { return this.gpu.maxTex; }
  // 画框: 照片按长宽比等比放进视图 (没有照片时 null = 整个视图)
  frame() {
    if (!this.bgImg || !(this.cssW > 0)) return null;
    const f = Math.min(this.cssW / this.bgImg.iw, this.cssH / this.bgImg.ih);
    return { w: this.bgImg.iw * f, h: this.bgImg.ih * f };
  }
  // 照片四角的屏幕坐标 [左上, 右上, 左下, 右下] (屏幕 y 向下, 所以逆时针 rot 对应 −rot)
  _bgQuad(fr) {
    const b = this.bg;
    if (!fr || !b.show) return null;
    const cx = this.cssW / 2 + b.ox * fr.w / 2, cy = this.cssH / 2 - b.oy * fr.h / 2;
    const hw = fr.w * b.scale / 2, hh = fr.h * b.scale / 2, a = -b.rot * DEG, c = Math.cos(a), s = Math.sin(a);
    const q = new Float32Array(8);
    [[-hw, -hh], [hw, -hh], [-hw, hh], [hw, hh]].forEach(([x, y], i) => { q[2 * i] = cx + x * c - y * s; q[2 * i + 1] = cy + x * s + y * c; });
    return q;
  }
  setScene(S) {
    this.S = S;
    this.gpu.setMeshes(S ? S.meshes || [] : []);
    this.gpu.setClouds(S ? S.clouds || [] : []);
    this.draw();
  }

  _bind() {
    const cv = this.cv;
    let drag = null;
    const bgEdit = () => this.bg.edit && this.bgImg;
    cv.addEventListener('pointerdown', e => {
      drag = { x: e.clientX, y: e.clientY, az: this.V.az, el: this.V.el, bg: bgEdit() ? { ...this.bg } : null };
      try { cv.setPointerCapture(e.pointerId); } catch (err) { /* 合成事件没有指针 */ }
      this.dragging = true;
    });
    cv.addEventListener('pointermove', e => {
      const r = cv.getBoundingClientRect();
      if (drag && drag.bg) {
        // 拖图: 左键拖动 = 移动 (按画框的显示尺寸换算, 与画框同一基准), Shift + 拖动 = 旋转 (0.3°/px)
        const dx = e.clientX - drag.x, dy = e.clientY - drag.y, fr = this.frame();
        if (e.shiftKey) this.bg.rot = wrap180(drag.bg.rot + dx * 0.3);
        else if (fr) {
          this.bg.ox = Math.max(-3, Math.min(3, drag.bg.ox + dx / (fr.w / 2)));
          this.bg.oy = Math.max(-3, Math.min(3, drag.bg.oy - dy / (fr.h / 2)));
        }
        this.draw();
        if (this.onBg) this.onBg(this.bg);
        return;
      }
      if (drag) {
        const dx = e.clientX - drag.x, dy = e.clientY - drag.y;
        const k = this.V.fov / Math.max(200, Math.min(r.width, r.height));
        if (this.V.lens === 5) this.V.az = wrap180(drag.az - dx * k * 0.5);
        else {
          this.V.az = wrap180(drag.az - dx * k);
          this.V.el = Math.max(-90, Math.min(90, drag.el + dy * k));
        }
        this.draw();
        if (this.onView) this.onView(this.V);
        return;
      }
      this._hoverAt(e.clientX - r.left, e.clientY - r.top);
    });
    const end = () => { if (drag) { drag = null; this.dragging = false; this.draw(); } };
    cv.addEventListener('pointerup', end);
    cv.addEventListener('pointercancel', end);
    cv.addEventListener('pointerleave', () => { this.hover = null; this.drawOverlayOnly(); if (this.onHover) this.onHover(null); });
    cv.addEventListener('wheel', e => {
      e.preventDefault();
      if (bgEdit()) {                       // 拖图模式: 滚轮缩放照片 (乘性)
        this.bg.scale = Math.max(0.01, Math.min(20, this.bg.scale * (e.deltaY < 0 ? 1.08 : 1 / 1.08)));
        this.draw();
        if (this.onBg) this.onBg(this.bg);
        return;
      }
      const f = Math.exp(e.deltaY * 0.0012);
      this.V.fov = Math.max(2, Math.min(fovMax(this.V.lens), this.V.fov * f));
      this.draw();
      if (this.onView) this.onView(this.V);
    }, { passive: false });
  }

  _hoverAt(x, y) {
    const d = this.proj.unproject(x, y);
    this.hover = d ? { x, y, d } : null;
    if (this.onHover) {
      if (!d || !this.S) { this.onHover(null); }
      else {
        const [az, el] = azElOf(d);
        const lv = sub3(this.S.lamp, this.S.eye), ln = Math.hypot(...lv);
        const ang = ln > 0 ? Math.acos(Math.max(-1, Math.min(1, (d[0] * lv[0] + d[1] * lv[1] + d[2] * lv[2]) / ln))) / DEG : NaN;
        this.onHover({ az, el, ang });
      }
    }
    this.drawOverlayOnly();
  }

  resize() {
    const r = this.cv.getBoundingClientRect(), dpr = Math.min(2, window.devicePixelRatio || 1);
    const W = Math.max(1, Math.round(r.width * dpr)), H = Math.max(1, Math.round(r.height * dpr));
    for (const c of [this.cv, this.glc]) if (c.width !== W || c.height !== H) { c.width = W; c.height = H; }
    this.cssW = r.width; this.cssH = r.height; this.dpr = dpr;
  }

  draw() {
    if (this._raf) return;
    this._raf = requestAnimationFrame(() => { this._raf = 0; this._draw(); });
  }
  drawNow() { if (this._raf) { cancelAnimationFrame(this._raf); this._raf = 0; } this._draw(); }

  // 悬停十字只重画线层: 第一次悬停时把干净的线层存一份 (拖动 / 重画时不存, 省掉整幅 getImageData)
  drawOverlayOnly() {
    if (this._raf || !this.cssW) { this.draw(); return; }
    const ctx = this.ctx;
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    if (this._base) ctx.putImageData(this._base, 0, 0);
    else this._base = ctx.getImageData(0, 0, this.cv.width, this.cv.height);
    ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    this._drawHover();
  }

  _draw() {
    this.resize();
    const W = this.cssW, H = this.cssH, ctx = this.ctx, V = this.V, p = this.proj;
    if (!(W > 2 && H > 2)) return;
    const fr = this.frame();
    p.set(V.az, V.el, V.fov, V.lens, W, H, fr);
    const S = this.S;
    this.gpu.render(p, V, this.dpr, !!(S && S.mirror), this._bgQuad(fr), this.haloA);
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, this.cv.width, this.cv.height);
    ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    this._drawGrid();
    if (S) { this._drawLines(); this._drawMarks(); this._drawLampEtc(); }
    this._drawFrame();
    this._base = null;
    if (this.hover && !this.dragging) { this._base = ctx.getImageData(0, 0, this.cv.width, this.cv.height); this._drawHover(); }
  }

  // 一条折线 → 画面路径. 有精确视方向 (d) 就用它; 只有位置时退回 P−E. 方向跳变 / 不可见 / 全景接缝处断开.
  _pathOf(seg, ys) {
    const p = this.proj, E = this.S.eye, ctx = this.ctx;
    const P = seg.p || seg, Dd = seg.d || null, n = P.length / 3;
    const cosMax = Math.cos(12 * DEG);
    let pen = false, qx = 0, qy = 0, qz = 0;
    for (let i = 0; i < n; i++) {
      let ux, uy, uz;
      if (Dd) { ux = Dd[3 * i]; uy = ys * Dd[3 * i + 1]; uz = Dd[3 * i + 2]; }
      else {
        const dx = P[3 * i] - E[0], dy = ys * P[3 * i + 1] - E[1], dz = P[3 * i + 2] - E[2];
        const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (r < 1e-9) { pen = false; continue; }
        ux = dx / r; uy = dy / r; uz = dz / r;
      }
      if (!p.project(ux, uy, uz)) { pen = false; continue; }
      if (pen && p.connect(qx, qy, qz, ux, uy, uz, cosMax)) ctx.lineTo(p.x, p.y);
      else ctx.moveTo(p.x, p.y);
      pen = true; qx = ux; qy = uy; qz = uz;
    }
  }

  _drawLines() {
    const S = this.S, ctx = this.ctx;
    ctx.lineJoin = 'round'; ctx.lineCap = 'round';
    const sides = S.mirror ? [1, -1] : [1];
    for (const ln of S.lines) {
      if (ln.only3d) continue;
      const al = (ln.alpha ?? 1) * (ln.body ? this.haloA : 1);      // body: 晕本身 (二自由度层线), 跟「晕透明」走
      if (!(al > 0.003)) continue;
      ctx.lineWidth = ln.width || 1.6;
      ctx.setLineDash(ln.dash || []);
      if (ln.segCols) {
        ln.segs.forEach((seg, k) => {
          ctx.strokeStyle = css(ln.segCols[k], al);
          ctx.beginPath();
          for (const ys of sides) this._pathOf(seg, ys);
          ctx.stroke();
        });
        continue;
      }
      ctx.strokeStyle = css(ln.col, al);
      ctx.beginPath();
      for (const ys of (ln.noMirror ? [1] : sides)) for (const seg of ln.segs) this._pathOf(seg, ys);
      ctx.stroke();
    }
    ctx.setLineDash([]);
  }

  _drawMarks() {
    const S = this.S, ctx = this.ctx, p = this.proj, E = S.eye;
    const sides = S.mirror ? [1, -1] : [1];
    for (const mk of S.marks) {
      const P = mk.pts, Dd = mk.dirs || null, n = P.length / 3;
      for (const ys of sides)
        for (let i = 0; i < n; i++) {
          let ux, uy, uz;
          if (Dd) { ux = Dd[3 * i]; uy = ys * Dd[3 * i + 1]; uz = Dd[3 * i + 2]; }
          else {
            const dx = P[3 * i] - E[0], dy = ys * P[3 * i + 1] - E[1], dz = P[3 * i + 2] - E[2];
            const r = Math.hypot(dx, dy, dz);
            if (r < 1e-9) continue;
            ux = dx / r; uy = dy / r; uz = dz / r;
          }
          if (!p.project(ux, uy, uz)) continue;
          glyph(ctx, p.x, p.y, mk.kind, mk.col, mk.size || 4);
        }
    }
  }

  _drawLampEtc() {
    const S = this.S, ctx = this.ctx, p = this.proj;
    const lv = sub3(S.lamp, S.eye), ln = Math.hypot(...lv);
    if (ln < 1e-9) return;
    const d = [lv[0] / ln, lv[1] / ln, lv[2] / ln];
    if (p.project(-d[0], -d[1], -d[2])) {
      ctx.strokeStyle = 'rgba(0,0,0,.55)'; ctx.lineWidth = 1.2;
      const x = p.x, y = p.y;
      ctx.beginPath(); ctx.moveTo(x - 5, y - 5); ctx.lineTo(x + 5, y + 5); ctx.moveTo(x - 5, y + 5); ctx.lineTo(x + 5, y - 5); ctx.stroke();
    }
    if (p.project(d[0], d[1], d[2])) {
      const x = p.x, y = p.y;
      ctx.fillStyle = '#f29a1a'; ctx.strokeStyle = '#8a4b00'; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.arc(x, y, 5.5, 0, 2 * Math.PI); ctx.fill(); ctx.stroke();
      ctx.strokeStyle = '#f29a1a'; ctx.lineWidth = 1.4;
      ctx.beginPath();
      for (let k = 0; k < 8; k++) {
        const a = k * Math.PI / 4;
        ctx.moveTo(x + 8 * Math.cos(a), y + 8 * Math.sin(a)); ctx.lineTo(x + 11 * Math.cos(a), y + 11 * Math.sin(a));
      }
      ctx.stroke();
    }
  }

  _drawGrid() {
    const ctx = this.ctx, p = this.proj, S = this.S;
    ctx.lineWidth = 1;
    ctx.strokeStyle = 'rgba(0,0,0,.10)';
    ctx.beginPath();
    for (const el of [-60, -30, 30, 60]) this._circleLat(el);
    for (let az = 0; az < 360; az += 30) this._meridian(az);
    ctx.stroke();
    ctx.strokeStyle = 'rgb(46,139,87)'; ctx.lineWidth = 1.5;
    ctx.beginPath(); this._circleLat(0); ctx.stroke();
    ctx.fillStyle = 'rgba(46,139,87,.95)'; ctx.font = `11px ${FONT}`; ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    for (let az = 0; az < 360; az += 30) {
      const a = -az * DEG;
      if (!p.project(Math.cos(a), Math.sin(a), 0)) continue;
      ctx.fillText(az === 0 ? (this.lampLabel || '0°') : `${az > 180 ? az - 360 : az}°`, p.x, p.y + 3);
    }
    if (p.project(0, 0, 1)) {
      ctx.fillStyle = 'rgba(0,0,0,.4)';
      ctx.beginPath(); ctx.arc(p.x, p.y, 2, 0, 2 * Math.PI); ctx.fill();
    }
    if (S && this.refCircle) {
      const lv = sub3(S.lamp, S.eye), ln = Math.hypot(...lv);
      if (ln > 1e-9) {
        ctx.strokeStyle = 'rgba(0,0,0,.28)'; ctx.setLineDash([3, 4]); ctx.lineWidth = 1;
        ctx.beginPath(); this._smallCircle([lv[0] / ln, lv[1] / ln, lv[2] / ln], this.refCircle); ctx.stroke();
        ctx.setLineDash([]);
      }
    }
  }

  _polyDirs(fn, n) {
    const p = this.proj, ctx = this.ctx;
    let pen = false, q = null;
    for (let i = 0; i <= n; i++) {
      const d = fn(i / n);
      if (!p.project(d[0], d[1], d[2])) { pen = false; continue; }
      if (pen && p.connect(q[0], q[1], q[2], d[0], d[1], d[2], Math.cos(20 * DEG))) ctx.lineTo(p.x, p.y);
      else ctx.moveTo(p.x, p.y);
      pen = true; q = d;
    }
  }
  _circleLat(el) {
    const ce = Math.cos(el * DEG), se = Math.sin(el * DEG);
    this._polyDirs(t => [ce * Math.cos(2 * Math.PI * t), ce * Math.sin(2 * Math.PI * t), se], 360);
  }
  _meridian(az) {
    const a = -az * DEG;
    this._polyDirs(t => { const e = (-89 + 178 * t) * DEG; return [Math.cos(e) * Math.cos(a), Math.cos(e) * Math.sin(a), Math.sin(e)]; }, 120);
  }
  _smallCircle(c, rDeg) {
    const t0 = Math.abs(c[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0];
    let e1 = cross3(c, t0); const l1 = Math.hypot(...e1); e1 = e1.map(v => v / l1);
    const e2 = cross3(c, e1);
    const cr = Math.cos(rDeg * DEG), sr = Math.sin(rDeg * DEG);
    this._polyDirs(t => {
      const a = 2 * Math.PI * t, ca = Math.cos(a), sa = Math.sin(a);
      return [cr * c[0] + sr * (ca * e1[0] + sa * e2[0]), cr * c[1] + sr * (ca * e1[1] + sa * e2[1]), cr * c[2] + sr * (ca * e1[2] + sa * e2[2])];
    }, 240);
  }

  _drawFrame() {
    const ctx = this.ctx, p = this.proj, V = this.V, fr = this.frame();
    if (fr) {                                    // 底图画框 (照片按比例放进视图的那块)
      ctx.strokeStyle = 'rgba(0,0,0,.28)'; ctx.lineWidth = 1; ctx.setLineDash([4, 4]);
      ctx.strokeRect(p.cx - fr.w / 2 + 0.5, p.cy - fr.h / 2 + 0.5, fr.w - 1, fr.h - 1);
      ctx.setLineDash([]);
    }
    if (V.lens === 0 || V.lens === 5) return;
    if (V.fov >= 359) return;
    ctx.strokeStyle = 'rgba(0,0,0,.35)'; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.arc(p.cx, p.cy, p.rpx, 0, 2 * Math.PI); ctx.stroke();
  }

  _drawHover() {
    const h = this.hover;
    if (!h) return;
    const ctx = this.ctx;
    ctx.strokeStyle = 'rgba(0,0,0,.55)'; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(h.x - 7, h.y); ctx.lineTo(h.x + 7, h.y); ctx.moveTo(h.x, h.y - 7); ctx.lineTo(h.x, h.y + 7); ctx.stroke();
  }

  /** 导出用: GPU 层 + 线层合成 */
  snapshot() {
    const hv = this.hover;
    this.hover = null;                 // 导出图里不要悬停十字
    this.drawNow();
    const c = document.createElement('canvas');
    c.width = this.cv.width; c.height = this.cv.height;
    const g = c.getContext('2d');
    g.fillStyle = '#fff'; g.fillRect(0, 0, c.width, c.height);
    g.drawImage(this.glc, 0, 0); g.drawImage(this.cv, 0, 0);
    this.hover = hv;
    return c;
  }
}

export function glyph(ctx, x, y, kind, col, size) {
  ctx.fillStyle = css(col); ctx.strokeStyle = '#000'; ctx.lineWidth = 0.9;
  ctx.beginPath();
  if (kind === 'tri') {
    ctx.moveTo(x, y - size * 1.15); ctx.lineTo(x + size, y + size * 0.75); ctx.lineTo(x - size, y + size * 0.75); ctx.closePath();
  } else if (kind === 'star') {
    for (let k = 0; k < 10; k++) {
      const r = k % 2 ? size * 0.45 : size * 1.2, a = -Math.PI / 2 + k * Math.PI / 5;
      if (k) ctx.lineTo(x + r * Math.cos(a), y + r * Math.sin(a)); else ctx.moveTo(x + r * Math.cos(a), y + r * Math.sin(a));
    }
    ctx.closePath();
  } else if (kind === 'ring') {
    ctx.arc(x, y, size, 0, 2 * Math.PI); ctx.stroke(); return;
  } else if (kind === 'pt') {
    ctx.arc(x, y, size, 0, 2 * Math.PI); ctx.fill(); return;
  } else ctx.arc(x, y, size, 0, 2 * Math.PI);
  ctx.fill(); ctx.stroke();
}

function sub3(a, b) { return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]; }
function cross3(a, b) { return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]; }
function wrap180(a) { a = ((a + 180) % 360 + 360) % 360 - 180; return a; }
