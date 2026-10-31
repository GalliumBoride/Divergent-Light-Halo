// 观测者天空视图: 站在眼睛处看出去. 晕面上每个晶体位置 P 映成方向 (P−E)/|P−E|, 再按镜头投影到画面.
// 镜头 (SkyProj) 与软件点栅 (Raster) 取自 LightPointTrace (与原窗体 / PreviewGLRenderer 一致):
//   方位 az 为配置方位角 (俯视顺时针增), az=0 朝灯; 画面右手 = 面朝视线的右手.
// 白底 + matplotlib 风格配色. 线与标记走 Canvas 2D, 点云走整数混色点栅.

const DEG = Math.PI / 180;
const FONT = '"PingFang SC", "Microsoft YaHei UI", "Microsoft YaHei", "Noto Sans SC", system-ui, sans-serif';

// ============================================================ 软件点栅
class Raster {
  constructor() { this.W = 0; this.H = 0; this.s = 1; this.img = null; this.buf = null; }
  resize(ctx, W, H, s) {
    if (W !== this.W || H !== this.H || !this.img) {
      this.img = ctx.createImageData(W, H);
      this.buf = new Uint32Array(this.img.data.buffer);
      this.W = W; this.H = H;
    }
    this.s = s;
  }
  clear(c) { this.buf.fill((0xFF000000 | (c[2] << 16) | (c[1] << 8) | c[0]) >>> 0); }
  // 定点整数混色 (sa: 0..256); x,y,size 为 CSS 像素
  dot(x, y, size, c, a) {
    let sa = (a * 256 + 0.5) | 0;
    if (sa <= 0) return; if (sa > 256) sa = 256;
    const s = this.s, W = this.W, H = this.H;
    const sz = Math.max(1, Math.round(size * s));
    let x0 = Math.round(x * s) - ((sz - 1) >> 1), y0 = Math.round(y * s) - ((sz - 1) >> 1);
    let x1 = x0 + sz - 1, y1 = y0 + sz - 1;
    if (x0 < 0) x0 = 0; if (y0 < 0) y0 = 0;
    if (x1 >= W) x1 = W - 1; if (y1 >= H) y1 = H - 1;
    if (x1 < x0 || y1 < y0) return;
    const sr = c[0] * sa, sg = c[1] * sa, sb = c[2] * sa, ia = 256 - sa, buf = this.buf;
    for (let yy = y0; yy <= y1; yy++) {
      const base = yy * W;
      for (let xx = x0; xx <= x1; xx++) {
        const d = buf[base + xx];
        const r = (sr + (d & 255) * ia + 128) >> 8;
        const g = (sg + ((d >> 8) & 255) * ia + 128) >> 8;
        const b = (sb + ((d >> 16) & 255) * ia + 128) >> 8;
        buf[base + xx] = 0xFF000000 | (b << 16) | (g << 8) | r;
      }
    }
  }
}

// ============================================================ 天球镜头
// lens: 0 直线 (透视) / 1 鱼眼等积 / 2 鱼眼等距 / 3 鱼眼体视 / 4 鱼眼正交 / 5 等距圆柱 (全景)
export const LENS_COUNT = 6;
export function fovMax(lens) {
  switch (lens) { case 0: return 170; case 3: return 355; case 4: return 180; default: return 360; }
}

export class SkyProj {
  set(az, el, fov, lens, W, H) {
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
    this.linFocal = 0.5 * Math.sqrt(W * W + H * H) / Math.tan(hfLin);
    this.eqScale = Math.min(W / 2.0, H) / Math.PI;
    this.cx = W / 2; this.cy = H / 2;
    this.rpx = Math.min(W, H) / 2 - 18;
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
    if (pn < 1e-12) pn = 1;
    const rt = this.rt, up = this.up;
    const tr = (px * rt[0] + py * rt[1] + pz * rt[2]) / pn;
    const tu = (px * up[0] + py * up[1] + pz * up[2]) / pn;
    const th = Math.acos(cc), hf = this.hf;
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

export const css = (c, a = 1) => a >= 1 ? `rgb(${c[0]},${c[1]},${c[2]})` : `rgba(${c[0]},${c[1]},${c[2]},${a})`;

// ============================================================ 视图
/**
 * 数据 (setScene) 由 main 组装: { eye:[x,y,z], lamp:[x,y,z], mirror:bool,
 *   clouds: [{ pts: Float32Array(3D 位置), col: [r,g,b] | null, cols: Uint8Array(rgb) | null, size, alpha }],
 *   lines:  [{ segs: [Float32Array(3D 位置)…], col, width, dash, alpha }],
 *   marks:  [{ pts: Float32Array, col, kind: 'dot'|'tri'|'star'|'ring', size }] }
 */
export class SkyView {
  constructor(canvas) {
    this.cv = canvas; this.ctx = canvas.getContext('2d');
    this.proj = new SkyProj(); this.ras = new Raster();
    this.V = { az: 0, el: 10, fov: 180, lens: 2 };
    this.S = null;
    this.hover = null;
    this.onView = null;       // 视角被拖动/缩放后回调 (同步到界面)
    this.onHover = null;      // 鼠标处读数回调
    this.dragging = false;
    this._bind();
  }

  setView(v) { Object.assign(this.V, v); this.V.fov = Math.min(this.V.fov, fovMax(this.V.lens)); this.draw(); }
  setScene(S) { this.S = S; this.draw(); }

  _bind() {
    const cv = this.cv;
    let drag = null;
    cv.addEventListener('pointerdown', e => {
      drag = { x: e.clientX, y: e.clientY, az: this.V.az, el: this.V.el };
      cv.setPointerCapture(e.pointerId);
      this.dragging = true;
    });
    cv.addEventListener('pointermove', e => {
      const r = cv.getBoundingClientRect();
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
    cv.addEventListener('pointerleave', () => { this.hover = null; this.draw(); if (this.onHover) this.onHover(null); });
    cv.addEventListener('wheel', e => {
      e.preventDefault();
      const f = Math.exp(e.deltaY * 0.0012);
      this.V.fov = Math.max(5, Math.min(fovMax(this.V.lens), this.V.fov * f));
      this.draw();
      if (this.onView) this.onView(this.V);
    }, { passive: false });
  }

  _hoverAt(x, y) {
    const d = this.proj.unproject(x, y);
    this.hover = d ? { x, y, d } : null;
    if (this.onHover) {
      if (!d || !this.S) { this.onHover(null); return; }
      const [az, el] = azElOf(d);
      const lv = sub3(this.S.lamp, this.S.eye), ln = Math.hypot(...lv);
      const ang = ln > 0 ? Math.acos(Math.max(-1, Math.min(1, (d[0] * lv[0] + d[1] * lv[1] + d[2] * lv[2]) / ln))) / DEG : NaN;
      this.onHover({ az, el, ang });
    }
    this.drawOverlayOnly();
  }

  resize() {
    const r = this.cv.getBoundingClientRect(), dpr = Math.min(2, window.devicePixelRatio || 1);
    const W = Math.max(1, Math.round(r.width * dpr)), H = Math.max(1, Math.round(r.height * dpr));
    if (this.cv.width !== W || this.cv.height !== H) { this.cv.width = W; this.cv.height = H; }
    this.cssW = r.width; this.cssH = r.height; this.dpr = dpr;
  }

  draw() {
    if (this._raf) return;
    this._raf = requestAnimationFrame(() => { this._raf = 0; this._draw(); });
  }
  drawNow() { if (this._raf) { cancelAnimationFrame(this._raf); this._raf = 0; } this._draw(); }

  drawOverlayOnly() {
    if (!this._base) { this.draw(); return; }
    const ctx = this.ctx;
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.putImageData(this._base, 0, 0);
    ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    this._drawHover();
  }

  _draw() {
    this.resize();
    const W = this.cssW, H = this.cssH, ctx = this.ctx, V = this.V, p = this.proj;
    if (!(W > 2 && H > 2)) return;
    p.set(V.az, V.el, V.fov, V.lens, W, H);
    // 1) 点栅: 白底 + 点云
    this.ras.resize(ctx, this.cv.width, this.cv.height, this.dpr);
    this.ras.clear([255, 255, 255]);
    const S = this.S;
    if (S) this._drawClouds();
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.putImageData(this.ras.img, 0, 0);
    ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    // 2) 网格 / 地平线 / 外框
    this._drawGrid();
    // 3) 线 + 标记
    if (S) { this._drawLines(); this._drawMarks(); this._drawLampEtc(); }
    this._drawFrame();
    this._base = ctx.getImageData(0, 0, this.cv.width, this.cv.height);
    this._drawHover();
  }

  _drawClouds() {
    const S = this.S, p = this.proj, ras = this.ras, E = S.eye;
    const thin = this.dragging ? 3 : 1;      // 拖动时抽稀
    const sides = S.mirror ? [1, -1] : [1];
    for (const cl of S.clouds) {
      const P = cl.pts, n = P.length / 3, cols = cl.cols;
      const col = cl.col || [70, 130, 180], size = cl.size || 1.6, alpha = cl.alpha ?? 0.5;
      for (const ys of sides) {
        for (let i = 0; i < n; i += thin) {
          const dx = P[3 * i] - E[0], dy = ys * P[3 * i + 1] - E[1], dz = P[3 * i + 2] - E[2];
          const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
          if (r < 1e-9) continue;
          if (!p.project(dx / r, dy / r, dz / r)) continue;
          if (cols) ras.dot(p.x, p.y, size, [cols[3 * i], cols[3 * i + 1], cols[3 * i + 2]], alpha);
          else ras.dot(p.x, p.y, size, col, alpha);
        }
      }
    }
  }

  // 一条 3D 折线 → 画面路径 (方向跳变 / 不可见 / 等距圆柱接缝处断开)
  _pathOf(seg, ys) {
    const p = this.proj, E = this.S.eye, ctx = this.ctx, n = seg.length / 3;
    const cosMax = Math.cos(25 * DEG);
    let pen = false, qx = 0, qy = 0, qz = 0;
    for (let i = 0; i < n; i++) {
      const dx = seg[3 * i] - E[0], dy = ys * seg[3 * i + 1] - E[1], dz = seg[3 * i + 2] - E[2];
      const r = Math.sqrt(dx * dx + dy * dy + dz * dz);
      if (r < 1e-9) { pen = false; continue; }
      const ux = dx / r, uy = dy / r, uz = dz / r;
      const ok = p.project(ux, uy, uz);
      if (!ok) { pen = false; continue; }
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
      ctx.strokeStyle = css(ln.col, ln.alpha ?? 1);
      ctx.lineWidth = ln.width || 1.6;
      ctx.setLineDash(ln.dash || []);
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
      const P = mk.pts, n = P.length / 3;
      for (const ys of sides)
        for (let i = 0; i < n; i++) {
          const dx = P[3 * i] - E[0], dy = ys * P[3 * i + 1] - E[1], dz = P[3 * i + 2] - E[2];
          const r = Math.hypot(dx, dy, dz);
          if (r < 1e-9 || !p.project(dx / r, dy / r, dz / r)) continue;
          glyph(ctx, p.x, p.y, mk.kind, mk.col, mk.size || 4);
        }
    }
  }

  _drawLampEtc() {
    const S = this.S, ctx = this.ctx, p = this.proj;
    const lv = sub3(S.lamp, S.eye), ln = Math.hypot(...lv);
    if (ln < 1e-9) return;
    const d = [lv[0] / ln, lv[1] / ln, lv[2] / ln];
    // 反灯点
    if (p.project(-d[0], -d[1], -d[2])) {
      ctx.strokeStyle = 'rgba(0,0,0,.55)'; ctx.lineWidth = 1.2;
      const x = p.x, y = p.y;
      ctx.beginPath(); ctx.moveTo(x - 5, y - 5); ctx.lineTo(x + 5, y + 5); ctx.moveTo(x - 5, y + 5); ctx.lineTo(x + 5, y - 5); ctx.stroke();
    }
    // 灯
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
    // 仰角圈 (每 30°) 与方位线 (每 30°)
    ctx.strokeStyle = 'rgba(0,0,0,.10)';
    ctx.beginPath();
    for (const el of [-60, -30, 30, 60]) this._circleLat(el);
    for (let az = 0; az < 360; az += 30) this._meridian(az);
    ctx.stroke();
    // 地平线
    ctx.strokeStyle = 'rgb(46,139,87)'; ctx.lineWidth = 1.5;
    ctx.beginPath(); this._circleLat(0); ctx.stroke();
    // 方位刻度 (地平线上): 0 = 朝灯
    ctx.fillStyle = 'rgba(46,139,87,.95)'; ctx.font = `11px ${FONT}`; ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    for (let az = 0; az < 360; az += 30) {
      const a = -az * DEG;
      if (!p.project(Math.cos(a), Math.sin(a), 0)) continue;
      ctx.fillText(az === 0 ? (this.lampLabel || '0°') : `${az > 180 ? az - 360 : az}°`, p.x, p.y + 3);
    }
    // 天顶
    if (p.project(0, 0, 1)) {
      ctx.fillStyle = 'rgba(0,0,0,.4)';
      ctx.beginPath(); ctx.arc(p.x, p.y, 2, 0, 2 * Math.PI); ctx.fill();
    }
    // 22° 参考圈 (以灯为心, 浅灰虚线)
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
    const ctx = this.ctx, p = this.proj, V = this.V;
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

  /** 导出用: 当前画面 (含线) 的离屏副本 */
  snapshot() { this.drawNow(); return this.cv; }
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
