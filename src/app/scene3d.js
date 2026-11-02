// 旁观者 3D 视图 (three.js): 灯、眼、地面网格, 晕面点云, 边界/焦散折线, 锥切片的入射锥, 发散度球壳.
// 世界系与引擎一致: Z 朝上, 灯在 +X 方向的 XZ 平面内, 眼在 (0,0,h). 按需渲染 (改动或拖动时才画).

import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import { LineSegments2 } from 'three/examples/jsm/lines/LineSegments2.js';
import { LineSegmentsGeometry } from 'three/examples/jsm/lines/LineSegmentsGeometry.js';
import { LineMaterial } from 'three/examples/jsm/lines/LineMaterial.js';

const rgbHex = c => (c[0] << 16) | (c[1] << 8) | c[2];
// sRGB 0–255 → 线性 (顶点色要线性值; 输出时 three 再转回 sRGB)
const SRGB2LIN = new Float32Array(256).map((_, i) => { const c = i / 255; return c <= 0.04045 ? c / 12.92 : Math.pow((c + 0.055) / 1.055, 2.4); });
const col3 = c => new THREE.Color(`rgb(${c[0]},${c[1]},${c[2]})`);

/** 折线集 ({ p, d } 或位置数组) → 线段端点数组 (LineSegmentsGeometry 用) */
const posOf = s => s.p || s;
function segPositions(segs) {
  let n = 0;
  for (const s of segs) n += Math.max(0, posOf(s).length / 3 - 1);
  const out = new Float32Array(n * 6);
  let k = 0;
  for (const sg of segs) {
    const s = posOf(sg), m = s.length / 3;
    for (let i = 0; i + 1 < m; i++) {
      out[k++] = s[3 * i]; out[k++] = s[3 * i + 1]; out[k++] = s[3 * i + 2];
      out[k++] = s[3 * i + 3]; out[k++] = s[3 * i + 4]; out[k++] = s[3 * i + 5];
    }
  }
  return out;
}

// 数据对象 → 编号 (显示对象按数据身份缓存: 数据没变就不重建几何体)
const OID = new WeakMap();
let oidN = 0;
const oid = o => {
  if (!o || typeof o !== 'object') return String(o);
  let k = OID.get(o);
  if (k === undefined) { k = ++oidN; OID.set(o, k); }
  return k;
};

// 圆点贴图 (点云 / 标记用)
function dotTexture() {
  const c = document.createElement('canvas'); c.width = c.height = 64;
  const g = c.getContext('2d');
  g.fillStyle = '#fff'; g.beginPath(); g.arc(32, 32, 28, 0, 2 * Math.PI); g.fill();
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = THREE.SRGBColorSpace;
  return t;
}

export class Scene3D {
  constructor(host) {
    this.host = host;
    this.renderer = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true });
    this.renderer.setClearColor(0xffffff, 1);
    host.appendChild(this.renderer.domElement);
    this.renderer.domElement.className = 'cv3d';
    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(42, 1, 0.05, 5000);
    this.camera.up.set(0, 0, 1);
    this.controls = new OrbitControls(this.camera, this.renderer.domElement);
    this.controls.enableDamping = true; this.controls.dampingFactor = 0.12;
    this.controls.screenSpacePanning = true;
    this.controls.addEventListener('change', () => this.dirty = true);
    this.dot = dotTexture();
    this.lineMats = [];
    this.mats = new Map();             // 材质缓存: 材质不随数据销毁 → 着色器程序常驻, 换数据时不重新编译 (编译会卡主线程)
    this.content = new THREE.Group(); this.scene.add(this.content);
    this.base = new THREE.Group(); this.scene.add(this.base);
    this.labels = host.querySelector('.labels3d');
    this.dirty = true;
    this.active = true;
    this.userMoved = false;
    this.renderer.domElement.addEventListener('pointerdown', () => { this.userMoved = true; });
    this.renderer.domElement.addEventListener('wheel', () => { this.userMoved = true; }, { passive: true });
    new ResizeObserver(() => { this._resize(); this.dirty = true; }).observe(host);
    this._resize();
    const loop = () => {
      requestAnimationFrame(loop);
      if (!this.active) return;
      this.controls.update();
      if (this.dirty) { this.dirty = false; this._render(); }
    };
    requestAnimationFrame(loop);
  }

  setActive(on) { this.active = on; if (on) { this._resize(); this.dirty = true; } }

  _resize() {
    const r = this.host.getBoundingClientRect();
    const W = Math.max(1, Math.floor(r.width)), H = Math.max(1, Math.floor(r.height));
    this.renderer.setPixelRatio(Math.min(2, window.devicePixelRatio || 1));
    this.renderer.setSize(W, H, false);
    this.renderer.domElement.style.width = W + 'px';
    this.renderer.domElement.style.height = H + 'px';
    this.camera.aspect = W / H; this.camera.updateProjectionMatrix();
    for (const m of this.lineMats) m.resolution.set(W, H);
    this.W = W; this.H = H;
  }

  // 共享材质: 同参数只建一次, 标记 keep (显示对象销毁时不连带销毁)
  _mat(key, make) {
    let m = this.mats.get(key);
    if (!m) { m = make(); m.userData.keep = true; this.mats.set(key, m); }
    return m;
  }

  _lineMat(col, width, o = {}) {
    const key = 'L' + JSON.stringify([col, width, !!o.dashed, o.dashSize, o.gapSize, o.alpha ?? 1, !!o.vc]);
    return this._mat(key, () => {
      const m = new LineMaterial({
        color: col3(col).getHex(), linewidth: width, worldUnits: false,
        dashed: !!o.dashed, dashSize: o.dashSize || 0.6, gapSize: o.gapSize || 0.45,
        transparent: (o.alpha ?? 1) < 1, opacity: o.alpha ?? 1,
        // 镜像副本用 scale.y = −1, three 会把它的三角形当成反面; 粗线是屏幕空间四边形, 方向不随镜像翻转 → 必须双面, 否则另一半的线全被剔掉
        side: THREE.DoubleSide,
      });
      if (o.vc) m.vertexColors = true;
      m.resolution.set(this.W || 1, this.H || 1);
      this.lineMats.push(m);
      return m;
    });
  }

  _clear(group) {
    for (const ch of [...group.children]) {
      group.remove(ch);
      ch.traverse(o => {
        if (o.geometry && !o.userData.shared) o.geometry.dispose();
        if (o.material && !o.userData.sharedMat) {
          const ms = Array.isArray(o.material) ? o.material : [o.material];
          for (const m of ms) if (!m.userData.keep) m.dispose();
        }
      });
    }
  }

  /** 灯 / 眼 / 地面 / 灯-眼连线 (场景尺度随灯距变) */
  setBase(B) {
    this._clear(this.base);
    const L = B.lamp, E = B.eye;
    const span = Math.max(20, B.extent || 60);
    this.span = span;
    const rS = Math.max(0.12, span * 0.008);
    if (B.ground) {
      const size = Math.ceil(span * 1.25 / 10) * 20;
      const grid = new THREE.GridHelper(size, Math.round(size / niceStep(size / 30)), 0xb9c0c9, 0xe3e6ea);
      grid.rotation.x = Math.PI / 2;
      grid.position.set(L[0] * 0.5, 0, 0);
      this.base.add(grid);
    }
    if (B.foci) {
      const sg = new THREE.SphereGeometry(rS, 24, 16);
      const lamp = new THREE.Mesh(sg, new THREE.MeshBasicMaterial({ color: 0xf29a1a }));
      lamp.position.set(L[0], L[1], L[2]);
      const eye = new THREE.Mesh(sg.clone(), new THREE.MeshBasicMaterial({ color: 0x2563eb }));
      eye.position.set(E[0], E[1], E[2]);
      this.base.add(lamp, eye);
      this.focusMeshes = [lamp, eye]; this.rS = rS;
      this._scaleFoci();
      const g = new LineSegmentsGeometry(); g.setPositions([L[0], L[1], L[2], E[0], E[1], E[2]]);
      const ln = new LineSegments2(g, this._lineMat([120, 120, 120], 1.2, { dashed: true, dashSize: span * 0.012, gapSize: span * 0.01 }));
      ln.computeLineDistances();
      this.base.add(ln);
      // 灯 / 眼 在地面上的投影竖线
      const g2 = new LineSegmentsGeometry(); g2.setPositions([L[0], L[1], 0, L[0], L[1], L[2], E[0], E[1], 0, E[0], E[1], E[2]]);
      this.base.add(new LineSegments2(g2, this._lineMat([150, 150, 150], 1, { alpha: 0.6 })));
    }
    this.lamp = L; this.eye = E;
    this.labelPts = B.foci ? [{ p: L, text: B.lampText || 'L', cls: 'lamp' }, { p: E, text: B.eyeText || 'E', cls: 'eye' }] : [];
    if (!this.userMoved) this.preset('oblique', true);
    this.dirty = true;
  }

  /** 晕面内容 (同 SkyView 的显示列表; 只取 3D 有意义的部分). 每个显示项按 (数据身份 + 样式) 缓存,
   *  计算途中频繁刷新时只建新到的部分 (旧的网格明暗、线几何不重算). */
  setScene(S) {
    const old = this.cache || new Map(), next = new Map();
    this.cache = next;
    const mir = S && S.mirror ? 1 : 0;
    const use = (key, make) => {
      key += '|' + mir;
      if (next.has(key)) return next.get(key);
      if (old.has(key)) { const o = old.get(key); next.set(key, o); old.delete(key); return o; }
      const obj = make();
      if (!obj) { next.set(key, null); return null; }
      if (mir && obj.userData.mirror !== false) {
        const g = new THREE.Group(), m = obj.clone();
        m.traverse(o => { o.userData.shared = true; o.userData.sharedMat = true; });
        m.scale.y = -1;
        g.add(obj, m);
        this.content.add(g);
        next.set(key, g);
        return g;
      }
      this.content.add(obj); next.set(key, obj);
      return obj;
    };
    if (S) {
      for (const cl of S.clouds) use(`C${oid(cl.pts)}|${cl.cols ? oid(cl.cols) : cl.col}|${cl.size}|${cl.alpha}`, () => this._cloud(cl));
      for (const m of S.meshes || []) use(`M${oid(m.pos)}|${m.cols ? oid(m.cols) : m.col}`, () => this._mesh(m));
      for (const ln of S.lines) {
        if (ln.skyOnly) continue;
        use(`L${ln.segs.map(oid).join(',')}|${ln.segCols ? oid(ln.segCols) : ln.col}|${ln.width}|${ln.dash ? 1 : 0}|${ln.alpha}|${this.span}`, () => this._line(ln));
      }
      for (const mk of S.marks) {
        if (mk.skyOnly) continue;
        use(`K${oid(mk.pts3d || mk.pts)}|${mk.col}|${mk.size}`, () => this._marks(mk));
      }
      if (S.cone) { const g = use(`cone${S.cone.h}|${S.cone.R}|${this.lamp}`, () => this._cone(S.cone)); if (g) g.userData.fill.opacity = S.cone.alpha; }
      if (S.shell && S.shell.sphere) { const g = use(`shell${S.shell.r}|${this.eye}`, () => this._shell(S.shell)); if (g) g.userData.fill.opacity = S.shell.alpha; }
    }
    for (const obj of old.values()) if (obj) { this.content.remove(obj); this._dispose(obj); }
    this.dirty = true;
  }

  _dispose(root) {
    root.traverse(o => {
      if (o.geometry && !o.userData.shared) o.geometry.dispose();
      if (o.material && !o.userData.sharedMat) {
        const ms = Array.isArray(o.material) ? o.material : [o.material];
        for (const m of ms) if (!m.userData.keep) m.dispose();
      }
    });
  }

  _cloud(cl) {
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(cl.pts, 3));
    const o = { size: (cl.size || 1.6) * 1.6, sizeAttenuation: false, map: this.dot, alphaTest: 0.4, transparent: true,
      opacity: Math.min(1, (cl.alpha ?? 0.5) + 0.3), depthWrite: false };
    if (cl.cols) { g.setAttribute('color', new THREE.BufferAttribute(cl.cols, 3, true)); o.vertexColors = true; }
    else o.color = col3(cl.col || [70, 130, 180]);
    return new THREE.Points(g, this._mat('P' + JSON.stringify([o.size, o.opacity, cl.cols ? 0 : cl.col]), () => new THREE.PointsMaterial(o)));
  }

  // 晕面网格: 顶点色 × 烘焙的双面明暗 (MeshBasic, 不依赖灯光单位). 不透明 (半透明时重叠的两层按三角形次序混合, 一片斑驳);
  // 明暗用引擎给的每顶点平滑法向; 光源都在 xz 平面 (另加 |n_y| 项), 镜像副本与原件明暗完全一样, 两者重合处看不出打架.
  // 网格略往后推 (polygonOffset), 压在面上的边界 / 焦散线不被挡.
  _mesh(m) {
    const P = m.pos, nV = P.length / 3;
    if (!nV) return null;
    const N = m.nrm || null;
    const col = new Float32Array(3 * nV), base = m.col || [70, 130, 180];
    const l1 = [-0.45, 0, 0.89], l2 = [0.8, 0, 0.6];
    const shade = (nx, ny, nz) => 0.42 + 0.33 * Math.abs(nx * l1[0] + nz * l1[2]) + 0.13 * Math.abs(nx * l2[0] + nz * l2[2]) + 0.12 * Math.abs(ny);
    for (let t = 0; t < nV / 3; t++) {
      const o = 9 * t;
      let fs = 0;
      if (!N) {
        const ux = P[o + 3] - P[o], uy = P[o + 4] - P[o + 1], uz = P[o + 5] - P[o + 2];
        const vx = P[o + 6] - P[o], vy = P[o + 7] - P[o + 1], vz = P[o + 8] - P[o + 2];
        let nx = uy * vz - uz * vy, ny = uz * vx - ux * vz, nz = ux * vy - uy * vx;
        const nl = Math.hypot(nx, ny, nz) || 1; nx /= nl; ny /= nl; nz /= nl;
        fs = shade(nx, ny, nz);
      }
      for (let k = 0; k < 3; k++) {
        const v = 3 * t + k;
        const sh = N ? shade(N[3 * v], N[3 * v + 1], N[3 * v + 2]) : fs;
        const r = m.cols ? m.cols[3 * v] : base[0], g = m.cols ? m.cols[3 * v + 1] : base[1], b = m.cols ? m.cols[3 * v + 2] : base[2];
        col[3 * v] = SRGB2LIN[r] * sh; col[3 * v + 1] = SRGB2LIN[g] * sh; col[3 * v + 2] = SRGB2LIN[b] * sh;
      }
    }
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(P, 3));
    g.setAttribute('color', new THREE.BufferAttribute(col, 3));
    return new THREE.Mesh(g, this._mat('mesh', () => new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.DoubleSide,
      polygonOffset: true, polygonOffsetFactor: 1.5, polygonOffsetUnits: 2 })));
  }

  _line(ln) {
    const pos = segPositions(ln.segs3d || ln.segs);
    if (pos.length === 0) return null;
    const g = new LineSegmentsGeometry(); g.setPositions(pos);
    if (ln.segCols) {                 // 每条一种颜色 (二自由度 h 层)
      const cols = new Float32Array(pos.length);
      let k = 0;
      ln.segs.forEach((s, i) => {
        const c = ln.segCols[i], r = SRGB2LIN[c[0] | 0], gg = SRGB2LIN[c[1] | 0], b = SRGB2LIN[c[2] | 0];
        for (let j = 0; j < 2 * Math.max(0, posOf(s).length / 3 - 1); j++) { cols[k++] = r; cols[k++] = gg; cols[k++] = b; }
      });
      g.setColors(cols);
      return new LineSegments2(g, this._lineMat([255, 255, 255], (ln.width || 1.3) * 1.15, { alpha: ln.alpha ?? 1, vc: true }));
    }
    const obj = new LineSegments2(g, this._lineMat(ln.col, (ln.width || 1.6) * 1.15, { dashed: !!ln.dash, alpha: ln.alpha ?? 1, dashSize: this.span * 0.012, gapSize: this.span * 0.008 }));
    if (ln.dash) obj.computeLineDistances();
    return obj;
  }

  _marks(mk) {
    const pts = mk.pts3d || mk.pts;
    if (!pts || pts.length < 3) return null;
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(pts, 3));
    return new THREE.Points(g, this._mat('K' + JSON.stringify([mk.size, mk.col]), () => new THREE.PointsMaterial({ size: (mk.size || 4) * 2.4,
      sizeAttenuation: false, color: col3(mk.col), map: this.dot, alphaTest: 0.4, transparent: true })));
  }

  // 锥切片: 灯处仰角 h 的入射锥 (光线方向 û 的仰角 = −h)
  _cone(c) {
    const { h, R, alpha } = c, L = this.lamp;
    const half = (90 - h) * Math.PI / 180, sA = Math.sin(half), cA = Math.cos(half);
    const seg = 72, rings = 16, pos = [], idx = [];
    for (let r = 0; r <= rings; r++) {
      const t = r / rings * R;
      for (let i = 0; i <= seg; i++) {
        const ph = 2 * Math.PI * i / seg;
        pos.push(L[0] + sA * Math.cos(ph) * t, L[1] + sA * Math.sin(ph) * t, L[2] - cA * t);
      }
    }
    for (let r = 0; r < rings; r++)
      for (let i = 0; i < seg; i++) {
        const a = r * (seg + 1) + i, b = a + 1, cc = a + seg + 1, d = cc + 1;
        idx.push(a, cc, b, b, cc, d);
      }
    const grp = new THREE.Group();
    grp.userData.mirror = false;                // 锥面本身左右对称, 不再镜像一份
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3)); g.setIndex(idx);
    grp.userData.fill = new THREE.MeshBasicMaterial({ color: 0xf2a33a, transparent: true, opacity: alpha, side: THREE.DoubleSide, depthWrite: false });
    grp.add(new THREE.Mesh(g, grp.userData.fill));
    const rim = [];
    for (let i = 0; i <= seg; i++) { const ph = 2 * Math.PI * i / seg; rim.push(L[0] + sA * Math.cos(ph) * R, L[1] + sA * Math.sin(ph) * R, L[2] - cA * R); }
    const rg = new LineSegmentsGeometry(); rg.setPositions(segPositions([new Float32Array(rim)]));
    grp.add(new LineSegments2(rg, this._lineMat([214, 130, 20], 1.2, { alpha: 0.8 })));
    return grp;
  }

  // 发散度球壳: 以眼为心、半径 r 的半透明球 + 三个大圆
  _shell(sh) {
    const { r, alpha } = sh, E = this.eye;
    const grp = new THREE.Group();
    grp.userData.mirror = false;
    grp.userData.fill = new THREE.MeshBasicMaterial({ color: 0x2bb39a, transparent: true, opacity: alpha, side: THREE.DoubleSide, depthWrite: false });
    const mesh = new THREE.Mesh(new THREE.SphereGeometry(1, 64, 32), grp.userData.fill);
    mesh.scale.setScalar(r); mesh.position.set(E[0], E[1], E[2]);
    grp.add(mesh);
    const rings = [];
    for (let k = 0; k < 3; k++) {
      const a = [];
      for (let i = 0; i <= 120; i++) {
        const t = 2 * Math.PI * i / 120, c = Math.cos(t) * r, s = Math.sin(t) * r;
        if (k === 0) a.push(E[0] + c, E[1] + s, E[2]); else if (k === 1) a.push(E[0] + c, E[1], E[2] + s); else a.push(E[0], E[1] + c, E[2] + s);
      }
      rings.push(new Float32Array(a));
    }
    const rg = new LineSegmentsGeometry(); rg.setPositions(segPositions(rings));
    grp.add(new LineSegments2(rg, this._lineMat([30, 140, 120], 1, { alpha: 0.55 })));
    return grp;
  }

  /** 新结果的取景范围 (中心, 半径); 用户没动过镜头时自动取景 */
  frame(c, r) {
    this.focus = { c, r: Math.max(2, r) };
    this._scaleFoci();
    if (!this.userMoved) this.preset(this.lastPreset || 'oblique', true);
  }

  // 灯 / 眼小球: 半径 ≈ 取景半径的 1.6%
  _scaleFoci() {
    if (!this.focusMeshes || !this.focus) return;
    const s = clampN(this.focus.r * 0.016, 0.05, 1e9) / this.rS;
    for (const m of this.focusMeshes) m.scale.setScalar(s);
    this.dirty = true;
  }

  /** 镜头预设: oblique 斜看 / side 侧面 / top 俯视 / eye 从眼睛看灯 */
  preset(kind, quiet) {
    const L = this.lamp || [15, 0, 3], E = this.eye || [0, 0, 0];
    const span = this.span || 40;
    this.lastPreset = kind;
    const f = this.focus;
    const tgt = f ? new THREE.Vector3(f.c[0], f.c[1], f.c[2]) : new THREE.Vector3((L[0] + E[0]) / 2, 0, (L[2] + E[2]) / 2);
    const tv = Math.tan(this.camera.fov * Math.PI / 360), th = tv * Math.max(0.2, this.camera.aspect);
    const d = f ? f.r / Math.min(tv, th) * 1.05 : span * 1.35;
    const cam = this.camera;
    if (kind === 'side') cam.position.set(tgt.x, tgt.y - d, tgt.z + d * 0.05);
    else if (kind === 'top') cam.position.set(tgt.x, tgt.y - d * 0.02, tgt.z + d * 1.1);
    else if (kind === 'eye') {
      cam.position.set(E[0] - span * 0.02, E[1], E[2] + 0.001);
      tgt.set(L[0], L[1], L[2]);
    } else cam.position.set(tgt.x - d * 0.55, tgt.y - d * 0.75, tgt.z + d * 0.55);
    this.controls.target.copy(tgt);
    cam.lookAt(tgt);
    this.controls.update();
    if (!quiet) this.userMoved = true;
    this.dirty = true;
  }

  _render() {
    this.renderer.render(this.scene, this.camera);
    // 灯 / 眼 文字标签
    if (this.labels) {
      let html = '';
      const v = new THREE.Vector3();
      for (const lb of this.labelPts || []) {
        v.set(lb.p[0], lb.p[1], lb.p[2]).project(this.camera);
        if (v.z > 1 || v.z < -1) continue;
        const x = (v.x + 1) / 2 * this.W, y = (1 - v.y) / 2 * this.H;
        html += `<span class="lb3 ${lb.cls}" style="left:${x.toFixed(1)}px;top:${y.toFixed(1)}px">${lb.text}</span>`;
      }
      this.labels.innerHTML = html;
    }
  }

  snapshot() { this._render(); return this.renderer.domElement; }
}

const clampN = (x, a, b) => Math.max(a, Math.min(b, x));
function niceStep(x) {
  const p = Math.pow(10, Math.floor(Math.log10(x))), m = x / p;
  return (m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10) * p;
}
export { rgbHex };
