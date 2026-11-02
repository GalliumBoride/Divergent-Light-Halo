// 术语解释: 先说它管什么事, 再说怎么算, 符号放最后. 键 → { zh: [标题, HTML], en: [...] }.
export default {
  scene: {
    zh: ['灯与眼 (发散光)', '<p>太阳在无穷远, 所有冰晶看到的太阳都在同一个方向. 路灯不是: 每个晶体看到的"太阳"(灯) 方向随它所在的位置变. 所以路灯下的晕不是天上的一条弧, 而是空间里的一张<b>晕面</b>.</p><p>灯在眼睛正前方 (方位 0°), 仰角 Σ、距离 D; 眼高是眼睛离地面的高度. 坐标: Z 朝上, 灯 L = (D cosΣ, 0, D sinΣ), 眼 E = (0, 0, 眼高).</p>'],
    en: ['Lamp and eye (divergent light)', '<p>The Sun is infinitely far away, so every ice crystal sees it in the same direction. A street lamp is not: the direction in which a crystal sees the "sun" (the lamp) depends on where the crystal is. So a lamp halo is not an arc in the sky but a <b>halo surface</b> in space.</p><p>The lamp is straight ahead of the eye (azimuth 0°) at elevation Σ and distance D; eye height is the eye above the ground. Coordinates: Z up, lamp L = (D cosΣ, 0, D sinΣ), eye E = (0, 0, eye height).</p>'],
  },
  path: {
    zh: ['光路', '<p>光在冰晶里依次经过哪些面: 第一个面是入射面, 最后一个是出射面, 中间的是内部反射面. 例: 3-5 = 从柱面 3 进、柱面 5 出 (22° 晕、幻日、切弧都是它).</p><p>面号: 底面 1 (顶) / 2 (底); 柱面 3–8 逆时针每 60°; 锥面 13–18 / 23–28 在柱面之上 / 之下. 多条光路用 ; 隔开, 一起画.</p>'],
    en: ['Ray path', '<p>The faces the light meets inside the crystal, in order: the first is the entry face, the last the exit face, the ones between are internal reflections. Example: 3-5 = in through prism face 3, out through prism face 5 (the 22° halo, sundogs and tangent arcs all use it).</p><p>Faces: basal 1 (top) / 2 (bottom); prism faces 3–8 every 60° counter-clockwise; pyramid faces 13–18 / 23–28 above / below the prism faces. Separate several paths with ;.</p>'],
  },
  orientation: {
    zh: ['晶体取向与自由度', '<p>冰晶在空中怎么摆. 一个自由度 = 只能绕竖直轴转 (板晶、Parry、自定义 ζ θ): 晕面是一张曲面. 两个自由度 = 还能再转一个角 (柱晶、Lowitz): 能到达的位置是一块"体", 画的是它的焦散或选定的一条边界曲面.</p><p>R(α, ζ, θ) = Rz(α)·Ry(ζ)·Rz(θ): α 方位, ζ c 轴天顶角, θ 绕 c 轴自旋.</p>'],
    en: ['Orientation and degrees of freedom', '<p>How the crystal is held in the air. One degree of freedom = it can only turn about the vertical (plate, Parry, custom ζ θ): the halo surface is a surface. Two degrees of freedom = one more angle is free (column, Lowitz): the reachable positions fill a volume, and what is drawn is its caustic or a chosen boundary surface.</p><p>R(α, ζ, θ) = Rz(α)·Ry(ζ)·Rz(θ): α azimuth, ζ zenith angle of the c-axis, θ spin about the c-axis.</p>'],
  },
  surface: {
    zh: ['晕面 · 光锥法', '<p>空间里所有能把灯光按这条光路送进眼睛的晶体位置. 观测者看到的晕 = 从眼睛看出去这张曲面的样子 (天空视图).</p><p><b>光锥法</b>: 按灯处入射仰角 h 一层层 (每层是以灯为顶点的圆锥), 在锥上闭式解出方位 ψ 和晶体位置 P = L + t·û.</p><p><b>网格</b> (一自由度): 在 (h, 晶体方位 α) 参数平面上铺网格, 每个格点闭式定位; 跨出有效区的三角形沿精确判据二分到零 (入射掠射 / 镜面掠射 / 出射临界 / ψ 转折 / 捏合 / 截断), 所以曲面的边恰好落在边界线上. 二自由度时每个 h 层画成一条线; 第一轮算完后在有曲线的 h 区间里再补层, 线够密才看得出曲面.</p><p>着色: <b>h</b> = 晶体看到灯的仰角; <b>ρ = r/D</b> = 晶体离眼的距离 ÷ 灯距 (发散度, ρ→0 就是平行光).</p>'],
    en: ['Halo surface · light-cone method', '<p>Every crystal position that sends lamp light along this path into the eye. What the observer sees = this surface as seen from the eye (the sky view).</p><p><b>Light-cone method</b>: layer by layer in the incidence elevation h at the lamp (each layer is a cone with its apex at the lamp), solve the azimuth ψ and the position P = L + t·û on the cone in closed form.</p><p><b>Mesh</b> (one degree of freedom): a grid on the (h, crystal azimuth α) parameter plane, every node placed exactly; triangles that leave the valid region are bisected to the exact criterion (entry grazing / mirror grazing / exit critical / ψ turn / pinch / truncation), so the edge of the surface lies on the boundary lines. With two degrees of freedom each h layer is drawn as a line; after the first pass more layers are added where curves exist, so the lines are dense enough to show the surface.</p><p>Colours: <b>h</b> = elevation at which the crystal sees the lamp; <b>ρ = r/D</b> = crystal-to-eye distance ÷ lamp distance (divergency; ρ→0 is parallel light).</p>'],
  },
  pts: {
    zh: ['光锥采样点 (原窗体)', '<p>原窗体画的点: 光锥法在均匀的 (h, α) 网格与 (α, ψ) 网格上闭式求出的晶体位置, 每个点都精确, 但分布不均 — 靠近出射临界边界时位置随 h 变化极快 (平方根型), 均匀取样在那里会留下大片空白, 所以点离边界线有一段距离. 晕面网格没有这个问题.</p>'],
    en: ['Light-cone samples (original)', '<p>The points the original window drew: crystal positions solved in closed form on uniform (h, α) and (α, ψ) grids. Each point is exact, but they are unevenly spread: near the exit-critical boundary the position changes very fast with h (square-root behaviour), so uniform sampling leaves wide gaps and the points stop short of the boundary lines. The surface mesh does not have this problem.</p>'],
  },
  bEntry: {
    zh: ['入射掠射边界', '<p>光正好擦着入射面进晶体 (入射角 90°). 再往外光就进不去了, 晕面在这里终止.</p><p>一维闭式: 入射方向躺在入射面切平面的单位圆上, 沿圆参数 τ 走一圈, 每点闭式定位; 角点 / 转折 / 捏合处二分到机器精度.</p>'],
    en: ['Entry-grazing boundary', '<p>The light just skims the entry face (incidence 90°). Beyond it the light cannot get in, so the surface ends here.</p><p>Closed form in one parameter: the incident direction lies on the unit circle in the entry face’s tangent plane; walk the circle (τ), place each point exactly, and bisect corners / turns / pinches to machine precision.</p>'],
  },
  bMirror: {
    zh: ['镜面掠射边界', '<p>光在某个内部反射面上正好擦过 (到不了那个面). 只有带内部反射的光路才有. 一维闭式构造同入射掠射.</p>'],
    en: ['Mirror-grazing boundary', '<p>Inside the crystal the light just skims one of the reflecting faces (it stops reaching that face). Only paths with internal reflections have it. Same one-parameter closed form as entry grazing.</p>'],
  },
  bTir: {
    zh: ['出射临界边界', '<p>光到出射面时正好处在全反射临界角: 再斜一点就出不去 (全反射). 二自由度时也包括"蓝边" (部分反射的台阶).</p>'],
    en: ['Exit-critical boundary', '<p>At the exit face the light is exactly at the critical angle: a little more oblique and it is totally reflected and cannot leave. With two degrees of freedom this bucket also holds blue edges (partial-reflection steps).</p>'],
  },
  bSeam: {
    zh: ['接缝 (二自由度)', '<p>内部反射面上光的法向冲量正好为零 (q = 0) 的那条线, 晕面在这里与相邻光路接上. 由平行光末端点的身份沿 h 延拓得到.</p>'],
    en: ['Seam (two degrees of freedom)', '<p>Where the light’s normal momentum on an internal reflecting face is exactly zero (q = 0); the surface joins the neighbouring path there. Obtained by continuing parallel-light endpoints along h.</p>'],
  },
  bTurn: {
    zh: ['转折缝 (二自由度)', '<p>焦散曲面碰到"可达带"的边 (晶体姿态在这里转不过去) 留下的线: 平行光焦散的转折点随 h 走出来的轨迹.</p>'],
    en: ['Turn seam (two degrees of freedom)', '<p>Where the caustic surface meets the edge of the reachable band (the orientation cannot turn further): the track of the parallel-light caustic’s turning points as h changes.</p>'],
  },
  bCusp: {
    zh: ['尖点脊线 (二自由度)', '<p>平行光焦散在天上的尖点 (折叠的折叠) 随 h 走出来的轨迹, 是焦散曲面自身的一条棱.</p>'],
    en: ['Cusp ridge (two degrees of freedom)', '<p>The track of the parallel-light caustic’s cusps (a fold of the fold) as h changes; an edge of the caustic surface itself.</p>'],
  },
  bRim: {
    zh: ['晕面外缘 (二自由度)', '<p>把相邻光锥层的层曲线端点连成的线, 也就是画出来的晕面实际的边. 上面几种边界线来自平行光曲线的末端点 (碰到别的边界、转折、尖点) 随 h 的延拓, 但晕面还会在别处收口: 有的末端点映到灯后或眼后 (看不见), 有的曲线挨着光路走不通的区域结束. 这条线把这些都包括在内; 与精确边界重合的地方被精确边界盖住.</p>'],
    en: ['Surface rim (two degrees of freedom)', '<p>The ends of the layer curves on neighbouring light-cone layers, joined up: the actual edge of the drawn surface. The boundary lines above come from continuing the parallel-light curve’s end points (meeting another boundary, turning, cusps) along h, but the surface can also close off elsewhere: some end points map behind the lamp or behind the eye (invisible), and some curves stop next to orientations where the ray path fails. This line covers all of them; where it coincides with an exact boundary, the exact one is drawn on top.</p>'],
  },
  caustic: {
    zh: ['发散焦散 (观测者)', '<p>观测者天上最亮的那些线: 很多晶体的光挤到同一个方向. 在空间里它是晕面从眼睛看过去的<b>轮廓</b> (曲面在这里转过身).</p><p>直接解, 不经过点云: 先追一条与灯无关的普适曲线 χ(ξ) = 0, 再用锥上闭式配对方位角并裁掉看不见的部分. 实线 = 有界段, 虚线 = 一直延伸到截断距离的段. 若没有, 结果里会说明原因 (真的没有 / 在灯后方 / 在眼后方 …).</p>'],
    en: ['Divergent caustic (observer)', '<p>The brightest lines in the observer’s sky, where light from many crystals piles up in one direction. In space it is the <b>outline</b> of the halo surface as seen from the eye (where the surface turns away).</p><p>Solved directly, not from the cloud: trace a lamp-independent universal curve χ(ξ) = 0, then pair the azimuth on the cone in closed form and cut what cannot be seen. Solid = bounded part, dashed = parts running out to the truncation distance. If there is none, the results panel says why (really none / behind the lamp / behind the eye …).</p>'],
  },
  cSeam: {
    zh: ['边界的天空像', '<p>掠射 / 全反射边界在观测者天上的样子 (K5 缝像): 亮度在这里突然变化.</p>'],
    en: ['Sky images of boundaries', '<p>The grazing / total-reflection boundaries as they appear in the observer’s sky (K5 seam images): the brightness steps there.</p>'],
  },
  cAsym: {
    zh: ['渐近圆', '<p>晕面延伸到无穷远时, 在天上趋近的圆 (逆反射方向, K4). 只有"无界"的晕面才有: 它的平行光弧经过反灯点.</p>'],
    en: ['Asymptotic circles', '<p>The circles in the sky that the halo surface approaches as it runs off to infinity (retro-reflection directions, K4). Only unbounded surfaces have them: their parallel-light arc passes through the anti-lamp point.</p>'],
  },
  parArc: {
    zh: ['平行光弧', '<p>把灯换成同一仰角、无穷远的太阳时, 这条光路的经典晕弧. 离眼极近的晶体看灯的方向和眼睛一样, 所以发散光晕总有一边贴着它 (只在天空视图里画).</p>'],
    en: ['Parallel-light arc', '<p>The classic halo arc this path would give if the lamp were a sun at the same elevation, infinitely far away. Crystals right next to the eye see the lamp in the same direction as the eye, so the divergent halo always touches it on one side (drawn in the sky view only).</p>'],
  },
  ref22: {
    zh: ['22° 参考圈', '<p>以灯为心、半径 22° 的圆, 用来估角度.</p>'],
    en: ['22° reference circle', '<p>A circle of radius 22° around the lamp, for judging angles.</p>'],
  },
  mirror: {
    zh: ['左右镜像补全', '<p>关于灯–眼竖直面左右翻转再画一遍. 很多光路要和它的镜像光路 (如 3-1-5 与 3-1-7) 合起来才左右对称; 打开它就不必两条都输入.</p>'],
    en: ['Left–right mirror fill-in', '<p>Draws everything again mirrored in the vertical plane through the lamp and the eye. Many paths are only symmetric together with their mirror path (3-1-5 with 3-1-7, for example); with this on you do not need to enter both.</p>'],
  },
  surf2: {
    zh: ['二自由度曲面', '<p>两个自由度时, 能把光送进眼睛的晶体位置填满一块体. 选一条平行光曲线做"曲面": 焦散 (最亮的缘) 或某条边界 (入射掠射 / 出射临界 / 反弹接缝 / 蓝边). 每个锥仰角 h 一层, 用光点法在取向平面上追这条曲线, 再映射到空间.</p><p>光路上没有选中的曲线时自动改用焦散.</p>'],
    en: ['2-DOF surface', '<p>With two degrees of freedom the crystal positions that reach the eye fill a volume. Pick one parallel-light curve to draw as the "surface": the caustic (brightest rim) or a boundary (entry grazing / exit critical / bounce seam / blue edge). Each cone elevation h is one layer: the light-point method traces the curve on the orientation plane, then it is mapped into space.</p><p>If the path has no such curve, the caustic is used.</p>'],
  },
  maxDist: {
    zh: ['截断距离', '<p>只画离灯这么远以内的晶体. 有的晕面 (Parry 3-5 等) 一直延伸到无穷远, 必须截断.</p>'],
    en: ['Truncation distance', '<p>Only crystals within this distance of the lamp are drawn. Some halo surfaces (Parry 3-5, for example) run off to infinity and have to be cut.</p>'],
  },
  cone: {
    zh: ['锥切片', '<p>灯处仰角为 h 的入射光锥 (光以仰角 −h 离开灯, 打到的晶体看灯的仰角正好是 h) 与晕面的交线: 光锥法里的一层. 一自由度时画锥上的晶体位置 (红点); 二自由度时是这一层的焦散 / 选定曲线 (红线). 天空视图里显示它们的视方向.</p>'],
    en: ['Cone slice', '<p>The incidence cone with elevation h at the lamp (light leaves the lamp downward at h, so the crystals it hits see the lamp at elevation h) cut with the halo surface: one layer of the light-cone method. One degree of freedom: the crystal positions on the cone (red dots); two degrees: that layer’s caustic / chosen curve (red lines). The sky view shows their directions.</p>'],
  },
  shell: {
    zh: ['发散度球壳', '<p>以眼为心、半径 r 的球面切晕面得到的晕线. 发散度 ρ = r/D (D = 眼到灯的距离), κ = lg ρ. ρ→0 是平行光晕, ρ→∞ 是反灯点剖面; 发散光照片大致是 ρ ≈ 1/4…4 这些层叠起来.</p><p>闭式: 一层壳上的晕线是晶体系太阳球上同一个二次型 Φ<sub>ρ</sub>(w) = 0, 被光路的边界剪成弧; 弧的末端 (白点) 是曲线与边界的交点, 折叠 (橙三角) 是两种落位重合处, 尖点 (品红星). κ 与 r 互相回填, 最后改的那个保持不变.</p>'],
    en: ['Divergency shell', '<p>The halo curve cut out of the surface by a sphere of radius r around the eye. Divergency ρ = r/D (D = eye-to-lamp distance), κ = lg ρ. ρ→0 is the parallel-light halo, ρ→∞ the anti-lamp profile; a photo of a lamp halo is roughly the layers ρ ≈ 1/4 … 4 stacked together.</p><p>Closed form: on one shell the curve is the zero set of one quadratic Φ<sub>ρ</sub>(w) = 0 on the crystal-frame sun sphere, cut into arcs by the path’s boundaries; arc ends (white dots) are where the curve meets a boundary, folds (orange triangles) are where the two placements coincide, cusps (magenta stars). κ and r fill each other in; the one you changed last stays fixed.</p>'],
  },
};
