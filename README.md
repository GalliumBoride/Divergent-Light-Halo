# Divergent Light Halos 3D

**English** · [中文](README.zh-CN.md)

Ice halos around a street lamp, in three dimensions. Sunlight arrives from one direction, so every crystal sees the Sun in the same place and a halo is an arc in the sky. A lamp is close by: each crystal sees it in a different direction, so the crystals that send lamp light into your eye fill a **halo surface** in space. This tool computes that surface, its boundaries and its caustics exactly, and shows them both from outside (3D) and as the observer sees them (sky view).

**Try it online: <https://galliumboride.github.io/Divergent-Light-Halo/>**

The whole tool is one HTML file. You can also download [`dist/DivergentLightHalo.html`](dist/DivergentLightHalo.html) and double-click it; it runs offline in any modern browser.

![Column crystals, path 3-5 (tangent arcs) under a lamp at 15°: the halo surface in 3D with its boundaries, the observer caustic and a divergency shell](docs/screenshot-column.png)
![Plate crystals, path 1-3 (circumzenithal arc) under a lamp at 25°: the surface mesh coloured by divergency ρ, ending exactly on its entry and exit boundaries, with a cone slice](docs/screenshot-plate.png)

## What it shows

- **The halo surface.** Every crystal position from which lamp light, following the chosen ray path, reaches the eye. With one degree of freedom it is drawn as a triangle mesh whose edges are bisected onto the exact boundaries; with two degrees of freedom each cone layer is drawn as a curve, and after the first pass extra layers are added where curves exist, so the family is dense enough to read as a surface. Colour by the elevation h at which the crystal sees the lamp, by the divergency ρ = r/D (distance from the eye ÷ lamp distance), by path, or in one colour. The original window's light-cone sample points can be shown too.
- **Boundaries.** Where the surface ends because a condition becomes critical: the light grazes the entry face, grazes an internal mirror face, or reaches the exit face exactly at the critical angle. With two degrees of freedom: seams, turn seams and cusp ridges of the caustic surface.
- **The observer caustic.** The brightest lines in the observer's sky, where light from many crystals piles up. In space it is the outline of the surface seen from the eye. Optional: sky images of the boundaries and the asymptotic circles of unbounded surfaces. When there is no caustic, the results panel says why (none at all, behind the lamp, behind the eye …).
- **Parallel-light arc.** The classic arc the same path gives if the lamp were a sun at the same elevation. The divergent halo always touches it.
- **Cone slice.** One layer of the computation: the cone of light leaving the lamp at elevation h, and the crystals on it.
- **Divergency shell.** The halo curve on a sphere of radius r around the eye (ρ = r/D, κ = lg ρ), with arc ends, folds and cusps. ρ → 0 is the parallel-light halo, ρ → ∞ the anti-lamp profile.

Two views, side by side or one at a time: a **3D scene** (orbit, zoom, pan; oblique / side / top / from-the-eye presets) and the **observer sky** (rectilinear, four fisheye projections or a panorama; drag to look around; hover for azimuth, elevation and distance from the lamp).

You can choose:
- **Lamp and eye**: lamp elevation (−90° to 90°), lamp distance (1–300 m), eye height.
- **Ray paths**: type them (`3-5; 3-1-5`) or click faces on the crystal; several paths at once. Pyramid faces 13–18 / 23–28 included.
- **Orientation**: plate, Parry, any fixed ζ/θ (one degree of freedom), column with fixed c-axis tilt or Lowitz with fixed spin (two degrees of freedom). For two degrees of freedom you can pick which parallel-light curve to draw as the surface (caustic or one of the boundaries).
- **More**: quality, truncation distance, refractive index, pyramid apex angles.

Exports: a PNG of the current views with a title, and a CSV of every surface point (position, azimuth/elevation of the exact viewing direction, distances, h, lg ρ) plus all boundary, caustic and shell curves. A glossary explains every term in plain language. Chinese and English interface.

## How it works

Coordinates: Z up, lamp L = (D cos Σ, 0, D sin Σ), eye E = (0, 0, h_eye), crystal orientation R = Rz(α)·Ry(ζ)·Rz(θ).

- **One degree of freedom — light-cone method.** For every incidence elevation h at the lamp, the incoming direction is û = Rz(ψ)·(−cos h, 0, −sin h). Deviation is symmetric about the vertical, so the outgoing direction is Rz(ψ)·d(α, h). The two-focus condition (E − L lies in the plane of û and the outgoing ray) becomes P sin ψ + Q cos ψ = R, solved in closed form; the crystal then sits at P = L + t·û with t from a 2×2 solve. A second sweep (for each α and ψ, bisect for h) fills the thin bands the first one misses.
- **Surface mesh.** Built row by row in h. On each row the valid α-intervals of each ψ branch are found exactly (sign changes of the validity margin bisected to 1e-13; intervals thinner than the sampling are found from local maxima of the margin and tracked from the neighbouring rows). Inside an interval the nodes follow α = α_a + (α_b − α_a)(1 − cos πτ)/2. Next to the exit-critical boundary and where the two ψ branches fold together the position varies like √(α − α_a); with α − α_a ∝ τ² it becomes a smooth function of τ, so the square-root singularity is absorbed by the coordinate rather than by refinement, and the boundary nodes are the exact roots (the two branches share them and join without a seam). Whole rows are inserted where neighbouring rows are far apart in 3D or in the sky, where the boundary bends, or where the intervals do not match one to one (an interval is born, dies, splits or merges); intervals are matched with an allowance of 20·Δh, because a thin interval can move further than its own width. Rows also continue past the coarse scan while intervals remain, so thin wedge tips are not cut off. Neighbouring rows are joined as quad strips: no T-junctions, and the exact boundary points lie on the mesh to 1 mm. Line intersections use |û × ê|² rather than 1 − (û·ê)², which keeps full precision when the two rays are nearly parallel (crystals far away). Regions thinner than the coarse scan are found from local maxima of the validity margin along each column and row (plate crystals under a low lamp exist only in a band of h just below Σ, 0.14° wide at Σ = 2°), and the minimum row spacing scales with |h| near horizontal rays. When the lamp is at eye height and every face on the path is vertical (plate prism paths), only horizontal rays connect lamp and eye and ψ becomes a free parameter: the surface is a flat region in the lamp–eye plane, meshed in (ψ, β = α − ψ), drawn in the sky as an arc on the horizon, with its boundary lines taken from the exact interval ends.
- **Boundaries (engine v4).** Each boundary family pins one point of the ray chain to a fixed circle, so it is a one-parameter closed form. Corners, turns and pinches are bisected to machine precision.
- **Drawn curves.** Every displayed line (boundaries, cone layers, edge tracks, caustics, shell arcs) is refined for display: points are added until neighbours are close in 3D, in the eye's view and in the lamp's view; the ends are bisected onto the exact place where the curve stops; pieces are joined only along the same ψ branch. The sky view uses each point's exact viewing direction (the reversed outgoing ray), never P − E, which loses all precision next to the eye. The surface mesh is drawn on the GPU through the chosen lens.
- **Two degrees of freedom.** At each h the lamp acts as a sun of elevation h, and each point of the chosen curve on the orientation plane is mapped onto the cone in closed form. The caustic Ĵ = J·√(exit margin) = 0 is traced by numerical continuation (predictor–corrector with step control), not assembled from a marching-squares grid: the grid only supplies seeds (sign changes on cell edges, and local minima of |J| for loops smaller than a cell). The tracer goes straight through crossings, steps over singular points where the gradient vanishes, bisects onto the exact boundary where it leaves the valid region and then solves the end exactly, so layer curves no longer stop part-way, retrace themselves or zigzag. The surface's edges come from continuing the parallel-light endpoints along h. Extra layers are placed by equidistributing the 3D spacing measured between the first-pass layers (more where layers are far apart, none where they crowd, e.g. near h = ±90° where the cone degenerates), and layers crowded closer than 0.3× the median spacing are not drawn.
- **Observer caustic.** χ(ξ) = det(ξ, F(ξ), DF⁻¹(ŵ×F) − ŵ×ξ) = 0 is a curve that does not depend on the lamp; the lamp only pairs the azimuth on the cone and decides what is visible.
- **Divergency shell.** On a shell the halo curve is the zero set of one quadratic Φ_ρ(w) = (S + ρe_z)² − 2ρ p u_z (S + ρe_z) + (ρ² − 1) u_z² on the crystal-frame sun sphere, cut by the path's boundaries; arcs are continued from their ends.

Nothing is estimated from random ray-tracing samples.

The engines in `src/engine/` are line-by-line JavaScript ports of the author's C# engines (the "Divergent Light Halos 3D" window of a private WinForms simulator, its divergency-shell and observer-caustic engines, and the light-point endpoint engine). They were cross-checked against the C# originals on 17 + 56 + 44 test cases (with the display refinement switched off): counts of points, curves, ends and flags are identical; coordinates agree to rounding (the remaining 1e-13 … 1e-6 m come from last-bit differences between the .NET and V8 sin/cos/atan2).

## Build from source

```bash
npm install
```

```bash
npm run build
```

This writes `dist/DivergentLightHalo.html` (one file with CSS, JavaScript, the Web Worker and three.js inlined) and a copy at `docs/index.html` for GitHub Pages. `npm run watch` rebuilds on change; `npm run serve` serves it at <http://127.0.0.1:5179>.

```
src/engine/   crystal.js (faces, paths, orientation)  cone.js (light-cone method, parallel arc, cone slice)  mesh1.js (surface mesh)
              bdry4.js (boundary engine v4)  epdiv.js (two degrees of freedom)  obscaustic.js (observer caustic)
              divshell.js (divergency shell)  lp.js / ep.js / aperture.js (light-point engine, shared with Light-Point-Trace)
src/worker.js computation worker (tasks are split across a worker pool)
src/app/      main.js (state, compute orchestration, panels), scene3d.js (three.js), sky.js (observer sky),
              crystal.js (clickable crystal), strings.js + terms.js (zh / en text and glossary)
```

Related: [Light-Point Halo Tracer](https://github.com/GalliumBoride/Light-Point-Trace) (parallel light, exact arcs, endpoints).

## License

[MIT](LICENSE)
