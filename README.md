# wormholesim

A symbolic tensor calculus engine written in OCaml that computes wormhole spacetime geometry from first principles. Derives Christoffel symbols, Riemann curvature, Ricci tensors, Einstein field equations, and geodesics — all from the metric tensor, with no hardcoded results. Includes a real-time 3D visualization in the browser via js_of_ocaml.

**[Live Demo →](https://graphpilled.github.io/wormholesim/wormhole-viz/index.html)**

---

## How it works

You give it a metric tensor. It does everything else.

```
metric.ml → christoffel.ml → riemann.ml → geodesic.ml → wormhole-viz/
 (g_μν)        (Γ^σ_μν)       (R^ρ_σμν)     (RK4)        (Three.js)
```

1. **Define the metric** — Morris-Thorne wormhole: `ds² = -e^{2Φ(r)}dt² + dr²/(1-b(r)/r) + r²dΩ²`
2. **Differentiate symbolically** — `expr.ml` performs exact symbolic differentiation with product rule, chain rule, quotient rule on an algebraic expression tree
3. **Compute Christoffels** — `Γ^σ_μν = ½ g^{σρ}(∂_μ g_{νρ} + ∂_ν g_{μρ} - ∂_ρ g_{μν})` by iterating over all index combinations, differentiating metric components, contracting with the inverse metric
4. **Compute Riemann tensor** — `R^ρ_σμν = ∂_μ Γ^ρ_νσ - ∂_ν Γ^ρ_μσ + Γ^ρ_μλ Γ^λ_νσ - Γ^ρ_νλ Γ^λ_μσ` by differentiating Christoffels and contracting products
5. **Contract to Ricci & Einstein** — `R_μν = R^ρ_μρν`, `R = g^μν R_μν`, `G_μν = R_μν - ½g_μν R`
6. **Integrate geodesics** — 4th-order Runge-Kutta on `d²x^μ/dτ² + Γ^μ_αβ u^α u^β = 0` for timelike and null paths
7. **Render in 3D** — The OCaml engine compiles to JavaScript via js_of_ocaml. Three.js renders the embedding diagram, geodesic trajectories, and KaTeX displays the live symbolic math.

## Architecture

### OCaml symbolic engine (`src/`)

| File | Lines | Description |
|---|---|---|
| `expr.ml` | 259 | Expression AST with symbolic differentiation — coordinates, variables, functions, trig, chain rule, product rule, LaTeX output |
| `simplify.ml` | 367 | CAS-style simplification via polynomial canonical forms, rational function cancellation, GCD, numerical zero detection |
| `derivation.ml` | 501 | Step-by-step derivation tracer — shows every algebraic manipulation with substitution, contraction, and simplification steps. JSON export for web |
| `metric.ml` | 159 | Metric tensor `g_μν` / `g^μν` with diagonal inverse. Morris-Thorne, Ellis, and constant-b wormhole constructors |
| `christoffel.ml` | 109 | Christoffel symbols from definition — symbolic metric derivatives contracted with inverse metric |
| `riemann.ml` | 221 | Riemann tensor `R^ρ_σμν`, Ricci tensor `R_μν`, Ricci scalar `R`, Einstein tensor `G_μν` — all computed by differentiating and contracting |
| `geodesic.ml` | 278 | RK4 geodesic integrator — timelike and null geodesics, 4-velocity normalization, impact parameter setup, embedding coordinates |
| `numeric.ml` | 255 | Numerical evaluation of symbolic expressions, wormhole parameter types (constant-b, Ellis), Christoffel numerical evaluation |
| `js_bridge.ml` | 263 | js_of_ocaml bridge — exposes `TensorCalc` API to JavaScript: `computeChristoffels`, `computeRiemann`, `traceGeodesic`, `getEmbeddingPoints`, `computeDerivation` |
| `main.ml` | 265 | CLI demo — computes full tensor hierarchy and prints results |

### Web visualization (`wormhole-viz/`)

| File | Lines | Description |
|---|---|---|
| `tensor_engine.js` | 24,582 | js_of_ocaml compiled output — the entire OCaml engine running in the browser |
| `wormhole-sim.js` | 1,050 | Three.js simulation — wormhole embedding surface, geodesic particle system, accretion disk, photon sphere |
| `wormhole-viz.js` | 502 | UI and rendering — parameter controls, KaTeX math display, embedding diagram generation, orbit camera |
| `index.html` | — | Application shell with Three.js, KaTeX CDN imports |

### Expression types

The symbolic engine represents all of general relativity's tensor algebra as an OCaml algebraic type:

```ocaml
type expr =
  | Num of float              (* Constants *)
  | Var of string             (* Variables: r, θ, φ *)
  | Func of string * expr     (* Functions: Φ(r), b(r) *)
  | Add | Sub | Mul | Div | Pow  (* Arithmetic *)
  | Exp | Log | Sin | Cos | Sqrt (* Transcendentals *)
  | Delta of coord * coord    (* Kronecker δ^μ_ν *)
  | Partial of expr * coord   (* Symbolic partial derivative *)
```

Smart constructors simplify at construction time (`0 * x → 0`, `x^0 → 1`, `e^(ln x) → x`). Full simplification uses polynomial canonical forms with rational function cancellation.

## Wormhole configurations

**Morris-Thorne** (general): `b(r)` and `Φ(r)` left as abstract symbolic functions — all results contain `b'(r)`, `Φ'(r)` as symbolic derivatives.

**Zero-tidal-force** (`Φ = 0`): Eliminates redshift terms. The metric simplifies to `ds² = -dt² + dr²/(1-b(r)/r) + r²dΩ²`.

**Constant-b**: `b(r) = b₀`. The throat radius is `b₀` and the shape function has zero derivative.

**Ellis drainhole**: `b(r) = b₀²/r`. Produces a different embedding curvature and geodesic behavior than constant-b.

## Building

### Native (CLI)

Requires OCaml and dune:

```bash
./build.sh
./tensorcalc
```

### JavaScript (browser)

Requires js_of_ocaml:

```bash
./build_js.sh
# Open wormhole-viz/index.html in a browser
```

## Project structure

```
src/
├── expr.ml            # Symbolic expression AST + differentiation
├── simplify.ml        # CAS simplification engine
├── derivation.ml      # Step-by-step derivation tracer
├── metric.ml          # Metric tensor definitions
├── christoffel.ml     # Γ^σ_μν computation
├── riemann.ml         # R^ρ_σμν → R_μν → R → G_μν
├── geodesic.ml        # RK4 geodesic integration
├── numeric.ml         # Numerical evaluation
├── js_bridge.ml       # OCaml → JavaScript bridge
├── main.ml            # CLI entry point
└── dune               # Build configuration

wormhole-viz/
├── index.html         # Application shell
├── tensor_engine.js   # Compiled OCaml engine (js_of_ocaml)
├── wormhole-sim.js    # Three.js wormhole simulation
└── wormhole-viz.js    # UI, controls, KaTeX rendering
```

## License

MIT
