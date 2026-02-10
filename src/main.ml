(* main.ml - Complete demonstration of TensorCalc v2 *)

open Expr
open Simplify
open Metric
open Christoffel
open Riemann
open Numeric
open Geodesic
open Derivation

let () =
  print_endline "";
  print_endline "╔══════════════════════════════════════════════════════════════════╗";
  print_endline "║     TensorCalc v2 - ACTUAL COMPUTATION FROM FIRST PRINCIPLES     ║";
  print_endline "╚══════════════════════════════════════════════════════════════════╝";
  print_endline "";
  
  (* ====================================================================== *)
  print_endline "═══════════════════════════════════════════════════════════════════";
  print_endline "                    PART 1: SYMBOLIC COMPUTATION";
  print_endline "═══════════════════════════════════════════════════════════════════";
  print_endline "";
  
  print_endline "Input: Morris-Thorne metric (Φ = 0)";
  print_endline "  ds² = -dt² + dr²/(1-b(r)/r) + r²dθ² + r²sin²θ dφ²";
  print_endline "";
  
  let metric = simple_wormhole_metric ~c:one () in
  
  print_endline "Step 1: Metric tensor g_μν";
  print_endline "──────────────────────────";
  print_metric metric;
  print_endline "";
  
  print_endline "Step 2: Compute Christoffels from Γ = (1/2)g⁻¹(∂g + ∂g - ∂g)";
  print_endline "─────────────────────────────────────────────────────────────";
  let cs = compute_christoffels metric in
  print_endline "Non-zero Christoffel symbols:";
  let coord_str = [|"t"; "r"; "θ"; "φ"|] in
  for s = 0 to 3 do
    for m = 0 to 3 do
      for n = m to 3 do
        let g = cs.gamma.(s).(m).(n) in
        let gs = simplify_complete g in
        if not (is_zero gs) then
          Printf.printf "  Γ^%s_{%s%s} = %s\n"
            coord_str.(s) coord_str.(m) coord_str.(n) (to_string gs)
      done
    done
  done;
  print_endline "";
  
  print_endline "Step 3: Compute Riemann from R = ∂Γ - ∂Γ + ΓΓ - ΓΓ";
  print_endline "────────────────────────────────────────────────────";
  let rt = compute_riemann cs in
  print_endline "Non-zero Riemann components (showing key ones):";
  (* Just show a few key components to avoid overwhelming output *)
  let r_rthrth = get_riemann rt R Theta R Theta in
  let r_rthrth_s = simplify_complete r_rthrth in
  Printf.printf "  R^r_{θrθ} = %s\n" (to_string r_rthrth_s);
  
  let r_thrthr = get_riemann rt Theta R Theta R in
  let r_thrthr_s = simplify_complete r_thrthr in
  Printf.printf "  R^θ_{rθr} = %s\n" (to_string r_thrthr_s);
  print_endline "";
  
  print_endline "Step 4: Ricci tensor R_μν = R^ρ_μρν";
  print_endline "───────────────────────────────────";
  let rc = compute_ricci rt in
  for mi = 0 to 3 do
    for ni = mi to 3 do
      let component = rc.ric.(mi).(ni) in
      let cs = simplify_complete component in
      if not (is_zero cs) then
        Printf.printf "  R_{%s%s} = %s\n"
          coord_str.(mi) coord_str.(ni) (to_string cs)
    done
  done;
  print_endline "";
  
  print_endline "Step 5: Ricci scalar R = g^μν R_μν";
  print_endline "──────────────────────────────────";
  let scalar = ricci_scalar metric rc in
  let scalar_s = simplify_complete scalar in
  Printf.printf "  R = %s\n" (to_string scalar_s);
  
  (* Numerical verification *)
  print_endline "";
  print_endline "  Numerical verification (should be ≈ 0):";
  let env = Numeric.wormhole_env ~b0:1.0 ~r:2.0 ~theta:(Float.pi /. 2.0) ~phi:0.0 in
  let params = Numeric.default_wormhole 1.0 in
  let r_numeric = Numeric.eval env params scalar in
  Printf.printf "    R at (r=2, θ=π/2) = %.10f\n" r_numeric;
  print_endline "";
  
  (* ====================================================================== *)
  print_endline "═══════════════════════════════════════════════════════════════════";
  print_endline "                   PART 2: NUMERICAL EVALUATION";
  print_endline "═══════════════════════════════════════════════════════════════════";
  print_endline "";
  
  let params = default_wormhole 1.0 in
  
  print_endline "Wormhole parameters: b₀ = 1.0, Φ = 0";
  print_endline "";
  
  print_endline "Christoffel symbols at r = 2.0, θ = π/2:";
  print_endline "─────────────────────────────────────────";
  let (g_r_rr, g_r_thth, g_r_phph, g_th_rth, g_th_phph, g_ph_rph, g_ph_thph) = 
    gamma_num params ~r:2.0 ~theta:(Float.pi /. 2.0) in
  Printf.printf "  Γ^r_{rr} = %.6f\n" g_r_rr;
  Printf.printf "  Γ^r_{θθ} = %.6f\n" g_r_thth;
  Printf.printf "  Γ^r_{φφ} = %.6f\n" g_r_phph;
  Printf.printf "  Γ^θ_{rθ} = %.6f\n" g_th_rth;
  Printf.printf "  Γ^θ_{φφ} = %.6f\n" g_th_phph;
  Printf.printf "  Γ^φ_{rφ} = %.6f\n" g_ph_rph;
  Printf.printf "  Γ^φ_{θφ} = %.6f\n" g_ph_thph;
  print_endline "";
  
  print_endline "Embedding diagram z(r):";
  print_endline "───────────────────────";
  let test_r = [1.0; 1.5; 2.0; 3.0; 5.0] in
  List.iter (fun r ->
    let z = embedding_z params r in
    Printf.printf "  r = %.1f → z = %.4f\n" r z
  ) test_r;
  print_endline "";
  
  print_endline "Proper distance from throat:";
  print_endline "────────────────────────────";
  List.iter (fun r ->
    let d = proper_distance params 1.001 r 1000 in
    Printf.printf "  ℓ(1 → %.1f) = %.4f\n" r d
  ) [1.5; 2.0; 3.0; 5.0];
  print_endline "";
  
  (* ====================================================================== *)
  print_endline "═══════════════════════════════════════════════════════════════════";
  print_endline "                   PART 3: GEODESIC INTEGRATION";
  print_endline "═══════════════════════════════════════════════════════════════════";
  print_endline "";
  
  print_endline "Radial infall from r = 5.0:";
  print_endline "───────────────────────────";
  let initial = radial_infall params ~r0:5.0 ~ur0:(-0.1) in
  Printf.printf "  Initial: r = %.2f, u^r = %.4f, u^t = %.4f\n" 
    initial.r initial.ur initial.ut;
  
  let trajectory = integrate_geodesic params ~dt:0.1 ~n_steps:100 initial in
  let final = trajectory.(100) in
  Printf.printf "  After 100 steps: r = %.4f, θ = %.4f\n" final.r final.theta;
  
  let turning = find_turning_points trajectory in
  if List.length turning > 0 then
    Printf.printf "  Turning points at r = %s\n" 
      (String.concat ", " (List.map (Printf.sprintf "%.3f") turning))
  else
    print_endline "  No turning points (passes through wormhole)";
  print_endline "";
  
  print_endline "Radial light ray (outgoing from r = 2.0):";
  print_endline "─────────────────────────────────────────";
  let light = radial_light params ~r0:2.0 ~outgoing:true in
  Printf.printf "  Initial: r = %.2f, k^r = %.4f, k^t = %.4f\n"
    light.r light.ur light.ut;
  
  let light_traj = integrate_geodesic params ~dt:0.05 ~n_steps:100 light in
  let light_final = light_traj.(100) in
  Printf.printf "  After 100 steps: r = %.4f\n" light_final.r;
  print_endline "";
  
  print_endline "Light ray with impact parameter (b = 3.0):";
  print_endline "──────────────────────────────────────────";
  let bent_light = light_with_impact params ~r0:10.0 ~impact_param:3.0 in
  Printf.printf "  Initial: r = %.2f, k^r = %.4f, k^φ = %.6f\n"
    bent_light.r bent_light.ur bent_light.uph;
  
  let bent_traj = integrate_geodesic params ~dt:0.1 ~n_steps:200 bent_light in
  let bent_final = bent_traj.(200) in
  Printf.printf "  After 200 steps: r = %.4f, φ = %.4f rad\n" 
    bent_final.r bent_final.phi;
  
  let bent_turning = find_turning_points bent_traj in
  if List.length bent_turning > 0 then
    Printf.printf "  Closest approach: r = %.4f\n" (List.hd bent_turning);
  print_endline "";
  
  (* ====================================================================== *)
  print_endline "═══════════════════════════════════════════════════════════════════";
  print_endline "                         VERIFICATION";
  print_endline "═══════════════════════════════════════════════════════════════════";
  print_endline "";
  print_endline "This system ACTUALLY COMPUTES (not looks up):";
  print_endline "";
  print_endline "  Symbolic:";
  print_endline "    ∂g/∂x^μ  →  via symbolic differentiation (expr.ml)";
  print_endline "    Γ^σ_μν   →  via definition + contraction (christoffel.ml)";
  print_endline "    R^ρ_σμν  →  via ∂Γ - ∂Γ + ΓΓ - ΓΓ (riemann.ml)";
  print_endline "    R_μν, R  →  via contraction (riemann.ml)";
  print_endline "";
  print_endline "  Numerical:";
  print_endline "    Metric components  →  eval at (r, θ, φ)";
  print_endline "    Christoffels       →  direct formulas for efficiency";
  print_endline "    Geodesics          →  RK4 integration";
  print_endline "    Proper distance    →  numerical integration";
  print_endline "    Embedding          →  analytic for constant/Ellis b(r)";
  print_endline "";
  print_endline "No lookup tables. No pre-encoded tensor components.";
  print_endline "The differentiate function in expr.ml is the engine.";
  print_endline "";
  
  (* ====================================================================== *)
  print_endline "═══════════════════════════════════════════════════════════════════";
  print_endline "                     COMPARISON: CONSTANT vs ELLIS";
  print_endline "═══════════════════════════════════════════════════════════════════";
  print_endline "";
  
  let const_params = default_wormhole 1.0 in
  let ellis_params = Numeric.ellis_wormhole 1.0 in
  
  print_endline "                    Constant b(r)=b₀    Ellis b(r)=b₀²/r";
  print_endline "                    ─────────────────   ─────────────────";
  
  List.iter (fun r ->
    let z_const = embedding_z const_params r in
    let z_ellis = embedding_z ellis_params r in
    Printf.printf "  z(%.1f):          %.4f              %.4f\n" r z_const z_ellis
  ) [1.5; 2.0; 3.0];
  
  print_endline "";
  
  List.iter (fun r ->
    let d_const = proper_distance const_params 1.001 r 1000 in
    let d_ellis = proper_distance ellis_params 1.001 r 1000 in
    Printf.printf "  ℓ(1→%.1f):         %.4f              %.4f\n" r d_const d_ellis
  ) [2.0; 3.0; 5.0];
  
  print_endline "";
  
  (* ====================================================================== *)
  print_endline "═══════════════════════════════════════════════════════════════════";
  print_endline "              PART 4: STEP-BY-STEP DERIVATION TRACING";
  print_endline "═══════════════════════════════════════════════════════════════════";
  print_endline "";
  
  print_endline "Showing how Γ^r_{rr} is computed step by step:";
  print_endline "";
  
  let deriv_christoffel = Derivation.trace_christoffel_component metric R R R in
  Derivation.print_derivation deriv_christoffel;
  
  (* ====================================================================== *)
  print_endline "═══════════════════════════════════════════════════════════════════";
  print_endline "              PART 5: PARAMETER CHANGE ANALYSIS";
  print_endline "═══════════════════════════════════════════════════════════════════";
  
  let from_config = Derivation.constant_wormhole 1.0 in
  let to_config = Derivation.ellis_wormhole 1.0 in
  
  Derivation.trace_parameter_change ~from_config ~to_config metric cs;
  
  print_endline "═══════════════════════════════════════════════════════════════════";
  print_endline "                            COMPLETE";
  print_endline "═══════════════════════════════════════════════════════════════════"
