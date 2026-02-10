(* geodesic.ml - Geodesic integration via RK4 *)

open Expr [@@warning "-33"]
open Numeric

(* ========================================================================
   GEODESIC STATE
   
   State vector: (t, r, θ, φ, u^t, u^r, u^θ, u^φ)
   where u^μ = dx^μ/dτ (or dx^μ/dλ for null geodesics)
   ======================================================================== *)

type geodesic_state = {
  t : float;
  r : float;
  theta : float;
  phi : float;
  ut : float;   (* dt/dτ *)
  ur : float;   (* dr/dτ *)
  uth : float;  (* dθ/dτ *)
  uph : float;  (* dφ/dτ *)
}

let state_to_array s =
  [| s.t; s.r; s.theta; s.phi; s.ut; s.ur; s.uth; s.uph |]

let array_to_state a =
  { t = a.(0); r = a.(1); theta = a.(2); phi = a.(3);
    ut = a.(4); ur = a.(5); uth = a.(6); uph = a.(7) }

(* ========================================================================
   GEODESIC EQUATION
   
   d²x^μ/dτ² + Γ^μ_αβ (dx^α/dτ)(dx^β/dτ) = 0
   
   Or as first-order system:
   dx^μ/dτ = u^μ
   du^μ/dτ = -Γ^μ_αβ u^α u^β
   ======================================================================== *)

let geodesic_rhs params state =
  let r = state.r in
  let theta = state.theta in
  let ut = state.ut in
  let ur = state.ur in
  let uth = state.uth in
  let uph = state.uph in
  
  (* Get Christoffel symbols *)
  let (g_r_rr, g_r_thth, g_r_phph, g_th_rth, g_th_phph, g_ph_rph, g_ph_thph) =
    gamma_num params ~r ~theta in
  
  (* For Φ = 0: Γ^t_αβ = 0, so du^t/dτ = 0 *)
  let dut = 0.0 in
  
  (* du^r/dτ = -Γ^r_αβ u^α u^β *)
  let dur = 
    -. g_r_rr *. ur *. ur
    -. g_r_thth *. uth *. uth
    -. g_r_phph *. uph *. uph in
  
  (* du^θ/dτ = -Γ^θ_αβ u^α u^β *)
  let duth =
    -. 2.0 *. g_th_rth *. ur *. uth
    -. g_th_phph *. uph *. uph in
  
  (* du^φ/dτ = -Γ^φ_αβ u^α u^β *)
  let duph =
    -. 2.0 *. g_ph_rph *. ur *. uph
    -. 2.0 *. g_ph_thph *. uth *. uph in
  
  (* dx^μ/dτ = u^μ *)
  { t = ut; r = ur; theta = uth; phi = uph;
    ut = dut; ur = dur; uth = duth; uph = duph }

(* ========================================================================
   RK4 INTEGRATOR
   ======================================================================== *)

let add_state s1 s2 =
  { t = s1.t +. s2.t; r = s1.r +. s2.r; 
    theta = s1.theta +. s2.theta; phi = s1.phi +. s2.phi;
    ut = s1.ut +. s2.ut; ur = s1.ur +. s2.ur;
    uth = s1.uth +. s2.uth; uph = s1.uph +. s2.uph }

let scale_state c s =
  { t = c *. s.t; r = c *. s.r; 
    theta = c *. s.theta; phi = c *. s.phi;
    ut = c *. s.ut; ur = c *. s.ur;
    uth = c *. s.uth; uph = c *. s.uph }

let rk4_step params dt state =
  let k1 = geodesic_rhs params state in
  let k2 = geodesic_rhs params (add_state state (scale_state (dt /. 2.0) k1)) in
  let k3 = geodesic_rhs params (add_state state (scale_state (dt /. 2.0) k2)) in
  let k4 = geodesic_rhs params (add_state state (scale_state dt k3)) in
  
  let weighted = add_state (scale_state (1.0 /. 6.0) k1)
                (add_state (scale_state (2.0 /. 6.0) k2)
                (add_state (scale_state (2.0 /. 6.0) k3)
                           (scale_state (1.0 /. 6.0) k4))) in
  
  add_state state (scale_state dt weighted)

(* ========================================================================
   INTEGRATION WITH BOUNDARY CHECKING
   ======================================================================== *)

let integrate_geodesic params ~dt ~n_steps initial =
  let b0 = match params.shape with ConstantB b -> b | EllisB b -> b in
  let trajectory = Array.make (n_steps + 1) initial in
  let current = ref initial in
  
  for i = 1 to n_steps do
    (* Check if we're too close to throat *)
    if !current.r > b0 *. 1.001 then begin
      current := rk4_step params dt !current;
      (* Keep θ in [0, π] *)
      if !current.theta < 0.0 then
        current := { !current with theta = -. !current.theta; uph = -. !current.uph };
      if !current.theta > Float.pi then
        current := { !current with theta = 2.0 *. Float.pi -. !current.theta };
      (* Keep φ in [0, 2π] *)
      let phi' = mod_float !current.phi (2.0 *. Float.pi) in
      current := { !current with phi = if phi' < 0.0 then phi' +. 2.0 *. Float.pi else phi' }
    end;
    trajectory.(i) <- !current
  done;
  
  trajectory

(* ========================================================================
   INITIAL CONDITIONS
   ======================================================================== *)

(* Normalize 4-velocity for timelike geodesic: g_μν u^μ u^ν = -1 *)
let normalize_timelike params state =
  let r = state.r in
  let theta = state.theta in
  let b = eval_shape params r in
  
  (* g_tt = -1, g_rr = r/(r-b), g_θθ = r², g_φφ = r²sin²θ *)
  let g_rr = r /. (r -. b) in
  let g_thth = r *. r in
  let g_phph = r *. r *. (sin theta) ** 2.0 in
  
  (* -1 = g_tt (u^t)² + g_rr (u^r)² + g_θθ (u^θ)² + g_φφ (u^φ)² *)
  (* -1 = -(u^t)² + g_rr (u^r)² + g_θθ (u^θ)² + g_φφ (u^φ)² *)
  let spatial = g_rr *. state.ur *. state.ur 
              +. g_thth *. state.uth *. state.uth 
              +. g_phph *. state.uph *. state.uph in
  let ut_sq = 1.0 +. spatial in
  if ut_sq < 0.0 then failwith "Cannot normalize: spacelike trajectory"
  else { state with ut = sqrt ut_sq }

(* Create initial state for radial infall *)
let radial_infall params ~r0 ~ur0 =
  let state = { t = 0.0; r = r0; theta = Float.pi /. 2.0; phi = 0.0;
                ut = 0.0; ur = ur0; uth = 0.0; uph = 0.0 } in
  normalize_timelike params state

(* Create initial state for circular orbit (if exists) *)
let circular_orbit params ~r0 =
  let b = eval_shape params r0 in
  let _b' = eval_shape_deriv params r0 in
  
  (* For circular orbit: Γ^r_tt (u^t)² + Γ^r_φφ (u^φ)² = 0 *)
  (* With Φ = 0: Γ^r_tt = 0, so we need Γ^r_φφ = 0, which means r = b *)
  (* No stable circular orbits in simple wormhole! *)
  
  (* Instead, give tangential velocity *)
  let theta = Float.pi /. 2.0 in
  let _g_phph = r0 *. r0 in
  
  (* Pick u^φ such that orbit is bound *)
  let uph = 0.5 /. r0 in  (* Arbitrary choice *)
  
  let state = { t = 0.0; r = r0; theta; phi = 0.0;
                ut = 0.0; ur = 0.0; uth = 0.0; uph } in
  normalize_timelike params state

(* ========================================================================
   NULL GEODESICS (LIGHT RAYS)
   
   For null geodesics: g_μν k^μ k^ν = 0
   ======================================================================== *)

(* Normalize for null geodesic: g_μν k^μ k^ν = 0 *)
let normalize_null params state =
  let r = state.r in
  let theta = state.theta in
  let b = eval_shape params r in
  
  let g_rr = r /. (r -. b) in
  let g_thth = r *. r in
  let g_phph = r *. r *. (sin theta) ** 2.0 in
  
  (* 0 = -(k^t)² + g_rr (k^r)² + g_θθ (k^θ)² + g_φφ (k^φ)² *)
  let spatial = g_rr *. state.ur *. state.ur 
              +. g_thth *. state.uth *. state.uth 
              +. g_phph *. state.uph *. state.uph in
  { state with ut = sqrt spatial }

(* Radial light ray *)
let radial_light params ~r0 ~outgoing =
  let theta = Float.pi /. 2.0 in
  let b = eval_shape params r0 in
  
  (* For radial null geodesic: -(k^t)² + g_rr (k^r)² = 0 *)
  (* k^t = √(g_rr) |k^r| = √(r/(r-b)) |k^r| *)
  let kr = if outgoing then 1.0 else -1.0 in
  let kt = sqrt (r0 /. (r0 -. b)) *. abs_float kr in
  
  { t = 0.0; r = r0; theta; phi = 0.0;
    ut = kt; ur = kr; uth = 0.0; uph = 0.0 }

(* Light ray with impact parameter *)
let light_with_impact params ~r0 ~impact_param =
  let theta = Float.pi /. 2.0 in
  let b = eval_shape params r0 in
  
  (* Impact parameter ξ = L/E = r² (dφ/dt) / (dt/dτ) *)
  (* For distant ray: ξ ≈ r * sin(angle from radial) *)
  let kph = impact_param /. (r0 *. r0) in
  let kr = -1.0 in  (* Ingoing *)
  
  let state = { t = 0.0; r = r0; theta; phi = 0.0;
                ut = 0.0; ur = kr; uth = 0.0; uph = kph } in
  normalize_null params state

(* ========================================================================
   TRAJECTORY ANALYSIS
   ======================================================================== *)

(* Find turning points (where ur changes sign) *)
let find_turning_points trajectory =
  let n = Array.length trajectory in
  let points = ref [] in
  for i = 1 to n - 1 do
    let prev = trajectory.(i-1) in
    let curr = trajectory.(i) in
    if prev.ur *. curr.ur < 0.0 then
      points := curr.r :: !points
  done;
  List.rev !points

(* Check if trajectory passes through wormhole *)
let passes_through_wormhole params trajectory =
  let b0 = match params.shape with ConstantB b -> b | EllisB b -> b in
  Array.exists (fun s -> s.r < b0 *. 1.1) trajectory

(* Compute proper time elapsed *)
let proper_time_elapsed trajectory =
  if Array.length trajectory < 2 then 0.0
  else
    let last = trajectory.(Array.length trajectory - 1) in
    let first = trajectory.(0) in
    last.t -. first.t  (* For Φ=0, coordinate time ≈ proper time at infinity *)

(* ========================================================================
   OUTPUT FOR VISUALIZATION
   ======================================================================== *)

(* Convert trajectory to (x, y, z) for 3D plotting *)
let trajectory_to_xyz trajectory =
  Array.map (fun s ->
    let x = s.r *. sin s.theta *. cos s.phi in
    let y = s.r *. sin s.theta *. sin s.phi in
    let z = s.r *. cos s.theta in
    (x, y, z)
  ) trajectory

(* Convert to (r, z_embedding) for 2D embedding diagram *)
let trajectory_to_embedding params trajectory =
  Array.map (fun s ->
    let z = embedding_z params s.r in
    (s.r, z)
  ) trajectory
