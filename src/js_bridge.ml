(* js_bridge.ml - Bridge between OCaml tensor engine and JavaScript *)

open Js_of_ocaml

(* ========================================================================
   DEFORMATIONS STATE
   ======================================================================== *)

type deformation = { def_r: float; strength: float; radius: float }
let deformations : deformation list ref = ref []

(* ========================================================================
   SHAPE FUNCTION with deformations
   ======================================================================== *)

let shape_with_deformations b0 r =
  let base = b0 in
  let deform_contrib = List.fold_left (fun acc def ->
    let dist = abs_float (r -. def.def_r) in
    let influence = def.strength *. exp (-. dist *. dist /. (def.radius *. def.radius)) in
    acc +. influence
  ) 0.0 !deformations in
  max 0.2 (base +. deform_contrib)

let shape_deriv_with_deformations b0 r =
  let dr = 0.001 in
  (shape_with_deformations b0 (r +. dr) -. shape_with_deformations b0 (r -. dr)) /. (2.0 *. dr)

(* ========================================================================
   CHRISTOFFEL SYMBOLS
   ======================================================================== *)

let compute_christoffel_at ~r ~theta ~b0 =
  let b = shape_with_deformations b0 r in
  let bp = shape_deriv_with_deformations b0 r in
  let sin_th = sin theta in
  let cos_th = cos theta in
  let rmb = r -. b in
  
  let g_r_rr = 
    if abs_float rmb < 1e-10 then 0.0
    else (b -. r *. bp) /. (2.0 *. r *. rmb) in
  let g_r_thth = -. rmb in
  let g_r_phph = -. rmb *. sin_th *. sin_th in
  let g_th_rth = 1.0 /. r in
  let g_th_phph = -. sin_th *. cos_th in
  let g_ph_rph = 1.0 /. r in
  let g_ph_thph = 
    if abs_float sin_th < 1e-10 then 0.0 
    else cos_th /. sin_th in
  
  object%js
    val rRr = g_r_rr
    val rThth = g_r_thth
    val rPhph = g_r_phph
    val thRth = g_th_rth
    val thPhph = g_th_phph
    val phRph = g_ph_rph
    val phThph = g_ph_thph
  end

(* ========================================================================
   GEODESIC INTEGRATION (RK4)
   ======================================================================== *)

let geodesic_step_internal ~b0 ~dt (t, r, theta, phi, ut, ur, uth, uph) =
  let rhs (_, r, theta, _, ut, ur, uth, uph) =
    let r_safe = max r (b0 *. 1.001) in
    let b = shape_with_deformations b0 r_safe in
    let bp = shape_deriv_with_deformations b0 r_safe in
    let sin_th = sin theta in
    let cos_th = cos theta in
    let rmb = r_safe -. b in
    
    let g_r_rr = if abs_float rmb < 1e-10 then 0.0 else (b -. r_safe *. bp) /. (2.0 *. r_safe *. rmb) in
    let g_r_thth = -. rmb in
    let g_r_phph = -. rmb *. sin_th *. sin_th in
    let g_th_rth = 1.0 /. r_safe in
    let g_th_phph = -. sin_th *. cos_th in
    let g_ph_rph = 1.0 /. r_safe in
    let g_ph_thph = if abs_float sin_th < 1e-10 then 0.0 else cos_th /. sin_th in
    
    let dut = 0.0 in
    let dur = -. g_r_rr *. ur *. ur -. g_r_thth *. uth *. uth -. g_r_phph *. uph *. uph in
    let duth = -. 2.0 *. g_th_rth *. ur *. uth -. g_th_phph *. uph *. uph in
    let duph = -. 2.0 *. g_ph_rph *. ur *. uph -. 2.0 *. g_ph_thph *. uth *. uph in
    
    (ut, ur, uth, uph, dut, dur, duth, duph)
  in
  
  let add8 (a1,a2,a3,a4,a5,a6,a7,a8) (b1,b2,b3,b4,b5,b6,b7,b8) c =
    (a1+.c*.b1, a2+.c*.b2, a3+.c*.b3, a4+.c*.b4, a5+.c*.b5, a6+.c*.b6, a7+.c*.b7, a8+.c*.b8)
  in
  
  let state = (t, r, theta, phi, ut, ur, uth, uph) in
  let k1 = rhs state in
  let k2 = rhs (add8 state k1 (dt /. 2.0)) in
  let k3 = rhs (add8 state k2 (dt /. 2.0)) in
  let k4 = rhs (add8 state k3 dt) in
  
  let (dt1,dr1,dth1,dph1,dut1,dur1,duth1,duph1) = k1 in
  let (dt2,dr2,dth2,dph2,dut2,dur2,duth2,duph2) = k2 in
  let (dt3,dr3,dth3,dph3,dut3,dur3,duth3,duph3) = k3 in
  let (dt4,dr4,dth4,dph4,dut4,dur4,duth4,duph4) = k4 in
  
  let h = dt /. 6.0 in
  (
    t +. h *. (dt1 +. 2.0*.dt2 +. 2.0*.dt3 +. dt4),
    r +. h *. (dr1 +. 2.0*.dr2 +. 2.0*.dr3 +. dr4),
    theta +. h *. (dth1 +. 2.0*.dth2 +. 2.0*.dth3 +. dth4),
    phi +. h *. (dph1 +. 2.0*.dph2 +. 2.0*.dph3 +. dph4),
    ut +. h *. (dut1 +. 2.0*.dut2 +. 2.0*.dut3 +. dut4),
    ur +. h *. (dur1 +. 2.0*.dur2 +. 2.0*.dur3 +. dur4),
    uth +. h *. (duth1 +. 2.0*.duth2 +. 2.0*.duth3 +. duth4),
    uph +. h *. (duph1 +. 2.0*.duph2 +. 2.0*.duph3 +. duph4)
  )

let normalize_null_internal ~b0 (t, r, theta, phi, _, kr, kth, kph) =
  let b = shape_with_deformations b0 r in
  let sin_th = sin theta in
  
  let g_rr = r /. (r -. b) in
  let g_thth = r *. r in
  let g_phph = r *. r *. sin_th *. sin_th in
  
  let spatial = g_rr *. kr *. kr +. g_thth *. kth *. kth +. g_phph *. kph *. kph in
  let kt = sqrt (max 0.0 spatial) in
  (t, r, theta, phi, kt, kr, kth, kph)

let integrate_geodesic_js ~b0 ~dt ~n_steps ~t0 ~r0 ~theta0 ~phi0 ~kr ~kth ~kph =
  let initial = normalize_null_internal ~b0 (t0, r0, theta0, phi0, 0.0, kr, kth, kph) in
  
  let results = Array.make (n_steps + 1) initial in
  let current = ref initial in
  
  for i = 1 to n_steps do
    let (t, r, theta, phi, kt, kr, kth, kph) = !current in
    
    let r' = max r (b0 *. 1.001) in
    let theta' = max 0.01 (min (Float.pi -. 0.01) theta) in
    let phi' = mod_float phi (2.0 *. Float.pi) in
    let phi'' = if phi' < 0.0 then phi' +. 2.0 *. Float.pi else phi' in
    
    current := geodesic_step_internal ~b0 ~dt (t, r', theta', phi'', kt, kr, kth, kph);
    results.(i) <- !current
  done;
  
  let js_results = Array.map (fun (ti, ri, thetai, phii, kti, kri, kthi, kphi) ->
    object%js
      val t = ti
      val r = ri
      val theta = thetai
      val phi = phii
      val kt = kti
      val kr = kri
      val kth = kthi
      val kph = kphi
    end
  ) results in
  Js.array js_results

(* ========================================================================
   EMBEDDING DIAGRAM
   ======================================================================== *)

let embedding_z_internal ~b0 r =
  if r <= b0 then 0.0
  else
    let n_steps = 80 in
    let dr = (r -. b0 *. 1.001) /. float_of_int n_steps in
    let sum = ref 0.0 in
    for i = 0 to n_steps - 1 do
      let ri = b0 *. 1.001 +. (float_of_int i +. 0.5) *. dr in
      let bi = shape_with_deformations b0 ri in
      if ri > bi then
        sum := !sum +. sqrt (bi /. (ri -. bi)) *. dr
    done;
    !sum

(* ========================================================================
   DERIVATION STEPS
   ======================================================================== *)

let get_derivation_steps ~r ~theta ~b0 =
  let b = shape_with_deformations b0 r in
  let bp = shape_deriv_with_deformations b0 r in
  let sin_th = sin theta in
  let rmb = r -. b in
  
  let g_rr = if abs_float rmb < 1e-10 then infinity else r /. rmb in
  let g_thth = r *. r in
  let g_phph = r *. r *. sin_th *. sin_th in
  
  let gamma_r_rr = if abs_float rmb < 1e-10 then 0.0 else (b -. r *. bp) /. (2.0 *. r *. rmb) in
  let gamma_r_thth = -. rmb in
  let gamma_r_phph = -. rmb *. sin_th *. sin_th in
  let gamma_th_rth = 1.0 /. r in
  let gamma_ph_rph = 1.0 /. r in
  
  let gaussian_K = -. bp *. b /. (2.0 *. r *. r *. r) in
  let ricci_R = 2.0 *. bp /. (r *. r) in
  let time_dilation = if rmb <= 0.0 then infinity else sqrt (r /. rmb) in
  
  object%js
    val r = r
    val theta = theta
    val b = b
    val bp = bp
    val gRr = g_rr
    val gThth = g_thth
    val gPhph = g_phph
    val gammaRRr = gamma_r_rr
    val gammaRThth = gamma_r_thth
    val gammaRPhph = gamma_r_phph
    val gammaThRth = gamma_th_rth
    val gammaPhRph = gamma_ph_rph
    val gaussianK = gaussian_K
    val ricciR = ricci_R
    val timeDilation = time_dilation
  end

(* ========================================================================
   DEFORMATION MANAGEMENT
   ======================================================================== *)

let add_deformation ~r ~strength ~radius =
  deformations := { def_r = r; strength; radius } :: !deformations;
  if List.length !deformations > 10 then
    deformations := List.rev (List.tl (List.rev !deformations))

let clear_deformations () =
  deformations := []

(* ========================================================================
   REGISTER WITH JAVASCRIPT
   ======================================================================== *)

let () =
  Js.export "TensorEngine" (object%js
    method computeChristoffel r theta b0 =
      compute_christoffel_at ~r ~theta ~b0
    
    method integrateGeodesic b0 dt nSteps t0 r0 theta0 phi0 kr kth kph =
      integrate_geodesic_js
        ~b0 ~dt
        ~n_steps:(int_of_float nSteps)
        ~t0 ~r0 ~theta0 ~phi0 ~kr ~kth ~kph
    
    method getDerivationSteps r theta b0 =
      get_derivation_steps ~r ~theta ~b0
    
    method embeddingZ b0 r =
      embedding_z_internal ~b0 r
    
    method addDeformation r strength radius =
      add_deformation ~r ~strength ~radius
    
    method clearDeformations =
      clear_deformations ()
    
    method getShapeAt b0 r =
      shape_with_deformations b0 r
  end)
