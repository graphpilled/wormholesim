(* numeric.ml - Numerical evaluation of symbolic expressions *)

open Expr

(* ========================================================================
   ENVIRONMENT: Maps variable names to values
   ======================================================================== *)

type env = (string * float) list

let empty_env : env = []

let bind name value (env : env) : env = (name, value) :: env

let lookup name (env : env) =
  try Some (List.assoc name env)
  with Not_found -> None

(* Standard wormhole environment *)
let wormhole_env ~b0 ~r ~theta ~phi =
  [("b₀", b0); ("r", r); ("θ", theta); ("φ", phi); ("t", 0.0)]

(* ========================================================================
   FUNCTION EVALUATION
   
   For b(r), Φ(r) etc, we need concrete implementations
   ======================================================================== *)

type shape_fn = 
  | ConstantB of float      (* b(r) = b₀ *)
  | EllisB of float         (* b(r) = b₀²/r *)

type redshift_fn =
  | ZeroPhi                 (* Φ(r) = 0 *)
  | ConstantPhi of float    (* Φ(r) = Φ₀ *)

type wormhole_params = {
  shape : shape_fn;
  redshift : redshift_fn;
}

let default_wormhole b0 = {
  shape = ConstantB b0;
  redshift = ZeroPhi;
}

let ellis_wormhole b0 = {
  shape = EllisB b0;
  redshift = ZeroPhi;
}

(* Evaluate shape function and derivatives *)
let eval_shape params r =
  match params.shape with
  | ConstantB b0 -> b0
  | EllisB b0 -> b0 *. b0 /. r

let eval_shape_deriv params r =
  match params.shape with
  | ConstantB _ -> 0.0
  | EllisB b0 -> -. b0 *. b0 /. (r *. r)

(* Evaluate redshift function and derivatives *)
let eval_redshift params _r =
  match params.redshift with
  | ZeroPhi -> 0.0
  | ConstantPhi phi0 -> phi0

let eval_redshift_deriv params _ =
  match params.redshift with
  | ZeroPhi -> 0.0
  | ConstantPhi _ -> 0.0

(* ========================================================================
   EXPRESSION EVALUATION
   ======================================================================== *)

let rec eval (env : env) (params : wormhole_params) expr : float =
  match expr with
  | Num n -> n
  | Var name ->
      (match lookup name env with
       | Some v -> v
       | None -> failwith (Printf.sprintf "Unbound variable: %s" name))
  
  | Func (name, arg) ->
      let arg_val = eval env params arg in
      (match name with
       | "b" -> eval_shape params arg_val
       | "b'" -> eval_shape_deriv params arg_val
       | "Φ" -> eval_redshift params arg_val
       | "Φ'" -> eval_redshift_deriv params arg_val
       | "sin" -> sin arg_val
       | "cos" -> cos arg_val
       | "exp" -> exp arg_val
       | "ln" -> log arg_val
       | "sqrt" -> sqrt arg_val
       | _ -> failwith (Printf.sprintf "Unknown function: %s" name))
  
  | Neg e -> -. (eval env params e)
  | Add (a, b) -> eval env params a +. eval env params b
  | Sub (a, b) -> eval env params a -. eval env params b
  | Mul (a, b) -> eval env params a *. eval env params b
  | Div (a, b) -> 
      let denom = eval env params b in
      if abs_float denom < 1e-15 then nan
      else eval env params a /. denom
  | Pow (a, b) -> (eval env params a) ** (eval env params b)
  | Exp e -> exp (eval env params e)
  | Log e -> log (eval env params e)
  | Sin e -> sin (eval env params e)
  | Cos e -> cos (eval env params e)
  | Sqrt e -> sqrt (eval env params e)
  | Delta (a, b) -> if a = b then 1.0 else 0.0
  | Partial _ -> failwith "Cannot numerically evaluate unevaluated Partial"

(* Convenience: evaluate with standard wormhole params *)
let eval_at ~b0 ~r ~theta expr =
  let env = wormhole_env ~b0 ~r ~theta ~phi:0.0 in
  let params = default_wormhole b0 in
  eval env params expr

(* ========================================================================
   METRIC NUMERICAL EVALUATION
   ======================================================================== *)

open Metric

let eval_metric_component metric env params mu nu =
  let expr = get_lower metric mu nu in
  eval env params expr

let eval_inverse_metric metric env params mu nu =
  let expr = get_upper metric mu nu in
  eval env params expr

(* Get full metric as 4x4 float array *)
let eval_metric_matrix metric env params =
  Array.init 4 (fun i ->
    Array.init 4 (fun j ->
      let mu = idx_to_coord i in
      let nu = idx_to_coord j in
      eval_metric_component metric env params mu nu
    )
  )

(* ========================================================================
   CHRISTOFFEL NUMERICAL EVALUATION
   ======================================================================== *)

open Christoffel

let eval_christoffel cs env params sigma mu nu =
  let expr = get_christoffel cs sigma mu nu in
  eval env params expr

(* Get all Christoffels as float array *)
let eval_christoffel_array cs env params =
  Array.init 4 (fun s ->
    Array.init 4 (fun m ->
      Array.init 4 (fun n ->
        let sigma = idx_to_coord s in
        let mu = idx_to_coord m in
        let nu = idx_to_coord n in
        eval_christoffel cs env params sigma mu nu
      )
    )
  )

(* ========================================================================
   DIRECT NUMERICAL CHRISTOFFELS (for efficiency in geodesic integration)
   
   These compute Γ directly from parameters without symbolic evaluation
   ======================================================================== *)

let gamma_num params ~r ~theta =
  let b = eval_shape params r in
  let b' = eval_shape_deriv params r in
  let sin_th = sin theta in
  let cos_th = cos theta in
  
  (* Γ^r_rr = (b - r*b') / (2*r*(r-b)) *)
  let g_r_rr = 
    if abs_float (r -. b) < 1e-10 then nan
    else (b -. r *. b') /. (2.0 *. r *. (r -. b)) in
  
  (* Γ^r_θθ = -(r - b) *)
  let g_r_thth = -. (r -. b) in
  
  (* Γ^r_φφ = -(r - b) * sin²θ *)
  let g_r_phph = -. (r -. b) *. sin_th *. sin_th in
  
  (* Γ^θ_rθ = 1/r *)
  let g_th_rth = 1.0 /. r in
  
  (* Γ^θ_φφ = -sinθ cosθ *)
  let g_th_phph = -. sin_th *. cos_th in
  
  (* Γ^φ_rφ = 1/r *)
  let g_ph_rph = 1.0 /. r in
  
  (* Γ^φ_θφ = cosθ/sinθ *)
  let g_ph_thph = if abs_float sin_th < 1e-10 then nan else cos_th /. sin_th in
  
  (* Return as record for clarity *)
  (g_r_rr, g_r_thth, g_r_phph, g_th_rth, g_th_phph, g_ph_rph, g_ph_thph)

(* ========================================================================
   EMBEDDING DIAGRAM
   
   z(r) for visualization in 3D
   ======================================================================== *)

let embedding_z params r =
  match params.shape with
  | ConstantB b0 ->
      if r < b0 then 0.0
      else 2.0 *. sqrt (b0 *. (r -. b0))
  | EllisB b0 ->
      if r < b0 then 0.0
      else b0 *. acosh (r /. b0)

(* Generate embedding curve data *)
let embedding_curve params r_max n_points =
  let b0 = match params.shape with ConstantB b -> b | EllisB b -> b in
  let dr = (r_max -. b0) /. float_of_int (n_points - 1) in
  Array.init n_points (fun i ->
    let r = b0 +. float_of_int i *. dr in
    let z = embedding_z params r in
    (r, z)
  )

(* ========================================================================
   PROPER DISTANCE (numerical integration)
   ======================================================================== *)

let proper_distance params r1 r2 n_steps =
  (* ℓ = ∫ dr / √(1 - b(r)/r) *)
  let dr = (r2 -. r1) /. float_of_int n_steps in
  let sum = ref 0.0 in
  for i = 0 to n_steps - 1 do
    let r = r1 +. (float_of_int i +. 0.5) *. dr in
    let b = eval_shape params r in
    let integrand = 1.0 /. sqrt (1.0 -. b /. r) in
    sum := !sum +. integrand *. dr
  done;
  !sum

(* ========================================================================
   LIGHT TRAVEL TIME
   ======================================================================== *)

let light_travel_time params r1 r2 n_steps =
  (* Δt = ∫ dr / (c * √(1 - b/r)), with c = 1 *)
  proper_distance params r1 r2 n_steps  (* Same integral for Φ = 0 *)
