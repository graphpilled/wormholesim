(* metric.ml - Metric tensor with actual computation *)

open Expr

(* A metric is a 4x4 matrix of expressions, indexed by coordinates *)
type metric = {
  components : expr array array;  (* g_μν - covariant, lower indices *)
  inverse : expr array array;     (* g^μν - contravariant, upper indices *)
}

let coord_to_idx = function
  | T -> 0 | R -> 1 | Theta -> 2 | Phi -> 3

let idx_to_coord = function
  | 0 -> T | 1 -> R | 2 -> Theta | 3 -> Phi
  | _ -> failwith "Invalid coordinate index"

(* Get metric component g_μν *)
let get_lower metric mu nu =
  metric.components.(coord_to_idx mu).(coord_to_idx nu)

(* Get inverse metric component g^μν *)  
let get_upper metric mu nu =
  metric.inverse.(coord_to_idx mu).(coord_to_idx nu)

(* ========================================================================
   MATRIX INVERSE FOR DIAGONAL METRICS
   
   For diagonal metric: g^μμ = 1/g_μμ, g^μν = 0 for μ≠ν
   
   This is exact for Morris-Thorne (and Schwarzschild, FLRW, etc.)
   General case would need symbolic determinant/cofactors.
   ======================================================================== *)

let invert_diagonal components =
  let inv = Array.make_matrix 4 4 zero in
  for i = 0 to 3 do
    for j = 0 to 3 do
      if i = j then
        inv.(i).(j) <- div one components.(i).(i)
      else
        inv.(i).(j) <- zero
    done
  done;
  inv

let is_diagonal components =
  let off_diag_zero = ref true in
  for i = 0 to 3 do
    for j = 0 to 3 do
      if i <> j && not (is_zero components.(i).(j)) then
        off_diag_zero := false
    done
  done;
  !off_diag_zero

(* Create metric from components, computing inverse *)
let make_metric components =
  if not (is_diagonal components) then
    failwith "Non-diagonal metrics not yet supported"
  else
    { components; inverse = invert_diagonal components }

(* ========================================================================
   MORRIS-THORNE WORMHOLE METRIC
   
   ds² = -e^{2Φ(r)} c² dt² + dr²/(1-b(r)/r) + r²dθ² + r²sin²θ dφ²
   
   We represent Φ(r) and b(r) as abstract functions that will be
   differentiated symbolically.
   ======================================================================== *)

let morris_thorne_metric ?(c=one) () =
  let r = Var "r" in
  let theta = Var "θ" in
  let phi_r = Func ("Φ", r) in    (* Redshift function *)
  let b_r = Func ("b", r) in      (* Shape function *)
  
  (* g_tt = -e^{2Φ(r)} c² *)
  let g_tt = neg (mul (Exp (mul two phi_r)) (mul c c)) in
  
  (* g_rr = 1/(1 - b(r)/r) = r/(r - b(r)) *)
  let g_rr = div r (sub r b_r) in
  
  (* g_θθ = r² *)
  let g_thth = mul r r in
  
  (* g_φφ = r² sin²θ *)
  let g_phph = mul (mul r r) (pow (Sin theta) two) in
  
  let components = Array.make_matrix 4 4 zero in
  components.(0).(0) <- g_tt;
  components.(1).(1) <- g_rr;
  components.(2).(2) <- g_thth;
  components.(3).(3) <- g_phph;
  
  make_metric components

(* Simplified: Φ = 0 (zero-tidal-force wormhole) *)
let simple_wormhole_metric ?(c=one) () =
  let r = Var "r" in
  let theta = Var "θ" in
  let b_r = Func ("b", r) in
  
  let g_tt = neg (mul c c) in
  let g_rr = div r (sub r b_r) in
  let g_thth = mul r r in
  let g_phph = mul (mul r r) (pow (Sin theta) two) in
  
  let components = Array.make_matrix 4 4 zero in
  components.(0).(0) <- g_tt;
  components.(1).(1) <- g_rr;
  components.(2).(2) <- g_thth;
  components.(3).(3) <- g_phph;
  
  make_metric components

(* Even simpler: constant b(r) = b₀ *)
let constant_b_wormhole ?(c=one) ?(b0=Var "b₀") () =
  let r = Var "r" in
  let theta = Var "θ" in
  
  let g_tt = neg (mul c c) in
  let g_rr = div r (sub r b0) in
  let g_thth = mul r r in
  let g_phph = mul (mul r r) (pow (Sin theta) two) in
  
  let components = Array.make_matrix 4 4 zero in
  components.(0).(0) <- g_tt;
  components.(1).(1) <- g_rr;
  components.(2).(2) <- g_thth;
  components.(3).(3) <- g_phph;
  
  make_metric components

(* ========================================================================
   PRINTING
   ======================================================================== *)

let print_metric metric =
  print_endline "Metric g_μν:";
  let coords = ["t"; "r"; "θ"; "φ"] in
  for i = 0 to 3 do
    for j = 0 to 3 do
      let g = metric.components.(i).(j) in
      if not (is_zero g) then
        Printf.printf "  g_{%s%s} = %s\n" 
          (List.nth coords i) (List.nth coords j) (to_string g)
    done
  done;
  print_endline "\nInverse metric g^μν:";
  for i = 0 to 3 do
    for j = 0 to 3 do
      let g = metric.inverse.(i).(j) in
      if not (is_zero g) then
        Printf.printf "  g^{%s%s} = %s\n"
          (List.nth coords i) (List.nth coords j) (to_string g)
    done
  done
