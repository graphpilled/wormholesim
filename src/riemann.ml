(* riemann.ml - Riemann tensor computed from Christoffel symbols *)

open Expr
open Metric
open Christoffel

(* ========================================================================
   RIEMANN TENSOR - ACTUALLY COMPUTED
   
   R^ρ_σμν = ∂_μ Γ^ρ_νσ - ∂_ν Γ^ρ_μσ + Γ^ρ_μλ Γ^λ_νσ - Γ^ρ_νλ Γ^λ_μσ
   
   This requires:
   1. Differentiating Christoffel symbols
   2. Contracting products of Christoffels
   ======================================================================== *)

(* Compute ∂_μ Γ^ρ_νσ *)
let christoffel_derivative cs mu rho nu sigma =
  let gamma = get_christoffel cs rho nu sigma in
  differentiate mu gamma |> full_simplify

(* Compute Γ^ρ_μλ Γ^λ_νσ (sum over λ) *)
let christoffel_product cs rho mu nu sigma =
  let sum = ref zero in
  List.iter (fun lambda ->
    let g1 = get_christoffel cs rho mu lambda in
    let g2 = get_christoffel cs lambda nu sigma in
    if not (is_zero g1) && not (is_zero g2) then
      sum := add !sum (mul g1 g2)
  ) all_coords;
  full_simplify !sum

(* Compute R^ρ_σμν from definition *)
let riemann cs rho sigma mu nu =
  (* R^ρ_σμν = ∂_μ Γ^ρ_νσ - ∂_ν Γ^ρ_μσ + Γ^ρ_μλ Γ^λ_νσ - Γ^ρ_νλ Γ^λ_μσ *)
  
  let term1 = christoffel_derivative cs mu rho nu sigma in   (* ∂_μ Γ^ρ_νσ *)
  let term2 = christoffel_derivative cs nu rho mu sigma in   (* ∂_ν Γ^ρ_μσ *)
  let term3 = christoffel_product cs rho mu nu sigma in      (* Γ^ρ_μλ Γ^λ_νσ *)
  let term4 = christoffel_product cs rho nu mu sigma in      (* Γ^ρ_νλ Γ^λ_μσ *)
  
  sub (add (sub term1 term2) term3) term4 |> full_simplify

(* Store computed Riemann tensor *)
type riemann_tensor = {
  r : expr array array array array;  (* R^ρ_σμν indexed [ρ][σ][μ][ν] *)
}

(* Compute all Riemann components *)
let compute_riemann cs =
  let r = Array.init 4 (fun _ ->
    Array.init 4 (fun _ ->
      Array.init 4 (fun _ ->
        Array.make 4 zero
      )
    )
  ) in
  
  (* Use symmetries to reduce computation:
     R^ρ_σμν = -R^ρ_σνμ  (antisymmetric in last two indices)
     R^ρ_σμν = -R^ρ_μσν + ... (more complex symmetries)
     
     For now, compute all and let simplification handle zeros *)
  List.iter (fun rho ->
    List.iter (fun sigma ->
      List.iter (fun mu ->
        List.iter (fun nu ->
          let ri = coord_to_idx rho in
          let si = coord_to_idx sigma in
          let mi = coord_to_idx mu in
          let ni = coord_to_idx nu in
          
          (* Use antisymmetry: R^ρ_σμν = -R^ρ_σνμ *)
          if mi < ni then begin
            let component = riemann cs rho sigma mu nu in
            r.(ri).(si).(mi).(ni) <- component;
            r.(ri).(si).(ni).(mi) <- neg component
          end else if mi = ni then
            r.(ri).(si).(mi).(ni) <- zero  (* Antisymmetric → diagonal = 0 *)
        ) all_coords
      ) all_coords
    ) all_coords
  ) all_coords;
  
  { r }

let get_riemann rt rho sigma mu nu =
  rt.r.(coord_to_idx rho).(coord_to_idx sigma).(coord_to_idx mu).(coord_to_idx nu)

(* ========================================================================
   RICCI TENSOR - Contract Riemann
   
   R_μν = R^ρ_μρν (sum over ρ)
   ======================================================================== *)

let ricci rt mu nu =
  let sum = ref zero in
  List.iter (fun rho ->
    let component = get_riemann rt rho mu rho nu in
    sum := add !sum component
  ) all_coords;
  full_simplify !sum

type ricci_tensor = {
  ric : expr array array;  (* R_μν *)
}

let compute_ricci rt =
  let ric = Array.make_matrix 4 4 zero in
  List.iter (fun mu ->
    List.iter (fun nu ->
      let mi = coord_to_idx mu in
      let ni = coord_to_idx nu in
      if mi <= ni then begin
        let component = ricci rt mu nu in
        ric.(mi).(ni) <- component;
        ric.(ni).(mi) <- component  (* Symmetric *)
      end
    ) all_coords
  ) all_coords;
  { ric }

let get_ricci rc mu nu =
  rc.ric.(coord_to_idx mu).(coord_to_idx nu)

(* ========================================================================
   RICCI SCALAR - Contract Ricci with metric
   
   R = g^μν R_μν
   ======================================================================== *)

let ricci_scalar metric rc =
  let sum = ref zero in
  List.iter (fun mu ->
    List.iter (fun nu ->
      let g_inv = get_upper metric mu nu in
      let r_mn = get_ricci rc mu nu in
      if not (is_zero g_inv) && not (is_zero r_mn) then
        sum := add !sum (mul g_inv r_mn)
    ) all_coords
  ) all_coords;
  full_simplify !sum

(* ========================================================================
   EINSTEIN TENSOR
   
   G_μν = R_μν - (1/2) g_μν R
   ======================================================================== *)

let einstein metric rc scalar mu nu =
  let r_mn = get_ricci rc mu nu in
  let g_mn = get_lower metric mu nu in
  sub r_mn (mul (div one two) (mul g_mn scalar)) |> full_simplify

type einstein_tensor = {
  g : expr array array;  (* G_μν *)
}

let compute_einstein metric rc =
  let scalar = ricci_scalar metric rc in
  let g = Array.make_matrix 4 4 zero in
  List.iter (fun mu ->
    List.iter (fun nu ->
      let mi = coord_to_idx mu in
      let ni = coord_to_idx nu in
      if mi <= ni then begin
        let component = einstein metric rc scalar mu nu in
        g.(mi).(ni) <- component;
        g.(ni).(mi) <- component
      end
    ) all_coords
  ) all_coords;
  { g }

let get_einstein et mu nu =
  et.g.(coord_to_idx mu).(coord_to_idx nu)

(* ========================================================================
   PRINTING
   ======================================================================== *)

let print_riemann rt =
  print_endline "Riemann tensor R^ρ_σμν (non-zero, independent components):";
  let coord_str = [|"t"; "r"; "θ"; "φ"|] in
  for ri = 0 to 3 do
    for si = 0 to 3 do
      for mi = 0 to 3 do
        for ni = mi + 1 to 3 do  (* Use antisymmetry *)
          let component = rt.r.(ri).(si).(mi).(ni) in
          if not (is_zero component) then
            Printf.printf "  R^%s_{%s%s%s} = %s\n"
              coord_str.(ri) coord_str.(si) coord_str.(mi) coord_str.(ni)
              (to_string component)
        done
      done
    done
  done

let print_ricci rc =
  print_endline "Ricci tensor R_μν (non-zero):";
  let coord_str = [|"t"; "r"; "θ"; "φ"|] in
  for mi = 0 to 3 do
    for ni = mi to 3 do
      let component = rc.ric.(mi).(ni) in
      if not (is_zero component) then
        Printf.printf "  R_{%s%s} = %s\n"
          coord_str.(mi) coord_str.(ni) (to_string component)
    done
  done

let print_einstein et =
  print_endline "Einstein tensor G_μν (non-zero):";
  let coord_str = [|"t"; "r"; "θ"; "φ"|] in
  for mi = 0 to 3 do
    for ni = mi to 3 do
      let component = et.g.(mi).(ni) in
      if not (is_zero component) then
        Printf.printf "  G_{%s%s} = %s\n"
          coord_str.(mi) coord_str.(ni) (to_string component)
    done
  done
