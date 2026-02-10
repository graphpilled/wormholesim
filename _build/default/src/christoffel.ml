(* christoffel.ml - Christoffel symbols computed from the metric *)

open Expr
open Metric

(* ========================================================================
   CHRISTOFFEL SYMBOLS - ACTUALLY COMPUTED
   
   Γ^σ_μν = (1/2) g^{σρ} (∂_μ g_{νρ} + ∂_ν g_{μρ} - ∂_ρ g_{μν})
   
   This is REAL COMPUTATION:
   1. Differentiate metric components symbolically
   2. Contract with inverse metric
   3. Simplify the result
   ======================================================================== *)

(* Compute ∂_ρ g_μν *)
let metric_derivative metric rho mu nu =
  let g_mu_nu = get_lower metric mu nu in
  differentiate rho g_mu_nu |> full_simplify

(* Compute Γ^σ_μν from definition *)
let christoffel metric sigma mu nu =
  (* Γ^σ_μν = (1/2) g^{σρ} (∂_μ g_{νρ} + ∂_ν g_{μρ} - ∂_ρ g_{μν}) *)
  
  (* Sum over ρ *)
  let sum = ref zero in
  List.iter (fun rho ->
    let g_inv_sigma_rho = get_upper metric sigma rho in
    
    (* Skip if g^{σρ} = 0 (common for diagonal metrics) *)
    if not (is_zero g_inv_sigma_rho) then begin
      let dmu_g_nu_rho = metric_derivative metric mu nu rho in   (* ∂_μ g_{νρ} *)
      let dnu_g_mu_rho = metric_derivative metric nu mu rho in   (* ∂_ν g_{μρ} *)
      let drho_g_mu_nu = metric_derivative metric rho mu nu in   (* ∂_ρ g_{μν} *)
      
      let bracket = sub (add dmu_g_nu_rho dnu_g_mu_rho) drho_g_mu_nu in
      let term = mul g_inv_sigma_rho bracket in
      sum := add !sum term
    end
  ) all_coords;
  
  mul (div one two) !sum |> full_simplify

(* Compute all Christoffel symbols *)
type christoffel_symbols = {
  gamma : expr array array array;  (* Γ^σ_μν indexed [σ][μ][ν] *)
}

let compute_christoffels metric =
  let gamma = Array.init 4 (fun _ -> 
    Array.init 4 (fun _ -> 
      Array.make 4 zero
    )
  ) in
  
  List.iter (fun sigma ->
    List.iter (fun mu ->
      List.iter (fun nu ->
        let s = coord_to_idx sigma in
        let m = coord_to_idx mu in
        let n = coord_to_idx nu in
        (* Use symmetry: Γ^σ_μν = Γ^σ_νμ *)
        if m <= n then begin
          let g = christoffel metric sigma mu nu in
          gamma.(s).(m).(n) <- g;
          gamma.(s).(n).(m) <- g  (* Symmetric *)
        end
      ) all_coords
    ) all_coords
  ) all_coords;
  
  { gamma }

(* Get Christoffel symbol Γ^σ_μν *)
let get_christoffel cs sigma mu nu =
  cs.gamma.(coord_to_idx sigma).(coord_to_idx mu).(coord_to_idx nu)

(* ========================================================================
   PRINTING
   ======================================================================== *)

let print_christoffels cs =
  print_endline "Christoffel symbols Γ^σ_μν (non-zero only):";
  let coord_str = [|"t"; "r"; "θ"; "φ"|] in
  for s = 0 to 3 do
    for m = 0 to 3 do
      for n = m to 3 do  (* Use symmetry *)
        let g = cs.gamma.(s).(m).(n) in
        if not (is_zero g) then
          Printf.printf "  Γ^%s_{%s%s} = %s\n"
            coord_str.(s) coord_str.(m) coord_str.(n) (to_string g)
      done
    done
  done

let print_christoffels_latex cs =
  print_endline "Christoffel symbols (LaTeX):";
  let coord_str = [|"t"; "r"; "\\theta"; "\\phi"|] in
  for s = 0 to 3 do
    for m = 0 to 3 do
      for n = m to 3 do
        let g = cs.gamma.(s).(m).(n) in
        if not (is_zero g) then
          Printf.printf "  \\Gamma^{%s}_{%s%s} &= %s \\\\\n"
            coord_str.(s) coord_str.(m) coord_str.(n) (to_latex g)
      done
    done
  done
