(* derivation.ml - Step-by-step derivation tracer for wormhole tensor computations
   
   Shows all algebraic manipulations when computing tensors from a metric.
   Useful for understanding HOW the computation proceeds, not just the result.
*)

open Expr
open Simplify

(* ========================================================================
   DERIVATION STEP TYPES
   ======================================================================== *)

type step_type =
  | Definition of string      (* Defining a quantity *)
  | Differentiation of string (* Computing a derivative *)
  | Substitution of string    (* Substituting a value *)
  | Simplification            (* Simplifying an expression *)
  | Contraction of string     (* Index contraction *)
  | Addition                  (* Adding terms *)
  | Result of string          (* Final result *)

type derivation_step = {
  step_num: int;
  step_type: step_type;
  description: string;
  expression: expr;
  simplified: expr option;
}

type derivation = {
  title: string;
  steps: derivation_step list;
}

(* Global step counter *)
let step_counter = ref 0

let reset_steps () = step_counter := 0

let make_step stype desc expr simp =
  incr step_counter;
  { step_num = !step_counter; step_type = stype; 
    description = desc; expression = expr; simplified = simp }

(* ========================================================================
   PRETTY PRINTING
   ======================================================================== *)

let step_type_str = function
  | Definition s -> Printf.sprintf "DEFINE %s" s
  | Differentiation s -> Printf.sprintf "∂/∂%s" s
  | Substitution s -> Printf.sprintf "SUBSTITUTE %s" s
  | Simplification -> "SIMPLIFY"
  | Contraction s -> Printf.sprintf "CONTRACT %s" s
  | Addition -> "ADD"
  | Result s -> Printf.sprintf "RESULT: %s" s

let print_step step =
  Printf.printf "  [%d] %s\n" step.step_num (step_type_str step.step_type);
  Printf.printf "      %s\n" step.description;
  Printf.printf "      = %s\n" (to_string step.expression);
  match step.simplified with
  | Some s when not (expr_equal step.expression s) ->
      Printf.printf "      → %s\n" (to_string s)
  | _ -> ()

let print_derivation deriv =
  Printf.printf "\n┌─────────────────────────────────────────────────────────────────┐\n";
  Printf.printf "│ DERIVATION: %-51s │\n" deriv.title;
  Printf.printf "└─────────────────────────────────────────────────────────────────┘\n\n";
  List.iter (fun step ->
    print_step step;
    print_endline ""
  ) deriv.steps

(* ========================================================================
   TRACED METRIC DERIVATIVE
   ======================================================================== *)

let trace_metric_derivative metric mu nu rho =
  let coord_names = [|"t"; "r"; "θ"; "φ"|] in
  let coord_to_int = function T -> 0 | R -> 1 | Theta -> 2 | Phi -> 3 in
  let mi, ni, ri = coord_to_int mu, coord_to_int nu, coord_to_int rho in
  
  let g_mn = Metric.get_lower metric mu nu in
  let steps = ref [] in
  
  (* Step 1: Starting metric component *)
  steps := [make_step 
    (Definition (Printf.sprintf "g_{%s%s}" coord_names.(mi) coord_names.(ni)))
    (Printf.sprintf "Metric component g_{%s%s}" coord_names.(mi) coord_names.(ni))
    g_mn None];
  
  (* Step 2: Compute derivative *)
  let coord_var = match rho with
    | T -> Var "t" | R -> Var "r" | Theta -> Var "θ" | Phi -> Var "φ" in
  let dg = differentiate rho g_mn in
  steps := !steps @ [make_step
    (Differentiation coord_names.(ri))
    (Printf.sprintf "∂g_{%s%s}/∂%s" coord_names.(mi) coord_names.(ni) coord_names.(ri))
    dg None];
  
  (* Step 3: Simplify *)
  let dg_simp = simplify_complete dg in
  if not (expr_equal dg dg_simp) then
    steps := !steps @ [make_step
      Simplification
      "Apply algebraic simplification"
      dg (Some dg_simp)];
  
  (* Result *)
  steps := !steps @ [make_step
    (Result (Printf.sprintf "∂g_{%s%s}/∂%s" coord_names.(mi) coord_names.(ni) coord_names.(ri)))
    "Final derivative"
    dg_simp None];
  
  { title = Printf.sprintf "Metric Derivative ∂g_{%s%s}/∂%s" 
      coord_names.(mi) coord_names.(ni) coord_names.(ri);
    steps = !steps }

(* ========================================================================
   TRACED CHRISTOFFEL COMPUTATION
   ======================================================================== *)

let trace_christoffel_component metric sigma mu nu =
  let coord_names = [|"t"; "r"; "θ"; "φ"|] in
  let coord_to_int = function T -> 0 | R -> 1 | Theta -> 2 | Phi -> 3 in
  let all_coords = [T; R; Theta; Phi] in
  let si, mi, ni = coord_to_int sigma, coord_to_int mu, coord_to_int nu in
  
  reset_steps ();
  let steps = ref [] in
  
  (* Title step *)
  steps := [make_step
    (Definition (Printf.sprintf "Γ^%s_{%s%s}" coord_names.(si) coord_names.(mi) coord_names.(ni)))
    (Printf.sprintf "Computing Christoffel symbol Γ^%s_{%s%s}" 
      coord_names.(si) coord_names.(mi) coord_names.(ni))
    zero None];
  
  (* Show formula *)
  steps := !steps @ [make_step
    (Definition "formula")
    "Using Γ^σ_μν = (1/2) g^σρ (∂_μ g_νρ + ∂_ν g_μρ - ∂_ρ g_μν)"
    zero None];
  
  let total = ref zero in
  
  (* Sum over ρ *)
  List.iter (fun rho ->
    let ri = coord_to_int rho in
    
    (* Get inverse metric *)
    let g_inv = Metric.get_upper metric sigma rho in
    if not (is_zero g_inv) then begin
      steps := !steps @ [make_step
        (Definition (Printf.sprintf "g^{%s%s}" coord_names.(si) coord_names.(ri)))
        (Printf.sprintf "Inverse metric component")
        g_inv None];
      
      (* Three derivative terms *)
      let d1 = differentiate mu (Metric.get_lower metric nu rho) in
      let d2 = differentiate nu (Metric.get_lower metric mu rho) in
      let d3 = differentiate rho (Metric.get_lower metric mu nu) in
      
      if not (is_zero d1) then
        steps := !steps @ [make_step
          (Differentiation coord_names.(mi))
          (Printf.sprintf "∂_{%s} g_{%s%s}" coord_names.(mi) coord_names.(ni) coord_names.(ri))
          d1 (Some (simplify_complete d1))];
      
      if not (is_zero d2) then
        steps := !steps @ [make_step
          (Differentiation coord_names.(ni))
          (Printf.sprintf "∂_{%s} g_{%s%s}" coord_names.(ni) coord_names.(mi) coord_names.(ri))
          d2 (Some (simplify_complete d2))];
      
      if not (is_zero d3) then
        steps := !steps @ [make_step
          (Differentiation coord_names.(ri))
          (Printf.sprintf "∂_{%s} g_{%s%s}" coord_names.(ri) coord_names.(mi) coord_names.(ni))
          d3 (Some (simplify_complete d3))];
      
      (* Combine: (d1 + d2 - d3) *)
      let bracket = sub (add d1 d2) d3 in
      let bracket_simp = simplify_complete bracket in
      
      if not (is_zero bracket_simp) then begin
        steps := !steps @ [make_step
          Addition
          (Printf.sprintf "(∂_%s g_{%s%s} + ∂_%s g_{%s%s} - ∂_%s g_{%s%s})" 
            coord_names.(mi) coord_names.(ni) coord_names.(ri)
            coord_names.(ni) coord_names.(mi) coord_names.(ri)
            coord_names.(ri) coord_names.(mi) coord_names.(ni))
          bracket (Some bracket_simp)];
        
        (* Multiply by (1/2) g^σρ *)
        let term = mul (mul (Num 0.5) g_inv) bracket_simp in
        let term_simp = simplify_complete term in
        
        steps := !steps @ [make_step
          (Contraction (Printf.sprintf "ρ=%s" coord_names.(ri)))
          (Printf.sprintf "(1/2) g^{%s%s} × bracket" coord_names.(si) coord_names.(ri))
          term (Some term_simp)];
        
        total := add !total term_simp
      end
    end
  ) all_coords;
  
  (* Final simplification *)
  let final = simplify_complete !total in
  steps := !steps @ [make_step
    (Result (Printf.sprintf "Γ^%s_{%s%s}" coord_names.(si) coord_names.(mi) coord_names.(ni)))
    "Sum over all ρ and simplify"
    !total (Some final)];
  
  { title = Printf.sprintf "Christoffel Symbol Γ^%s_{%s%s}" 
      coord_names.(si) coord_names.(mi) coord_names.(ni);
    steps = !steps }

(* ========================================================================
   TRACED RIEMANN COMPONENT
   ======================================================================== *)

let trace_riemann_component christoffels rho sigma mu nu =
  let coord_names = [|"t"; "r"; "θ"; "φ"|] in
  let coord_to_int = function T -> 0 | R -> 1 | Theta -> 2 | Phi -> 3 in
  let all_coords = [T; R; Theta; Phi] in
  let rh, si, mi, ni = coord_to_int rho, coord_to_int sigma, 
                        coord_to_int mu, coord_to_int nu in
  
  reset_steps ();
  let steps = ref [] in
  
  steps := [make_step
    (Definition (Printf.sprintf "R^%s_{%s%s%s}" 
      coord_names.(rh) coord_names.(si) coord_names.(mi) coord_names.(ni)))
    "Using R^ρ_σμν = ∂_μ Γ^ρ_νσ - ∂_ν Γ^ρ_μσ + Γ^ρ_μλ Γ^λ_νσ - Γ^ρ_νλ Γ^λ_μσ"
    zero None];
  
  let gamma = christoffels.Christoffel.gamma in
  
  (* Term 1: ∂_μ Γ^ρ_νσ *)
  let g_rho_nu_sigma = gamma.(rh).(ni).(si) in
  let term1 = differentiate mu g_rho_nu_sigma in
  let term1_simp = simplify_complete term1 in
  
  steps := !steps @ [make_step
    (Differentiation coord_names.(mi))
    (Printf.sprintf "∂_%s Γ^%s_{%s%s}" 
      coord_names.(mi) coord_names.(rh) coord_names.(ni) coord_names.(si))
    term1 (Some term1_simp)];
  
  (* Term 2: -∂_ν Γ^ρ_μσ *)
  let g_rho_mu_sigma = gamma.(rh).(mi).(si) in
  let term2 = neg (differentiate nu g_rho_mu_sigma) in
  let term2_simp = simplify_complete term2 in
  
  steps := !steps @ [make_step
    (Differentiation coord_names.(ni))
    (Printf.sprintf "-∂_%s Γ^%s_{%s%s}" 
      coord_names.(ni) coord_names.(rh) coord_names.(mi) coord_names.(si))
    term2 (Some term2_simp)];
  
  (* Term 3: Γ^ρ_μλ Γ^λ_νσ (sum over λ) *)
  let term3 = ref zero in
  List.iter (fun lambda ->
    let li = coord_to_int lambda in
    let g1 = gamma.(rh).(mi).(li) in
    let g2 = gamma.(li).(ni).(si) in
    let prod = mul g1 g2 in
    if not (is_zero (simplify_complete prod)) then begin
      steps := !steps @ [make_step
        (Contraction (Printf.sprintf "λ=%s" coord_names.(li)))
        (Printf.sprintf "Γ^%s_{%s%s} Γ^%s_{%s%s}"
          coord_names.(rh) coord_names.(mi) coord_names.(li)
          coord_names.(li) coord_names.(ni) coord_names.(si))
        prod (Some (simplify_complete prod))];
      term3 := add !term3 prod
    end
  ) all_coords;
  let term3_simp = simplify_complete !term3 in
  
  (* Term 4: -Γ^ρ_νλ Γ^λ_μσ (sum over λ) *)
  let term4 = ref zero in
  List.iter (fun lambda ->
    let li = coord_to_int lambda in
    let g1 = gamma.(rh).(ni).(li) in
    let g2 = gamma.(li).(mi).(si) in
    let prod = neg (mul g1 g2) in
    if not (is_zero (simplify_complete prod)) then begin
      steps := !steps @ [make_step
        (Contraction (Printf.sprintf "λ=%s" coord_names.(li)))
        (Printf.sprintf "-Γ^%s_{%s%s} Γ^%s_{%s%s}"
          coord_names.(rh) coord_names.(ni) coord_names.(li)
          coord_names.(li) coord_names.(mi) coord_names.(si))
        prod (Some (simplify_complete prod))];
      term4 := add !term4 prod
    end
  ) all_coords;
  let term4_simp = simplify_complete !term4 in
  
  (* Combine all terms *)
  let total = add (add (add term1_simp term2_simp) term3_simp) term4_simp in
  let final = simplify_complete total in
  
  steps := !steps @ [make_step
    (Result (Printf.sprintf "R^%s_{%s%s%s}" 
      coord_names.(rh) coord_names.(si) coord_names.(mi) coord_names.(ni)))
    "Sum all four terms"
    total (Some final)];
  
  { title = Printf.sprintf "Riemann Component R^%s_{%s%s%s}"
      coord_names.(rh) coord_names.(si) coord_names.(mi) coord_names.(ni);
    steps = !steps }

(* ========================================================================
   WORMHOLE PARAMETER CHANGE TRACER
   ======================================================================== *)

type wormhole_config = {
  name: string;
  shape_function: expr;     (* b(r) *)
  shape_derivative: expr;   (* b'(r) *)
  redshift_function: expr;  (* Φ(r) *)
  throat_radius: float;
}

let constant_wormhole b0 = {
  name = Printf.sprintf "Constant b(r) = %.2f" b0;
  shape_function = Num b0;
  shape_derivative = zero;
  redshift_function = zero;
  throat_radius = b0;
}

let ellis_wormhole b0 = {
  name = Printf.sprintf "Ellis b(r) = %.2f²/r" b0;
  shape_function = Div (Num (b0 *. b0), Var "r");
  shape_derivative = Neg (Div (Num (b0 *. b0), Pow (Var "r", Num 2.0)));
  redshift_function = zero;
  throat_radius = b0;
}

let custom_wormhole ~name ~b_expr ~b_deriv ~phi_expr ~b0 = {
  name;
  shape_function = b_expr;
  shape_derivative = b_deriv;
  redshift_function = phi_expr;
  throat_radius = b0;
}

(* Substitute b(r) and b'(r) with specific functions *)
let rec substitute_shape config expr =
  match expr with
  | Func ("b", Var "r") -> config.shape_function
  | Func ("b'", Var "r") -> config.shape_derivative
  | Func ("Φ", Var "r") -> config.redshift_function
  | Func ("Φ'", Var "r") -> differentiate R config.redshift_function
  | Neg e -> Neg (substitute_shape config e)
  | Add (a, b) -> Add (substitute_shape config a, substitute_shape config b)
  | Sub (a, b) -> Sub (substitute_shape config a, substitute_shape config b)
  | Mul (a, b) -> Mul (substitute_shape config a, substitute_shape config b)
  | Div (a, b) -> Div (substitute_shape config a, substitute_shape config b)
  | Pow (a, b) -> Pow (substitute_shape config a, substitute_shape config b)
  | Sin e -> Sin (substitute_shape config e)
  | Cos e -> Cos (substitute_shape config e)
  | Exp e -> Exp (substitute_shape config e)
  | Log e -> Log (substitute_shape config e)
  | Sqrt e -> Sqrt (substitute_shape config e)
  | _ -> expr

let trace_parameter_change ~from_config ~to_config metric christoffels =
  let coord_names = [|"t"; "r"; "θ"; "φ"|] in
  
  Printf.printf "\n";
  Printf.printf "╔═══════════════════════════════════════════════════════════════════╗\n";
  Printf.printf "║  PARAMETER CHANGE: %s\n" from_config.name;
  Printf.printf "║                 → %s\n" to_config.name;
  Printf.printf "╚═══════════════════════════════════════════════════════════════════╝\n\n";
  
  (* Show shape function change *)
  Printf.printf "Shape function change:\n";
  Printf.printf "  b(r): %s → %s\n" 
    (to_string from_config.shape_function)
    (to_string to_config.shape_function);
  Printf.printf "  b'(r): %s → %s\n"
    (to_string from_config.shape_derivative)
    (to_string to_config.shape_derivative);
  Printf.printf "\n";
  
  (* Show how key Christoffels change *)
  Printf.printf "Effect on Christoffel symbols:\n";
  Printf.printf "─────────────────────────────────\n\n";
  
  let gamma = christoffels.Christoffel.gamma in
  
  (* Show Γ^r_rr as example *)
  let g_r_rr = gamma.(1).(1).(1) in
  let g_r_rr_from = simplify_complete (substitute_shape from_config g_r_rr) in
  let g_r_rr_to = simplify_complete (substitute_shape to_config g_r_rr) in
  
  Printf.printf "  Γ^r_{rr}:\n";
  Printf.printf "    General form: %s\n" (to_string (simplify_complete g_r_rr));
  Printf.printf "    With %s:\n" from_config.name;
  Printf.printf "      = %s\n" (to_string g_r_rr_from);
  Printf.printf "    With %s:\n" to_config.name;
  Printf.printf "      = %s\n" (to_string g_r_rr_to);
  Printf.printf "\n";
  
  (* Show Γ^r_θθ *)
  let g_r_thth = gamma.(1).(2).(2) in
  let g_r_thth_from = simplify_complete (substitute_shape from_config g_r_thth) in
  let g_r_thth_to = simplify_complete (substitute_shape to_config g_r_thth) in
  
  Printf.printf "  Γ^r_{θθ}:\n";
  Printf.printf "    General form: %s\n" (to_string (simplify_complete g_r_thth));
  Printf.printf "    With %s:\n" from_config.name;
  Printf.printf "      = %s\n" (to_string g_r_thth_from);
  Printf.printf "    With %s:\n" to_config.name;
  Printf.printf "      = %s\n" (to_string g_r_thth_to);
  Printf.printf "\n";
  
  (* Numerical comparison at sample point *)
  Printf.printf "Numerical comparison at r = 2.0, θ = π/2:\n";
  Printf.printf "───────────────────────────────────────────\n";
  
  let eval_at config expr =
    let substituted = substitute_shape config expr in
    let env = [("r", 2.0); ("θ", Float.pi /. 2.0); ("φ", 0.0)] in
    Simplify.eval_at_point substituted env
  in
  
  Printf.printf "  Symbol          %15s  %15s\n" from_config.name to_config.name;
  Printf.printf "  ──────────────  ───────────────  ───────────────\n";
  
  for s = 0 to 3 do
    for m = 0 to 3 do
      for n = m to 3 do
        let g = gamma.(s).(m).(n) in
        let gs = simplify_complete g in
        if not (is_zero gs) then begin
          let v_from = eval_at from_config g in
          let v_to = eval_at to_config g in
          Printf.printf "  Γ^%s_{%s%s}         %15.6f  %15.6f\n"
            coord_names.(s) coord_names.(m) coord_names.(n) v_from v_to
        end
      done
    done
  done;
  Printf.printf "\n"

(* ========================================================================
   EXPORT DERIVATION TO JSON (for web visualization)
   ======================================================================== *)

let escape_json_string s =
  let b = Buffer.create (String.length s) in
  String.iter (fun c ->
    match c with
    | '"' -> Buffer.add_string b "\\\""
    | '\\' -> Buffer.add_string b "\\\\"
    | '\n' -> Buffer.add_string b "\\n"
    | '\t' -> Buffer.add_string b "\\t"
    | c -> Buffer.add_char b c
  ) s;
  Buffer.contents b

let derivation_to_json deriv =
  let step_to_json step =
    let type_str = match step.step_type with
      | Definition s -> Printf.sprintf "\"definition\", \"target\": \"%s\"" (escape_json_string s)
      | Differentiation s -> Printf.sprintf "\"differentiation\", \"variable\": \"%s\"" (escape_json_string s)
      | Substitution s -> Printf.sprintf "\"substitution\", \"value\": \"%s\"" (escape_json_string s)
      | Simplification -> "\"simplification\""
      | Contraction s -> Printf.sprintf "\"contraction\", \"index\": \"%s\"" (escape_json_string s)
      | Addition -> "\"addition\""
      | Result s -> Printf.sprintf "\"result\", \"name\": \"%s\"" (escape_json_string s)
    in
    let simp_str = match step.simplified with
      | Some s -> Printf.sprintf ", \"simplified\": \"%s\"" (escape_json_string (to_string s))
      | None -> ""
    in
    Printf.sprintf "    {\"step\": %d, \"type\": %s, \"description\": \"%s\", \"expression\": \"%s\"%s}"
      step.step_num type_str 
      (escape_json_string step.description)
      (escape_json_string (to_string step.expression))
      simp_str
  in
  let steps_json = String.concat ",\n" (List.map step_to_json deriv.steps) in
  Printf.sprintf "{\n  \"title\": \"%s\",\n  \"steps\": [\n%s\n  ]\n}"
    (escape_json_string deriv.title) steps_json

let save_derivation_json filename deriv =
  let json = derivation_to_json deriv in
  let oc = open_out filename in
  output_string oc json;
  close_out oc;
  Printf.printf "Saved derivation to %s\n" filename
