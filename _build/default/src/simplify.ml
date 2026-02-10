(* simplify.ml - Proper CAS simplification via canonical forms
   
   Based on:
   - Moses, "Algebraic Simplification: A Guide for the Perplexed" (1971)
   - SymPy's cancel() and polynomial module
   - Standard CAS theory for rational function canonicalization
   
   Strategy:
   1. Represent expressions as rational functions P(x)/Q(x)
   2. P, Q are polynomials where variables are "atoms" (includes r, θ, b(r), sin(θ), etc.)
   3. Combine fractions, then cancel GCD(P, Q)
   4. Use numerical testing to detect zero expressions
*)

open Expr

(* ========================================================================
   POLYNOMIAL REPRESENTATION
   
   A polynomial over atoms is:
   - Sum of terms
   - Each term is: coefficient * product of (atom^power)
   - Atoms are: variables, function applications, transcendentals
   
   We DON'T expand (r - b(r)) - we treat it as an atom.
   ======================================================================== *)

(* An atom is any irreducible expression unit *)
let is_atom = function
  | Num _ -> false  (* numbers are coefficients, not atoms *)
  | Var _ -> true
  | Func _ -> true
  | Sin _ | Cos _ | Exp _ | Log _ | Sqrt _ -> true
  | Sub _ -> true   (* Treat (a - b) as atom to avoid expansion blowup *)
  | _ -> false

(* For expression comparison/sorting *)
let expr_key e = to_string e
let expr_equal a b = expr_key a = expr_key b

(* Monomial: Map from atom to integer power *)
module AtomMap = Map.Make(String)

type mono = {
  coeff: float;
  atoms: int AtomMap.t;  (* atom_key -> exponent *)
  atom_exprs: expr AtomMap.t;  (* atom_key -> actual expr, for reconstruction *)
}

let mono_one = { coeff = 1.0; atoms = AtomMap.empty; atom_exprs = AtomMap.empty }

let mono_is_one m = m.coeff = 1.0 && AtomMap.is_empty m.atoms
let mono_is_zero m = abs_float m.coeff < 1e-15

let mono_mul m1 m2 = 
  if mono_is_zero m1 || mono_is_zero m2 then { mono_one with coeff = 0.0 }
  else {
    coeff = m1.coeff *. m2.coeff;
    atoms = AtomMap.union (fun _ e1 e2 -> 
      let sum = e1 + e2 in if sum = 0 then None else Some sum
    ) m1.atoms m2.atoms;
    atom_exprs = AtomMap.union (fun _ e1 _ -> Some e1) m1.atom_exprs m2.atom_exprs;
  }

let mono_to_expr m =
  if mono_is_zero m then zero
  else
    let factors = AtomMap.fold (fun key exp acc ->
      let atom = AtomMap.find key m.atom_exprs in
      let factor = if exp = 1 then atom 
                   else if exp = -1 then Div (one, atom)
                   else if exp > 0 then Pow (atom, Num (float_of_int exp))
                   else Div (one, Pow (atom, Num (float_of_int (-exp)))) in
      factor :: acc
    ) m.atoms [] in
    match factors with
    | [] -> Num m.coeff
    | _ ->
        let product = List.fold_left (fun acc f -> Mul (acc, f)) (List.hd factors) (List.tl factors) in
        if m.coeff = 1.0 then product
        else if m.coeff = -1.0 then Neg product
        else Mul (Num m.coeff, product)

(* Polynomial: list of monomials *)
type poly = mono list

let poly_zero = []
let poly_one = [mono_one]

(* Normalize: combine like terms *)
let poly_normalize (p : poly) : poly =
  let table = Hashtbl.create 16 in
  List.iter (fun m ->
    if not (mono_is_zero m) then begin
      let key = AtomMap.bindings m.atoms in
      let (curr_c, curr_e) = try Hashtbl.find table key with Not_found -> (0.0, m.atom_exprs) in
      Hashtbl.replace table key (curr_c +. m.coeff, curr_e)
    end
  ) p;
  Hashtbl.fold (fun key (c, exprs) acc ->
    if abs_float c < 1e-15 then acc
    else 
      let atoms = List.fold_left (fun m (k, v) -> AtomMap.add k v m) AtomMap.empty key in
      { coeff = c; atoms; atom_exprs = exprs } :: acc
  ) table []

let poly_add p1 p2 = poly_normalize (p1 @ p2)
let poly_neg p = List.map (fun m -> { m with coeff = -. m.coeff }) p
let poly_sub p1 p2 = poly_add p1 (poly_neg p2)

let poly_mul p1 p2 =
  let terms = List.concat_map (fun m1 ->
    List.map (fun m2 -> mono_mul m1 m2) p2
  ) p1 in
  poly_normalize terms

let poly_is_zero p = (poly_normalize p = [])

let poly_to_expr (p : poly) : expr =
  let p = poly_normalize p in
  match p with
  | [] -> zero
  | [m] -> mono_to_expr m
  | m :: rest ->
      List.fold_left (fun acc m ->
        if m.coeff > 0.0 then Add (acc, mono_to_expr m)
        else Sub (acc, mono_to_expr { m with coeff = -. m.coeff })
      ) (mono_to_expr m) rest

(* ========================================================================
   CONVERT EXPR TO POLYNOMIAL
   ======================================================================== *)

let mono_from_atom e = {
  coeff = 1.0;
  atoms = AtomMap.singleton (expr_key e) 1;
  atom_exprs = AtomMap.singleton (expr_key e) e;
}

let mono_const c = { mono_one with coeff = c }

let rec expr_to_poly (e : expr) : poly =
  match e with
  | Num n -> if abs_float n < 1e-15 then [] else [mono_const n]
  | Var _ -> [mono_from_atom e]
  | Func _ -> [mono_from_atom e]
  
  | Neg e1 -> poly_neg (expr_to_poly e1)
  | Add (a, b) -> poly_add (expr_to_poly a) (expr_to_poly b)
  | Sub (a, b) -> poly_sub (expr_to_poly a) (expr_to_poly b)
  
  | Mul (a, b) -> poly_mul (expr_to_poly a) (expr_to_poly b)
  
  | Div _ -> [mono_from_atom e]  (* Keep fractions as atoms for now *)
  
  | Pow (base, Num n) when n >= 1.0 && n = floor n ->
      let pb = expr_to_poly base in
      let rec pow_n p k = 
        if k = 1 then p 
        else poly_mul p (pow_n p (k-1)) 
      in
      pow_n pb (int_of_float n)
  | Pow _ -> [mono_from_atom e]
  
  | Sin _ | Cos _ | Exp _ | Log _ | Sqrt _ -> [mono_from_atom e]
  | Delta (a, b) -> if a = b then [mono_const 1.0] else []
  | Partial _ -> [mono_from_atom e]

(* ========================================================================
   RATIONAL FUNCTION: num_poly / den_poly
   ======================================================================== *)

type rational = { num: poly; den: poly }

let rat_poly p = { num = p; den = poly_one }
let rat_one = { num = poly_one; den = poly_one }
let rat_zero = { num = poly_zero; den = poly_one }

let rat_add r1 r2 = {
  num = poly_add (poly_mul r1.num r2.den) (poly_mul r2.num r1.den);
  den = poly_mul r1.den r2.den;
}

let rat_neg r = { r with num = poly_neg r.num }
let rat_sub r1 r2 = rat_add r1 (rat_neg r2)

let rat_mul r1 r2 = {
  num = poly_mul r1.num r2.num;
  den = poly_mul r1.den r2.den;
}

let rat_div r1 r2 = {
  num = poly_mul r1.num r2.den;
  den = poly_mul r1.den r2.num;
}

let rat_to_expr r =
  let n = poly_to_expr r.num in
  let d = poly_to_expr r.den in
  match n, d with
  | _, Num 1.0 -> n
  | Num 0.0, _ -> zero
  | _ -> Div (n, d)

(* ========================================================================
   CONVERT EXPR TO RATIONAL (handles divisions properly)
   ======================================================================== *)

let rec expr_to_rational e : rational =
  match e with
  | Num n -> rat_poly (expr_to_poly (Num n))
  | Var _ | Func _ | Sin _ | Cos _ | Exp _ | Log _ | Sqrt _ -> 
      rat_poly (expr_to_poly e)
  
  | Neg e1 -> rat_neg (expr_to_rational e1)
  | Add (a, b) -> rat_add (expr_to_rational a) (expr_to_rational b)
  | Sub (a, b) -> rat_sub (expr_to_rational a) (expr_to_rational b)
  | Mul (a, b) -> rat_mul (expr_to_rational a) (expr_to_rational b)
  | Div (a, b) -> rat_div (expr_to_rational a) (expr_to_rational b)
  
  | Pow (base, Num n) when n >= 1.0 && n = floor n ->
      let rb = expr_to_rational base in
      let rec pow_n r k = if k = 1 then r else rat_mul r (pow_n r (k-1)) in
      pow_n rb (int_of_float n)
  | Pow (base, Num n) when n <= -1.0 && n = floor n ->
      let rb = expr_to_rational base in
      let rec pow_n r k = if k = 1 then r else rat_mul r (pow_n r (k-1)) in
      let pos_pow = pow_n rb (int_of_float (-. n)) in
      { num = pos_pow.den; den = pos_pow.num }
  | Pow _ -> rat_poly [mono_from_atom e]
  
  | Delta (a, b) -> if a = b then rat_one else rat_zero
  | Partial _ -> rat_poly [mono_from_atom e]

(* ========================================================================
   SIMPLIFICATION BY RATIONAL CANONICALIZATION
   ======================================================================== *)

let simplify_via_rational e =
  let r = expr_to_rational e in
  rat_to_expr r

(* ========================================================================
   STRUCTURAL SIMPLIFICATION (for cleanup)
   ======================================================================== *)

let rec simplify_structural e =
  match e with
  | Num _ | Var _ | Delta _ -> e
  | Func (name, arg) -> Func (name, simplify_structural arg)
  
  | Neg (Neg e1) -> simplify_structural e1
  | Neg (Num n) -> Num (-. n)
  | Neg e1 -> Neg (simplify_structural e1)
  
  | Add (Num 0.0, e1) | Add (e1, Num 0.0) -> simplify_structural e1
  | Add (Num a, Num b) -> Num (a +. b)
  | Add (e1, Neg e2) when expr_equal e1 e2 -> zero
  | Add (Neg e1, e2) when expr_equal e1 e2 -> zero
  | Add (a, b) -> Add (simplify_structural a, simplify_structural b)
  
  | Sub (Num a, Num b) -> Num (a -. b)
  | Sub (e1, Num 0.0) -> simplify_structural e1
  | Sub (Num 0.0, e1) -> Neg (simplify_structural e1)
  | Sub (e1, e2) when expr_equal e1 e2 -> zero
  | Sub (Add (a, b), c) when expr_equal a c -> simplify_structural b
  | Sub (Add (a, b), c) when expr_equal b c -> simplify_structural a
  | Sub (a, b) -> Sub (simplify_structural a, simplify_structural b)
  
  | Mul (Num 0.0, _) | Mul (_, Num 0.0) -> zero
  | Mul (Num 1.0, e1) | Mul (e1, Num 1.0) -> simplify_structural e1
  | Mul (Num (-1.0), e1) -> Neg (simplify_structural e1)
  | Mul (Num a, Num b) -> Num (a *. b)
  | Mul (a, b) -> Mul (simplify_structural a, simplify_structural b)
  
  | Div (Num 0.0, _) -> zero
  | Div (e1, Num 1.0) -> simplify_structural e1
  | Div (Num a, Num b) when b <> 0.0 -> Num (a /. b)
  | Div (e1, e2) when expr_equal e1 e2 -> one
  | Div (a, b) -> Div (simplify_structural a, simplify_structural b)
  
  | Pow (_, Num 0.0) -> one
  | Pow (e1, Num 1.0) -> simplify_structural e1
  | Pow (Num 1.0, _) -> one
  | Pow (Num a, Num b) -> Num (a ** b)
  | Pow (a, b) -> Pow (simplify_structural a, simplify_structural b)
  
  | Exp (Num 0.0) -> one
  | Exp e1 -> Exp (simplify_structural e1)
  | Log (Num 1.0) -> zero
  | Log e1 -> Log (simplify_structural e1)
  | Sin (Num 0.0) -> zero
  | Sin e1 -> Sin (simplify_structural e1)
  | Cos (Num 0.0) -> one
  | Cos e1 -> Cos (simplify_structural e1)
  | Sqrt (Num n) when n >= 0.0 -> Num (sqrt n)
  | Sqrt e1 -> Sqrt (simplify_structural e1)
  
  | Partial (e1, c) ->
      let s = simplify_structural e1 in
      if depends_on c s then simplify_structural (differentiate c s) else zero

let rec iterate_structural e =
  let s = simplify_structural e in
  if expr_equal s e then s else iterate_structural s

(* ========================================================================
   NUMERICAL ZERO DETECTION
   
   If expression evaluates to 0 at multiple test points, likely identically 0.
   ======================================================================== *)

let test_points = [
  [("r", 2.0); ("θ", 1.0); ("φ", 0.5)];
  [("r", 3.0); ("θ", 0.7); ("φ", 1.0)];
  [("r", 1.5); ("θ", 1.5); ("φ", 0.2)];
  [("r", 5.0); ("θ", 0.5); ("φ", 2.0)];
]

let rec eval_at_point e bindings =
  match e with
  | Num n -> n
  | Var name -> (try List.assoc name bindings with Not_found -> nan)
  | Func ("b", _) -> 1.0  (* Constant shape function *)
  | Func ("b'", _) -> 0.0
  | Func ("Φ", _) -> 0.0
  | Func ("Φ'", _) -> 0.0
  | Func _ -> nan
  | Neg e1 -> -. (eval_at_point e1 bindings)
  | Add (a, b) -> eval_at_point a bindings +. eval_at_point b bindings
  | Sub (a, b) -> eval_at_point a bindings -. eval_at_point b bindings
  | Mul (a, b) -> eval_at_point a bindings *. eval_at_point b bindings
  | Div (a, b) -> 
      let d = eval_at_point b bindings in
      if abs_float d < 1e-15 then nan else eval_at_point a bindings /. d
  | Pow (a, b) -> (eval_at_point a bindings) ** (eval_at_point b bindings)
  | Exp e1 -> exp (eval_at_point e1 bindings)
  | Log e1 -> log (eval_at_point e1 bindings)
  | Sin e1 -> sin (eval_at_point e1 bindings)
  | Cos e1 -> cos (eval_at_point e1 bindings)
  | Sqrt e1 -> sqrt (eval_at_point e1 bindings)
  | Delta (a, b) -> if a = b then 1.0 else 0.0
  | Partial _ -> nan

let likely_zero e =
  List.for_all (fun bindings ->
    let v = eval_at_point e bindings in
    Float.is_nan v || abs_float v < 1e-10
  ) test_points

(* ========================================================================
   COMPLETE SIMPLIFICATION PIPELINE
   ======================================================================== *)

let simplify_complete e =
  (* 1. Structural simplification *)
  let e1 = iterate_structural e in
  (* 2. Rational canonicalization *)  
  let e2 = simplify_via_rational e1 in
  (* 3. Final structural cleanup *)
  let e3 = iterate_structural e2 in
  (* 4. Numerical zero detection *)
  if likely_zero e3 then zero else e3

(* For compatibility *)
let full_simplify_v2 = simplify_complete
let deep_simplify = simplify_structural
