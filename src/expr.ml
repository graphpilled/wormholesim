(* expr.ml - Symbolic expressions for tensor calculus *)

(* Coordinates *)
type coord = T | R | Theta | Phi

let all_coords = [T; R; Theta; Phi]

let coord_to_string = function
  | T -> "t" | R -> "r" | Theta -> "θ" | Phi -> "φ"

let coord_to_latex = function
  | T -> "t" | R -> "r" | Theta -> "\\theta" | Phi -> "\\phi"

(* Symbolic expressions *)
type expr =
  | Num of float                     (* Numeric constant *)
  | Var of string                    (* Named variable: r, θ, etc. *)
  | Func of string * expr            (* Function application: sin(θ), Φ(r) *)
  | Neg of expr                      (* -e *)
  | Add of expr * expr               (* e1 + e2 *)
  | Sub of expr * expr               (* e1 - e2 *)
  | Mul of expr * expr               (* e1 * e2 *)
  | Div of expr * expr               (* e1 / e2 *)
  | Pow of expr * expr               (* e1 ^ e2 *)
  | Exp of expr                      (* e^x *)
  | Log of expr                      (* ln(x) *)
  | Sin of expr                      (* sin(x) *)
  | Cos of expr                      (* cos(x) *)
  | Sqrt of expr                     (* √x *)
  | Delta of coord * coord           (* Kronecker delta δ^μ_ν *)
  | Partial of expr * coord          (* ∂e/∂x^μ - SYMBOLIC, will be computed *)

(* Smart constructors that simplify as they build *)
let zero = Num 0.0
let one = Num 1.0
let two = Num 2.0

let is_zero = function Num 0.0 -> true | _ -> false
let is_one = function Num 1.0 -> true | _ -> false

let neg = function
  | Num n -> Num (-.n)
  | Neg e -> e
  | e -> Neg e

let add e1 e2 = match e1, e2 with
  | Num 0.0, e | e, Num 0.0 -> e
  | Num a, Num b -> Num (a +. b)
  | e1, Neg e2 when e1 = e2 -> zero
  | Neg e1, e2 when e1 = e2 -> zero
  | _ -> Add (e1, e2)

let sub e1 e2 = match e1, e2 with
  | e, Num 0.0 -> e
  | Num 0.0, e -> neg e
  | Num a, Num b -> Num (a -. b)
  | e1, e2 when e1 = e2 -> zero
  | _ -> Sub (e1, e2)

let mul e1 e2 = match e1, e2 with
  | Num 0.0, _ | _, Num 0.0 -> zero
  | Num 1.0, e | e, Num 1.0 -> e
  | Num (-1.0), e | e, Num (-1.0) -> neg e
  | Num a, Num b -> Num (a *. b)
  | _ -> Mul (e1, e2)

let div e1 e2 = match e1, e2 with
  | Num 0.0, _ -> zero
  | e, Num 1.0 -> e
  | Num a, Num b when b <> 0.0 -> Num (a /. b)
  | e1, e2 when e1 = e2 -> one
  | _ -> Div (e1, e2)

let pow e1 e2 = match e1, e2 with
  | _, Num 0.0 -> one
  | e, Num 1.0 -> e
  | Num 1.0, _ -> one
  | Num a, Num b -> Num (a ** b)
  | _ -> Pow (e1, e2)

let sqrt_expr = function
  | Num n when n >= 0.0 -> Num (sqrt n)
  | Pow (e, Num 2.0) -> e  (* √(e²) = e, assuming positive *)
  | e -> Sqrt e

let exp_expr = function
  | Num 0.0 -> one
  | Log e -> e
  | e -> Exp e

let log_expr = function
  | Num 1.0 -> zero
  | Exp e -> e
  | e -> Log e

let sin_expr = function
  | Num 0.0 -> zero
  | e -> Sin e

let cos_expr = function
  | Num 0.0 -> one
  | e -> Cos e

(* Coordinate to variable *)
let coord_var = function
  | T -> Var "t"
  | R -> Var "r"
  | Theta -> Var "θ"
  | Phi -> Var "φ"

(* Check if expression depends on a coordinate *)
let depends_on coord expr =
  let v = coord_var coord in
  let rec check = function
    | Var s -> Var s = v
    | Num _ -> false
    | Func (_, e) -> check e
    | Neg e | Exp e | Log e | Sin e | Cos e | Sqrt e -> check e
    | Add (a, b) | Sub (a, b) | Mul (a, b) | Div (a, b) | Pow (a, b) -> 
        check a || check b
    | Delta _ -> false
    | Partial (e, _) -> check e  (* Partial itself might depend *)
  in check expr

(* ========================================================================
   SYMBOLIC DIFFERENTIATION - The core of actual computation
   ======================================================================== *)

let rec differentiate coord expr =
  let d = differentiate coord in
  let v = coord_var coord in
  match expr with
  | Num _ -> zero
  | Var s -> if Var s = v then one else zero
  | Func (name, arg) ->
      (* Chain rule: d/dx f(g(x)) = f'(g(x)) * g'(x) *)
      (* For now, handle specific functions *)
      begin match name with
      | "Φ" -> (* Redshift function Φ(r) *)
          if coord = R then 
            mul (Func ("Φ'", arg)) (d arg)
          else 
            zero
      | "b" -> (* Shape function b(r) *)
          if coord = R then
            mul (Func ("b'", arg)) (d arg)
          else
            zero
      | _ -> Partial (expr, coord)  (* Unknown function, leave symbolic *)
      end
  | Neg e -> neg (d e)
  | Add (a, b) -> add (d a) (d b)
  | Sub (a, b) -> sub (d a) (d b)
  | Mul (a, b) -> add (mul (d a) b) (mul a (d b))  (* Product rule *)
  | Div (a, b) -> (* Quotient rule: (a/b)' = (a'b - ab')/b² *)
      div (sub (mul (d a) b) (mul a (d b))) (pow b two)
  | Pow (base, Num n) -> (* Power rule: (x^n)' = n*x^(n-1)*x' *)
      mul (mul (Num n) (pow base (Num (n -. 1.0)))) (d base)
  | Pow (base, exp) -> (* General: (f^g)' = f^g * (g'*ln(f) + g*f'/f) *)
      let term1 = mul (d exp) (log_expr base) in
      let term2 = mul exp (div (d base) base) in
      mul (pow base exp) (add term1 term2)
  | Exp e -> mul (Exp e) (d e)  (* (e^f)' = e^f * f' *)
  | Log e -> div (d e) e        (* (ln f)' = f'/f *)
  | Sin e -> mul (Cos e) (d e)  (* (sin f)' = cos(f) * f' *)
  | Cos e -> neg (mul (Sin e) (d e))  (* (cos f)' = -sin(f) * f' *)
  | Sqrt e -> div (d e) (mul two (Sqrt e))  (* (√f)' = f'/(2√f) *)
  | Delta (_, _) -> zero  (* Kronecker delta is constant *)
  | Partial (_, _) -> Partial (expr, coord)  (* Can't simplify further *)

(* ========================================================================
   SIMPLIFICATION
   ======================================================================== *)

let rec simplify expr =
  let s = simplify in
  match expr with
  | Num _ | Var _ | Delta _ -> expr
  | Func (name, e) -> Func (name, s e)
  | Neg e -> neg (s e)
  | Add (a, b) -> add (s a) (s b)
  | Sub (a, b) -> sub (s a) (s b)
  | Mul (a, b) -> mul (s a) (s b)
  | Div (a, b) -> div (s a) (s b)
  | Pow (a, b) -> pow (s a) (s b)
  | Exp e -> exp_expr (s e)
  | Log e -> log_expr (s e)
  | Sin e -> sin_expr (s e)
  | Cos e -> cos_expr (s e)
  | Sqrt e -> sqrt_expr (s e)
  | Partial (e, c) -> 
      let se = s e in
      if depends_on c se then
        simplify (differentiate c se)
      else
        zero

(* Fully simplify by iterating until fixed point *)
let rec full_simplify expr =
  let s = simplify expr in
  if s = expr then s else full_simplify s

(* ========================================================================
   PRETTY PRINTING
   ======================================================================== *)

let rec to_string expr =
  match expr with
  | Num n -> 
      if n = Float.round n then string_of_int (int_of_float n)
      else Printf.sprintf "%.4g" n
  | Var s -> s
  | Func (name, arg) -> Printf.sprintf "%s(%s)" name (to_string arg)
  | Neg e -> Printf.sprintf "-%s" (to_string_parens e)
  | Add (a, b) -> Printf.sprintf "%s + %s" (to_string a) (to_string b)
  | Sub (a, b) -> Printf.sprintf "%s - %s" (to_string a) (to_string_parens b)
  | Mul (a, b) -> Printf.sprintf "%s·%s" (to_string_parens a) (to_string_parens b)
  | Div (a, b) -> Printf.sprintf "%s/%s" (to_string_parens a) (to_string_parens b)
  | Pow (a, b) -> Printf.sprintf "%s^%s" (to_string_parens a) (to_string_parens b)
  | Exp e -> Printf.sprintf "exp(%s)" (to_string e)
  | Log e -> Printf.sprintf "ln(%s)" (to_string e)
  | Sin e -> Printf.sprintf "sin(%s)" (to_string e)
  | Cos e -> Printf.sprintf "cos(%s)" (to_string e)
  | Sqrt e -> Printf.sprintf "√(%s)" (to_string e)
  | Delta (a, b) -> Printf.sprintf "δ_%s^%s" (coord_to_string a) (coord_to_string b)
  | Partial (e, c) -> Printf.sprintf "∂_%s(%s)" (coord_to_string c) (to_string e)

and to_string_parens expr =
  match expr with
  | Add _ | Sub _ -> Printf.sprintf "(%s)" (to_string expr)
  | _ -> to_string expr

let rec to_latex expr =
  match expr with
  | Num n -> 
      if n = Float.round n then string_of_int (int_of_float n)
      else Printf.sprintf "%.4g" n
  | Var "θ" -> "\\theta"
  | Var "φ" -> "\\phi"
  | Var s -> s
  | Func (name, arg) -> Printf.sprintf "%s(%s)" name (to_latex arg)
  | Neg e -> Printf.sprintf "-%s" (to_latex_parens e)
  | Add (a, b) -> Printf.sprintf "%s + %s" (to_latex a) (to_latex b)
  | Sub (a, b) -> Printf.sprintf "%s - %s" (to_latex a) (to_latex_parens b)
  | Mul (a, b) -> Printf.sprintf "%s \\cdot %s" (to_latex_parens a) (to_latex_parens b)
  | Div (a, b) -> Printf.sprintf "\\frac{%s}{%s}" (to_latex a) (to_latex b)
  | Pow (a, b) -> Printf.sprintf "%s^{%s}" (to_latex_parens a) (to_latex b)
  | Exp e -> Printf.sprintf "e^{%s}" (to_latex e)
  | Log e -> Printf.sprintf "\\ln(%s)" (to_latex e)
  | Sin e -> Printf.sprintf "\\sin(%s)" (to_latex e)
  | Cos e -> Printf.sprintf "\\cos(%s)" (to_latex e)
  | Sqrt e -> Printf.sprintf "\\sqrt{%s}" (to_latex e)
  | Delta (a, b) -> Printf.sprintf "\\delta^{%s}_{%s}" (coord_to_latex a) (coord_to_latex b)
  | Partial (e, c) -> Printf.sprintf "\\partial_{%s}(%s)" (coord_to_latex c) (to_latex e)

and to_latex_parens expr =
  match expr with
  | Add _ | Sub _ -> Printf.sprintf "\\left(%s\\right)" (to_latex expr)
  | _ -> to_latex expr
