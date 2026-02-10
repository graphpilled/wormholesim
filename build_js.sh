#!/bin/bash
# build_js.sh
cd "$(dirname "$0")"

echo "Building TensorCalc JS..."

ocamlfind ocamlc -package js_of_ocaml -package js_of_ocaml-ppx \
  -linkpkg -ppx "js_of_ocaml-ppx" \
  -I src \
  src/expr.ml \
  src/simplify.ml \
  src/metric.ml \
  src/christoffel.ml \
  src/riemann.ml \
  src/numeric.ml \
  src/geodesic.ml \
  src/derivation.ml \
  src/js_bridge.ml \
  -o tensor_engine.bc

js_of_ocaml tensor_engine.bc -o tensor_engine.js

echo "Done. Output: tensor_engine.js"
