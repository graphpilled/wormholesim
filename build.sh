#!/bin/bash
# Build TensorCalc v2

cd "$(dirname "$0")"

echo "Building TensorCalc v2..."

# Compile in dependency order
ocamlopt -c src/expr.ml -I src
ocamlopt -c src/simplify.ml -I src
ocamlopt -c src/metric.ml -I src
ocamlopt -c src/christoffel.ml -I src
ocamlopt -c src/riemann.ml -I src
ocamlopt -c src/numeric.ml -I src
ocamlopt -c src/geodesic.ml -I src
ocamlopt -c src/derivation.ml -I src
ocamlopt -c src/main.ml -I src

# Link
ocamlopt -o tensorcalc \
  src/expr.cmx \
  src/simplify.cmx \
  src/metric.cmx \
  src/christoffel.cmx \
  src/riemann.cmx \
  src/numeric.cmx \
  src/geodesic.cmx \
  src/derivation.cmx \
  src/main.cmx

echo "Done. Run with: ./tensorcalc"
