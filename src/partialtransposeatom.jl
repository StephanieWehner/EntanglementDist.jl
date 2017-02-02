#############################################################################
# ptranspose.jl
# Adapted from Convex.jl's transpose.jl for use with Convex itself. 
# Returns the partial transpose of a matrix rhoAB, where A and B are 
# of equal dimension. 
#
# All expressions and atoms are subtpyes of AbstractExpr.
# Please read expressions.jl first.
#############################################################################

using Convex
using SCS

import Base.sign
import Base.transpose
import Convex.monotonicity
import Convex.conic_form!
import Convex.curvature

export ptranspose, PTransposeAtom
export sign, curvature, monotonicity, evaluate, conic_form!


type PTransposeAtom <: AbstractExpr
  head::Symbol
  id_hash::UInt64
  children::Tuple{AbstractExpr}
  size::Tuple{Int, Int}

  function PTransposeAtom(x::AbstractExpr)
    children = (x,)
    return new(:ptranspose, hash(children), children, (x.size[2], x.size[1]))
  end
end

function sign(x::PTransposeAtom)
  return sign(x.children[1])
end

function monotonicity(x::PTransposeAtom)
  return (Nondecreasing(),)
end

function curvature(x::PTransposeAtom)
  return ConstVexity()
end

# function evaluate(x::PTransposeAtom)
  # return evaluate(x.children[1])'
# end

# Since everything is vectorized, we simply need to multiply x by a permutation
# matrix such that coeff * vectorized(x) - vectorized(x') = 0
function conic_form!(x::PTransposeAtom, unique_conic_forms::UniqueConicForms)
  if !has_conic_form(unique_conic_forms, x)
    objective = conic_form!(x.children[1], unique_conic_forms)

    sz = get_vectorized_size(x)

    num_rows = x.size[1]
    num_cols = x.size[2]

    I = Array(Int, sz)
    J = Array(Int, sz)

    matrixCoordToVecCoord(i::Integer, j::Integer, n::Integer) = i + (j-1)*n

    """
    Determine the new coordinates `(x, y)` of an entry
    with coordinates `(i, j)` after partial transposition.
    """
    function partialtranseposecoord(i::Integer, j::Integer, n::Integer, l::Integer)
      # calculate which block (i, j) is in
      a = floor(Int, (i-1)/l)+1
      b = floor(Int, (j-1)/l)+1
      xj = j + (a-b)*l
      yj = i + (b-a)*l
      return (xj, yj)
    end

    matrixCoordToVecCoord(ij::Tuple{Int, Int}, n::Integer) = matrixCoordToVecCoord(ij[1], ij[2], n)
    partialtranseposecoord(ij::Tuple{Int, Int}, n::Integer, l::Integer) = partialtranseposecoord(ij[1], ij[2], n, l)

    index = 1
    n = num_rows
    l = round(Int, sqrt(num_rows))
    for i in 1:n
      for j in 1:n
        (a, b) = partialtranseposecoord(i, j, n, l)
        I[index] = matrixCoordToVecCoord(i, j, n)
        J[index] = matrixCoordToVecCoord(a, b, n)
        index += 1
      end
    end

    transpose_matrix = sparse(I, J, 1)

    objective = transpose_matrix * objective
    cache_conic_form!(unique_conic_forms, x, objective)
  end
  return get_conic_form(unique_conic_forms, x)
end

ptranspose(x::AbstractExpr)= PTransposeAtom(x);
