# Kron to work with Convex.jl

using Convex

import Base.kron

export kron

function kron(a::Union{AbstractArray, Convex.Constant, Convex.AbstractExpr}, b::Union{AbstractArray, Convex.Constant, Convex.AbstractExpr})
  Rs = AbstractExpr[]
  for i in 1:size(a)[1]
    Vs = Convex.AbstractExpr[]
    for j in 1:size(a)[2]
      push!(Vs, a[i, j] * b)
    end
    push!(Rs, foldl(hcat, Vs))
  end
  return foldl(vcat, Rs)
end

function kron(a::AbstractMatrix, b::Convex.AbstractExpr)
  return kron(Constant(a), b)
end

function kron(a::AbstractExpr, b::AbstractMatrix)
  return kron(a, Constant(b))
end
