# use Convex to get access to its types and functions
using Convex

# import reshape to extend it.
import Base.reshape

export partialtransposematrix, partialtranspose, partialtransposenaive

# switch between coordinates of a matrix and its vectorised form
vecCoordToMatrixCoord(k::Integer, n::Integer) = (1 + (k-1)%n, round(Int, 1+floor((k-1)/n)%n))
matrixCoordToVecCoord(i::Integer, j::Integer, n::Integer) = i + (j-1)*n


#Determine the new coordinates `(x, y)` of an entry
#with coordinates `(i, j)` after partial transposition.

function partialtranseposecoord(i::Integer, j::Integer, n::Integer, l::Integer)
  @assert n % l == 0 "l should divide n"
  @assert n ≥ l "n should be greater than l"

  # calculate which block (i, j) is in
  a = floor(Int, (i-1)/l)+1
  b = floor(Int, (j-1)/l)+1

  x = j + (a-b)*l
  y = i + (b-a)*l

  @assert x > 0
  @assert y > 0
  @assert n ≥ x
  @assert n ≥ y

  return (x, y)
end


#Construct a permutation matrix for a `n × n` matrix
#in vector form where we transpose in `l × l` blocks.

function partialtransposematrix(n::Integer, l::Integer)
  @assert n % l == 0 "l should divide n"
  @assert n ≥ l "n should be greater than l"

  # Julia will iterate over I and J simulatiously
  # and set (I[k], J[k]) to 1 in the returned matrix
  # The value in I corresponds to the vectorised coordinates
  # before partial transposition, whereas the value in J
  # corresponds to the coordinates after partial transposition.
  IM = Array(Int, n^2)
  JM = Array(Int, n^2)
  index = 1
  for i in 1:n
    for j in 1:n
      (x, y) = partialtranseposecoord(i, j, n, l)
      IM[index] = matrixCoordToVecCoord(i, j, n)
      JM[index] = matrixCoordToVecCoord(x, y, n)
      index += 1
    end
  end

  return sparse(IM, JM, 1)
end

""" `rhoOut = partialtranspose(rho)` or `rhoOut = partialtranspose(rho, sys)` or `rhoOut = partialtranspose(rho, sys, dim)`

Applies the partial tanspose to the matrix *rho*.

Inputs:
- *rho* input matrix
- *systems* (optional) single system or a vector of systems to be transposed. By default it is equal to 2. If *dim* is not provided, the only two possible values of *systems* are 1 or 2.
- *dim* (optional) vector with dimensions of the systems. By default it assumes two systems of equal dimensions.

Outputs:
- *rhoOut* matrix *rho* with partial transpose applied
"""
function partialtranspose(ρ::Union{AbstractArray, AbstractExpr}, systems::Vector, dim::Vector)
  result = ρ
  for sys in systems
    result = partialtranspose(result, sys, dim)
  end
  return result
end

#Apply the partial tanspose to a single sytem. If no
#dimensions are specified, it will assume there are
#two systems of identical dimension. If no system is
#specified, it will transpose the last system. This
#will always be the second system, because it's not
#possible to leave out the system while providing
#the dimensions.

function partialtransposesystem(ρ::Union{AbstractArray, AbstractExpr}, sys::Integer = 2, dim::Vector = [])
  # do naive partial transposition
  if dim == []
    l = round(Int, sqrt(size(ρ)[1]))
    if sys == 2
      return partialtranspose(ρ, l, naive = true)
    else
      dim = [l; l]
    end
  end


  @assert 1 ≤ sys ≤ length(dim) "System doesnt exist, expected sys to be between 1 and ", length(dim), ". Is your dim vector correct?"

  # permute systems such that the system to apply PT is last
  perm = collect(1:length(dim))
  deleteat!(perm, sys)
  push!(perm, sys)

  l = dim[sys]
  x = permutesystems(ρ, perm, dim)
  y = partialtransposenaive(x, l)
  z = permutesystems(y, perm, dim)
  return z
end


#Apply partial transpose on systems if naive is false.
#Otherwise do the naive partial transpose in blocks
#of `l × l`, where `l` takes the value `systems`.
#This is equivalent to applying the partial
#transepose to the last system if that particular system
#has dimension `l`.

function partialtranspose(ρ::Union{AbstractArray, AbstractExpr}, systems::Integer = 2, dim::Vector = []; naive::Bool = false)
  if naive
    return partialtransposenaive(ρ, systems)
  end

  return partialtransposesystem(ρ, systems, dim)
end



#This does the naive version of the partial transpose
#that I wrote first it simply applies a transpose in blocks
#of `l × l`. This is equivalent to applying the partial
#transepose to the last system if that particular system
#has dimension `l`.

function partialtransposenaive(ρ::Union{AbstractExpr, AbstractArray, Convex.Variable}, l::Integer)
  # if the input is sparse, keep the output sparse as well
  if !isa(ρ, Convex.Variable) && !isa(ρ, AbstractExpr) && issparse(ρ)
    vectorised = sparsevec(ρ)
  else
    vectorised = vec(ρ)
  end
  n = size(ρ)[1]

  @assert n % l == 0 "l should divide n"
  @assert n ≥ l "n should be greater than l"

  return reshape(partialtransposematrix(n, l) * vectorised, n, n)
end

# some aliases for convenience
matrixCoordToVecCoord(ij::Tuple{Int, Int}, n::Integer) = matrixCoordToVecCoord(ij[1], ij[2], n)
partialtranseposecoord(ij::Tuple{Int, Int}, n::Integer, l::Integer) = partialtranseposecoord(ij[1], ij[2], n, l)
reshape(x::AbstractExpr, mn::Tuple{Int, Int}) = reshape(x, mn[1], mn[2])
