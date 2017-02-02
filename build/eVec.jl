# Makes a vector in the standard basis

export eVec

""" `vec = eVec(d, j)`

Returns a sparse vector with all 0's except at the position j.
"""

function eVec(dim::Int, j::Int)
	@assert dim >= 0 "Dimension must be positive."
	@assert j >= 0 "Index must be positive."
	@assert j <= dim "Index must be less than dimension."
	return sparsevec([j], [1], dim)
end
