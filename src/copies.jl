# Produces n copies of a density matrix

export copies

""" `rhoOut = copies(rho,n)`

Returns *rho*^(tensor *n*), i.e., *n* copies of the input *rho*
"""

function copies(rho::AbstractMatrix, n::Int)
	return reduce(kron, repeated(rho, n))
end

function copies(rho::AbstractMatrix, n::Float64)
	return copies(rho, convert(Int, n))
end
