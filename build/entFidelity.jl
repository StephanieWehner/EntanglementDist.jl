# Computes the fidelity with an EPR pair

export entFidelity, eprFidelity

""" `F = entFidelity(rho)`

Returns the fidelity with a maximally entangled state of equal dimension.

"""

function entFidelity(rho::AbstractMatrix)

	@assert isQuantumState(rho) "Input is not a valid quantum state."

	# Construct the maximally entangled state of equal local dimension
	(d1, d2) = size(rho);
	ent = maxEnt(convert(Int,sqrt(d1)));

	# Compute the fidelity with the pure state ent
	out = trace(ent * rho);

	return out;

end

""" `F = eprFidelity(rho)`

For locally two-dimensional states, returns the fidelity with the closest bell state.

"""

function eprFidelity(rho::AbstractMatrix)

	(d, db) = size(rho);
	@assert isQuantumState(rho) "Input is not a valid quantum state."
    	@assert d == 4 "The state is not locally two-dimensional" 

    	M = toBellBasis(rho)

    	# Compute the fidelity with the closest bell state
    	out = max(diag(M)...);

    	return out;
end
