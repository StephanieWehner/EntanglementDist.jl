# Computes the fidelity with an EPR pair

export entFidelity

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
