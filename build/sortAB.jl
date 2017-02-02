# Given rhoAB^(tensor n), sort the state so all A parts come first

export sortAB

""" `rhoSorted = sortAB(rhoBig, d, n)`

Given a state of the form rhoBig = rhoAB^(tensor n), rearrange, so all A parts are first, followed by all B parts. *d* is the dimension of A, assumed to be the same as for B. *n* is the number of copies. 
"""

function sortAB(rhoBig::AbstractMatrix, d::Int, n::Int)

	# Check whether the input is a valid state and the dimensions match
	@assert isQuantumState(rhoBig) "The input is not valid quantum state."

	(d1,d2) = size(rhoBig);
	@assert d1 == (d^2)^n "Input dimensions don't match."


	rhoOut = permutesystems(rhoBig, [collect(0:2:2(n - 1)) + 1; collect(2:2:2n)], dim = fill(d, 2n)) 	

	return rhoOut;
end
