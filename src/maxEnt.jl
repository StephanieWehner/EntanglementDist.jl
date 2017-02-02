# Produces a maximally entangled state as vector or density matrix

export maxEnt, maxEntVec;

""" `vec = maxEnt(d)`

Returns a projector (density matrix) onto the maximally entangled state of a given dimension *d*.

"""
function maxEnt(d::Number)
	v = maxEntVec(d);
	return (v*v')/d;
end

""" `vec = maxEntVec(d)`

Returns a the maximally entangled state of a given dimension *d* as a vector.

"""
function maxEntVec(d::Number)
	v = zeros(d^2);
	for j = 1:d
		v += kron(eVec(d,j),eVec(d,j));
	end
	return v;
end
