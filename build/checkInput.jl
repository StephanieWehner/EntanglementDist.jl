# Functions for checking the validity of quantum elements such as states, unitarties,...

export isQuantumState, isHermitian, isUnitary, isPPT;

""" `b = isQuantumState(rho)` or `b = isQuantumState(rho,prec)`

Checks whether the input matrix rho is a valid quantum state, that is, a matrix that is positive semidefinite and has trace 1. *prec* is the precision that small eigenvalues are dealt with: If the input matrix has a negative eigenvalue that is no smaller than *-prec*, then this eigenvalue is considered zero.

Returns true/false.

"""
function isQuantumState(rho::AbstractMatrix, prec::Number=0.00001)

	# Check whether it's hermitian
	if !isHermitian(rho, prec)
		print("isQuantumState: Matrix is not Hermitian.\n")
		return false;
	end

	# check for positivity, real is to deal with small numerical errors
	l = real(eig(rho)[1]);
	if (minimum(l) < - prec)
		print("isQuantumState: Matrix is not positive semidefinite.\n")
		return false;
	end

	# Check the normalization
	if (real(trace(rho)) < 1 - prec) && (real(trace(rho)) > 1 + prec)
		print("isQuantumState: Matrix is does not have trace one.\n");
		return false;
	end

	return true;
end

""" 'b = isHermitian(rho)' or 'b = is Hermitian(rho, prec)'

Checks whether the input matrix rho is a Hermitian matrix. Entries are compares up to precision *prec* on average.

Returns true/false.
"""

function isHermitian(rho::AbstractMatrix, prec::Number=0.00001)

	# check whether the matrix is equal to its conjugate transpose
	(d1,d2) = size(rho);
	diff = sum(rho - rho');
	if(real(diff) + imag(diff) > prec * d1)
		return false;
	end
	return true;
end

""" 'b = isUnitary(U)'

Checks whether the input matrix is unitary.

Returns true/false.

"""

function isUnitary(U::AbstractMatrix)

	# Check dimensions
	(d,da) = size(U);
	if (d != da) 
		print("isUnitary: Input is not a square matrix.\n");
		return false;
	end

	# Check unitarity
	if (U*U' != eye(d))
		return false;
	end
	if (U'*U != eye(d))
		return false;
	end

	return true;
end
	
""" `b = isPPT(rhoAB, nA, nB)`

Checks whether the input *rhoAB* is PPT across dimension *nA* and *nB*.

Returns true/false
"""

function isPPT(rho::AbstractMatrix,nA::Number,nB::Number)

	# defined the cutoff when we consider something to be positive
	epsilon = 10.0^(-15);

	# Check whether it's hermitian
	@assert isQuantumState(rho) "Input is not a quantum state."

	# Check whether dimensions match
	(da, db) = size(rho)
	@assert (da ==  nA * nB) "Input does not match given dimensions."

	rhoPT = partialTranspose(rho, nA, nB, 1);

	# check for positivity
	l = eig(rhoPT)[1]
	if (minimum(l) < - epsilon)
		return false
	end

	return true
end
