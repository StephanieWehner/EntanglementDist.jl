# Generate random unitaries and density matrices for testing

export random_unitary, random_densitymatrix

""" `U = random_unitary(dim, real=true/false, complex=true/false)`

Generate a random unitary of dimension *dim*. Depending on *real* and *complex* the matrix has real and complex parts.

"""

function random_unitary(dim::Integer; real::Bool = true, complex::Bool = false)

  	@assert dim > 1 "Dimension must be at least 2."

  	generated = randn(dim, dim)

  	if complex || !real
    		generated += im * randn(dim, dim)
  	end

	# Compute the QR decomposition
  	(Q, R) = qr(generated)

	# Get the upper diagonal part of R
  	F = diagm(diag(R))

  	R = sign(diagm(diag(R)))
  	R[find((x)->  x == 0, R)] = 1
  	F = diagm(diag(R))

  	return Q*F
end

""" `rho = random_densitymatrix(dim)`

Generate a random density matrix in dimension *dim*.

"""

function random_densitymatrix(dim::Integer)
 
	@assert dim > 1 "Dimension must be at least 2."

        # Draw random numbers
	rList = rand(dim);

	# Normalize to form a distribution
	s = sum(rList);
	nList = rList/s;
	
	# Construct a matrix with the numbers on the diagonal
	M = diagm(nList);

	# Draw a random unitary
   	U = random_unitary(dim, real = true, complex = true)

	# Roate the diagonal matrix
	return(U*M*U');
end
