# Implement the more stringent relaxation of PPT plus a number of symmetric extensions

using Convex
using SCS


export pptRelax2Ext, pptRelax2ExtCopies
export pptRelax1Ext, pptRelax1ExtCopies

""" `(problem, F, psucc) = pptRelax1Ext(rhoAB, n, k, delta, verbose=true/false, eps=1e^-4)`

Implements the PPT relaxation plus 1 symmetric extension for distillable entanglement with fixed desired success probability.

Inputs:
- *rhoAB* quantum state to be distilled
- *n* number of qubits on one side, assuming dimensions nA=nB =2^n in the input.
- *k* local dimension of the maximally entangled output state
- *delta* desired success probability

Outputs:
- *problem* problem object given by Convex,
- *F* fidelity bound
- *psucc* success probability actually attained (should equal *delta* up to numerical imprecisions)
"""

function pptRelax1Ext(rho::AbstractMatrix, n::Number, k::Number, delta::Number; verbose::Bool = true, eps::Number = 1e-4)
    # Presume the input dimensions are equally divided and we're using qubits
  nA = 2^n;
  nB = nA;

  return pptRelax1Ext(rho, nA, nB, k, delta; verbose=verbose, eps=eps)
end

""" `(problem, F, psucc) = pptRelax1Ext(rhoAB, nA, nB, k, delta, verbose=true/false, eps=1e^-4)`

Implements the PPT relaxation plus 1 symmetric extension for distillable entanglement with fixed desired success probability.

Inputs:
- *rhoAB* quantum state to be distilled
- *nA*  dimension of the A system
- *nB*  dimension of the B system
- *k*   local dimension of the maximally entangled output state
- *delta* desired success probability

Outputs:
- *problem* problem object given by Convex,
- *F* fidelity bound
- *psucc* success probability actually attained (should equal *delta* up to numerical imprecisions)
"""

function pptRelax1Ext(rho::AbstractMatrix, nA::Number, nB::Number, k::Number, delta::Number; verbose::Bool = true, eps::Number = 1e-4)

  # Check whether rho is a quantum state
  @assert isQuantumState(rho) "Input is not a quantum state."

  # Check whether dimensions match
  (d, db) = size(rho)
  @assert d == nA*nB "Input dimensions don't match."

  #dimensions of the symmetric and antisymmetric subspaces:
  dimSym = Int( (1/2) * k^2 * nA * (k * nA + 1) * nB )
  dimAsym = Int( (1/2) * k^2 * nA * (k * nA - 1) * nB )

  # define the variable on the symmetric subspace Ws
  if(real(rho) == rho)
      Ws = Semidefinite(dimSym);
  else
      Ws= ComplexVariable(dimSym);
  end

  #Since we extend in the symmetric part, the antisymmetric part is zero
  Wa = spzeros(dimAsym, dimAsym)

  #Construct the change of basis matrix
  Pcb = getPcb(nA, nB, k)

  # output state
  epr = maxEnt(k);

  # dimensions of W
    dims = [k, nA, # Ahat_1 A'_1
            k, nA, # Ahat_2 A'_2
            k, nB] # Bhat B'

  #Define direct sum
  dirsum(x, y) = [x zeros(size(x)[1], size(y)[2]); zeros(size(y)[1], size(x)[2]) y]

  #transform the variable from the symmetric-antisymmetric basis to the standard basis
  W_A1A2B = Pcb * dirsum(Ws, Wa) * Pcb'

  # Choi with first A system
  W_A1B = partialtrace(W_A1A2B, [3, 4], dims)

  # Choi with second A system
  W_A2B = partialtrace(W_A1A2B, [1, 2], dims)

  # define the objective
  X =  Pcb' * nA * nB * kron( eye(k * nA), permutesystems(kron(epr, transpose(rho)), [1, 3, 2, 4], dim = [k, k, nA, nB]) ) * Pcb
  problem = maximize(trace(X[1:dimSym, 1:dimSym] * Ws))
  problem.constraints += Ws in :SDP;

  #Define probability of success
  Y = Pcb' * nA * nB * kron( eye(k * nA), permutesystems(kron(eye(k^2), transpose(rho)), [1, 3, 2, 4], dim = [k, k, nA, nB]) ) * Pcb
  p_succ = trace(Y[1:dimSym, 1:dimSym] * Ws)

  # define constraints
  problem.constraints += [
    eye(d)/d - partialtrace(W_A2B , [1, 3], [k, nA, k, nB]) ⪰ 0

    # constrain probability of success
    p_succ == delta

    #PPT condition
    ptranspose(W_A2B) ⪰ 0
  ]

  # Maximize P
  solve!(problem, SCSSolver(verbose = verbose, eps = eps))

  # Output
  Psuccess = nA * nB * trace(permutesystems(kron(eye(k^2), transpose(rho)), [1, 3, 2, 4], dim = [k, k, nA, nB]) * partialtrace(Pcb * dirsum(Ws.value, Wa) * Pcb', [1, 2], dims))
  F = problem.optval/Psuccess
  return (problem, F, Psuccess)
end

""" `(problem, F, p) = pptRelax1ExtCopies(rho, n, nA, nB, k, delta, verbose, eps, max_iters)`

Calls pptRelax1Ext for an input state of the form *rho&*^(tensor *n*), reordered appropriately.

"""
function pptRelax1ExtCopies(rho::AbstractMatrix,n::Number, nA::Number, nB::Number, k::Number, delta::Number; verbose::Bool = true, eps::Number = 1e-4)

@assert n > 1 "Need at least two copies."

  # copy rho n times
  rhoNew = copies(rho, n);

    # sort all the entries, such that all the systems on A are first and all
  # the systems on B are second
  rhoSorted = sortAB(rhoNew, nA, n)

  return pptRelax1Ext(rhoSorted, nA^n, nB^n, k, delta; verbose=verbose, eps=eps)
end


""" `(problem, F, psucc) = pptRelax2Ext(rhoAB, n, k, delta, verbose=true/false, eps=1e^-4)`

Implements the PPT relaxation plus 2 extensions for distillable entanglement with fixed desired success probability.

Inputs:
- *rhoAB* quantum state to be distilled
- *n* number of qubits on one side, assuming dimensions nA=nB =2^n in the input.
- *k* local dimension of the maximally entangled output state
- *delta* desired success probability

Outputs:
- *problem* problem object given by Convex,
- *F* fidelity bound
- *psucc* success probability actually attained (should equal *delta* up to numerical imprecisions)
"""

function pptRelax2Ext(rho::AbstractMatrix, n::Number, k::Number, delta::Number; verbose::Bool = true, eps::Number = 1e-4)
    # Presume the input dimensions are equally divided and we're using qubits
  nA = 2^n;
  nB = nA;

  return pptRelax2Ext(rho, nA, nB, k, delta; verbose=verbose, eps=eps)
end


""" `(problem, F, psucc) = pptRelax2Ext(rhoAB, nA, nB, k, delta, verbose=true/false, eps=1e^-4)`

Implements the PPT relaxation plus 2 extensions for distillable entanglement with fixed desired success probability.

Inputs:
- *rhoAB* quantum state to be distilled
- *nA*  dimension of the A system
- *nB*  dimension of the B system
- *k*   local dimension of the maximally entangled output state
- *delta* desired success probability

Outputs:
- *problem* problem object given by Convex,
- *F* fidelity bound
- *psucc* success probability actually attained (should equal *delta* up to numerical imprecisions)
"""

function pptRelax2Ext(rho::AbstractMatrix, nA::Number, nB::Number, k::Number, delta::Number; verbose::Bool = true, eps::Number = 1e-4)

  	# Check whether rho is a quantum state
  	@assert isQuantumState(rho) "Input is not a quantum state."

  	# Check whether dimensions match
  	(d, db) = size(rho);
  	@assert d == nA*nB "Input dimensions don't match."

  	# define the variables W
	if(real(rho) == rho)
  		W_total = Semidefinite(k^4 * nA^3 * nB);;
	else
		ms = k^4 * nA^3 * nB;
  		W_total= ComplexVariable(ms,ms);
	end
	

	# output state
	epr = maxEnt(k);

  	# dimensions of W_total
  	dims = [k, nA, # Ahat_1 A'_1
          	k, nA, # Ahat_2 A'_2
          	k, nA, # Ahat_3 A'_3
          	k, nB]; # B_hat B_'

  	# Choi with first A system
  	W_A1B = partialtrace(W_total, [3,4, 5, 6], dims);

  	# Choi with second A system
  	W_A2B = partialtrace(W_total, [1, 2, 5, 6], dims);

	# Choi with third A system
  	W_A3B = partialtrace(W_total, [1, 2, 3, 4], dims);

  	# define the objective
	tv = Variable(1);
  	problem = maximize(tv);
	problem.constraints += tv == nA * nB * trace(permutesystems(kron(epr, transpose(rho)), [1, 3, 2, 4], dim = [k, k, nA, nB]) * W_A1B);
	problem.constraints += W_total in :SDP;

  	# define constraints
  	problem.constraints += [
    		eye(d)/d - partialtrace(W_A1B , [1, 3], [k, nA, k, nB]) in :SDP

    		# constrain probability of success
    		nA * nB * trace(permutesystems(kron(eye(k^2), transpose(rho)), [1, 3, 2, 4], dim = [k, k, nA, nB]) * W_A1B) == delta
    		

		# PPT condition
    		ptranspose(W_A1B) in :SDP

    		# Choi should be equal if we trace out one of the extensions
    		W_A1B == W_A2B
    		W_A2B == W_A3B
  	]

  	# Maximize P
  	solve!(problem, SCSSolver(verbose = verbose, eps = eps))

  	# Output
  	Psuccess = nA * nB * trace(permutesystems(kron(eye(k^2),transpose(rho)), [1, 3, 2, 4], dim = [k, k, nA, nB]) * partialtrace(W_total.value, [3,4,5,6], dims))
  	F = problem.optval/Psuccess
  	return (problem, F, Psuccess)
end

""" `(problem, F, p) = pptRelax2ExtCopies(rho, n, nA, nB, k, delta, verbose, eps, max_iters)`

Calls pptRelax2Ext for an input state of the form *rho&*^(tensor *n*), reordered appropriately.

"""
function pptRelax2ExtCopies(rho::AbstractMatrix,n::Number, nA::Number, nB::Number, k::Number, delta::Number; verbose::Bool = true, eps::Number = 1e-4)

@assert n > 1 "Need at least two copies."

  # copy rho n times
  rhoNew = copies(rho, n);

    # sort all the entries, such that all the systems on A are first and all
  # the systems on B are second
  rhoSorted = sortAB(rhoNew, nA, n)

  return pptRelax2Ext(rhoSorted, nA^n, nB^n, k, delta; verbose=verbose, eps=eps)
end

""" `Pcb = getPcb(nA,nB,k)`

Generates a change of basis matrix that takes us from the basis of the symmetric-antisymmetric subspace to the standard basis

Inputs:
- *nA*  dimension of the A system
- *nB*  dimension of the B system
- *k*   local dimension of the maximally entangled output state

Outputs:
- *Pcb* the change of basis matrix
"""
function getPcb(nA::Number, nB::Number, k::Number)
  
  #generate basis vectors on the k * nA dimensional system Ahat_i A'_i
  v0 = spzeros(k * nA, 1)
  v = Any[]
  for i = 1:k * nA
    x = copy(v0)
    x[i] = 1
    push!(v, x)
  end

  #generate basis vectors on the k * nB dimensional system Bhat B'
  w0 = spzeros(k * nB, 1)
  w = Any[]
  for i = 1:k * nB
    x = copy(w0)
    x[i] = 1
    push!(w, x)
  end

  #Generate basis vectors of the symmetric subspace
  #Product vectors
  u = Any[]
  for i = 1:k * nA
    x = kron( v[i] , v[i] )
    push!(u, x)
  end

  # Entangled vectors
  for j = 1:k * nA - 1
    for i = 1:k * nA-j
      x = ( 1/sqrt( 2 ) ) * ( kron( v[j] , v[i+j] ) + kron( v[i+j] , v[j] ) )
      push!(u, x)
    end
  end

  #Generate basis vectors of the antisymmetric subspace
  # Entangled vectors
  for j = 1: k * nA - 1
    for i = 1:k * nA-j
      x = ( 1/sqrt( 2 ) ) * ( kron( v[j] , v[i+j] ) - kron( v[i+j] , v[j] ))
      push!(u, x)
    end
  end

  #Generate the change of basis matrix
  Pcb = spzeros(k^3 * nA^2 * nB, k^3 * nA^2 * nB)
  for i = 1:k^2 * nA^2
    for j = 1:k * nB
      Pcb[:, k * nB * (i - 1) + j] = kron(u[i], w[j] )
    end
  end

  return Pcb
end
