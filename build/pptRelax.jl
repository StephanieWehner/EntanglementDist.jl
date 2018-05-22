# Implements the method of PPT relaxations of MX operations

using Convex
using SCS

export pptRelax, pptRelaxCopies;

""" `(problem, F, p_succ) = pptRelax(rho,n,k,delta, verbose, eps, max_iters)`

Implements the PPT relaxation computing an upper bound on the fidelity achievable using measure and exchange operations for a given input state and fixed succes probability.

Inputs:
- *rho* quantum state to be distilled
- *n* number of qubits on one side, assuming dimensions nA=nB =2^n in the input.
- *k* desired output dimension of the maximally entangled state
- *delta* desired success probability
- *verbose* flag indicating whether the SCS solver should be verbose (default false)
- *eps* desired precision (default 1e-4) in SCS solver
- *max_iters* maximum number of iterations in SCS solver (default 2e4)

Outputs:
- *problem* SDP object from Convex.jl allowing further processing.
- *F* Fidelity achieved
- *p_succ* success probability attained in SDP (should equal *delta*)
"""

function pptRelax(rho::AbstractMatrix, n::Int, k::Int, delta::Real; verbose::Bool = false, eps::Number = 1e-4, max_iters = 2e4)
	# Presume the input dimensions are equally divided and we're using qubits
	nA = 2^n;
	nB = nA;

	return pptRelax(rho, nA, nB, k, delta, verbose = verbose, eps = eps, max_iters = max_iters)
end

""" `(problem, F, p_succ) = pptRelax(rho,nA, nB, k,delta, verbose, eps, max_iters)`

Implements the PPT relaxation computing an upper bound on the fidelity achievable using measure and exchange operations for a given input state and fixed succes probability.

Inputs:
- *rho* quantum state to be distilled on A and B
- *nA* dimension of the A system
- *nB* dimension of the B system
- *k* desired output dimension of the maximally entangled state
- *delta* desired success probability
- *verbose* flag indicating whether the SCS solver should be verbose (default true)
- *eps* desired precision (default 1e-4) in SCS solver
- *max_iters* maximum number of iterations in SCS solver (default 2e4)

Outputs:
- *problem* SDP object from Convex.jl allowing further processing.
- *F* Fidelity achieved
- *p_succ* success probability attained in SDP (should equal *delta*)
"""

function pptRelax(rho::AbstractArray, nA::Int, nB::Int, k::Int, delta::Number; verbose::Bool = true, eps::Number = 1e-4, max_iters::Number = 2e4)

	# Check whether rho is a quantum state
	@assert isQuantumState(rho) "The input state must be a valid quantum state."

	# Check whether dimensions match
	(d, db) = size(rho);
	@assert d == nA * nB "Input state doesn't match given dimensions." , d , " != ", nA*nB

	# Check whether delta is a valid probability
	@assert delta <= 1 "Success probability must be less than one."
	@assert 0 <= delta "Success probability must be positive."

	# define the identity matrix on the whole space
	id = eye(d);

	# define the variables - if the input is real it is sufficient to optimize over reals
	if real(rho) == rho
		D = Semidefinite(d);
		E = Semidefinite(d);
	else
		# complex variable support advertises HermitianSemidefinite
		# but it's buggy
		D = ComplexVariable(d,d);
		E = ComplexVariable(d,d);
	end

	# define the objective
	tv = Variable(1);
	problem = maximize(tv);
	problem.constraints += tv == trace(nA * nB * D * rho.');

	problem.constraints += D in :SDP;
	problem.constraints += E in :SDP;
	problem.constraints += (id/d - (D+E)) in :SDP;

	# Constraints relating to the PPT Condition
	EPT = ptranspose(E);
	DPT = ptranspose(D);
	problem.constraints += (DPT + EPT/(k+1)) in :SDP;
	problem.constraints += (- DPT + EPT/(k-1)) in :SDP;
	problem.constraints += (id/d - (DPT+EPT)) in :SDP;

	# Constraint coming from the success probability
	problem.constraints += nA * nB * trace(rho.'*(D+E)) == delta;

	solve!(problem, SCSSolver(verbose = verbose, eps = eps, max_iters = max_iters),verbose=verbose);
	p_succ = nA * nB * trace(rho.'*(D.value + E.value));
	F = problem.optval / p_succ;

	return (problem, F, p_succ)
end

""" `(problem, F, p) = pptRelaxCopies(rho, n, nA, nB, k, delta, verbose, eps, max_iters)`

Calls pptRelax for an input state of the form *rho&*^(tensor *n*), reordered appropriately.

"""

function pptRelaxCopies(rho::AbstractMatrix,n::Int,nA::Int,nB::Int,k::Int,delta::Real; verbose::Bool = false, eps::Number = 1e-4, max_iters::Number = 2e4)

	@assert n > 1 "Need at least two copies."

	# copy rho n times
	rhoNew = copies(rho, n);

  	# sort all the entries, such that all the systems on A are first and all
	# the systems on B are second
	rhoSorted = sortAB(rhoNew, nA, n)

	return pptRelax(rhoSorted,nA^n,nB^n, k, delta, verbose=verbose, eps=eps,max_iters=max_iters);
end
