#
## Purity heuristic
#
#  Next to the PPT condition we use a purity enforcing heuristic to find choi states corresponding
#  to RO operation in which the size of the environment is known. 
# 
# Input:
#    rho   - the input state rhoAB to be distilled, assuming dim(A)=dim(B)
#    dout  - desired dimension of the maximally entangled output state (on one node).
#    nE_A  - size of the environment to be used at node A
#    nE_B  - size of the environment to be used at node B
#    maxIter - maximum number of iterations
#

using Convex
using SCS

export purity

function purity(rho::AbstractArray, dout::Int, nE_A::Int, nE_B::Int, maxIter::Int)

	# Define a tolerance parameter: eigenvalues below that number are considered 0
	tol = (10//1)^(-5);

	# In analogy to the notation used in the paper, set the output dimensions of Ahat and Bhat
        # to contain the (ideally) maximally entangled output state.
  	nAhat = dout;
  	nBhat = nAhat;

	# Dimension of the input state
	(d1,d2) = size(rho);
	din = convert(Int,sqrt(d1));
  	nAprime = din;
  	nBprime = nAprime;

	# Dimension of the flag registers, only 2 possible values: success or failure
  	nF_A = 2;
  	nF_B = 2;

	# Compute the total dimension of the A and B maps
  	nA = convert(Int,nAhat * nF_A * nAprime * nE_A);
  	nB = convert(Int,nBhat * nF_B * nBprime * nE_B);
  	nC = nA * nB;

	# vector of system dimensions: specifies how we order them
	sysDims = [nAhat, nBhat, nF_A, nF_B, nE_A, nE_B, nAprime, nBprime];

	# At the first level of the iteration, no additional constraints will be placed 
  	t_A = 1;
  	t_B = 1;

  	# instantiate direction matrices
  	W_Abar = zeros(nA, nA);
  	W_Bbar = zeros(nB, nB);

  	# do iterations until purity is reached or 100 iterations were done
  	for i = 1:maxIter
    		println("Starting iteration ", i)

		# Solve the optimization problem
    		@time (C_AbarBbar, F, p_succ) = findChoi(rho, W_Abar, W_Bbar, t_A, t_B, nE_A, nE_B, dout)

    		# sometimes Convex.jl fails.
    		# try again in that case
    		j = 1;
    		while C_AbarBbar == zeros((nC, nC))
      			println(" trying again ", j, "...")
      			j += 1
      			@time (C_AbarBbar, F, p_succ) = findChoi(rho, W_Abar, W_Bbar, t_A, t_B, nE_A, nE_B, dout);
    		end

		# Compute the marginal states. Ideally these are pure.
		C_Abar = partialtrace(CABarBbar,  [2, 4, 6, 8], sysDims);
		C_Bbar = partialtrace(CABarBbar,  [1, 3, 5, 7], sysDims);

		# Verify whether the A and B part are pure
    		eA = eigvals(C_Abar);
    		eA = -sort(-real(eA));
    		println(" eA = ", round(eA[1:10], 2));

		eB = eigvals(C_Bbar);
		eB = -sort(-real(eB));
		println(" eB = ", round(eB[1:10],2));

		# Print the fidelity and success probability attained
    		println(" F = ", F)
    		println(" p_succ = ", p_succ)

    		# check with some tolerance if the state is pure
    		# if so, return the found Choi state
    		if (eA <= tol) && (eB <= tol) 
      			return (C_Abar,C_Bbar);
    		end

    		@time (W_Abar, t_A) = findW(C_Abar);
    		@time (W_Bbar, t_B) = findW(C_Bbar);

  	end

  	return(C_Abar, C_Bbar);
end

#
# Solve the original SDP, including the additional constraint on the matrices enforcing the purity heuristic
#
# Input:
#   rho		- input state to be distilled
#   W_Abar	- matrix enforcing the purity heuristic on system A
#   W_Bbar      - matrix enforcing the purity heuristic on system B
#   t_A		- constraint t_A enforcing purity on system A
#   t_B		- constraint t_B enforcing purity on system B
#   nE_A 	- size of the environment on A
#   nE_B 	- size of the environment on B
#   dout	- dimension of the output system (on one side)
#
# Output:
#
#   C		- overall choi states
#   F		- fidelity achieved
#   psucc	- success probability achieved
#
	
function findChoi(rho::AbstractArray, W_Abar::AbstractArray, W_Bbar::AbstractArray, t_A::Number, t_B::Number, nE_A::Number, nE_B::Number, dout::Number)
 
	# verify the input 
	@assert isQuantumState(rho) "expected rho to be a quantum state"
  	@assert 0 <= t_A <= 1
  	@assert 0 <= t_B <= 1

	# In analogy to the notation used in the paper, set the output dimensions of Ahat and Bhat
        # to contain the (ideally) maximally entangled output state.
  	nAhat = dout;
  	nBhat = nAhat;

	# Dimension of the input state
	(d1,d2) = size(rho);
	din = convert(Int,sqrt(d1));
  	nAprime = din;
  	nBprime = nAprime;

	# Dimension of the flag registers, only 2 possible values: success or failure
  	nF_A = 2;
  	nF_B = 2;

	# Compute the total dimension of the A and B maps
  	nA = convert(Int,nAhat * nF_A * nAprime * nE_A);
  	nB = convert(Int,nBhat * nF_B * nBprime * nE_B);
  	nC = convert(Int,nA * nB);

	# vector of system dimensions: specifies how we order them
	sysDims = [nAhat, nBhat, nF_A, nF_B, nE_A, nE_B, nAprime, nBprime];

	# At the first level of the iteration, no additional constraints will be placed 
  	t_A = 1;
  	t_B = 1;

	# Compute projector onto the 11 flag indicating success
  	e0 = eVec(2, 1);
  	e1 = eVec(2, 2);

	e11 = kron(e1,e1);
  	e11 = e11 * e11';

	# Compute a projector onto the Epr pair
	epr = maxEnt(nAhat);

	# Set up the optimization problem. This SDP optimizes (For now XXX) the product of the 
	# fidelity and the success probability

	# The overall choi state
  	C_AbarBbar = Semidefinite(nC);

	# The partial trace components on A and B
  	C_Abar = partialtrace(C_AbarBbar, [2, 4, 6, 8], sysDims);
  	C_Bbar = partialtrace(C_AbarBbar, [1, 3, 5, 7], sysDims);

  	C_AprimeBprime = partialtrace(C_AbarBbar, [1, 2, 3, 4, 5, 6] , sysDims);

	# Define the objective function: maximize F * psucc
  	problem = maximize(trace(kron(kron(epr,eye(nE_A * nE_B)),kron(e11,rho')) * C_AbarBbar));

  	# constraints from purity
  	problem.constraints += [trace(C_Abar * W_Abar) <= t_A];
  	problem.constraints += [trace(C_Bbar * W_Bbar) <= t_B];

  	# constraints for Choi state
  	problem.constraints += [C_AprimeBprime == eye(din^2)/din^2];
  	problem.constraints += [trace(C_AbarBbar) == 1];

  	# PPT constraint on B
	print("Before PPT\n");
  	problem.constraints += partialtranspose(C_AbarBbar, [2, 4, 6, 8], sysDims) in :SDP;

	print("Running SDP\n");
  	solve!(problem, SCSSolver(verbose = true));

  	p_succ = nAprime * nBprime  * trace(kron(kron(eye(nAhat * nBhat * nE_A * nE_B),e11),rho') * CAbarBbar.value);
  	F = problem.optval * nAprime * nBprime/ p_succ;

  	return (C_AbarBbar.value, F, p_succ);
end

#
## findW
#
# Solves the iterative SDP of the purity heuristic, trying to force the relevant part of the input
# to be pure by searching for an updated W. The new t will be the new SDP value
#

function findW(C_opt::AbstractArray, verbose::Bool = false)

	# Make sure the input is actually a quantum state
  	@assert isQuantumState(C_opt) "expected C_opt to be a quantum state"

	# Recover the size of the input
  	nC = size(C_opt)[1];

	# Search for the new W
  	W_Xbar = Semidefinite(nC);

	# Define the objective 
  	problem = minimize(trace(C_opt' * W_Xbar));

	# Normalize W
  	problem.constraints += ([W_Xbar <= eye(nC)]);

	# If W was pure itself, then it'll have trace 1
  	problem.constraints += ([trace(W_Xbar) == nC - 1]);

  	solve!(problem, SCSSolver(verbose = verbose));
  	return (W_Xbar.value, problem.optval)
end

