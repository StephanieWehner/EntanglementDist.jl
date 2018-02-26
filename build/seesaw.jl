# Perform a seesaw optimization searching for a good distillation scheme
# in the vicinity of a known one

using Convex
using SCS

export seesaw, seesawAlice, seesawBob;

"""
`(CA,CB, F, psucc) = seesaw(rhoAB, nA, nB, k, startCA, startCB)` or
`(CA,CB, F, psucc) = seesaw(rhoAB, nA, nB, k, startCA, startCB, P)` or
`(CA,CB, F, psucc) = seesaw(rhoAB, nA, nB, k, startCA, startCB, P, max_iter)` or
`(CA,CB, F, psucc) = seesaw(rhoAB, nA, nB, k, startCA, startCB, P, max_iter, delta)`

Performs a seesaw optimization searching for a good scheme to distill *rhoAB*
into a maximally entangled state of local dimension *k*. Good here means that we want to maximize the fidelity of the output state to the maximally entangled state, where the success probability should be the same as the existing one - or, if supplied, it should achieve a specific success probability *delta*.

Alice's map is assumed to be from
A -> hatA, F_A, and Bob's from B -> hatB, F_B, where F_A and F_B are the flag registered that yield success or failure. After the map *P* is measured on the flag registeres to determine success where *P* is the projector determining success. Example P=|11><11| on the flag registers F_A and F_B.

Input:
- *rhoAB* state to be distilled
- *nA* dimension of A
- *nB* dimension of B
- *k*  local dimension of the desired maximally entangled output state
- *startCA* choi state of Alice's map to start from (ordering of systems: hatA (output state), F_A (flag register), A (input))
- *startCB* choi state of Bob's map to start from (ordering of systems: hatB (output state), F_B (flag register), B (input))
- *P* (optional) projector onto success: Default is |11><11|
- *max_iter* (optional) maximum number of iterations (default 10)
- *delta* (optional) desired success probability if different from the starting map.

Output
- *CA* new map for Alice as choi state
- *CB* new map for Bob as choi state
- *F* fidelity obtained
- *psucc* success probability attained
"""

function seesaw(
	rho::AbstractMatrix,
	nA::Int, nB::Int, k::Int,
	startCA::AbstractMatrix, startCB::AbstractMatrix,
	P::AbstractMatrix = [0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 1],
	max_iter::Number=10,
	delta::Number=0)

	# Validate the essential elements of the input
	@assert isQuantumState(rho) "Input is not a valid quantum state."
	nTot = nA * nB;
	(d1,d2) = size(rho);
	@assert d1 == nTot "Input state doesn't match indicated dimensions."
	@assert delta <= 1 "Success probability cannot exceed 1."
	@assert delta >= 0 "Success probability must be positive."

	# The output state we want to produce
	want = maxEnt(k);

	# Dimensions of the systems in ordering hatA, FA, A, hatB, FB, B
	dimSys = [k, k, 2, 2, nA, nB];

	# The current probability of success
				M = permutesystems(kron(kron(eye(k^2),P),rho'), [1,3,5,2,4,6], dimSys);
	psuccOrig = nA * nB * trace(M * kron(startCA, startCB));
	bestP = psuccOrig;

	# The current fidelity achieved
				M = permutesystems(kron(kron(want,P),rho'), [1,3,5,2,4,6], dimSys);
	FOrig = nA * nB * trace(M * kron(startCA, startCB))/psuccOrig;
	bestF = FOrig;

	print("\n\nSEESAW Optimization\n\n");
	print("Starting with fidelity F=",round(FOrig,3), " and p_succ=", round(psuccOrig,3),"\n\n");

	# Set the desired probability of success
	if delta > 0
		# Desired delta is supplied: make that the desired acceptable succes probability
		fixedSucc = delta;
	else
		# Delte is not supplied: we want to keep the success probability of the starting scheme
		fixedSucc = psuccOrig;
	end


	# Iteration to search for a better scheme
	currCA = startCA;
	currCB = startCB;

	for j = 1:max_iter

		# A flag whether we achieved an improvement on either side
		betterA = false;
		betterB = false;

		# Fix Bob and find a better map on Alice's side
		(newCA, FA, psuccA) = seesawAlice(rho, nA, nB, k, currCB, P, fixedSucc);

		if (FA > bestF)
			# We improved the fidelity while achieving the desired success prob.
			# Update Alice's map
			currCA = newCA;
			betterA = true;
			bestF = FA;
			bestP = psuccA;
		end

		# Fix Alice and find a better map on Bob's side
		(newCB, FB, psuccB) = seesawBob(rho, nA, nB, k, currCA, P, fixedSucc);

		if (FB > bestF)
			# We improved the fidelity while maintaining the desired success prob.
			# Update Bob's map
			currCB = newCB;
			betterB = true;
			bestF = FB;
			bestP = psuccB;
		end

		# Check whether any improvement was achieved. Otherwise just return
		if !(betterA || betterB)
			print("\nStopping early after ",j," iterations.\n");
			print("Final fidelity F=",round(bestF,3), " and psucc=", round(bestP,3),"\n");
			return(currCA, currCB, bestF, bestP);
		end
	end
	print("\nFinal fidelity F=",round(bestF,3), " and psucc=", round(bestP,3),"\n");

	return(currCA, currCB, bestF, bestP);
end

""" `(newCA, FA, psuccA) = seesawAlice(rhoAB, nA, nB, k, CB, P, fixedSucc)`

Searches for a better distillation map for Alice, optimizing the product of the fidelity and success probability, for a fixed desired success probability. Called from seesaw.

Input (assumed to be valid):
- *rhoAB* state to be distilled
- *nA* dimension of A
- *nB* dimension of B
- *k* local dimension of desired max. entangled state
- *CB* choi state of Bob's map
- *P* measurement indicating success
- *fixedSucc * desired probability of success

Output:
- *newCA* new map for Alice as choi state
- *F* fidelity of new scheme
- *psucc* success probability of new scheme
"""

function seesawAlice(
	rho::AbstractMatrix,
	nA::Int, nB::Int, k::Int,
	CB::AbstractMatrix,
	P::AbstractMatrix,
	fixedSucc::Number)

	# Dimension of the maps
	d = nA * 2 * k;

	# Desired output state
	want = maxEnt(k);

	# Check whether we want a real or complex optimization problem
	if (real(rho) == rho) && (real(CB) == CB)
		CA = Semidefinite(d);
	else
		CA = ComplexVariable(d,d);
	end

	# Dimensions in ordering hatA, Fa, A,...
	dimSys = [k, k, 2, 2, nA, nB];

	# Set up the optimization problem
	# The current probability of success
				Mwant = permutesystems(kron(kron(want,P),rho'), [1,3,5,2,4,6], dimSys);
	tv = Variable(1);
	problem = maximize(tv);
	problem.constraints += tv == nA * nB * trace(Mwant * simplekron(CA, CB));
	problem.constraints += CA in :SDP;

	# Constraints
	# CA should be a choi state
	problem.constraints += ([trace(CA) == 1]);

	CAin = partialtrace(CA,[1],[2 * k, nA]);
	problem.constraints += ([CAin == eye(nA)/nA]);

	# We want a fixed success probability
				Mid = permutesystems(kron(kron(eye(k^2),P),rho'), [1,3,5,2,4,6], dimSys);
	problem.constraints += ([nA * nB * trace(Mid * simplekron(CA, CB)) == fixedSucc]);

	# Solve the SDP
	solve!(problem, SCSSolver(verbose=false));

	# Compute the new success probability
	psuccNew = nA * nB * trace(Mid * kron(CA.value, CB));

	# The current fidelity achieved
	Fnew = nA * nB * trace(Mwant * kron(CA.value, CB))/psuccNew;


	return(CA.value, Fnew, psuccNew);
end



""" `(newCB, FB, psuccB) = seesawBob(rhoAB, nA, nB, k, CA, P, fixedSucc)`

Searches for a better distillation map for Bob, optimizing the product of the fidelity and success probability, for a fixed desired success probability. Called from seesaw.

Input (assumed to be valid):
- *rhoAB* state to be distilled
- *nA* dimension of A
- *nB* dimension of B
- *k* local dimension of desired max. entangled state
- *CA* choi state of Alice's map
- *P* measurement indicating success
- *fixedSucc* desired probability of success

Output:
- *newCB* new map for Bob as choi state
- *F* fidelity of new scheme
- *psucc* success probability of new scheme
"""

function seesawBob(
	rho::AbstractMatrix,
	nA::Int, nB::Int, k::Int,
	CA::AbstractMatrix,
	P::AbstractMatrix,
	fixedSucc::Number)

	# Dimension of the maps
	d = nB * 2 * k;

	# Desired output state
	want = maxEnt(k);

	# Check whether we want a real or complex optimization problem
	if (real(rho) == rho) && (real(CA) == CA)
		CB = Semidefinite(d);
	else
		CB = ComplexVariable(d,d);
	end

	# Dimensions in ordering hatA, Fa, A,...
	dimSys = [k, k, 2, 2, nA, nB];

	# Set up the optimization problem
	# The current probability of success
	Mwant = permutesystems(kron(kron(want,P),rho'), [1,3,5,2,4,6], dimSys);
	tv = Variable(1);
	problem = maximize(tv);
	problem.constraints += tv == nA * nB * trace(Mwant * kron(CA, CB));
	problem.constraints += CB in :SDP;

	# Constraints
	# CA should be a choi state
	problem.constraints += ([trace(CB) == 1]);

	CBin = partialtrace(CB,[1],[2 * k, nB]);
	problem.constraints += ([CBin == eye(nA)/nA]);

	# We want a fixed success probability
	Mid = permutesystems(kron(kron(eye(k^2),P),rho'), [1,3,5,2,4,6], dimSys);
	problem.constraints += ([nA * nB * trace(Mid * kron(CA, CB)) == fixedSucc]);

	# Solve the SDP
	solve!(problem, SCSSolver(verbose=false));

	# Compute the new success probability
	psuccNew = nA * nB * trace(Mid * kron(CA, CB.value));

	# The current fidelity achieved
	Fnew = nA * nB * trace(Mwant * kron(CA, CB.value))/psuccNew;

	return(CB.value, Fnew, psuccNew);
end
