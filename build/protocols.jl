# Implements some well known distillation/filtering protocols
#
# It will create a QC state of the form
# rho_hatA, hatB, F_A, F_B where hatA, and hatB contain the outputs, 
# and FA, FB are the flag registers
#


export filtering, filteringMakeChoi, measureScheme, measureSchemeMakeChoi;

""" 
``(rhoQC,P,F,psucc) = filtering(rhoAB,eps)`

Implements the filtering protocol on input state *rho* with filtering parameter *eps*. This filtering protocol is designed
to work well on states of the form p * EPR + (1-p) |10><10|. 

Input:
- *rhoAB* 2x2 state to be filtered
- *eps* filtering parameter, smaller eps typically means higher fidelity but lower probability of success

Output:
- *rhoQC* filtered state with flags: Q=hatA, hatB contains the output state, C  F_A,F_B contains the flags (1 success, 0 failure)
- *P* projector specifying how Alice and Bob determine success |11><11|, both succeed
- *F* fidelity achieved for the given input state
- *psucc* probability of success
"""

function filtering(rho::AbstractMatrix, eps::Number; makeChoi::Bool=false)

	(d1,d2) = size(rho);
	@assert d1 == 4 "Filtering takes a 2x2 input state."
	@assert eps <= 1 "Filtering parameter eps needs to be less than 1."
	@assert eps >= 0 "Filtering parameter eps must be positive."

	# Measurement operator for Alice and Bob 
	MokAlice = sqrt(eps) * eVec(2,1)*eVec(2,1)' + eVec(2,2)*eVec(2,2)';
	MfailAlice = sqrt(1-eps) * eVec(2,1) * eVec(2,1)';
	MokBob = sqrt(eps) * eVec(2,2)*eVec(2,2)' + eVec(2,1)*eVec(2,1)';
	MfailBob= sqrt(1-eps) * eVec(2,2) * eVec(2,2)';

	# How success and failure are indicated on the local flag registers
	# Here the same for both Alice and Bob
	# success = |1><1|, failure = |0><0|
	fok = eVec(2,2)*eVec(2,2)';
	ffail = eVec(2,1) * eVec(2,1)';

	# Overall success only when both succeed
	P = kron(fok,fok);

	# Compute the state rhoQC
	(rhoQC, F, psucc) = measureScheme(rho, 2, MokAlice, MfailAlice, MokBob, MfailBob, fok, ffail, fok, ffail, P);

	return(rhoQC, P, F, psucc);
end

""" 
``(choiA, choiB) = filteringMakeChoi(eps)`

Computes the choi states for Alice and Bob in the filtering protocol
with filtering parameter *eps*. This filtering protocol is designed
to work well on states of the form p * EPR + (1-p) |10><10|. 

Input:
- *eps* filtering parameter, smaller eps typically means higher fidelity but lower probability of success

Output:
- *choi* state for Alice hatA,FA,A
- *choi* state for Bob hatB,FB,B
"""

function filteringMakeChoi(eps::Number)

	@assert eps <= 1 "Filtering parameter eps needs to be less than 1."
	@assert eps >= 0 "Filtering parameter eps must be positive."

	# Measurement operator for Alice and Bob 
	MokAlice = sqrt(eps) * eVec(2,1)*eVec(2,1)' + eVec(2,2)*eVec(2,2)';
	MfailAlice = sqrt(1-eps) * eVec(2,1) * eVec(2,1)';
	MokBob = sqrt(eps) * eVec(2,2)*eVec(2,2)' + eVec(2,1)*eVec(2,1)';
	MfailBob= sqrt(1-eps) * eVec(2,2) * eVec(2,2)';

	# How success and failure are indicated on the local flag registers
	# Here the same for both Alice and Bob
	# success = |1><1|, failure = |0><0|
	fok = eVec(2,2)*eVec(2,2)';
	ffail = eVec(2,1) * eVec(2,1)';

	choiA = measureSchemeMakeChoi(MokAlice,MfailAlice,fok, ffail);
	choiB = measureSchemeMakeChoi(MokBob,MfailBob,fok, ffail);

	return(choiA, choiB);
end

""" `(rhoQC, F, psucc) = measureScheme(rhoAB, k, MokAlice, MfailAlice, MokBob, MfailBob, fokAlice, ffailAlice, fokBob, ffailBob, P)`

Implements a distillation scheme in which Alice and Bob perform the indicated measurement and then output a success or failure flag individually. Global success is determined by the projector P which may compare the to local flag registers.

Input:
- *rhoAB* input state to be distilled
- *k* local dimension of the max entangled state to be produced
- *MokAlice* Kraus operator corresponding to Alice success
- *MfailAlice* Kraus operator corresponding to Alice failure
- *MokBob* Kraus operator corresponding to Bob success
- *MfailBob* Kraus operator corresponding to Bob failure
- *fokAlice* flag Alice outputs for success
- *ffailAlice* flag Alice outputs for failure
- *fokBob* flag Bob outputs for success
- *ffailBob* flag Bob outputs for failure
- *P* projector determining for which joint flags they declare global success

Output:
- *rhoQC* output state where Q=hatA, hatB is the output state and C=F_A,F_B the flag register
- *F* fidelity achieved
- *psucc* success probability achieved

"""

function measureScheme(rho::AbstractMatrix, k::Number, MokAlice::AbstractMatrix, MfailAlice::AbstractMatrix, MokBob::AbstractMatrix, MfailBob::AbstractMatrix, fokAlice::AbstractMatrix, ffailAlice::AbstractMatrix, fokBob::AbstractMatrix, ffailBob::AbstractMatrix, P::AbstractMatrix)

	@assert isQuantumState(rho) "Input is not a quantum state."
	@assert k >= 2 "Output dimension is at least 2."

	# Desired output state
	want = maxEnt(k);

	# Construct the QC state of the output
	rhoQC = kron(kron(MokAlice,MokBob) * rho * kron(MokAlice,MokBob)', kron(fokAlice, fokBob));
	rhoQC += kron(kron(MokAlice,MfailBob) * rho * kron(MokAlice,MfailBob)', kron(fokAlice, ffailBob));
	rhoQC += kron(kron(MfailAlice,MokBob) * rho * kron(MfailAlice,MokBob)', kron(ffailAlice, fokBob));
	rhoQC += kron(kron(MfailAlice,MfailBob) * rho * kron(MfailAlice,MfailBob)', kron(ffailAlice, ffailBob));

	# Compute prob. of success and the post-measurement state conditioned on success, includes flag register
	rhoSuccUN = kron(eye(k^2),P) * rhoQC * kron(eye(k^2),P)';
	psucc = trace(rhoSuccUN);
	rhoSucc  = rhoSuccUN/psucc;

	# Tracing out the flag register
	(d1,d2) = size(P);
	rhoSuccNoFlag = partialtrace(rhoSucc,[2], [k^2, d1]);

	# Compute the fidelity with max entangled state
	F = entFidelity(rhoSuccNoFlag);

	return(rhoQC, F, psucc);
end


""" `choiState = measureSchemeMakeChoi(MokAlice, MfailAlice, fokAlice, ffailAlice)`

Computes the choi state of a measuring distillation/filtering protocol of one party: Outputs choi_QA where A is a copy of the input and Q is the output of the distillation map A->hatA, fA
 
Input:
- *MokAlice* Kraus operator corresponding to Alice success
- *MfailAlice* Kraus operator corresponding to Alice failure
- *fokAlice* flag Alice outputs for success
- *ffailAlice* flag Alice outputs for failure

"""

function measureSchemeMakeChoi(MokAlice::AbstractMatrix, MfailAlice::AbstractMatrix, fokAlice::AbstractMatrix, ffailAlice::AbstractMatrix)

	@assert isQuantumState(fokAlice) "Success tag is not a valid quantum state."
	@assert isQuantumState(ffailAlice) "Failure tag is not a valid quantum state."

	# Compute input dimensions
	(outA, inA) = size(MokAlice);
	(fA, fA2) = size(fokAlice);

	# We'll produce the choi state starting with the maximally entangled one.
	rho = maxEnt(inA);

	# Construct the QC state of the output
	id = eye(inA);
	choi = kron(kron(MokAlice,id) * rho * kron(MokAlice,id)', fokAlice);
	choi += kron(kron(MfailAlice,id) * rho * kron(MfailAlice,id)', ffailAlice);

	# Reorder to hatA, F_A, A
	choi = permutesystems(choi, [1, 3, 2], dim=[outA, inA, fA]);

	return choi;
end
