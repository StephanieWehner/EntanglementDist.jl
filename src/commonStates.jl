# Defines some commonly used states in quantum information

export bell, bellDiagState, isotropicState, rState, sState, sStateQutrit, rStatePhase, rStateCorrPhase, rStateCorrPhaseCopies;

""" `bell`

The Bell states (as vectors), bell[j,1:4] gives the j-th bell state.
"""
bell = zeros(4,4);
bell[1,1:4] = maxEntVec(2);
bell[2,1:4] = (kron(eVec(2,1),eVec(2,1)) -  kron(eVec(2,2),eVec(2,2)))/sqrt(2);
bell[3,1:4] = (kron(eVec(2,1),eVec(2,2)) +  kron(eVec(2,2),eVec(2,1)))/sqrt(2);
bell[4,1:4] = (kron(eVec(2,1),eVec(2,2)) -  kron(eVec(2,2),eVec(2,1)))/sqrt(2);

""" `epr` and `singlet`

Two of these have special names: we will use epr for the EPR pair, and singlet for the singlet.
"""
epr = bell[1,1:4];
singlet = bell[4,1:4];

""" `bellS`

An array of density matrices corresponding to the 4 Bell states.
"""
bS1 = epr*epr';
bS2 = bell[2,1:4]* bell[2,1:4]';
bS3 = bell[3,1:4]* bell[3,1:4]';
bS4 = bell[4,1:4]*bell[4,1:4]';
bellS = [bS1; bS2; bS3; bS4];

""" `PS`
Projector onto the symmetric subspace of 2 qubits.
"""
PS = bS1 + bS2 + bS3;

""" `PA`
Projector onto the antisymmetric subspace of 2 qubits.
"""
PA = bS4;

""" `rho = bellDiagState(p1,p2,p3)`

Returns a bell diagonal state that is a mixture between the 4 Bell states: p1 * phi^+ + p2 * psi^+ + p3 * phi^- + (1 - p1 - p2 - p3) * psi^-.
"""

function bellDiagState(p1::Number, p2::Number, p3::Number)

	@assert 0 <= p1 "Probilities must be positive."
	@assert 0 <= p2 "Probilities must be positive."
	@assert 0 <= p3 "Probilities must be positive."

	@assert p1 + p2 + p3 <= 1 "Probabilities cannot exceed 1."

	# Produce the desired mixture
	out = p1 * bS1 + p2 * bS3 + p3 * bS2 + (1-p1-p2-p3) * bS4;
	return out;
end

""" `rho = sState(p)`

Returns a state that is an a mixture between the EPR pair (with probability *p*), and the state |11><11|.
"""

function sState(p::Number)

	@assert 0 <= p "Probilities must be positive."
	@assert p <= 1 "Probabilities cannot exceed 1."

	# Generate the maximally entangled state
	epr = maxEnt(2);

	# Generate |11>
	e0 = [1 0];
	e1 = [0 1];
	v11 = kron(e1,e1);

	# Produce the desired mixture
	out = p * epr + (1-p) * v11'*v11;
	return out;
end

""" `rho = sStateQutrit(p)`

Returns a state that is a mixture between a 3 dimensional maximally entangled state (with probability *p*) and the state |00><00|
"""

function sStateQutrit(p::Number)

	@assert 0 <= p "Probilities must be positive."
	@assert p <= 1 "Probabilities cannot exceed 1."

	# Maximally entangled state
	epr = maxEnt(3);

	# Generate |00>
	e0 = [1 0 0];
	e1 = [0 0 1];
	v00 = kron(e0,e0);

	# Produce the desired mixture
	out = p * epr + (1-p) * v00'*v00;

  	return out;
end

""" `rho = rState(p)`

Returns a mixture between a state proportional to |01> + |10> (with probability *p*) and the state |11><11|.
"""

rState(p::Number) = rStatePhase(p, 0)

""" `rho = rStatePhase(p, phi)` or `rho = rStatePhase(p)`

Returns a mixture between a state proportional to |01> + e^(i *phi*) |10> (with probability *p*) and the state |11><11|. No value for *phi* defaults to *phi*=0.

"""

function rStatePhase(p::Number, phi::Number = 0.0)

	@assert 0 <= p "Probilities must be positive."
	@assert p <= 1 "Probabilities cannot exceed 1."

	# Produce the state |01> + e^(i phi) |10>
	e0 = eVec(2, 1);
	e1 = eVec(2, 2);
	vec = (kron(e0, e1) + e^(im*phi) * kron(e1,e0))/sqrt(2);

	if phi == 0 || phi == pi
		vec = real(vec);
	end

	# Produce a state orthogonal to the one above
	v11 = kron(e1,e1);

	# Construct the desire mixture
	out = p * vec*vec' + (1-p) * v11*v11';

  	return out
end

""" `rho = rStateCorrPhase(p, pd)`

Returns a state of the form integral phi [pd * r(p,phi) + (1-pd) * r(p,phi+pi)] tensor r(p,phi), where
r(p,phi) = rStatePhase(p,phi). Default for pd is 1.
"""
function rStateCorrPhase(p::Number, pd::Number = 1)

	@assert 0 <= p "Probilities must be positive."
	@assert p <= 1 "Probabilities cannot exceed 1."


    rStatePhaseDep(phi) = pd * rStatePhase(p, phi) + (1-pd) * rStatePhase(p, phi + pi)

	integrand(phi) = ( 1/(2*pi) ) * kron(rStatePhaseDep(phi), rStatePhase(p, phi));

	eps = 10.0^(-4)
	out = real( quadgk( integrand, 0, 2 * pi; reltol=sqrt(eps), abstol=0, maxevals=10^7, order=7, norm=vecnorm)[1])
  	return out
end


""" `rho = rStateCorrPhaseCopies(p)`

Returns a state of the form integral phi r(p,phi)^(tensor n) , where
r(p,phi) = rStatePhase(p,phi) and n is the number of desired copies.
"""
function rStateCorrPhaseCopies(p::Number, n::Int)

	@assert 0 <= p "Probilities must be positive."
	@assert p <= 1 "Probabilities cannot exceed 1."
	@assert n > 1 "Number of copies must be at least 1."

	integrand(phi) = ( 1/(2*pi) ) * copies(rStatePhase(p, phi), n);

	eps = 10.0^(-4)
	out = real( quadgk( integrand, 0, 2 * pi; reltol=sqrt(eps), abstol=0, maxevals=10^7, order=7, norm=vecnorm)[1])
  	return out
end

""" `rho = isotropicState(p)` or `rho = isotropicState(p,d)`

Returns an isotropic state, i.e., a mixture of a maximally entangled pair (with probability *p*) and the maximally mixed state in local dimension *d*. If no argument *d* is given the default is *d*=2, that is, we take the mixture of the EPR pair with the maximally mixed state of local dimension 2.
"""

function isotropicState(p::Number; d::Int = 2)

	@assert 0 <= p "Probilities must be positive."
	@assert p <= 1 "Probabilities cannot exceed 1."

        epr = maxEnt(d);
	out = p * epr*epr' + (1-p) * eye(d^2)/d^2;

	return out;
end
