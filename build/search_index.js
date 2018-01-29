var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Optimizing Entanglement Distillation",
    "title": "Optimizing Entanglement Distillation",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#Optimizing-Entanglement-Distillation-1",
    "page": "Optimizing Entanglement Distillation",
    "title": "Optimizing Entanglement Distillation",
    "category": "section",
    "text": "This Julia package is a companion to the arXiv paper ... It implements the numerical methods for computing bounds on the fidelity of the best possible entanglement distillation schemes using realistic operations (RO), when we desire to achieve a particular probability of success. It can also be used to compare the performance of various existing schemes to such bounds, and search for new schemes improving existing methods for particular quantum states.CurrentModule = EntanglementDist"
},

{
    "location": "index.html#Contents-1",
    "page": "Optimizing Entanglement Distillation",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"intro.md\", \"quantum.md\", \"states.md\", \"schemes.md\", \"optimize.md\"]\nDepth = 1"
},

{
    "location": "index.html#Examples-1",
    "page": "Optimizing Entanglement Distillation",
    "title": "Examples",
    "category": "section",
    "text": "Pages = [\"examplePPT.md\", \"exampleKEXT.md\", \"exampleSeesaw.md\"]\nDepth = 1"
},

{
    "location": "index.html#Index-1",
    "page": "Optimizing Entanglement Distillation",
    "title": "Index",
    "category": "section",
    "text": "Pages= [\"quantum.md\", \"states.md\", \"schemes.md\", \"optimize.md\"]"
},

{
    "location": "intro.html#",
    "page": "Introduction",
    "title": "Introduction",
    "category": "page",
    "text": ""
},

{
    "location": "intro.html#Introduction-1",
    "page": "Introduction",
    "title": "Introduction",
    "category": "section",
    "text": ""
},

{
    "location": "intro.html#Installation-1",
    "page": "Introduction",
    "title": "Installation",
    "category": "section",
    "text": "First, you need the latest version of Convex which is not yet available  via Pkg.add. The new version of Convex.jl can handle complex matrices, which we will need. If you have a previously installed version remove it (Pkg.rm) and then install the most recent one by doing:Pkg.clone(\"https://github.com/JuliaOpt/Convex.jl.git\")If you have loaded Convex, a restart of Julia is required before you can use the new package.Second, you want to install the code of this package itself. Pkg.add(\"EntanglementDist\")In the usual way, you can now start using the package:using EntanglementDist;\nimportall EntanglementDist;You are now ready to go!  We recommend trying out the example IJulia sheet found in the example  directory."
},

{
    "location": "intro.html#Basic-functions-1",
    "page": "Introduction",
    "title": "Basic functions",
    "category": "section",
    "text": "Let's see how we can compute a bound on a basic filtering operation. First let us produce a state to be filtered:	rho = isotropicState(0.9);\n	F = pptRelax(rho, 0.8);	\nendWe then run the PPT relaxation to find an upper bound on the fidelity achievable using realistic operations, when trying to obtain an EPR pair with probability of success 0.9.\n# Dimensions of the input state (i.e. what is Alice and Bob in rho)\nnA = 2;\nnB = 2;\n\n# Dimension of the maximally entangled state to produce\nk = 2;\n\n# Run the PPT relaxation, asking for success probability 0.9\n(problem, F, psucc) = pptRelax(rho,2,2,2,0.9)"
},

{
    "location": "optimize.html#",
    "page": "Tools for optimizing distillation",
    "title": "Tools for optimizing distillation",
    "category": "page",
    "text": ""
},

{
    "location": "optimize.html#EntanglementDist.pptRelax",
    "page": "Tools for optimizing distillation",
    "title": "EntanglementDist.pptRelax",
    "category": "Function",
    "text": "(problem, F, p_succ) = pptRelax(rho,n,k,delta, verbose, eps, max_iters)\n\nImplements the PPT relaxation computing an upper bound on the fidelity achievable using realistic operations for a given input state and fixed succes probability.\n\nInputs:\n\nrho quantum state to be distilled\nn number of qubits on one side, assuming dimensions nA=nB =2^n in the input.\nk desired output dimension of the maximally entangled state\ndelta desired success probability\nverbose flag indicating whether the SCS solver should be verbose (default false)\neps desired precision (default 1e-4) in SCS solver\nmax_iters maximum number of iterations in SCS solver (default 2e4)\n\nOutputs:\n\nproblem SDP object from Convex.jl allowing further processing.\nF Fidelity achieved\np_succ success probability attained in SDP (should equal delta)\n\n\n\n(problem, F, p_succ) = pptRelax(rho,nA, nB, k,delta, verbose, eps, max_iters)\n\nImplements the PPT relaxation computing an upper bound on the fidelity achievable using realistic operations for a given input state and fixed succes probability.\n\nInputs:\n\nrho quantum state to be distilled on A and B\nnA dimension of the A system\nnB dimension of the B system\nk desired output dimension of the maximally entangled state\ndelta desired success probability\nverbose flag indicating whether the SCS solver should be verbose (default true)\neps desired precision (default 1e-4) in SCS solver\nmax_iters maximum number of iterations in SCS solver (default 2e4)\n\nOutputs:\n\nproblem SDP object from Convex.jl allowing further processing.\nF Fidelity achieved\np_succ success probability attained in SDP (should equal delta)\n\n\n\n"
},

{
    "location": "optimize.html#EntanglementDist.pptRelaxCopies",
    "page": "Tools for optimizing distillation",
    "title": "EntanglementDist.pptRelaxCopies",
    "category": "Function",
    "text": "(problem, F, p) = pptRelaxCopies(rho, n, nA, nB, k, delta, verbose, eps, max_iters)\n\nCalls pptRelax for an input state of the form rho&^(tensor n), reordered appropriately.\n\n\n\n"
},

{
    "location": "optimize.html#EntanglementDist.pptRelax1Ext",
    "page": "Tools for optimizing distillation",
    "title": "EntanglementDist.pptRelax1Ext",
    "category": "Function",
    "text": "(problem, F, psucc) = pptRelax1Ext(rhoAB, n, k, delta, verbose=true/false, eps=1e^-4)\n\nImplements the PPT relaxation plus 1 symmetric extension for distillable entanglement with fixed desired success probability.\n\nInputs:\n\nrhoAB quantum state to be distilled\nn number of qubits on one side, assuming dimensions nA=nB =2^n in the input.\nk local dimension of the maximally entangled output state\ndelta desired success probability\n\nOutputs:\n\nproblem problem object given by Convex,\nF fidelity bound\npsucc success probability actually attained (should equal delta up to numerical imprecisions)\n\n\n\n(problem, F, psucc) = pptRelax1Ext(rhoAB, nA, nB, k, delta, verbose=true/false, eps=1e^-4)\n\nImplements the PPT relaxation plus 1 symmetric extension for distillable entanglement with fixed desired success probability.\n\nInputs:\n\nrhoAB quantum state to be distilled\nnA  dimension of the A system\nnB  dimension of the B system\nk   local dimension of the maximally entangled output state\ndelta desired success probability\n\nOutputs:\n\nproblem problem object given by Convex,\nF fidelity bound\npsucc success probability actually attained (should equal delta up to numerical imprecisions)\n\n\n\n"
},

{
    "location": "optimize.html#EntanglementDist.pptRelax1ExtCopies",
    "page": "Tools for optimizing distillation",
    "title": "EntanglementDist.pptRelax1ExtCopies",
    "category": "Function",
    "text": "(problem, F, p) = pptRelax1ExtCopies(rho, n, nA, nB, k, delta, verbose, eps, max_iters)\n\nCalls pptRelax1Ext for an input state of the form rho&^(tensor n), reordered appropriately.\n\n\n\n"
},

{
    "location": "optimize.html#EntanglementDist.pptRelax2Ext",
    "page": "Tools for optimizing distillation",
    "title": "EntanglementDist.pptRelax2Ext",
    "category": "Function",
    "text": "(problem, F, psucc) = pptRelax2Ext(rhoAB, n, k, delta, verbose=true/false, eps=1e^-4)\n\nImplements the PPT relaxation plus 2 extensions for distillable entanglement with fixed desired success probability.\n\nInputs:\n\nrhoAB quantum state to be distilled\nn number of qubits on one side, assuming dimensions nA=nB =2^n in the input.\nk local dimension of the maximally entangled output state\ndelta desired success probability\n\nOutputs:\n\nproblem problem object given by Convex,\nF fidelity bound\npsucc success probability actually attained (should equal delta up to numerical imprecisions)\n\n\n\n(problem, F, psucc) = pptRelax2Ext(rhoAB, nA, nB, k, delta, verbose=true/false, eps=1e^-4)\n\nImplements the PPT relaxation plus 2 extensions for distillable entanglement with fixed desired success probability.\n\nInputs:\n\nrhoAB quantum state to be distilled\nnA  dimension of the A system\nnB  dimension of the B system\nk   local dimension of the maximally entangled output state\ndelta desired success probability\n\nOutputs:\n\nproblem problem object given by Convex,\nF fidelity bound\npsucc success probability actually attained (should equal delta up to numerical imprecisions)\n\n\n\n"
},

{
    "location": "optimize.html#EntanglementDist.pptRelax2ExtCopies",
    "page": "Tools for optimizing distillation",
    "title": "EntanglementDist.pptRelax2ExtCopies",
    "category": "Function",
    "text": "(problem, F, p) = pptRelax2ExtCopies(rho, n, nA, nB, k, delta, verbose, eps, max_iters)\n\nCalls pptRelax2Ext for an input state of the form rho&^(tensor n), reordered appropriately.\n\n\n\n"
},

{
    "location": "optimize.html#EntanglementDist.seesaw",
    "page": "Tools for optimizing distillation",
    "title": "EntanglementDist.seesaw",
    "category": "Function",
    "text": "(CA,CB, F, psucc) = seesaw(rhoAB, nA, nB, k, startCA, startCB) or (CA,CB, F, psucc) = seesaw(rhoAB, nA, nB, k, startCA, startCB, P) or (CA,CB, F, psucc) = seesaw(rhoAB, nA, nB, k, startCA, startCB, P, max_iter) or (CA,CB, F, psucc) = seesaw(rhoAB, nA, nB, k, startCA, startCB, P, max_iter, delta)\n\nPerforms a seesaw optimization searching for a good scheme to distill rhoAB into a maximally entangled state of local dimension k. Good here means that we want to maximize the fidelity of the output state to the maximally entangled state, where the success probability should not fall below the existing one - or, if supplied, not below the minimum success probability delta.\n\nAlice's map is assumed to be from A -> hatA, F_A, and Bob's from B -> hatB, F_B, where F_A and F_B are the flag registered that yield success or failure. After the map P is measured on the flag registeres to determine success where P is the projector determining success. Example P=|11><11| on the flag registers F_A and F_B.\n\nInput:\n\nrhoAB state to be distilled\nnA dimension of A\nnB dimension of B\nk  local dimension of the desired maximally entangled output state\nstartCA choi state of Alice's map to start from (ordering of systems: hatA (output state), F_A (flag register), A (input))\nstartCB choi state of Bob's map to start from (ordering of systems: hatB (output state), F_B (flag register), B (input))\nP (optional) projector onto success: Default is |11><11|\nmax_iter (optional) maximum number of iterations (default 10)\ndelta (optional) minimum acceptable success probability.\n\nOutput\n\nCA new map for Alice as choi state\nCB new map for Bob as choi state\nF fidelity obtained\npsucc success probability attained\n\n\n\n"
},

{
    "location": "optimize.html#Tools-for-optimizing-distillation-1",
    "page": "Tools for optimizing distillation",
    "title": "Tools for optimizing distillation",
    "category": "section",
    "text": "The following functions implement the various methods to optimize entanglement distillation described in arXiv:XXX. pptRelax\npptRelaxCopies\npptRelax1Ext\npptRelax1ExtCopies\npptRelax2Ext\npptRelax2ExtCopies\nseesaw"
},

{
    "location": "quantum.html#",
    "page": "Tools for creating and manipulating quantum states",
    "title": "Tools for creating and manipulating quantum states",
    "category": "page",
    "text": ""
},

{
    "location": "quantum.html#EntanglementDist.isQuantumState",
    "page": "Tools for creating and manipulating quantum states",
    "title": "EntanglementDist.isQuantumState",
    "category": "Function",
    "text": "b = isQuantumState(rho) or b = isQuantumState(rho,prec)\n\nChecks whether the input matrix rho is a valid quantum state, that is, a matrix that is positive semidefinite and has trace 1. prec is the precision that small eigenvalues are dealt with: If the input matrix has a negative eigenvalue that is no smaller than -prec, then this eigenvalue is considered zero.\n\nReturns true/false.\n\n\n\n"
},

{
    "location": "quantum.html#EntanglementDist.isHermitian",
    "page": "Tools for creating and manipulating quantum states",
    "title": "EntanglementDist.isHermitian",
    "category": "Function",
    "text": "b = isHermitian(rho) or b = is Hermitian(rho, prec)\n\nChecks whether the input matrix rho is a Hermitian matrix. Entries are compares up to precision prec on average.\n\nReturns true/false.\n\n\n\n"
},

{
    "location": "quantum.html#EntanglementDist.isUnitary",
    "page": "Tools for creating and manipulating quantum states",
    "title": "EntanglementDist.isUnitary",
    "category": "Function",
    "text": "b = isUnitary(U)\n\nChecks whether the input matrix is unitary.\n\nReturns true/false.\n\n\n\n"
},

{
    "location": "quantum.html#EntanglementDist.eVec",
    "page": "Tools for creating and manipulating quantum states",
    "title": "EntanglementDist.eVec",
    "category": "Function",
    "text": "vec = eVec(d, j)\n\nReturns a sparse vector with all 0's except at the position j.\n\n\n\n"
},

{
    "location": "quantum.html#EntanglementDist.maxEnt",
    "page": "Tools for creating and manipulating quantum states",
    "title": "EntanglementDist.maxEnt",
    "category": "Function",
    "text": "vec = maxEnt(d)\n\nReturns a projector (density matrix) onto the maximally entangled state of a given dimension d.\n\n\n\n"
},

{
    "location": "quantum.html#EntanglementDist.maxEntVec",
    "page": "Tools for creating and manipulating quantum states",
    "title": "EntanglementDist.maxEntVec",
    "category": "Function",
    "text": "vec = maxEntVec(d)\n\nReturns the maximally entangled state of a given dimension d as a vector.\n\n\n\n"
},

{
    "location": "quantum.html#EntanglementDist.copies",
    "page": "Tools for creating and manipulating quantum states",
    "title": "EntanglementDist.copies",
    "category": "Function",
    "text": "rhoOut = copies(rho,n)\n\nReturns rho^(tensor n), i.e., n copies of the input rho\n\n\n\n"
},

{
    "location": "quantum.html#EntanglementDist.sortAB",
    "page": "Tools for creating and manipulating quantum states",
    "title": "EntanglementDist.sortAB",
    "category": "Function",
    "text": "rhoSorted = sortAB(rhoBig, d, n)\n\nGiven a state of the form rhoBig = rhoAB^(tensor n), rearrange, so all A parts are first, followed by all B parts. d is the dimension of A, assumed to be the same as for B. n is the number of copies.\n\n\n\n"
},

{
    "location": "quantum.html#EntanglementDist.entFidelity",
    "page": "Tools for creating and manipulating quantum states",
    "title": "EntanglementDist.entFidelity",
    "category": "Function",
    "text": "F = entFidelity(rho)\n\nReturns the fidelity with a maximally entangled state of equal dimension.\n\n\n\n"
},

{
    "location": "quantum.html#EntanglementDist.eprFidelity",
    "page": "Tools for creating and manipulating quantum states",
    "title": "EntanglementDist.eprFidelity",
    "category": "Function",
    "text": "F = eprFidelity(rho)\n\nFor locally two-dimensional states, returns the fidelity with the closest bell state.\n\n\n\n"
},

{
    "location": "quantum.html#EntanglementDist.permutesystems",
    "page": "Tools for creating and manipulating quantum states",
    "title": "EntanglementDist.permutesystems",
    "category": "Function",
    "text": "rhoOut = permutesystems(rho, perm) or rhoOut = permutesystems(rho, perm, dim) or rhoOut = permutesystems(rho,perm, dim, row_only) or rhoOut = permutesystems(rho,perm, dim, row_only, inv_perm)\n\nPermutes the systems in the matrix rho.\n\nInputs:\n\nrho input matrix\nperm a permutation vector on all the systems\ndim (optional keyword argument) vector with dimensions of all the systems. By default all systems are of equal dimension.\nrow_only (optional keyword argument) if set to true, only the rows of rho are permuted.\ninv_perm (optional keyword argument) if set to true, the subsystems are permuted according to the inverse of perm rather than perm itself.\n\nOutputs:\n\nrhoOut matrix rho with permuted systems\n\n\n\n"
},

{
    "location": "quantum.html#EntanglementDist.partialtrace",
    "page": "Tools for creating and manipulating quantum states",
    "title": "EntanglementDist.partialtrace",
    "category": "Function",
    "text": "rhoOut = partialtrace(rhoIn, sysNum, dimVec)\n\nComputes the partial trace of rhoIn, tracing out system number sysNum, where the dimensions of the systems are specified in the vector dimVec. \n\nExample: rhoOut = partialtrace(epr, 2, [2 2])\n\nTraces out the second qubit of a 2x2 state, here the EPR pair.\n\n\n\n"
},

{
    "location": "quantum.html#EntanglementDist.partialtranspose",
    "page": "Tools for creating and manipulating quantum states",
    "title": "EntanglementDist.partialtranspose",
    "category": "Function",
    "text": "rhoOut = partialtranspose(rho) or rhoOut = partialtranspose(rho, sys) or rhoOut = partialtranspose(rho, sys, dim)\n\nApplies the partial tanspose to the matrix rho.\n\nInputs:\n\nrho input matrix\nsystems (optional) single system or a vector of systems to be transposed. By default it is equal to 2. If dim is not provided, the only two possible values of systems are 1 or 2.\ndim (optional) vector with dimensions of the systems. By default it assumes two systems of equal dimensions.\n\nOutputs:\n\nrhoOut matrix rho with partial transpose applied\n\n\n\n"
},

{
    "location": "quantum.html#Tools-for-creating-and-manipulating-quantum-states-1",
    "page": "Tools for creating and manipulating quantum states",
    "title": "Tools for creating and manipulating quantum states",
    "category": "section",
    "text": "Within our package we provide a number of useful functions for dealing with quantum states. isQuantumState\nisHermitian\nisUnitary\neVec\nmaxEnt\nmaxEntVec\ncopies\nsortAB\nentFidelity\neprFidelity\npermutesystems\npartialtrace\npartialtranspose"
},

{
    "location": "schemes.html#",
    "page": "Known distillation/filtering schemes",
    "title": "Known distillation/filtering schemes",
    "category": "page",
    "text": ""
},

{
    "location": "schemes.html#EntanglementDist.DEJMPSParam",
    "page": "Known distillation/filtering schemes",
    "title": "EntanglementDist.DEJMPSParam",
    "category": "Function",
    "text": "(F, p_succ) = DEJMPSParam(rho)\n\nDetermines the performance of DEJMPS for two copies of the input state rho (consisting of two qubits) and returns the fidelity and success probability achieved. The procedure optimises over local single qubit rotations that permute the bell diagonal coefficients before the CNOTs to achieve the highest output fidelity. returns (F_out, p_succ)\n\n\n\n"
},

{
    "location": "schemes.html#EntanglementDist.BBPSSWParam",
    "page": "Known distillation/filtering schemes",
    "title": "EntanglementDist.BBPSSWParam",
    "category": "Function",
    "text": "(F, p_succ) = BBPSSWParam(F)\n\nDetermines the performance of BBPSSW for two copies of the input state (consisting of two qubits) with fidelity F each. returns (F_out, p_succ)\n\n\n\n(F, p_succ) = BBPSSWParam(rho)\n\nDetermines the performance of BBPSSW for two copies of the input state rho (consisting of two qubits). returns (F_out, p_succ)\n\n\n\n"
},

{
    "location": "schemes.html#EntanglementDist.EPLDParam",
    "page": "Known distillation/filtering schemes",
    "title": "EntanglementDist.EPLDParam",
    "category": "Function",
    "text": "(F, p_succ) = EPLDParam(p, pd)\n\nDetermines the performance of EPL-D for rStateCorrPhase. The inputs are the p and pd parameters of rStateCorrPhase. returns (F_out, p_succ)\n\n\n\n(F, p_succ) = EPLDParam(rho)\n\nDetermines the performance of EPL-D for an arbitrary bipartite state of dimensions nA=nB=4. returns (F_out, p_succ)\n\n\n\n"
},

{
    "location": "schemes.html#EntanglementDist.filtering",
    "page": "Known distillation/filtering schemes",
    "title": "EntanglementDist.filtering",
    "category": "Function",
    "text": "`(rhoQC,P,F,psucc) = filtering(rhoAB,eps)\n\nImplements the filtering protocol on input state rho with filtering parameter eps. This filtering protocol is designed to work well on states of the form p * EPR + (1-p) |10><10|.\n\nInput:\n\nrhoAB 2x2 state to be filtered\neps filtering parameter, smaller eps typically means higher fidelity but lower probability of success\n\nOutput:\n\nrhoQC filtered state with flags: Q=hatA, hatB contains the output state, C  F_A,F_B contains the flags (1 success, 0 failure)\nP projector specifying how Alice and Bob determine success |11><11|, both succeed\nF fidelity achieved for the given input state\npsucc probability of success\n\n\n\n"
},

{
    "location": "schemes.html#EntanglementDist.filteringMakeChoi",
    "page": "Known distillation/filtering schemes",
    "title": "EntanglementDist.filteringMakeChoi",
    "category": "Function",
    "text": "`(choiA, choiB) = filteringMakeChoi(eps)\n\nComputes the choi states for Alice and Bob in the filtering protocol with filtering parameter eps. This filtering protocol is designed to work well on states of the form p * EPR + (1-p) |10><10|.\n\nInput:\n\neps filtering parameter, smaller eps typically means higher fidelity but lower probability of success\n\nOutput:\n\nchoi state for Alice hatA,FA,A\nchoi state for Bob hatB,FB,B\n\n\n\n"
},

{
    "location": "schemes.html#EntanglementDist.measureScheme",
    "page": "Known distillation/filtering schemes",
    "title": "EntanglementDist.measureScheme",
    "category": "Function",
    "text": "(rhoQC, F, psucc) = measureScheme(rhoAB, k, MokAlice, MfailAlice, MokBob, MfailBob, fokAlice, ffailAlice, fokBob, ffailBob, P)\n\nImplements a distillation scheme in which Alice and Bob perform the indicated measurement and then output a success or failure flag individually. Global success is determined by the projector P which may compare the to local flag registers.\n\nInput:\n\nrhoAB input state to be distilled\nk local dimension of the max entangled state to be produced\nMokAlice Kraus operator corresponding to Alice success\nMfailAlice Kraus operator corresponding to Alice failure\nMokBob Kraus operator corresponding to Bob success\nMfailBob Kraus operator corresponding to Bob failure\nfokAlice flag Alice outputs for success\nffailAlice flag Alice outputs for failure\nfokBob flag Bob outputs for success\nffailBob flag Bob outputs for failure\nP projector determining for which joint flags they declare global success\n\nOutput:\n\nrhoQC output state where Q=hatA, hatB is the output state and C=F_A,F_B the flag register\nF fidelity achieved\npsucc success probability achieved\n\n\n\n"
},

{
    "location": "schemes.html#EntanglementDist.measureSchemeMakeChoi",
    "page": "Known distillation/filtering schemes",
    "title": "EntanglementDist.measureSchemeMakeChoi",
    "category": "Function",
    "text": "choiState = measureSchemeMakeChoi(MokAlice, MfailAlice, fokAlice, ffailAlice)\n\nComputes the choi state of a measuring distillation/filtering protocol of one party: Outputs choi_QA where A is a copy of the input and Q is the output of the distillation map A->hatA, fA\n\nInput:\n\nMokAlice Kraus operator corresponding to Alice success\nMfailAlice Kraus operator corresponding to Alice failure\nfokAlice flag Alice outputs for success\nffailAlice flag Alice outputs for failure\n\n\n\n"
},

{
    "location": "schemes.html#Known-distillation/filtering-schemes-1",
    "page": "Known distillation/filtering schemes",
    "title": "Known distillation/filtering schemes",
    "category": "section",
    "text": "Compute how well several known entanglement distillation schemes perform on specific input states. We also implement a simple filtering scheme and provide function for computing the performance of general filtering schemes given the relevant measurement operators. The functions computing the choi states of said filtering protocols can be used to test the performance of the seesaw method.DEJMPSParam\nBBPSSWParam\nEPLDParam\nfiltering\nfilteringMakeChoi\nmeasureScheme\nmeasureSchemeMakeChoi"
},

{
    "location": "states.html#",
    "page": "Special quantum states",
    "title": "Special quantum states",
    "category": "page",
    "text": ""
},

{
    "location": "states.html#EntanglementDist.isotropicState",
    "page": "Special quantum states",
    "title": "EntanglementDist.isotropicState",
    "category": "Function",
    "text": "rho = isotropicState(p) or rho = isotropicState(p,d)\n\nReturns an isotropic state, i.e., a mixture of a maximally entangled pair (with probability p) and the maximally mixed state in local dimension d. If no argument d is given the default is d=2, that is, we take the mixture of the EPR pair with the maximally mixed state of local dimension 2.\n\n\n\n"
},

{
    "location": "states.html#EntanglementDist.bellDiagState",
    "page": "Special quantum states",
    "title": "EntanglementDist.bellDiagState",
    "category": "Function",
    "text": "rho = bellDiagState(p1,p2,p3)\n\nReturns a bell diagonal state that is a mixture between the 4 Bell states: p1 * phi^+ + p2 * psi^+ + p3 * phi^- + (1 - p1 - p2 - p3) * psi^-.\n\n\n\n"
},

{
    "location": "states.html#EntanglementDist.sState",
    "page": "Special quantum states",
    "title": "EntanglementDist.sState",
    "category": "Function",
    "text": "rho = sState(p)\n\nReturns a state that is an a mixture between the EPR pair (with probability p), and the state |11><11|.\n\n\n\n"
},

{
    "location": "states.html#EntanglementDist.rState",
    "page": "Special quantum states",
    "title": "EntanglementDist.rState",
    "category": "Function",
    "text": "rho = rState(p)\n\nReturns a mixture between a state proportional to |01> + |10> (with probability p) and the state |11><11|.\n\n\n\n"
},

{
    "location": "states.html#EntanglementDist.sStateQutrit",
    "page": "Special quantum states",
    "title": "EntanglementDist.sStateQutrit",
    "category": "Function",
    "text": "rho = sStateQutrit(p)\n\nReturns a state that is a mixture between a 3 dimensional maximally entangled state (with probability p) and the state |00><00|\n\n\n\n"
},

{
    "location": "states.html#EntanglementDist.rStatePhase",
    "page": "Special quantum states",
    "title": "EntanglementDist.rStatePhase",
    "category": "Function",
    "text": "rho = rStatePhase(p, phi) or rho = rStatePhase(p)\n\nReturns a mixture between a state proportional to |01> + e^(i phi) |10> (with probability p) and the state |11><11|. No value for phi defaults to phi=0.\n\n\n\n"
},

{
    "location": "states.html#EntanglementDist.rStateCorrPhase",
    "page": "Special quantum states",
    "title": "EntanglementDist.rStateCorrPhase",
    "category": "Function",
    "text": "rho = rStateCorrPhase(p, pd)\n\nReturns a state of the form integral phi [pd * r(p,phi) + (1-pd) * r(p,phi+pi)] tensor r(p,phi), where r(p,phi) = rStatePhase(p,phi). Default for pd is 1.\n\n\n\n"
},

{
    "location": "states.html#EntanglementDist.rStateCorrPhaseCopies",
    "page": "Special quantum states",
    "title": "EntanglementDist.rStateCorrPhaseCopies",
    "category": "Function",
    "text": "rho = rStateCorrPhaseCopies(p)\n\nReturns a state of the form integral phi r(p,phi)^(tensor n) , where r(p,phi) = rStatePhase(p,phi) and n is the number of desired copies.\n\n\n\n"
},

{
    "location": "states.html#Special-quantum-states-1",
    "page": "Special quantum states",
    "title": "Special quantum states",
    "category": "section",
    "text": "The following functions generate specific sets of quantum states which commonly occur either in physical implemtations or in mathematical studies of entanglement distillation.isotropicState\nbellDiagState\nsState\nrState\nsStateQutrit\nrStatePhase\nrStateCorrPhase\nrStateCorrPhaseCopies"
},

{
    "location": "examplePPT.html#",
    "page": "Using the PPT Relaxation method",
    "title": "Using the PPT Relaxation method",
    "category": "page",
    "text": ""
},

{
    "location": "examplePPT.html#Using-the-PPT-Relaxation-method-1",
    "page": "Using the PPT Relaxation method",
    "title": "Using the PPT Relaxation method",
    "category": "section",
    "text": "Let us now see how we can use the method of PPT relaxations to compute an upper bound on the maximum achievable fidelity of any entanglement distillation scheme using realistic operations.Let us first define a state to be distilled. In this example, let's consider the isotropic state. This state will be a mixture between the EPR pair (with probability 0.9), and the maximally mixed state (with probability 0.1).\nrho = isotropicState(0.9);\nLet's now go and compute the maximum fidelity of distilling a few of these Werner states.\n# Local dimension of what we'll call Alice in the input state rho\nnA = 2;\n\n# Similarly for Bob\nnB = 2;\n\n# Number of copies to distill\nn = 3;\n\n# Dimension of the maximally entangled state with want to produce\nk = 2;\n\n# Desired probability of success\ndelta = 0.8;\n\n# Compute the maximum fidelity F\n(problem, F, psucc) = pptRelaxCopies(rho, n, nA, nB, k, delta); "
},

{
    "location": "exampleKEXT.html#",
    "page": "Using the method of k extensions",
    "title": "Using the method of k extensions",
    "category": "page",
    "text": ""
},

{
    "location": "exampleKEXT.html#Using-the-method-of-k-extensions-1",
    "page": "Using the method of k extensions",
    "title": "Using the method of k extensions",
    "category": "section",
    "text": "This example is almost identical to the one of PPT relaxations, except that we will use instead the method of extensions to try and get a sharper bound. In the example below, we will allow one extra symmetric extension.Let us first define a state to be distilled. In this example, let's consider the isotropic state. This state will be a mixture between the EPR pair (with probability 0.9), and the maximally mixed state (with probability 0.1).\nrho = isotropicState(0.9);\nLet's now go and compute the maximum fidelity of distilling a few of these Werner states.\n# Local dimension of what we'll call Alice in the input state rho\nnA = 2;\n\n# Similarly for Bob\nnB = 2;\n\n# Number of copies to distill\nn = 3;\n\n# Dimension of the maximally entangled state with want to produce\nk = 2;\n\n# Desired probability of success\ndelta = 0.8;\n\n# Compute the maximum fidelity F\n(problem, F, psucc) = pptRelax1ExtCopies(rho, n, nA, nB, k, delta); "
},

{
    "location": "exampleSeesaw.html#",
    "page": "Using the seesaw method",
    "title": "Using the seesaw method",
    "category": "page",
    "text": ""
},

{
    "location": "exampleSeesaw.html#Using-the-seesaw-method-1",
    "page": "Using the seesaw method",
    "title": "Using the seesaw method",
    "category": "section",
    "text": "Remember from arXiv:... that the seesaw method attempts to find a better entanglement distillation scheme for a given state, starting with an existing protocol. Below, we will apply this method to improve a simple filtering scheme that takes one 2x2 quantum state, and returns another 2x2 quantum state with higher fidelity to the EPR pair - at least with some probability of success. To use the seesaw method, the existing entanglement distillation scheme must be supplied in terms of the relevant choi states for Alice and Bob. The scheme is therby formulated as a map that takes Alice's system A to an output hat A and a flag register F_A, and similarly for Bob. A projector P must be supplied which is applied to the flag register that Alice and Bob will use to decide success. For example, P=|11><11| would correspond to a protocol in which Alice and Bob declare success if each of them locally declared success corresponding to a flag of 1. \n# Epsilon parameter for the filtering protocol (see arXiv paper)\neps = 0.2;\n\n# Let's generate a test state to filter\np = 0.7\nrho = sState(p);\n\n# Output initial fidelity\nFinit = entFidelity(rho);\nprint(\"Initial fidelity is F=\", round(Finit,3), \"\\n\");\n\n# Test the seesaw optimization for the filtering scheme\n# Adapted to work well for states of the form p EPR + (1-p) |01><01|\nprint(\"Testing the seesaw method\\n\");\nprint(\"Epsilon is set to \", eps, \"\\n\");\n\n# To use the seesaw method we will need choi states of the existing distillation protocol\nprint(\"Computing Choi states\\n\");\n(CA,CB) = filteringMakeChoi(eps);\nprint(\"Compute the filtering performance without Choi state to double check:\\n\");\n(rhoQC, P, Fwo, pwo) = filtering(rho, eps);\nprint(\"Checked F=\",round(Fwo,3), \" and psucc=\",round(pwo,3),\"\\n\");\n\n# We now run the seesaw optimization. \n# If no success probability is given, then we try to improve the \n# fidelity without hurting the success probability.\nprint(\"Run the SEESAW Optimization\\n\");\n(newCA, newCB, Fnew, pnew) = seesaw(rho,2,2,2,CA, CB, P);"
},

]}
