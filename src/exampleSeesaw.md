# Using the seesaw method

Remember from arXiv:... that the seesaw method attempts to find a better entanglement distillation scheme for a given state, starting with an existing protocol. Below, we will apply this method to improve a simple filtering scheme that takes one 2x2 quantum state, and returns another 2x2 quantum state with higher fidelity to the EPR pair - at least with some probability of success. 

To use the seesaw method, the existing entanglement distillation scheme must be supplied in terms of the relevant choi states for Alice and Bob. The scheme is therby formulated as a map that takes Alice's system A to an output hat A and a flag register F_A, and similarly for Bob. A projector P must be supplied which is applied to the flag register that Alice and Bob will use to decide success. For example, P=|11><11| would correspond to a protocol in which Alice and Bob declare success if each of them locally declared success corresponding to a flag of 1. 

```julia

# Epsilon parameter for the filtering protocol (see arXiv paper)
eps = 0.2;

# Let's generate a test state to filter
p = 0.7
rho = p * maxEnt(2) + (1-p) * eVec(4,1)*eVec(4,1)';

# Output initial fidelity
Finit = entFidelity(rho);
print("Initial fidelity is F=", round(Finit,3), "\n");

# Test the seesaw optimization for the filtering scheme
# Adapted to work well for states of the form p EPR + (1-p) |01><01|
print("Testing the seesaw method\n");
print("Epsilon is set to ", eps, "\n");

# To use the seesaw method we will need choi states of the existing distillation protocol
print("Computing Choi states\n");
(CA,CB) = filteringMakeChoi(eps);
print("Compute the filtering performance without Choi state to double check:\n");
(rhoQC, P, Fwo, pwo) = filtering(rho, eps);
print("Checked F=",round(Fwo,3), " and psucc=",round(pwo,3),"\n");

# We now run the seesaw optimization. 
# If no success probability is given, then we try to improve the 
# fidelity without hurting the success probability.
print("Run the SEESAW Optimization\n");
(newCA, newCB, Fnew, pnew) = seesaw(rho,2,2,2,CA, CB, P);
```
