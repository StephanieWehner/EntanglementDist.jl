using EntanglementDist

# Change the following to perform further tests:
eps = 0.2;

# rho = wernerState(0.9);
# rho = rState(0.7);

# Customized state:
p = 0.7
rho = sState(p);

# Output initial fidelity
Finit = entFidelity(rho);
print("Initial fidelity is F=", round(Finit,3), "\n");

# Test the seesaw optimization for the filtering scheme
# Adapted to work well for states of the form p EPR + (1-p) |01><01|
print("Testing the seesaw method\n");
print("Epsilon is set to ", eps, "\n");

print("Computing Choi states\n");
(CA,CB) = filteringMakeChoi(eps);
print("Compute the filtering performance without Choi state to double check:\n");
(rhoQC, P, Fwo, pwo) = filtering(rho, eps);
print("Checked F=",round(Fwo,3), " and psucc=",round(pwo,3),"\n");

print("Run the SEESAW Optimization\n");
(newCA, newCB, Fnew, pnew) = seesaw(rho,2,2,2,CA, CB, P);
 


