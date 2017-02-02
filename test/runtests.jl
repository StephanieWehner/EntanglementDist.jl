using EntanglementDist
using Base.Test

# Test PPT Relax - no improvements should be possible when demanding success probability 1
rho = wernerState(0.9);
initF = entFidelity(rho);
(problem,F,psucc) = pptRelax(rho,2,2,2,1);
@test round(F,3) == round(initF,3)

# Similarly, no improvement for 2 copies
rho = maxEnt(2);
initF = entFidelity(rho);
(problem,F,psucc) = pptRelaxCopies(rho,2,2,2,2,1);
@test round(F,3) == round(initF,3)

# Same for 1 extension
rho = rState(0.9);
initF = entFidelity(rho);
(problem,F,psucc) = pptRelax1Ext(rho,2,2,2,1);
@test round(F,3) == round(initF,3)

# or even 2 extension
(problem,F,psucc) = pptRelax2Ext(rho,2,2,2,1);
@test round(F,3) == round(initF,3)

