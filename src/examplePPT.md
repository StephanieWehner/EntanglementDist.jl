# Using the PPT Relaxation method

Let us now see how we can use the method of PPT relaxations to compute an upper bound on the maximum achievable fidelity of any entanglement distillation scheme using realistic operations.

Let us first define a state to be distilled. In this example, let's consider the isotropic state. This state will be a mixture between the EPR pair (with probability 0.9), and the maximally mixed state (with probability 0.1).

```julia

rho = isotropicState(0.9);

```

Let's now go and compute the maximum fidelity of distilling a few of these Werner states.

```julia

# Local dimension of what we'll call Alice in the input state rho
nA = 2;

# Similarly for Bob
nB = 2;

# Number of copies to distill
n = 3;

# Dimension of the maximally entangled state with want to produce
k = 2;

# Desired probability of success
delta = 0.8;

# Compute the maximum fidelity F
(problem, F, psucc) = pptRelaxCopies(rho, n, nA, nB, k, delta); 
```
