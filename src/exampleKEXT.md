# Using the method of k extensions

This example is almost identical to the one of PPT relaxations, except that we will use instead the method of symmetric extensions to try and get a sharper bound. In the example below, we will allow one extra extension.

Let us first define a state to be distilled. In this example, let's consider the Werner state. This state will be a mixture between the EPR pair (with probability 0.9), and the maximally mixed state (with probability 0.1).

```julia

rho = wernerState(0.9);

```

Let's now go and compute the maximum fidelity of distillting a few of these Werner states.

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
(problem, F, psucc) = pptRelax1ExtCopies(rho, n, nA, nB, k, delta); 
```
