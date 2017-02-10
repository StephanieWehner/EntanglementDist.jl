# Introduction
## Installation

First, you need the latest version of Convex which is not yet available 
via Pkg.add. The new version of Convex.jl can handle complex matrices, which we will need. If you have a previously installed version remove it (Pkg.rm) and then install the most recent one by doing:

```julia
Pkg.clone("https://github.com/JuliaOpt/Convex.jl.git")
```

If you have loaded Convex, a restart of Julia is required before you can use the new package.

Second, you want to install the code of this package itself. 
```julia
Pkg.add("EntanglementDist")
```

In the usual way, you can now start using the package:
```julia
using EntanglementDist;
importall EntanglementDist;
```
You are now ready to go! 
We recommend trying out the example IJulia sheet found in the example 
directory.

## Basic functions
Let's see how we can compute a bound on a basic filtering operation. First let us produce a state to be filtered:
```julia
	rho = wernerState(0.9);
	F = pptRelax(rho, 0.8);	
end
```
We then run the PPT relaxation to find an upper bound on the fidelity achievable using realistic operations, when trying to obtain an EPR pair with probability of success 0.9.
```julia

# Dimensions of the input state (i.e. what is Alice and Bob in rho)
nA = 2;
nB = 2;

# Dimension of the maximally entangled state to produce
k = 2;

# Run the PPT relaxation, asking for success probability 0.9
(problem, F, psucc) = pptRelax(rho,2,2,2,0.9)
```

