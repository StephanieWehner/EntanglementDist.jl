# EntanglementDist

[![Build Status](https://travis-ci.org/StephanieWehner/EntanglementDist.jl.svg?branch=master)](https://travis-ci.org/StephanieWehner/EntanglementDist.jl)

[![Coverage Status](https://coveralls.io/repos/StephanieWehner/EntanglementDist.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/StephanieWehner/EntanglementDist.jl?branch=master)

[![codecov.io](http://codecov.io/github/StephanieWehner/EntanglementDist.jl/coverage.svg?branch=master)](http://codecov.io/github/StephanieWehner/EntanglementDist.jl?branch=master)

This package implements the schemes for optimizing entanglement distillation from arXiv:...

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

## Documentation

For further documentation see build/index.html
