# Fix of kron to work with Convex.jl in seesawAlice
# This is a fix of the bug in the kron extension
# of the latest registered version of Convex for Julia v0.5.
# The bug does not allow to apply kron to a variable and fixed matrix.
# To solve that we define a new function "simplekron" which deals with this case.

using Convex

export simplekron

function simplekron(a::Union{AbstractArray, Convex.Constant, Convex.AbstractExpr}, b::Union{AbstractArray, Convex.Constant, Convex.AbstractExpr})
	Rs = AbstractExpr[]
	for i in 1:size(a)[1]
		Vs = Convex.AbstractExpr[]
		for j in 1:size(a)[2]
			push!(Vs, a[i, j] * b)
		end
		push!(Rs, foldl(hcat, Vs))
	end
	return foldl(vcat, Rs)
end

function simplekron(a::AbstractMatrix, b::Convex.AbstractExpr)
	return simplekron(Constant(a), b)
end

function simplekron(a::AbstractExpr, b::AbstractMatrix)
	return simplekron(a, Constant(b))
end
