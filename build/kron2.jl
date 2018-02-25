# Fix of kron to work with Convex.jl in seesawAlice

using Convex

export kron2

function kron2(a::Union{AbstractArray, Convex.Constant, Convex.AbstractExpr}, b::Union{AbstractArray, Convex.Constant, Convex.AbstractExpr})
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

function kron2(a::AbstractMatrix, b::Convex.AbstractExpr)
	return kron2(Constant(a), b)
end

function kron2(a::AbstractExpr, b::AbstractMatrix)
	return kron2(a, Constant(b))
end
