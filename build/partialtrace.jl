# Partial trace

using Convex
using SCS

import Base.sign
import Convex.monotonicity
import Convex.conic_form!
import Convex.curvature

export ptrace, partialtrace

"""
Type to be used internally by Convex.jl
"""

type PartialTraceAtom <: AbstractExpr
	head::Symbol
	id_hash::UInt64
	children::Tuple{AbstractExpr}
	size::Tuple{Integer, Integer}
	sys::Integer
	dims::Vector

	function PartialTraceAtom(x::AbstractExpr, sys::Integer, dims::Vector)
		if x.size[1] ≠ x.size[2]
			error("Only square matrices are supported")
		end

		if ! (1 ≤ sys ≤ length(dims))
			error("Invalid system, should between 1 and ", length(dims))
		end

		if x.size[1] ≠ prod(dims)
			error("Dimension of system doesnt correspond to dimension of subsystems")
		end

		children = (x, )
		newsize = (round(Int, x.size[2]/dims[sys]), round(Int, x.size[1]/dims[sys]))

		return new(:partialTrace,
		hash(children),
		children,
		newsize,
		sys,
		dims)
	end
end

function sign(x::PartialTraceAtom)
	 return NoSign()
	 return sign(x.children[0])
 end

function curvature(x::PartialTraceAtom)
  return ConstVexity()
end

function monotonicity(x::PartialTraceAtom)
	return (NoMonotonicity(),)
end

function conic_form!(x::PartialTraceAtom, unique_conic_forms::UniqueConicForms)
	if !has_conic_form(unique_conic_forms, x)
		sys = x.sys
		dims = x.dims
		function entry(ρ, j::Integer)
			bra = speye(1)
			ket = speye(1)
			i_sys = 1
			for dim in dims
				if i_sys == sys
					vO = sparsevec([j], [1], dim);
					bra = kron(bra, vO')
					ket = kron(ket, vO)
				else
					bra = kron(bra, speye(dim))
					ket = kron(ket, speye(dim))
				end
				i_sys += 1
			end
			return bra * ρ * ket
		end

		objective = conic_form!(sum([entry(x.children[1], j) for j in 1:dims[sys]]), unique_conic_forms)
		cache_conic_form!(unique_conic_forms, x, objective)
	end
	return get_conic_form(unique_conic_forms, x)
end

function partialtrace(x, systems::Vector, dims::Vector)
	dims = copy(dims)
	out = x
	for sys in systems
		out = partialtrace(out, sys, dims)
		dims[sys] = 1
	end
	return out
end

# General function, works in Convex.jl and with normal matrices

""" `rhoOut = partialtrace(rhoIn, sysNum, dimVec)`

Computes the partial trace of *rhoIn*, tracing out system number *sysNum*, where the dimensions of the systems are specified in the vector *dimVec*. 

Example: rhoOut = partialtrace(epr, 2, [2 2])

Traces out the second qubit of a 2x2 state, here the EPR pair.
"""

function partialtrace(x, sys::Integer, dims::Vector)

	# Define a helper function returning matrix entries
	function entry(rho, j::Integer)
		bra = speye(1)
		ket = speye(1)
		i_sys = 1
		for dim in dims
			if i_sys == sys
				vO = sparsevec([j], [1], dim)
				bra = kron(bra, vO')
				ket = kron(ket, vO)
			else
				bra = kron(bra, speye(dim))
				ket = kron(ket, speye(dim))
			end
			i_sys += 1
		end
		return bra * rho * ket
	end

	# Compute the partial trace by picking the appropriate entries
	return sum([entry(x, j) for j in 1:dims[sys]])
end


ptrace(x::AbstractExpr, sys::Integer, dim::Vector) = partialtrace(x, sys, dim)
ptrace(x::AbstractArray, sys::Integer, dim::Vector) = partialtrace(x, sys, dim)
