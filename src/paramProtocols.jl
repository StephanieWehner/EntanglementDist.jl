# Gives the fidelity and success probability of several well known
# distillation protocols for 2 copies

export bennettParam
export deutschParam
export eplParam

import Base.convert
convert(::Type{Float64}, x::Array{Float64,1}) = x[1]
convert(::Type{Float64}, x::Array{Float64,2}) = x[1]

toMatrix(x::VecOrMat) = x*x'

e0 = [1; 0]
e1 = [0; 1]

⊗(x, y) = kron(x, y)
Φ_plus = 1/√2 * (e0 ⊗ e0 + e1 ⊗ e1)
Φ_min  = 1/√2 * (e0 ⊗ e0 - e1 ⊗ e1)
Ψ_plus = 1/√2 * (e0 ⊗ e1 + e1 ⊗ e0)
Ψ_min  = 1/√2 * (e0 ⊗ e1 - e1 ⊗ e0)
basis = Any[Φ_plus, Φ_min, Ψ_plus, Ψ_min]

function toBellBasis(ρ::AbstractMatrix)
  output = zeros(size(ρ))
  i = 1
  for x in basis
    j = 1
    for y in basis
      output[i, j] = (x' * ρ * y)[1]
      j += 1
    end
    i += 1
  end
  return output
end

""" `(F, p_succ) = deutschParam(rho)`

Determines the performance of DEJMPS for an input state rho and returns the fidelity and success probability achieved.

"""
function deutschParam(rho::AbstractMatrix)
  M = toBellBasis(rho);
  m1 = M[1, 1] # phi_+
  m2 = M[2, 2] # phi_-
  m3 = M[3, 3] # psi_+
  m4 = M[4, 4] # psi_-

  p1, p2, p3, p4 = sort([m1; m2; m3; m4], rev = true)

  p = (p1 + p4)^2 + (p2 + p3)^2
  F = (p1^2 + p4^2)/p
  return (F, p)
end

""" `(F, p_succ) = bennettParam(F)
Determines the performance of BBPSSW for an input state with fidelity F
returns (F_out, p_succ)
"""
function bennettParam(F::Real)
  p_succ = F^2 + 2F*(1-F)/3 + 5((1-F)/3)^2
  F_out = (F^2 + ((1-F)/3)^2)/p_succ
  return (F_out, p_succ)
end

"""
Determines the performance of BBPSSW for an input state ρ
returns (F_out, p_succ)
"""
function bennettParam(ρ::AbstractMatrix)
  F = max(diag(toBellBasis(ρ))...)
  return bennettParam(F)
end

"""
Determines the performance of EPL for Ronald 2 states
returns (F_out, p_succ)
"""
function eplParam(p::Number, pd::Number)
  return (pd, 0.5*p^2)
end
