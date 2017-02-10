# Gives the fidelity and success probability of several well known
# distillation protocols for 2 copies of a bipartite state consisting of two qubits

export toBellBasis
export BBPSSWParam
export DEJMPSParam
export EPLParam

import Base.convert
convert(::Type{Float64}, x::Array{Float64,1}) = x[1]
convert(::Type{Float64}, x::Array{Float64,2}) = x[1]

toMatrix(x::VecOrMat) = x*x'

e0 = [1; 0]
e1 = [0; 1]

proj0 = toMatrix(e0)
proj1 = toMatrix(e1)

phi_plus = 1/√2 * (kron(e0,e0) + kron(e1,e1))
phi_min  = 1/√2 * (kron(e0,e0) - kron(e1,e1))
psi_plus = 1/√2 * (kron(e0,e1) + kron(e1,e0))
psi_min  = 1/√2 * (kron(e0,e1) - kron(e1,e0))
basis = Any[phi_plus, phi_min, psi_plus, psi_min]

function toBellBasis(rho::AbstractMatrix)
  output = zeros(size(rho))
  i = 1
  for x in basis
    j = 1
    for y in basis
      output[i, j] = (x' * rho * y)[1]
      j += 1
    end
    i += 1
  end
  return output
end

""" `(F, p_succ) = DEJMPSParam(rho)`

Determines the performance of DEJMPS for two copies of the input state *rho* (consisting of two qubits) and returns the fidelity and success probability achieved.
The procedure optimises over local single qubit rotations that permute the bell diagonal coefficients before the CNOTs to achieve the highest output fidelity.
returns (F_out, p_succ)
"""
function DEJMPSParam(rho::AbstractMatrix)
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

""" `(F, p_succ) = BBPSSWParam(F)`

Determines the performance of BBPSSW for two copies of the input state (consisting of two qubits) with fidelity *F* each.
returns (F_out, p_succ)
"""
function BBPSSWParam(F::Real)
  p_succ = F^2 + 2F*(1-F)/3 + 5((1-F)/3)^2
  F_out = (F^2 + ((1-F)/3)^2)/p_succ
  return (F_out, p_succ)
end

""" `(F, p_succ) = BBPSSWParam(rho)`

Determines the performance of BBPSSW for two copies of the input state *rho* (consisting of two qubits).
returns (F_out, p_succ)
"""
function BBPSSWParam(rho::AbstractMatrix)
  F = max(diag(toBellBasis(rho))...)
  return BBPSSWParam(F)
end

""" `(F, p_succ) = EPLParam(p, pd)`

Determines the performance of EPL for rStateCorrPhase.
The inputs are the *p* and *pd* parameters of rStateCorrPhase.
returns (F_out, p_succ)
"""
function EPLParam(p::Number, pd::Number)
  return (pd, 0.5*p^2)
end

""" `(F, p_succ) = EPLParam(rho)`

Determines the performance of EPL for an arbitrary bipartite state of dimensions nA=nB=4.
returns (F_out, p_succ)
"""
function EPLParam(rho::AbstractMatrix)

  #Define the CNOT operation:
  CNOT = [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]

  #Alice and Bob apply CNOT:
  rho2 = kron(CNOT, CNOT) * rho * kron(CNOT, CNOT)'

  #Alice and Bob measure the flags and postselect them to be 11.
  #Projector on the success space on one side
  projSucc = kron(eye(2), proj1)

  #Success probability:
  p_succ = trace(kron(projSucc,projSucc) * rho2)

  #State on 2 qubit space after postselection:
  rhoOut = partialtrace((kron(projSucc,projSucc) * rho2 * kron(projSucc,projSucc)')/p_succ, [2,4], [2,2,2,2])
  
  #Output fidelity
  M = toBellBasis(rhoOut)

  F = max(M[1, 1], M[2, 2], M[3, 3], M[4, 4])

  return (F, p_succ)
end
