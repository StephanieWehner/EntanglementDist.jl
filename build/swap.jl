function swap(X, swap; dim = false, row_only = false, inv_perm = false, stop = false)
  @assert dim !== false
  perm = collect(1:length(dim))
  perm[swap[1]] = swap[2]
  perm[swap[2]] = swap[1]
  return permutesystems(X, perm, dim = dim, row_only = row_only, inv_perm = inv_perm, stop = stop)
end
