# Permute systems in a state rho_ABCD....
#
# This is a port of QETLAB's permutesystems to Julia.

export permutesystems

""" `rhoOut = permutesystems(rho, perm)` or `rhoOut = permutesystems(rho, perm, dim)` or `rhoOut = permutesystems(rho,perm, dim, row_only)` or `rhoOut = permutesystems(rho,perm, dim, row_only, inv_perm)`

Permutes the systems in the matrix rho.

Inputs:
- *rho* input matrix
- *perm* a permutation vector on all the systems
- *dim* (optional keyword argument) vector with dimensions of all the systems. By default all systems are of equal dimension.
- *row_only* (optional keyword argument) if set to true, only the rows of *rho* are permuted.
- *inv_perm* (optional keyword argument) if set to true, the subsystems are permuted according to the inverse of *perm* rather than *perm* itself.

Outputs:
- *rhoOut* matrix *rho* with permuted systems
"""

function permutesystems(X, perm, dim = false, row_only = false, inv_perm = false)
	dX = collect(size(X));
  	is_vec = (minimum(dX) == 1) || (length(dX) == 1);
  	num_sys = length(perm);
  	if is_vec
    		if length(dX) == 1
      			vec_orien = 2;
    		else
      			vec_orien = 3 - findfirst(x -> x == 1, dX); # 1 if column vector, 2 if row vector;
    		end
  	end

  	# set optional argument defaults: dim=round(lX^(1/num_sys)), row_only=0, inv_perm=0
  	if dim == false
    		dim = [round(dX[1]^(1/num_sys))*ones(1, num_sys); round(dX[2]^(1/num_sys))*ones(1, num_sys)]
  	end

  	# allow the user to enter a vector for dim if X is square
  	if minimum(size(dim)) == 1 || length(size(dim)) == 1
    		dim_tmp = dim[:]' # force dim to be a row vector
    		if is_vec
      			dim = ones(2, length(dim))
      			dim[vec_orien, :] = dim_tmp
    		else
      			dim = [dim_tmp; dim_tmp]
    		end
  	end

  	prod_dimR = prod(dim[1, :])
  	prod_dimC = prod(dim[2, :])

  	# Do some basic input checking.
  	@assert length(perm) == num_sys "length(PERM) must equal length(DIM)."
  	@assert sort(perm) == collect(1:num_sys) "PERM must be a permutation vector."
  	# Permuting systems for pure states is easy enough, so just make the vector
  	# full and then perform the permutation (new-ish versions of MATLAB don't
  	# like sparse multidimensional arrays).
  	if is_vec
    		dim = round(Int, dim)
    		if inv_perm
      			PX = reshape(ipermutedims(reshape(full(X), tuple(dim[vec_orien, end:-1:1]...)), num_sys + 1 - perm[end:-1:1]), dX[1])
    		else
      			PX = reshape(permutedims(reshape(full(X), tuple(dim[vec_orien,end:-1:1]...)), num_sys + 1 - perm[end:-1:1]), dX[1])
    		end
    		#
    		# Preserve the sparsity of X.
    		if issparse(X)
        		PX = sparse(PX)
    		end
    		return PX
  	end

  	# If X is not a pure state, it's slightly trickier... do *not* just use the
  	# same pure state trick with repeated indices though, since that has an
  	# intermediate step of making the matrix a multidimensional array, which
  	# you can't do with sparse matrices in new-ish version of MATLAB. The trick
  	# used here reduces the problem to the pure state version of the problem in
  	# another way that plays nicely with both full and sparse matrices
  	dim = round(Int, dim)
  	row_perm = permutesystems(collect(1:dX[1]), perm, dim[1, :], false, inv_perm)
  	PX = X[row_perm, :]
  	if ! row_only
    		col_perm = permutesystems(collect(1:dX[2]), perm, dim[2, :], false, inv_perm)
    		PX = PX[:, col_perm]
  	end
  	return PX;
end
