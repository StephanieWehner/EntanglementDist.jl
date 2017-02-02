# Permute systems in a state rho_ABCD....
#
# This is a port of QETLAB's permutesystems to Julia.

export permutesystems

function permutesystems(X, perm; dim = false, row_only = false, inv_perm = false, stop = false)
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

  	if stop
    		return
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
  	row_perm = permutesystems(collect(1:dX[1]), perm, dim = dim[1, :], row_only = false, inv_perm = inv_perm, stop = false)
  	PX = X[row_perm, :]
  	if ! row_only
    		col_perm = permutesystems(collect(1:dX[2]), perm, dim = dim[2, :], row_only = false, inv_perm = inv_perm, stop = false)
    		PX = PX[:, col_perm]
  	end
  	return PX;
end
