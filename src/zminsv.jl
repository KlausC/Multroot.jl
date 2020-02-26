"""
  `zminsv(A::Array, tol) -> sval, svec`

  zminsv calculates the smallest singular value
  of matrix A and the associated right singular vector
  via an implicit inverse iteration in the form of the
  Gauss-Newton iteration

    input  A --- the matrix
           tol --- the error tolerence

   output  s --- the smallest singular value
           x --- the associated right singular vector
"""
function zminsv(A::AbstractArray{S,2}, tol::U) where {S<:Real,U<:AbstractFloat}
    E = one(S)
    m, n = size(A)           	# get the dimensions of A
    if m < n-1 throw(ArgumentError("zminsv only if m >= n-1 but $m < $n-1")) end
   
    Q, R = qr(A)        # QR decomp. of A, maybe input
    if m == n - 1				# Corner case when s == 0 is guaranteed 
        y = normalize([backsub(R[1:m,1:m], R[1:m,end]); -E])
        ss = zero(S)
	    # println("zminsv (m=$m, n=$n) returns ss= $ss y = $y")
        return ss, y
    end

    cr = U(0)
    scale = norm(A, Inf)     		    # get the magnitude of rows as scaler
	a = scale * normalize(randn(S,n))	# random initial vector (row)
    b = [scale; zeros(n)]

    T, trans = hessqr([a'; R])
    z = hqrt(trans, b)

    x = backsub(T[1:n,1:n], z[1:n])
    normalize!(x)
	
    r = [scale * x'; R] * x
    r = r - b

    for k = 1:5
		D = [ x' * scale * 2; R]		# Least squares problem min |D z - b|²
        T, trans = hessqr(D) 			# Hessenberg QR decomp. of stacked matrix
        z = hqrt(trans, r)        		# same QQ on b
        u = backsub(T[1:n,1:n], z[1:n])	# backward substitution to solve RR z = QQ' b
        y = normalize(x - u)			# new iteration
        r = [scale * y'; R] * y - b
        ss = norm(r[2:n])
        crk = norm(y-x)
		#println("k = $k ss = $ss crk = $crk")
      
        if (k == 1 && crk < tol ) || (k > 1 &&  (crk >= cr || crk^2 / (cr - crk) < tol)) 
            break
        end
        x = y
        cr = crk
    end
	# println("zminsv (m=$m, n=$n) returns ss= $ss x = $x")
    ss, y
end

