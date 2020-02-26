"""
Global choice of form of sylvester matrix.
"""
global const sylves = sylves2
global const extract_sylves = extract_sylves2

"""
  gcdroot(z::Vector{Number} [, tol=1e-10]) -> PolyZeros, bkerror

  Calculate the multiplicity structure and initial root 
  approximation of a real or complex polynomial p,
 
        p(x) = a_1 x^n + a_2 x^n-1 + ... + a_n x + a_n+1,

  given by an (n+1)-row vector p = (a_1, ..., a_n+1) 

  METHOD : details can be found in a paper of the author  Zhonggang Zeng

  INPUT :  p = polynomial coefficients. 
           tol = zero remainder threshold; 
                 default = 10^{-10}, set internally
                 if no value is specified.

  OUTPUT : z = roots of polynomial p
           l = corresponding multiplicities
           bkerr = backward error

  CALL :   It is best to call pzero with only one argument as   
               `z = gcdroot(p)`. 

"""
function gcdroot(p::AbstractVector{T}, tol::S = 1e-10) where {T<:Number,S<:AbstractFloat}
    gamma = S(100)		# residual growth factor, default = 100
    delta = S(100)	  # threshold growth factor, default = 100
    thresh = tol*100	# zero singular threshold, default = 100*tol
    drop = S(5.0e-5)	# try Gauß-Newton if sigma(k) < drop * sigma(k-1)

    E = one(T); RE = real(E); Z = zero(T); RZ = real(Z)

    # make the polynomial monic
    if p[1] == Z
        p = p[findfirst(p):end]
    end
    n = length(p) - 1
    p = p / p[1]
    q = deriv(p) / n #  q(x) = p'(x) / degree(p)

    f = float(p); g = float(q)	# working copies of the polynomials 
    nf = maximum(abs, f)	 			# the largest coefficient

    mx = n; wtol = tol; s0 = RZ; s = RZ; wthrh = thresh

    k = n;            # the degree of working polynomial
    z = Complex{T}[]
    l = Int64[]
    while k >= 1
        m = 1 # ensure this variable is used as loop index and value preserved after break
        if k == 1     # the polynomial is linear, GCD = 1
            h = E; u = f; v = E
        else
            for mm = 1:k
                m = mm
                A = sylves(f, g, m)
                scalerows!(A)
                s0 = s  # save previous sigma
            
                s, x = zminsv(A, tol)
                # @printf("s. value %g,%g,%g\n", m, k, s)
                
                # I don't understand why s is compared to nf.
                if s < wthrh * nf || m == mx || s < drop * s0
                    h, u, v, res0, res, sm = gcd_refinement(x, m, f, g, A)
                    if res < wtol || m == mx
                        wtol = max(wtol, res * gamma)	# increase tolerance by factor 
                        wthrh = max(wthrh, (s / nf) * delta)	# increase threshold by factor
                        break # for mm
                    end
                end
            end # for mm
        end
        # println("after m-loop: m = $m mx=$mx k=$k n=$n s=$s thr=$(wthrh*nf) drop=$(drop*s0) ")
        if k == n			# the root values of u contain all roots of f 
            z = roots(u)	# u has only simple roots
            # assert length(z) == m
            length(z) == m || error("number of roots of u $(length(z)) is not $m")
            l = ones(Int64, m) 
            if m == 1
                # only one root value with multiplicity n
                l[1] = k
                k = 0
            end
            # verify roots of u are roots of f: horner_analysis(z, f)
        else
            t = roots(u)	# u has only simple roots - they should be in z
            jj = 0
            for j = 1:m
                tj = t[j]
                tz, jj = findmin(abs.(z .- tj))	# find root closest to tj
                ljp = l[jj] + 1
                l[jj] = ljp
				# z[jj] += (tj - z[jj]) / ljp # store mean value with weights l[jj] and 1
                println("k = $k m = $m: z[$jj] = $(z[jj]), tz = $tz") 
                # tz <= 0.01 * abs(z[jj]) || error("inacceptable root, k = $k tz = $tz")
            end
            if m == 1
                l[jj] += k - 1
                k = 0
            end
        end
        k = max(k - m, 0)
        if k <= 0
            return finish(z, l, p)
        end
        f = h
        g = deriv(f) / k
        nf = norm(f)
        mx = m
    end # k-while loop
end

function finish(z, l, p)
    f = polymult(z, l)
    w = ones(eltype(p), length(p))
    for j = 1:length(p)
        apj = abs(p[j])
        if apj > 1.0
            w[j] = 1.0 / apj
        end
    end
    # println("f = $(size(f)) p = $(size(p)) w = $(size(w))")
    bkerr = maximum(abs, (f - p) .* w)
   
    PolyZeros(z, l), bkerr
end

function gcd_refinement(x, m, f, g, A)
    u0, v0 = extract_sylves(x)
    #println("refinement u0: $(norm(conv(u0,g) + conv(v0,f))/norm(f))")
    B = cauchymt(u0, length(f) - length(u0) + 1)
    g0 = scalsq(B, f)	# scaled least squares solver g0 = div(f, u0)
    g0 = g0 / g0[1]		# first approximation of gcd(p, p')
    #printresidue("g0", g0, u0, v0, f, g)
    h, u, v, res0, res = gcdgn(f, g, g0, u0, v0)	# refinement of g, u, v by G-N
    #println("refinement u: $(norm(conv(u,g) - conv(v,f))/norm(f))")
    #printresidue("h ", h, u, v, f, g)
    sm = norm(A * [u; v]) / norm([u;v])
    h, u, v, res0, res, sm
end

function printresidue(text, g0, u0, v0, f, g)
    uuu = conv(g0, u0) - f	# should be zero
    vvv = conv(g0, v0) - g	# should be zero

    # println("$text p(x)  $f")
    # println("$text p'(x)  $g")
    println("$text u(x)  $u0")
    println("$text v(x)  $v0")
    # println("$text h(x)  $g0")
    # println("$text Δp(x)  $uuu")
    # println("$text Δp'(x) $vvv")
end

function analysis_sv(k, m, A, x, s) 

    #println("svcheck: $k-$m norm(Ax) = $(norm(A*x)) s = $s")


end

function scalerows!(A::Array{Float64,2})
	s = map(i->2.0^-exponent(norm(A[i,:])), 1:size(A,1))
    scale!(s, A)
	s
end

scale!(s, A) = lmul!(Diagonal(s), A)

