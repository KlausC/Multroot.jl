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
    delta = S(100)	    # threshold growth factor, default = 100
    thresh = tol*100	# zero singular threshold, default = 100*tol
    drop = S(5.0e-5)	# try Gauß-Newton if sigma(k) < drop * sigma(k-1)

    E = one(T); RE = real(E); Z = zero(T); RZ = real(Z)

    a = findfirst(!iszero, p)
    a !== nothing || throw(ArgumentError("gcdroot of zero polynomial not allowed"))
    # make the polynomial monic
    k = length(p)
    n = k - a
    p = f = view(p, a:k) / p[a] # working copy of the monic polynomial 
    g = deriv(f) / n            #  g(x) = f'(x) / degree(f)
    nf = norm(f, Inf)           # the largest coefficient

    mx = n; wtol = tol; s0 = RZ; s = RZ; wthrh = thresh

    z = complex(T)[]        # all different roots
    l = Int64[]             # the multiplicities of those roots as known so far
    nroots = 0              # invariant: length(z) == length(l) == nroots

    k = n;                  # the degree of working polynomial
    while k >= 1
        m = 1 # ensure this variable is used as loop index and value preserved after break
        if k == 1     # the polynomial is linear, GCD = 1
            h = E; u = f; v = E
        else
            for mm = 1:k
                m = mm
                A = sylves(f, g, m)
                s0 = s  # save previous sigma
                
                s, x = zminsv(A, tol)
                if m > 1 && s > s0 + norm(f, Inf) * eps()
                    m -= 1
                    break
                end
                # @printf("s. value %g,%g,%g\n", m, k, s)
                if x[1] == 0
                    x[1] = eps()
                end
                # I don't understand why s is compared to nf.
                if s < wthrh * nf || m == mx || s < drop * s0
                    # h = gcd(f, g); u = f / h; v = g / h
                    h, u, v, res0, res, sm = gcd_refinement(x, m, f, g, A)
                    if m == length(u) - 1 && (res < wtol || m == mx)
                        wtol = max(wtol, res * gamma)	     # increase tolerance by factor 
                        wthrh = max(wthrh, (s / nf) * delta) # increase threshold by factor
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
            nroots = m
            l = ones(Int64, nroots) 
            if m == 1
                # only one root value with multiplicity n
                l[1] = k
                k = 0
            end
            # verify roots of u are roots of f: horner_analysis(z, f)
        else
            t = roots(u)	# u has only simple roots - they should be in z
            jj = 0
            ix = Vector(1:nroots)
            for j = 1:m
                tj = t[j]
                tz, jm = findmin(abs.(z[ix] .- tj))	# find root closest to tj
                jj = ix[jm]
                deleteat!(ix, jm)
                l[jj] += 1                     # increase multiplicity of that root
                # println("k = $k m = $m: z[$jj] = $(z[jj]) tz = $tz")
                tz <= 0.1 * (abs(z[jj])+0.01) || error("inacceptable root, k = $k tz = $tz")
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
        f = h               # f = gcd(f, f') -- all positive multiplicities decreased by 1 
        g = deriv(f) / k
        nf = norm(f, Inf)
        mx = m
    end # k-while loop
end

function finish(z, l, p)
    f = polymult(z, l)
    w = [ max(abs(pj), 1) for pj in p]
    # println("f = $(size(f)) p = $(size(p)) w = $(size(w))")
    bkerr = maximum(abs, (f - p) ./ w)

    # println("finish: bkerr = $bkerr")
    # display([f p (f-p) ./ w])

    PolyZeros(z, l), bkerr
end

function gcd_refinement(x, m, f, g, A)
    u0, v0 = extract_sylves(x)
    # println("after extract_sylves: x = $x\nu0 = $u0\nv0 = $v0")
    # println("refinement u0: $(norm(convolute(u0,g) - convolute(v0,f))/norm(f))")
    B = cauchymt(u0, length(f) - length(u0) + 1)
    g0 = scalsq(B, f)	# scaled least squares solver g0 = div(f, u0)
    g0 = g0 / g0[1]		# first approximation of gcd(p, p')
    # printresidue("g0", g0, u0, v0, f, g)
    h, u, v, res0, res = gcdgn(f, g, g0, u0, v0)	# refinement of g, u, v by G-N
    # println("refinement u:  $(norm(convolute(u,g) - convolute(v,f))/norm(f))")
    # printresidue("h ", h, u, v, f, g)
    sm = norm(A * [u; v]) / norm([u;v])
    h, u, v, res0, res, sm
end

function printresidue(text, g0, u0, v0, f, g)
    uuu = convolute(g0, u0) - f	# should be zero
    vvv = convolute(g0, v0) - g	# should be zero

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

function scalerows!(A::Array{<:Number,2})
    s = map(i->2.0^-exponent(norm(A[i,:])), 1:size(A,1))
    scale!(s, A)
    s
end

scale!(s, A) = lmul!(Diagonal(s), A)

