"""
    multroot(p; tol, scale=false)
  
  Finds all roots of a real or complex polynomial p,
 
        p(x) = a_1 x^n + a_2 x^n-1 + ... + a_n x + a_n+1,

  given by an (n+1)-row vector `p = (a_1, ..., a_n+1)` 

  METHOD : for details, contact the author 

             Zhonggang Zeng
             Department of Mathematics
             Northeastern Illinois University
             Chicago, IL 60625
          
             email: zzeng@neiu.edu

 This code is freely released for research exchange only. The author 
 is not responsible for any demage caused by using this code.

  INPUT :  `p = polynomial coefficients.` 
           
  OUTPUT : `z = distinct roots of polynomial p`

            l = corresponding multiplicities
  CALL :   
               z = multroot(p). 
"""
function multroot(p::AbstractVector{T}; tol::S = 1e-10, scaling::Bool=false) where {T<:Number,S<:Real}
    #
    # clear leading/trailing zeros
    #
    jj = findall(!iszero, p) 
    isempty(jj) && throw(ArgumentError("zero polynomial has no roots"))
	Z = zero(T)
	E = one(T)
    n = length(p)
    j1, j2 = extrema(jj)
    q = p[j1:j2]
    if j1 >= j2
        return PolyZeros(T[], Int[]), Z, 0
    end

    scale = scale!(q, scaling)
	z0, bke = gcdroot(q, tol)

    if bke < S(1.0e200)
println("before pjeroot")

        z1, bkerr, pjcnd, job = pejroot(q, z0.z, z0.mult)
        
        # z1, bkerr, pcjnd, job = z0, bke, Z, 1
        z1.z .*= scale
        z = z1.z
        l = z1.mult
        if j2 < n
            push!(z1.z,  Z)
            push!(z1.mult, n-j2)
        end
        if job == 1
            #
            # show off results
            #
			nz = length(z)
            @printf("\n");
            @printf("    !!!THE COMPUTATION IS SUCCESSFUL!!!\n")
            @printf("\n")
            @printf("THE PEJORATIVE CONDITION NUMBER:       %g \n", pjcnd)
            @printf("THE BACKWARD ERROR:                    %6.2e %6.2e \n", bkerr, bke)
            @printf("THE ESTIMATED FORWARD ROOT ERROR:      %6.2e \n", 2 * bkerr * pjcnd)
            @printf("\n");
            if norm(imag(z)) == Z 
                @printf("        computed roots         multiplicities\n")
                @printf("\n")
				for j = 1:nz @printf("%25.15f \t \t \t %3d \n", z[j], l[j]) end
            else
                @printf("                  computed roots ")
                @printf("                 multiplicities\n");
                @printf("\n");
                for j = 1:nz
				    @printf("%22.15f + %22.15f i \t %3d \n", real(z[j]), imag(z[j]), l[j])
				end
            end
        else
            z1 = roots(p)	# fallback to standard rootfinder
			l = ones(n - 1)
        end
    else
        z1 = roots(p)		# fallback to standard rootfinder
		l = ones(T, n - 1)
		job = 0
		bkerr = S(Inf)
		pjcnd = S(Inf)
    end
    z1, bkerr, pjcnd, job
end

function scale!(p::Vector, scaling::Bool)
    n = length(p) - 1
    a = p[1]
    c = p[1]
    p ./= a
    if scaling
        b = abs(p[n+1])
        p[n÷2+2:n+1] ./= b
        c = b ^ (1 / n)
        cc = c
        for k = 2:(n+1)÷2
            p[k] /= cc
            p[n+2-k] *= cc
            cc *= c
        end
        if iseven(n)
            p[n÷2+1] /= cc
        end
    end
    c
end
