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
    n = length(p)
    n == 0 && throw(ArgumentError("zero polynomial has no roots"))
	Z = zero(T)
	E = one(T)
    if p[1] == Z || p[n] == Z 
        jj = findall(!iszero, p) 
        isempty(jj) && throw(ArgumentError("zero polynomial has no roots"))
        j1 = minimum(jj)
	    j2 = maximum(jj)
        q = p[j1:j2]
    else
        j2 = n;
        q = p
    end
    q = q / q[1]
    #
    # scaling
    #
    m = length(q) - 1
    if m > 0 && scaling
        c = E / (abs(q[m+1])) ^ (1/m)
        # println("c = $c, m = $m, q = $q")
   	    q = q .* (c .^ (0:m))
    else
        c = E
    end
println("before gcdroot")
	z0, bke = gcdroot(q, tol)
  
    println(z0)

    if bke < S(1.0e-2)
println("before pjeroot")
        z1, bkerr, pjcnd, job = pejroot(q, z0.z, z0.mult)
        z1.z ./= c
        z = z1.z
        l = z1.mult
        if j2 < n
            push!(z1.z,  Z)
            append!(z1.mult, zeros(Int, n-j2))
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
            @printf("THE BACKWARD ERROR:                    %6.2e \n", bkerr)
            @printf("THE ESTIMATED FORWARD ROOT ERROR:      %6.2e \n", 2 * bkerr * pjcnd)
            @printf("\n");
            if norm(imag(z)) == Z 
                @printf("        computed roots         multiplicities\n")
                @printf("\n")
				for j = 1:nz @printf("%25.15f \t \t \t %3g \n", z[j], l[j]) end
            else
                @printf("        computed roots ")
                @printf("   \t\t\t\t\t\t     multiplicities\n");
                @printf("\n");
                for j = 1:nz
				    @printf("%22.15f + %22.15f i \t \t %3g \n", real(z[j]), imag(z[j]), l[j])
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

