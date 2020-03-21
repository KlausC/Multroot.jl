
export PolyZeros
import Polynomials: roots, degree

"""
    PolyZeros(z, mult)

represents a set of distinct (potential) roots z[j] of a polynomial with their
multiplicities. If the same root value appears repeatedly, the surplus occurrences
are discarded and the multiplicity addapted accordingly.
"""
struct PolyZeros{S<:Number}
    z::AbstractVector{S}		# the distinct root values
    mult::AbstractArray{Int}	# the multiplicities of the roots

    function PolyZeros(zz::AbstractVector{S}, ll::AbstractVector{<:Integer} = Int[]) where S<:Number
        zz, ll = filter_duplicates(zz, ll)
        new{S}(zz, ll)
    end
end

function PolyZeros(A::AbstractMatrix{T}) where T
    m, n = size(A)
    n == 2 || throw(ArgumentError("need exactly 2 colums"))
    zz = A[:,1]
    ll = Int.(A[:,2])
    PolyZeros(zz, ll)
end

Base.isreal(z::PolyZeros) = z.real_coefficients
multiplicities(z::PolyZeros) = z.mult
roots(z::PolyZeros) = z.z
degree(z::PolyZeros) = sum(abs.(z.mult))

"""
  `lessabs(a,b) -> Bool`

Order by absolute value of a and b.
For equal absolute values first negative real part, then positive imaginary part first.

For example 0 < -1 < 1 < -4+3im < -4-3im < -3-4im < 3+4im < 5  

"""
function lessabs(a::S, b::T) where {S<:Number,T<:Number}
    diff = abs2(a) - abs2(b)
    Z = zero(diff)
    if diff != Z || ( isreal(a) && isreal(b) )
      return diff < Z || ( diff == Z && real(a) < real(b) )
    end
    real(a) < real(b) || ( real(a) == real(b) && imag(a) > imag(b) )
end

function less_real_first(a::S, b::T) where {S<:Number, T<:Number}
    ra, rb = real(a), real(b)
    ia, ib = imag(a), imag(b)
    ra < rb ? true : ra > rb ? false : abs(ia) < abs(ib) ? true : abs(ia) > abs(ib) ? false : ia < ib
end

"""
    zerostructure(zz::PolyZero) -> Type, Vector{Int}
    zerostructure(z::Vector, mult::Vector)

Return the needed coefficients type (real or complex) of polynomial with given zeros,
and for each zero an integer indication with the following meaning:

```
== 0:  real root
 < 0:  complex root without corresponding conjugate root of same multiplicity
 > 0:  complex root with corresponding conjugate root at position given by indication
```
 """
function zerostructure(zz::PolyZeros{<:Number})
    zerostructure(zz.z, zz.mult)
end
function zerostructure(z::AbstractVector{S}, ll::AbstractVector{<:Integer}) where S<:Number
    n = length(z)
    res = zeros(Int, n)
    isreal(z) && return real(S), res
    for k = 1:n
        res[k] != 0 && continue
        zk = z[k]
        isreal(zk) && continue
        res[k] = -1
        lk = ll[k]
        for j = k+1:n
            zj = z[j]
            lj = ll[j]
            if lj == lk && zk == conj(zj)
                res[k] = j
                res[j] = k
                break
            end
        end
    end
    T = all(res .>= 0) ? real(S) : S
    return T, res
end

"""
    `filter_duplicates(roots, multiplicities)`

Remove roots with multiplicities <= 0 and join duplicates
(add multiplicities for equal root values).
"""
function filter_duplicates(z::AbstractVector{<:Number}, ll::AbstractVector{<:Integer})
    n = length(z)
    if length(ll) == 0
        ll = ones(Int, n)
    else
        length(ll) == n ||
            throw(ArgumentError("roots and multiplicities have differing sizes"))
        all( ll .> 0 ) || throw(ArgumentError("only positive multiplicities allowed"))
    end
    lorig = ll
    for k = 1:n
        zk = z[k]
        for j = k+1:n
            zj = z[j]
            lj = ll[j]
            if zk == z[j] && lj > 0
                if ll === lorig
                    ll = copy(ll)
                end
                ll[k] += lj
                ll[j] = 0
            end
        end
    end
    if ll !== lorig
        jj = findall(iszero, ll)
        z = deleteat!(copy(z), jj)
        deleteat!(ll, jj)
    end
    z, ll
end

#
# calculate product over j of (x - z[j])^ll[j]
#
function evaluate(z::PolyZeros{S}, x::T, ll::AbstractVector{Int}) where {S,T}

  l = ll
	s = one(promote_type(T,S))
	n = length(ll)
	i = 1
	zz = z.z
	while i <= n
	    zi = x - zz[i]
		li = l[i]
	    if i < n && zz[i+1] == conj(zz[i])
		    m = min(li, l[i+1])
			if m > 0
		        s = s * (real(zi)^2 + imag(zi)^2) ^ m
			end
			if li > m
			    s = s * zi ^ (li - m)
			else
			    s = s * conj(zi) ^ (l[i+1] - m)
			end
			i += 2
		else
	        s = s * zi ^ li
			i += 1
		end
    end
    isreal(s) ? real(s) : s
end
evaluate(z::PolyZeros{S}, x::T) where {S,T} = evaluate(z, x, z.mult)

"""
"""
function Base.sort(z::PolyZeros)
    m = length(z.z)
    flist = [
             a::Int->z.mult[a],
             a::Int->abs(z.z[a]),
             a::Int->real(z.z[a]),
             a::Int->imag(z.z[a])
            ]
    wlist = [
             1000,
             1.0,
             0.1,
             0.01
            ]
    ev(f, a) = f(a)
    function lt(a::Int, b::Int)
        dot(ev.(flist, a), wlist) < dot(ev.(flist, b), wlist)
    end
    jj = sortperm(1:m, lt=lt)
    PolyZeros(z.z[jj], z.mult[jj])
end 

function diffmult(za::PolyZeros, zb::PolyZeros)
    maximum(abs, za.mult - zb.mult)
end
function diffzeros(za::PolyZeros, zb::PolyZeros)
    norm(za.z - zb.z, Inf)
end

"""
    subsetapprox(p::Vector{S}, q::Vector{S}, dist::Function)

Let `length(p) <= length(q)` and `S` have a distance function `dist(::S,::S)`.
Find a vector `x` of indices into `q` with `length(jj) == length(p)`
which minimizes the term `maximum(dist(p[i], q[x[i]) for i in 1:length(p))`.
"""
function subsetapprox(a::AbstractVector{T}, b::AbstractVector{T}, dist::Function, ::Type{S}) where {T,S<:Real}

    branch_and_bound(a, b, dist, typemax(S))
end

function branch_and_bound(a, b, dist, inbound)
    x = branch_and_bound1(a, b, dist, inbound)
    if x[1]
        println("called b&b($inbound)")
        display(dist.(b, transpose(a)))
        println("returning $x")
    end
    x
end

function branch_and_bound1(a, b, dist, inbound::S) where S
    bound, ix, y, isopt = heuristics(a, b, dist, inbound)
    isopt == 1 && return true, bound, y
    isopt == 2 && return false, Inf, y
    m = length(a)
    n = length(b)
    ii = [1:ix-1;ix+1:m]
    ai = a[ix]
    ci = bound
    xi = Vector{Int}(undef, m-1)
    ji = 0
    for j = 1:n
        cij = dist(ai, b[j])
        # cij > ci && continue
        jj = [1:j-1;j+1:n]
        ok, bij, xx = branch_and_bound(view(a, ii), view(b, jj), dist, ci)
        ok || continue
        if bij < ci
            ci = bij
            xi = jj[xx]
            ji = j
        end
    end
    return true, ci, [xi[1:ix-1];ji;xi[ix:m-1]]
end

function heuristics(a::AbstractVector{T}, b::AbstractVector{T}, dist::Function, inbounds::S) where {T,S<:Real}

    m = length(a)
    n = length(b)
    m <= n || throw(ArgumentError("first vector longer than second vector"))
    function findmin2(ai, jj)
        vi = findmin([dist(ai, b[j]) for j in jj])
        vi !== nothing ? vi : (typemax(S), 0)
    end
    c = Vector{S}(undef, m)
    y = Vector{Int}(undef, m)
    for i = 1:m
        c[i], y[i] = findmin2(a[i], 1:n)
    end
    ii = collect(1:m)
    jj = collect(1:n)
    cc = zero(S)
    ix = 0
    isopt = 1
    while !isempty(ii)
        kc = 0
        ci, k = findmax(view(c, ii))
        ic = ii[k]
        jc = y[ic]
        if ci > cc
            if iszero(cc) && ci >= inbounds
                return ci, ic, y, 2
            end
            cc = ci
            ix = ic
        end
        deleteat!(ii, k)
        deleteat!(jj, findfirst(isequal(jc), jj))
        for i in ii
            if y[i] == jc
                ci, k = findmin2(a[i], jj)
                y[i] = jj[k]
                c[i] = ci
                if ci > cc
                    cc = ci
                    ix = i
                    isopt = 0
                end
            end
        end
    end
    cc, ix, y, isopt
end

