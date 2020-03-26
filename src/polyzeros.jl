
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
function subsetapprox(a::AbstractVector{T}, b::AbstractVector{T}, dist::Function) where T
    A = dist.(b, transpose(a))
    subsetapprox(A)
end

function subsetapprox(A::AbstractMatrix{S}) where S<:Real
    fi = zeros(Int, size(A, 2))
    branch_and_bound(A, fi, typemax(S))
end

function branch_and_bound2(A, fi, inbound)
    println("called b&b($fi <= $inbound)")
    res = branch_and_bound1(A, fi, inbound)
    v, x = res
    println("returning $v")
    if minimum(x) < 1 || maximum(x) > size(A,1)
        println("wrong x after b&b:")
        display(x')
    end
    if false
        vv, xx = result(A, x)
        display([x xx]')
    end
    res
end

function result(A, y)
    z = getindex.(Ref(A), y, axes(y)[1])
    target(z), z
end

function target(z)
    norm(z, 2)
end

function fix!(fj, fi, j, i)
    if !iszero(fi[i])
        fj[fi[i]] = 0
    end
    fi[i] = j
    fj[j] = i
end

function checkfi(f)
    g = f[f .> 0]
    length(g) == length(unique(g))
end

function copyfi(fi, n)
    fi = copy(fi)
    fj = zeros(eltype(fi), n)
    for k = 1:length(fi)
        if fi[k] != 0
            fj[fi[k]] = k
        end
    end
    fi, fj
end

function branch_and_bound(A::Matrix{S}, fi, inbound::S) where S<:Real
    isopt, bound, y, ix = heuristics(A, fi, inbound)
    isopt == 1 && return bound, y
    isopt == 2 && return Inf, y
    n, m = size(A)
    fi, fj = copyfi(fi, n)
    fix!(fj, fi, y[ix], ix)
    ci, xi = branch_and_bound(A, fi, bound)
    if ci >= bound
        ci, xi = bound, y
    end
    for j = sortperm(view(A, :, ix))
        iszero(fj[j]) || continue
        cij = A[j,ix]
        cij > ci && continue
        fix!(fj, fi, j, ix)
        bij, xx = branch_and_bound(A, fi, ci)
        isinf(bij) && continue
        if bij < ci
            ci = bij
            xi = xx
        end
    end
    ci, xi
end

function heuristics(A::AbstractMatrix{S}, fi::V, inbounds::S) where {S<:Real,V<:AbstractVector{Int}}

    n, m = size(A)
    m <= n || throw(ArgumentError("first vector longer than second vector"))
    function findmin2(i, fj)
        cm, jm = typemax(S), 0
        for j = axes(fj)[1]
            if iszero(fj[j]) && A[j,i] < cm
                cm = A[j,i]
                jm = j
            end
        end
        cm, jm
    end
    findmin2(i, fj, fi) = fi[i] == 0 ? findmin2(i, fj) : (A[fi[i],i], fi[i])
    function findmax2(c, fi)
        cm, ii = zero(S), 0
        for i = axes(fi)[1]
            if iszero(fi[i]) && c[i] > cm
                cm = c[i]
                ii = i
            end
        end
        cm, ii
    end

    fi, fj = copyfi(fi, n)
    c = similar(fi, S)
    y = similar(fi)

    for i = 1:m
        c[i], y[i] = findmin2(i, fj, fi)
    end
    ci = target(c)
    if ci > inbounds # the most optimistic estimation can't be better than inbounds
        return 2, ci, y, 0
    end
    ix = 0
    bb = zero(S)
    isopt = 1
    free = count(iszero, fi)
    while free > 0
        ci, ic = findmax2(c, fi)
        jc = y[ic]
        if ci > bb
            bb = ci
            ix = ic
        end
        fix!(fj, fi, jc, ic)
        free -= 1
        for i in 1:m
            if fi[i] == 0 && y[i] == jc
                c[i], y[i] = findmin2(i, fj)
                if c[i] > bb
                    bb = ci
                    ix = i
                end
                isopt = 0
            end
        end
    end
    cc = target(c)
    isopt, cc, y, ix
end

# brute-force calculation of minimum
function bf(A)
    c = typemax(eltype(A))
    x = Int[]
    for jj in allcomb(size(A,2), size(A,1))
        cjj = target(A[jj[i], i] for i in axes(A,2))
        if cjj < c
            c = cjj
            x = jj
        end
    end
    c, x
end
function allcomb(a::Integer, b::Integer)
    allcomb.(b-a+1, b, 0:binomial(b, a) * factorial(a)-1)
end
function allcomb(m::Integer, n::Integer, z::Integer)
    function op((a, b), k)
        x, y = divrem(b, k)
        [a; y], x
    end
    function op2((a,b), k)
        j = 0
        for l = 1:length(b)
            if b[l] == 0
                if j == k
                    b[l] = 1
                    return [a; l], b
                else
                    j += 1
                end
            end
        end
        (a, b)
    end
    foldl(op2, foldl(op, n:-1:m, init = (Int[], z))[1], init = (Int[], zeros(n)))[1]
end

