import Polynomials
using Polynomials

function polymult(z::AbstractVector{T}, l::AbstractVector{Int}) where {T<:Number}
#
# construct monic polynomial from its roots
# use sorting and a special order of processing to increase accuray.
# special implementation of poly()
# returns Real polynomal, if all conjugate roots with equal multiplicities.
#
#   INPUT:
#            z vector of distinct roots of polynomial
#            l vector of corresponding multiplicities
#   OUTPUT:
#            ff Poly representing the product (x - z[k]) ^ l[k]
#
    # println("polymult($(z), $(l)")
    # z = PolyZeros(z)
    products(z, l)
end

# Calculate the coefficient operator of the pejorative manifold
# PolyRoots is used to leverage calculation of (x-z[j]) - products
#
function polyG(z, mult)
    products(z, mult)[2:end]
end

# The derivative of G at point z
function polyDG(zz::AbstractVector{S}, ll) where S
    m = length(zz)
    s = products(zz, ll .- 1)
    n0 = length(s)
    n = sum(ll)
    Df = zeros(S, n, m)

    for j = 1:m
        Df[1:n0,j] .= s
        ix = n0
        for k = 1:m
            if k != j && ll[j] > 0
                ix = product1!(ix, view(Df,1:n,j), 1, -zz[k])
            end
        end
        Df[1:n,j] .*= -ll[j]
    end
    isreal(Df) ? real(Df) : Df
end

"""

Calculate polynomial product of `p`and `q` given as vectors.
"""
product

function product!(r::AbstractVector{U}, p::AbstractVector{S}, q::AbstractVector{T}) where {U<:Number,S<:Number,T<:Number}

  np = length(p)
  nq = length(q)
  nr = length(r)
  nr >= np + nq - 1 || throw(ArgumentError("size of result vector must be exactly (length(a) + length(b) - 1 but $np + $nq - 1 != $nr"))

  nq > np && return product(q, p)
  U == promote_type(S, T, U) || throw(ArgumentError("result type not compatible with input"))

  if nq == 0
    return S[]
  end
  E = one(U)
  copy!(r, p)
  fill!(view(r, np+1:nr), zero(U))
  q1 = U(q[1])
  if q1 != E
    for j = 1:np
      r[j] *= q1
    end
  end
  for k = 2:nq
    q1 = U(q[k])
    for j = 1:np 
      jk = j + k - 1
      r[jk] = muladd(q1, p[j], r[jk])
    end
  end
  r
end

function product(p::AbstractVector{S}, q::AbstractVector{T}) where {S<:Number,T<:Number}
  U = promote_type(S,T)
  product!(zeros(U, length(p) + length(q) - 1), p, q)
end



"""
    products(z::Vector{T}, ll::Vector{Int})

Calculate coefficients of polynomial given by its roots with multiplicities.
The z should be ordered by increasing absolute value for better accuracy. 
"""
function products(z::AbstractVector{T}, ll::AbstractVector{Int}) where T <:Number
    m = length(z)
    m == length(ll) || throw(ArgumentError("root vector and muliplicities have differnet length"))
    all(ll .>= 0) || throw(ArgumentError("multiplicities must not be negative"))
    n = sum(ll)
    S = coeffstype(z, ll)
    r = zeros(S, n+1)
    r[1] = 1
    Z = zero(S)
    ix = 1
    jj = sortperm(z, by=abs2, rev=false)
    z = view(z, jj)
    ll =  ll[jj]
        while sum(ll) > 0
            p = minpos(ll)
            jm = findall(x-> x >= p, ll)
            s = singprod(z[jm], S)
            s = power(s, p)
            ix = convolute!(ix, r, s)
            ll .-= p
        end
    r
end

"""
    singprod(z::Vector, S::Type))

Return coefficients of polynomial `prod((x - zz) for zz in z)`

"""
function singprod(z::AbstractVector{T}, ::Type{S}) where {T<:Number,S<:Number}
    m = length(z)
    r = Vector{S}(undef, m+1)
    r[1] = 1
    ix = 1
    if S<:Real && T<:Complex
        for k = 1:m
            zk = z[k]
            rr = -real(zk)
            ri = imag(zk)
            if iszero(ri) 
                ix = product1!(ix, r, 1, rr)
            elseif ri > 0
                ix = product2!(ix, r, 1, 2*rr, abs2(zk))
            end
        end
    else # T is complex, but coefficient type is real
        for k = 1:m
            ix = product1!(ix, r, 1, -z[k])
        end
    end
    r
end

"""
    minpos(v::Vector{<:Integer})

Return the minimal positive integer in v
"""
function minpos(x::AbstractVector{T}) where T<:Integer
    y = typemax(T)
    minimum(x .+ y) - y
end

"""
    convolute!(ix, r, p)

If `r` is a vector representing a polynomial in `r[1:ix]` and p is a vector
containing another polynomial, modify `r` to contain the product.
Assume `ix >= 0` and `length(r) >= ix + length(p) - 1`.
"""
function convolute!(ix::Integer, r::AbstractVector, p::AbstractVector)
    ix <= 0 && return ix
    n = length(p) - 1
    ri = r[ix]
    for j = 0:n
        r[ix+j] = ri * p[j+1]
    end
    for i = ix-1:-1:1
        ri = r[i]
        for j = n:-1:1
            r[i+j] += ri * p[j+1]
        end
        r[i] = ri * p[1]
    end
    ix + n
end

"""
    convolute(p, q)

Return the convolution (polynomial product) of vectors p and q.
The coefficient order does not matter. No checks for zero coefficients.
"""
function convolute(p::AbstractVector{T}, q::AbstractVector{S}) where {T,S}
    R = promote_type(T,S)
    m = length(p)
    n = length(q)
    (n == 0 || m == 0) && return R[]
    r = similar(p, R, n + m - 1)
    copyto!(r, p)
    convolute!(m, r, q)
    r
end

"""
    product1!(ix, r, ll, p1)

multiply by (x + p1)^ll. ix is the current used size of r.
"""
function product1!(ix::Int, r::AbstractVector{T}, ll::Int, zk::T) where T<:Number
    for j = 1:ll
        r[ix+1] = r[ix] * zk
        for i = ix:-1:2
            r[i] = muladd(r[i-1], zk, r[i])
        end
        ix += 1
    end
    ix  
end

"""
    product2!(ix, r, ll, p1, p2)

multiply by (x^2 + p1*x + p2)^ll. ix is the current used size of r.
"""
function product2!(ix::Int, r::AbstractVector{T}, ll::Int, p1::T, p2::T) where T<:Number
    for j = 1:ll
        r[ix+1] = r[ix+2] = 0  
        for i = ix:-1:1
            rr = muladd(r[i+1], p1, r[i+2])
            r[i+2] = muladd(r[i], p2, rr)
        end
        r[2] = muladd(r[1], p1, r[2])
        ix += 2
    end
    ix
end

function power(p::AbstractArray, n::Integer)
    (Poly(p)^n).a
end

function Polynomials.roots(pc::AbstractVector{T}) where T<:Number
# redirect to implementation in Polynomials
	p = Poly(reverse(pc))
	roots(p)
end

function deriv(p::AbstractVector{T}) where T<:Number

    n = length(p) - 1
    map(k -> (p[k] * (n-k+1)), 1:n)
end
