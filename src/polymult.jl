import Polynomials
using Polynomials

function polymult{T}(z::AbstractVector{T}, l::AbstractVector{Int})
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
    z = PolyZeros(z)
    products(z, l)
end

# Calculate the coefficient operator of the pejorative manifold
# PolyRoots is used to leverage calculation of (x-z[j]) - products
#
function polyG{S}(pz::PolyZeros{S})
  products(isreal(pz), pz.z, pz.mult)[2:end]
end

# The derivative of G at point z
function polyDG{S}(pz::PolyZeros{S})
    m = length(pz.z)
    ll = pz.mult
    s = products(isreal(pz), pz.z, ll-1)
    n = sum(ll)
    Df = zeros(S, n, m)
    p1 = zeros(S, 2)

    for j = 1:m
      copy!(view(Df,1:n,j), s)
      ix = length(s)
      for k = 1:m
        if k != j && ll[j] > 0
          zk = pz.z[k]
          ix = product1!(ix, view(Df,1:n,j), 1, zk)
          scale!(view(Df,1:n,j), -ll[k])
        end
      end
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
  `products(z::PolyZeros{S}`

Calculate coefficients of polynomial given by its roots with multiplicities.
The z should be ordered by increasing absolute value for better accuracy. 
"""
function products(creal::Bool, z::AbstractVector{T}, ll::AbstractVector{Int}) where T <:Number
  m = length(z)
  m == length(ll) || throw(ArgumentError("root vector and muliplicities have differnet length"))
  all(ll .>= 0) || throw(ArgumentError("multiplicities must not be negative"))
  n = sum(ll)
  S = coeffstype(z, ll)
  r = zeros(S, n+1)
  r[1] = one(S)
  Z = zero(S)
  ix = 1
  if !((S<:Real) && !(T<:Real))
    for k = 1:m
      zk = -z[k]
      ix = product1!(ix, r, ll[k], zk)
    end
  else # T is not real, but coefficient type is real
    for k = 1:m
      zk = -z[k]
      zi = imag(zk)
      zk = real(zk)
      if zi == Z
        ix = product1!(ix, r, ll[k], zk)
      elseif zi > Z
        zk = 2zk
        zi = abs2(z[k])
        product2!(ix, r, ll[k], zk, zi)
      end
    end
  end
  r
end

"""
multiply by (x + p1)^l. ix is the current used size of r.
"""
function product1!(ix::Int, r::AbstractVector{S}, ll::Int, zk::T) where {S<:Number,T<:Number}
  for j = 1:ll
    for i = ix:-1:1
      r[i+1] = muladd(r[i], zk, r[i+1])
    end
    ix += ll
  end
  ix  
end


"""
multiply by (x^2 + p1*x + p2)^ll. ix is the current used size of r.
"""
function product2!(ix::Int, r::AbstractVector{S}, ll::Int, p1::T, p2::T) where {S<:Number,T<:Number}
  for j = 1:ll
    for i = ix:-1:1
      rr = muladd(r[i+1], p1, r[i+2])
      r[i+2] = muladd(r[i], p2, rr)
    end
    r[2] = muladd(r[1], p1, r[2])
    ix += 2ll
  end
  ix
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
