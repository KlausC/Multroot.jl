
# PolyZeros represents a set of distinct (potential) roots z[j] of a polynomial
# The z-values are ordered by absolute value to obtain better numeric stability
# when calculating products of the form  Prod over j (x-z[j]) ^mult[j]
# where l defines a pejorative multiplicity structure of the polynomial.
# The stored values are obtaind by z[jj] and mult[jj] from the original arrays.
#
struct PolyZeros{S<:Number}
    z::AbstractVector{S}		  # the distinct root values, already ordered
    mult::AbstractArray{Int}	# the multiplicity structure of the zeros, reordered
    real_coefficients::Bool   # indicates, that polynomial coefficients are real
    indicator::BitArray
end

# Constructor for PolyZeros - do not re-order elements
function PolyZeros(zz::AbstractVector{S}, ll::AbstractVector{Int} = Int[]) where S
    if length(ll) == 0
        ll = ones(Int, length(zz))
    end
    zz, ll = filter_duplicates(zz, ll)
    m = length(zz)
    jj = sortperm(zz, lt = lessabs)
    zz = zz[jj]
    ll = ll[jj]
    CT = coeffstype(zz, ll)
    creal = CT <: Real
    ind = BitArray(zeros(Bool,2m))
    _variables!(creal, zz, ind)
    PolyZeros(zz, ll, creal, ind)
end

import Polynomials: roots, degree
Base.isreal(z::PolyZeros) = z.real_coefficients
multiplicities(z::PolyZeros) = z.mult
roots(z::PolyZeros) = z.z
degree(z::PolyZeros) = sum(abs.(z.mult))
function var_real(z::PolyZeros, j::Int)
  checkvar(z, j)
  z.indicator[j]
end
function var_imag(z::PolyZeros, j::Int)
  m = checkvar(z, j)
  z.indicator[j+m]
end
function var_simple(z::PolyZeros, j::Int)
  m = checkvar(z, j)
  ! (z.indicator[j] || z.indicator[j+m])
end

@inline function checkvar(z::PolyZeros, j::Int)
  m = length(z.z)
  0 < j <= m || throw(ArgumentError("variable number($j) not in 1..$m"))
  m
end

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

# Determine type of polynomial coefficients given the roots (with type)
function coeffstype(z::AbstractVector{S}, ll::AbstractVector{Int}) where S<:Number
    if isreal(z)
        return real(S)
    end
    jj = findall(x-> !isreal(x), z) # extract all complex elements
    for j in jj
        jc = findfirst(x-> x == conj(z[j]), z)
        if jc === nothing || ll[j] != ll[jc]
            return S
        end
    end
    real(S)
end
coeffstype(z::PolyZeros) = coeffstype(z.z, z.mult)

"""
    `filter_duplicates(roots, multiplicities)`

Remove roots with multiplicities <= 0 and join duplicates
(add multiplicities for equal root values).
"""
function filter_duplicates(zz::AbstractVector{S}, ll::AbstractVector{Int}) where S<:Number
  m = length(zz)
  length(ll) == m ||
    throw(ArgumentError("roots and multiplicities have differing sizes $m and $(length(ll))"))
  all( ll .> 0 ) || throw(ArgumentError("only positive multiplicities allowed"))

  z = S[]
  l = Int[]
  for k = 1:m
    zk = zz[k]
    j = findfirst(x -> x == zk, z)
    if j !== nothing
      l[j] += ll[k]
    else
      push!(z, zk)
      push!(l, ll[k])
    end
  end
  z, l
end

"""

Return vector of variables corresponding to roots.
Output is vector of real or complex variables, which can be used for G-N iterations.
If polynomial is real, return real variables, where a complex root pair `x + y*im` is
represented by two real variables with values `v[j] = 2y` and `v[j+1] = (x^2+y^2)`.
The latter case is indicated by negation of multiplicities in the structure.
(`z.mult[j] == z.mult[j+1] < 0`)
"""
function variables(z::PolyZeros)
  _variables!(isreal(z), z.z, z.indicator)
end

function _variables!(creal::Bool, z::Vector{S}, ind::BitArray) where S<:Number
  m = length(z)
  T = real(S)
  if creal && (! (S<:Real) || any(ind) )
    v = Array{T}(undef, m)
    k = 1
    while k <= m
      zk = z[k]
      zi = imag(zk)
      if isreal(zk)
        if !ind[k]
          v[k] = real(zk)
        else
          zi = z[k+1]
          v[k] = -(zk + zi)
          v[k+1] = zk * zi
          k += 1
        end
      else
        if z[k+1] != conj(zk)
          throw(ArgumentError("missing second of conjugate complex pair"))
        end
        v[k] = -2real(zk)
        v[k+1] = abs2(zk)
        ind[k] = true
        ind[k+m+1] = true
        k += 1
      end
      k += 1
    end
  else
    v = copy(z)
  end
  v
end

"""
Create a copy of the PolyZeros with root values replaced by variables.
In the real case, the variable characteristic is used - the variables must
contain the real coefficients of quadratic divisor at the appropriate positions.
"""
function Base.copy(pz::PolyZeros{S}, v::AbstractVector{T}) where {S<:Number,T<:Number}
  
  m = length(v)
  m == length(pz.z) || throw(ArgumentError("length mismatch"))
  

  newind = pz.indicator
  if isreal(pz) && (T <: Real) && any(newind)
    z = Array{complex(S)}(m)
    Z = zero(T)
    TWO = T(2)
    k = 1
    while k <= m
      zk = v[k]
      if var_real(pz, k)
        zr = -zk / TWO
        zi = zr^2 - v[k+1]
        if zi < Z
          zi = sqrt(-zi)
          z[k] = complex(zr, zi)
          z[k+1] = complex(zr, -zi)
        else
          zi = sqrt(zi)
          z[k] = zr - zi
          z[k+1] = zr + zi
        end
        k += 2
      else
        z[k] = zk
        k += 1
      end
    end
    if isreal(z)
      z = real(z)
    end
  else
    z = v
  end
  PolyZeros(z, pz.mult, pz.real_coefficients, pz.indicator)
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

