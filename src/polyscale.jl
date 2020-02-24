import Polynomials
using Polynomials

#
# given the coefficients of a polynomial p(x) = p(n)*x^n + + p(1)*x + p(0)
#
# determne a scaling factor for x and a linear factor t that the transformed
# polynomial q(z) = t * p(s*z) has coefficient p(n) = 1 and all |p(j)| <= 1
#
# Invariant:
# if q, shift, smax, u = polyscale(p) then q = polytransform(p, shift, smax) * u
# up to numerical precision

FloatOrComplex = Union{AbstractFloat, Complex, Complex{Float32}, Complex{BigFloat}}

function polyscale(P::Poly{T}) where {T<:FloatOrComplex}
	  n = length(P.a) -1
    shift = - P.a[n] / n	
    P = polytransform(P, shift)	
	  p = P.a
    p[n] = 0
    pni = 1 / p[n+1]
    smaxlog = maximum(log(abs(p[1:n] * pni)) ./ (n+1 - (1:n)))
    smax = exp(smaxlog)
    u = pni / exp(smaxlog * n)
    P = polytransform(Poly(p), 0.0, smax) * u
    P.a[n+1] = T(1)
    P, shift, smax, u
end


"""
given a monic polynomial (x^n + ...) determine a scaling factor s that the transformed
polynomial q(z) = p(s*z) is half of the time < 1, rest >= 1. Ignore zero coefficients.
"""
function polyscale2(P::Poly{T}) where T<:FloatOrComplex

  p = abs.(P.a)
  n = length(p) - 1
  pni = T(1) / p[n+1]
  scale!(p, pni)
  p = view(p,1:n)
  bigg = log(realmax(T))
  tiny = log(realmin(T))
  clip(x::T) = x == -Inf ? -x : x

  llim = exp(maximum(((log.(p) - bigg) ./ (n:-1:1))))
  ulim = exp(minimum(clip, ((log.(p) - tiny) ./ (n:-1:1)))) 
  med = exp(sort(log.(p) ./ (n:-1:1))[n÷2+1])




  med, llim, ulim

end

"""

Find scaling factor for a polynomial p_1 + p_2*x^1 + ... + p_(n+1)*x^n,
which has the minimal possible quotient of coefficients.
Coefficients which are zero ignored.
"""
function polyscale_mini_quotient(P::Poly{T}) where T<:FloatOrComplex
  p = abs.(P.a)
  n = length(p) - 1
  pni = log(abs(p[n+1]))
  
  plog = log.(p[1:n]) - pni

  f = Array{real(T)}(2n)
  copy!(f, plog)
  copy!(f, n+1, -f, 1, n)
  g = Array{real(T)}(2n)
  g =
  copy!(g, -(n:-1:1))
  copy!(g, n+1, -g, 1, n)

  xa, xb, val = minimize(f, g)
  exp(xa), exp(xb), exp(val)
end

"""

  `minimize(f::Array, g::Array of same size) -> xa, xb, vmin`


Solve the following special linear optimization problem:
Define `f(x) = max{ f1[k] + g1[k] * x / k = 1..n1}`.
Minimize `f(x)` for real x.
The function f is convex as maximum of linear (convex) functions is convex.
If f is bounded below, it has minmum value `vmin` and a solution interval `[xa, xb]`
with `f(x) = vmin for all x in [xa, xb]`.
If f is not bounded below, return `vmin = -Inf` and `xa = xb = ±Inf`.
"""
function minimize(f::AbstractVector{T}, g::AbstractVector{T}) where T<:Real

  n::Int = length(f)
  n == length(g) || error(ArgumentError("lengths of f1 and g1 are different"))

  Z = zero(T)
  inf = typemax(T)

  # remove unwanted infinities
  f, g = unzip(Iterators.filter(x -> x[1] !=  inf, zip(f, g)))::NTuple{2,AbstractArray{T}}
  n = length(f)
  
  if n == 0
    return -inf, inf, -inf
  elseif n == 1
    if g[1] == Z
      return -inf, inf, f[1]
    elseif g[1] < Z
      return inf, inf, -inf
    else
      return -inf, -inf, -inf
    end
  end

  eva(x::T)::Tuple{T,T,Int} = maximum(1:n) do k; ((f[k] + g[k] * x), g[k], k) end

  function nextto(x::T)::Tuple{T, T, T}
    v, gk, k = eva(x)
    cross(j::Int) = -(f[j] - f[k]) / (g[j] - gk)
    
    if isnan(v)
      xa = xb = x
    elseif gk > Z
      xa = xb = maximum(filter(y->y<x, cross.(filter(j->j!=k, 1:n))))
    elseif gk < Z
      xa = xb = minimum(filter(y->y>x, cross.(filter(j->j!=k, 1:n))))
    elseif gk == Z
      xa = maximum(filter(y->y<x, cross.(filter(j->j!=k, 1:n))))
      xb = minimum(filter(y->y>x, cross.(filter(j->j!=k, 1:n))))
    else
      xa, xb = inf, -inf
    end
    xa, xb, gk
  end

  # find extreme slopes
  gmax::T, fmax::T, kmax::Int = maximum(1:n) do k; (g[k], f[k], k) end
  gmin::T, fmin::T, kmin::Int = minimum(1:n) do k; (g[k], f[k], k) end
  
  xmax::T = maximum(filter(x->!isnan(x), -(f-f[kmax]) ./ (g - gmax)))
  xmin::T = minimum(filter(x->!isnan(x), -(f-f[kmin]) ./ (g - gmin)))
  xa, xb = xmax, xmin

  while xmin < xa || xb < xmax
    xa, xb, gab = nextto((xmin + xmax) / 2)
    gs::Int = cmp(gab, Z)
    if gs >= 0
      xmax = xb
    end
    if gs <= 0
      xmin = xa
    end
  end
  xa, xb, eva((xa + xb)/ 2)...
end

maximum1(f::Function, v0, itr) = mapreduce(f, Base.scalarmax, v0, itr)
maximum1(v0, itr) = mapreduce(identity, Base.scalarmax, v0, itr)
maximum1(f::Function, itr) = mapreduce(f, Base.scalarmax, lowerbound(f, itr), itr)
maximum1(itr) = mapreduce(identity, Base.scalarmax, lowerbound(itr), itr)
lowerbound(itr) = lowerbound(eltype(itr))
lowerbound(f::Function, itr) = lowerbound(promote_type(Base.return_types(f, (eltype(itr),))))
lowerbound(T::Type) = error("no lower bound for type $T")
lowerbound(::Type{T}) where T<:Union{AbstractFloat,Unsigned} = typemin(T)
lowerbound(S::Type{Rational{T}}) where T<:Integer = typemin(S)
lowerbound(::Type{String}) = ""

#Base.mr_empty(::typeof(identity),::typeof(Base.scalarmax),::Type{T}) where T<:Real = typemin(T)
#Base.mr_empty(::typeof(identity),::typeof(Base.scalarmin),::Type{T}) where T<:Real = typemax(T)

"""
  `unzip(itr)`

The iterable object itr must deliver N-tuples.
Convert an array of N-tuples into N-tuple of arrays.
"""
function unzip(aot)
  tuple(map(x->Array(collect(x)), zip(aot...))...)
end

