"""

  Construct the sylvester matrix of degree k from polynomials f and g = f'.
  The matrix has 2k+1 columns. Each column contains a shifted version of f or g.
  Order of columns: k+1 * g, k * f
"""
function sylves1(f::AbstractVector{T}, g::AbstractVector{T}, k::Int) where {T<:Number}
#
#  The k-th sylvester matrix of f(x), f'(x) = g(x)
#
    l = length(f)
    m = l + k - 1
    n = 2 * k + 1
    A = zeros(T, m, n)
	for j = 1:k+1
        A[j:j+l-2,j] = g
    end
    for j = 1:k
        A[j:j+l-1,j+k+1] = f
    end
    A
end

"""
Extract result polynomials u and v if x is a solution of sylves1(f,g) * x == 0.
  u * g + v * f == 0, degree of v is k.
"""
function extract_sylves1(x::AbstractVector{T}) where {T<:Number}
  k = length(x)
  k > 2 && k % 2 == 1 || error("vector too short or not odd length")
  k = (k + 1) ÷ 2
  x[1:k] / x[1], x[k+1:end] / x[k+1]  # u, v
end

"""

  Construct the sylvester matrix of degree k from polynomials f and g = f'.
  The matrix has 2k+1 columns. Each column contains a shifted version of f or g.
  Order of columns: g, k * [f g]
"""
function sylves2(f::AbstractVector{T}, g::AbstractVector{T}, k::Int) where {T<:Number}
#
#  The k-th sylvester matrix of f(x), f'(x) = g(x)
#
    l = length(f)
    m = l + k - 1
    n = 2 * k + 1
    A = zeros(T, m, n)
	for j = 1:k+1
        A[j:j+l-2,j*2-1] = g
    end
    for j = 1:k
        A[j:j+l-1,j*2] = f
    end
    A
end
   
"""
Extract result polynomials u and v if x is a solution of sylves2(f,g) * x == 0.
  u * g + v * f == 0, degree of v is k.
"""
function extract_sylves2(x::AbstractVector{T}) where {T<:Number}
  k = length(x)
  k > 2 && k % 2 == 1 || error("vector too short or not odd length")
  x[1:2:k] / x[1], x[2:2:end] / x[2]  # u, v
end

