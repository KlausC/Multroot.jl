function forsub(A::AbstractArray{T,2}, b::AbstractVector{T}) where {T<:Number}
#
# assume A is an lower triangular matrix
# forsub performs forward substitution to solve Ax=b
#
    n = size(A,1)
    LowerTriangular(A) \ b[1:n]
end
