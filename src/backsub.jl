function backsub(A::AbstractArray{T,2}, b::AbstractVector{T}) where {T<:Number}
#
# assume A is an upper triangular matrix
# backsub performs backward substitution to solve Ax=b
#
    n = size(A,1)
    N = zero(T)
    s = maximum(abs, diag(A)) * eps(real(T))
    for k = 1:n
        if iszero(A[k,k])
            A[k,k] = s
        end
    end
    UpperTriangular(A) \ b[1:n]
end

