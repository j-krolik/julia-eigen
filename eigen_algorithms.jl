using LinearAlgebra

function eigen_power(A, print_it=false)
    d=typeof(A)
    if (d!=Matrix{Float64})
        A=convert(Matrix{Float64},A)
    end
   return power_metod(A, print_it)
end

function power_metod(A::Matrix{T}, print_it=false, iteration=nothing, delta=5*10^-6) where {T<:Float64}
    dim::Int32 = size(A)[1]
    A=A[1:dim,1:dim]

    eigen_vector = ones(dim,1)
    eigen_value = norm(eigen_vector)
    if iteration==nothing iteration=10000 end
    last_value::T = 0
     @inbounds for k in 1:iteration
        eigen_vector = A*eigen_vector./eigen_value
        eigen_value = norm(eigen_vector)
        # if enought delta enought small, break
        if abs(last_value-eigen_value)<delta print_it&&k!=1&&println(k); break end
        last_value=eigen_value
    end
    return Eigen([eigen_value], eigen_vector)
end


function eigen_qr(A::Matrix, print_it=false)
    if (typeof(A)!=Matrix{Float64})
        A=convert(Matrix{Float64},A)
    end
   return qr_q_eigen(A, print_it)
end

function qr_q_eigen(A::Matrix{T}, print_iteration=false, delta=5*10^-6) where {T<:Float64}
    dim::Int32 = size(A)[1]
    i::Int32 = size(A)[2]
    if (i < dim)
        dim=i
    end
    A=A[1:dim,1:dim]
    A_cpy = A
    Q = Matrix{T}(undef,dim,dim)
    last_value::T = last(A)
    for k in 1:10000
        # QR factorization
        @inbounds for col in 1:dim
           Q[:,col] =(A[:,col])
           @inbounds for col_proj in 1:(col-1)
               Q[:,col] -= (Q[:,col_proj]'*A[:,col]) / (Q[:,col_proj]'*Q[:,col_proj]) * Q[:,col_proj]
           end
        end
        @inbounds for col in 1:dim
           Q[:,col] /= sqrt(sum(Q[:,col].^2))
        end
        # calc eigen
        A = Q'*A*Q
        # if enought delta enought small, break
        if abs(last(A)-last_value)<delta print_iteration&&println(k); break end
        last_value=last(A)

    end
    eigenvalue = Vector{T}(undef, dim)
    eigenvector = Matrix{Float64}(undef,dim,dim)
    @inbounds for k in 1:dim
        eigenvalue[k] = A[k,k]
        eigenvector[:,k] = (A_cpy - Array{Float64,2}(I,dim,dim)*A[k,k])*ones(dim,1)
        eigenvector[:,k] /= max(broadcast(abs, eigenvector[:,k])...)
    end
    return Eigen([eigenvalue], eigenvector)
end

function qr_q!(A::Matrix{T},U::Matrix{T}) where {T<:Float64}
    dim::Int32 = size(A)[1]
    #U = zeros(dim, dim)
     @inbounds for col in 1:dim
        U[:,col] =(A[:,col])
        @inbounds for col_proj in 1:(col-1)
            U[:,col] -= (U[:,col_proj]'*A[:,col]) / (U[:,col_proj]'*U[:,col_proj]) * U[:,col_proj]
        end
    end
    @inbounds for col in 1:dim
        U[:,col] /= sqrt(sum(U[:,col].^2))
    end
    nothing
end


function eigen_householder(A::Matrix, print_it=false)
    d=typeof(A)
    if (d!=Matrix{Float64})
        A=convert(Matrix{Float64},A)
    end
   return qr_householder_eigen!(A, print_it)
end

function qr_householder_eigen!(A::Matrix{T}, print_iteration::Bool=false, delta=5*10^-6) where {T<:Float64}
    dim::Int32 = size(A)[1]
    A_cpy = A
    # QR decomposition
    last_value::T = last(A)
    Q = Array{T,2}(I,dim,dim)
    for k in 1:100
        # Householder factorization
        for k in 1:dim
            z = A[k:dim, k]
            v = [ -sign(z[1])*norm(z) - z[1]; -z[2:end] ]
            v /= sqrt(v'*v)
            Q[k:dim,:] = Q[k:dim,:] - v*( 2*(v'*Q[k:dim,:]) )
        end
        A = Q*A*Q'
        # if enought delta enought small, break
        if abs(last(A)-last_value)<delta print_iteration&&println(k); break end
        last_value=last(A)
        Q = Array{Float64,2}(I,dim,dim)
    end
    # determine eigenvector
    eigenvalue = Vector{Float64}(undef, dim)
    eigenvector = A_cpy
    @inbounds for k in 1:dim
        eigenvalue[k] = A[k,k]
        eigenvector[:,k] = (A_cpy - Array{Float64,2}(I,dim,dim)*A[k,k])*ones(dim,1)
        eigenvector[:,k] /= max(broadcast(abs, eigenvector[:,k])...)
    end
    return Eigen([eigenvalue], eigenvector)
    end
