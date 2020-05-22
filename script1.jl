using LinearAlgebra
using Plots

function start_power(A, print_it=false)
    d=typeof(A)
    if (d!=Matrix{Float64})
        A=convert(Matrix{Float64},A)
    end
   return power_metod(A, print_it)
end

function power_metod(A::Matrix{T}, print_it=false, iteration=nothing, delta=10^-10) where {T<:Float64}
    dimeansions::Int32 = size(A)[1]
    i::Int32 = size(A)[2]
    if (i < dimeansions)
        dimeansions=i
    end
    A=A[1:dimeansions,1:dimeansions]

    eigen_vector = ones(dimeansions,1)
    eigen_value = norm(eigen_vector)
    if iteration==nothing iteration=10000 end
    last_value::T = 0
     @inbounds for k in 1:iteration
        eigen_vector = A*eigen_vector./eigen_value
        eigen_value = norm(eigen_vector)
        # chech if enought
        if abs(last_value-eigen_value)<delta print_it&&println(k); break end
        last_value=eigen_value
    end
    return Eigen([eigen_value], eigen_vector)
end


function start_qr(A::Matrix, print_it=false)
    if (typeof(A)!=Matrix{Float64})
        A=convert(Matrix{Float64},A)
    end
   #return E_qr=qr_q_eigen(A)
   return qr_q_eigen(A, print_it)
end

function qr_q_eigen(A::Matrix{T}, print_iteration=false, delta=10^-8) where {T<:Float64}
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
        # A-(QR)->Q
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
        # chech if enought
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


function start_house(A::Matrix, print_it=false)
    d=typeof(A)
    if (d!=Matrix{Float64})
        A=convert(Matrix{Float64},A)
    end
   return qr_house_eigen!(A, print_it)
end

function qr_house_eigen!(A::Matrix{T}, print_iteration::Bool=false, delta=10^-8) where {T<:Float64}
    m:: Int32,n:: Int32 = size(A)
    if (m>n)
     m=n
    end
    if (n>m)
        n=m
    end
    A=A[1:m,1:n]
    #Q = Typeof(A)(I, m,n)
    Q = Array{T,2}(I,m,n)
    for k in 1:n
        z = A[k:m, k]
        v = [ -sign(z[1])*norm(z) - z[1]; -z[2:end] ]
        v /= sqrt(v'*v)
        #for j in 1:n
        #    A[k:m, j] = A[k:m, j] - v*( 2*(v'*A[k:m,j]) )
        #end
        for j in 1:m
            Q[k:m, j] = Q[k:m, j] - v*( 2*(v'*Q[k:m,j]) )
        end
    end
    Q = Q'
    #R = triu(A)
    A_cpy = A
    last_value::T = last(A)
    for k in 1:100
        #QR = qr(A)
        #A = QR.Q'*A*QR.Q
        #Q = qr_q(A)
        Q = Array{Float64,2}(I,m,n)
        for k in 1:n
            z = A[k:m, k]
            v = [ -sign(z[1])*norm(z) - z[1]; -z[2:end] ]
            v /= sqrt(v'*v)

            #for j in 1:n
            #    A[k:m, j] = A[k:m, j] - v*( 2*(v'*A[k:m,j]) )
            #end
            for j in 1:m
                Q[k:m, j] = Q[k:m, j] - v*( 2*(v'*Q[k:m,j]) )
            end
        end

        #Q = Q'
        A = Q*A*Q'

        if abs(last(A)-last_value)<delta print_iteration&&println(k); break end
        last_value=last(A)
        #A = R*Q
        #?A = Q*A*Q'
    end
    eigenvalue = Vector{Float64}(undef, m)
    eigenvector = A_cpy
    @inbounds for k in 1:m
        eigenvalue[k] = A[k,k]
        eigenvector[:,k] = (A_cpy - Array{Float64,2}(I,m,n)*A[k,k])*ones(m,1)
        eigenvector[:,k] /= max(broadcast(abs, eigenvector[:,k])...)
    end
    return Eigen([eigenvalue], eigenvector)
    end

    function qr_house!(A,Q)
        n::Int64,m::Int64 = size(A)
        for k in 1:n
            z = A[k:m, k]
            v = [ -sign(z[1])*norm(z) - z[1]; -z[2:end] ]
            v /= sqrt(v'*v)

            #for j in 1:n
            #    A[k:m, j] = A[k:m, j] - v*( 2*(v'*A[k:m,j]) )
            #end
            for j in 1:m
                Q[k:m, j] = Q[k:m, j] - v*( 2*(v'*Q[k:m,j]) )
            end
        end
        Q = Q'
        nothing
    end

# A = [3 1; 1 2]
A = [0 1; 1 1]
foo = [12 -51 4; 6 167 -68; -4 24 -41] * 1.0
#qr_q(foo)
E = eigen(foo)
E_power = power_metod(foo)
E_qr = eigen_qr(foo)
