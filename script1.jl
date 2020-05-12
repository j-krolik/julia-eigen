using LinearAlgebra
using Plots

function power_metod(A::Matrix{Float64}, iteration=40)
    dimensions = size(A)[1]
    eigen_vector = ones(dimeansions,1)
    eigen_value = norm(eigen_vector)
    for k in 1:iteration
        eigen_vector = A*eigen_vector/eigen_value
        eigen_value = norm(eigen_vector)
    end
    return Eigen([eigen_value], eigen_vector)
end

function qr_q(A)
    dim = size(A)[1]
    U = zeros(dim, dim)
    for col in 1:dim
        U[:,col] = A[:,col]
        for col_proj in 1:(col-1)
            U[:,col] -= (U[:,col_proj]'*A[:,col]) / (U[:,col_proj]'*U[:,col_proj]) * U[:,col_proj]
        end
    end
    for col in 1:dim
        U[:,col] /= sqrt(sum(U[:,col].^2))
    end
    U
end

function eigen_qr(A, iteration=70)
    dim = size(A)[1]
    A_cpy = A
    for k in 1:iteration
        #QR = qr(A)
        #A = QR.Q'*A*QR.Q

        Q = qr_q(A)
        #Q = qr_house(A)
        A = Q'*A*Q
        #A = R*Q
        #?A = Q*A*Q'
    end
    eigenvalue = Vector{Float64}(undef, dim)
    eigenvector = A_cpy
    for k in 1:dim
        eigenvalue[k] = A[k,k]
        eigenvector[:,k] = (A_cpy - Array{Float64,2}(I,dim,dim)*A[k,k])*ones(dim,1)
        eigenvector[:,k] /= max(broadcast(abs, eigenvector[:,k])...)
    end

    return eigenvalue, eigenvector
end

function qr_house(A)
    m,n = size(A)
    #Q = Typeof(A)(I, m,n)
    Q = Array{Float64,2}(I,m,n)
    for k in 1:n
        z = A[k:m, k]
        v = [ -sign(z[1])*norm(z) - z[1]; -z[2:end] ]
        v /= sqrt(v'*v)

        for j in 1:n
            A[k:m, j] = A[k:m, j] - v*( 2*(v'*A[k:m,j]) )
        end
        for j in 1:m
            Q[k:m, j] = Q[k:m, j] - v*( 2*(v'*Q[k:m,j]) )
        end
    end
    Q = Q'
    R = triu(A)
    return Q, R
end

struct Eigen_foo{T,V,S<:AbstractMatrix,U<:AbstractVector} <: Factorization{T}
    values::U
    vectors::S
    Eigen_foo{T,V,S,U}(values::AbstractVector{V}, vectors::AbstractMatrix{T}) where {T,V,S,U} =
        new(values, vectors)
endq
Eigen_foo(values::AbstractVector{V}, vectors::AbstractMatrix{T}) where {T,V} =
    Eigen_foo{T,V,typeof(vectors),typeof(values)}(values, vectors)

# A = [3 1; 1 2]
A = [0 1; 1 1]
foo = [12 -51 4; 6 167 -68; -4 24 -41] * 1.0
#qr_q(foo)
E = eigen(foo)
E_power = power_metod(foo)
E_qr = eigen_qr(foo)

"""const ITERATION = 20
val = Vector{Float64}(undef, ITERATION)
vector = Matrix{}(undef, ITERATION, 2)
for k in 10:ITERATION
    foo = power_metod(A,k)
    vector[k,1] = foo[1][1]
    vector[k,2] = foo[1][2]
end"""

"""plot(vector[10:20,1])
plot!(vector[:,2])"""
