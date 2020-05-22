using LinearAlgebra
using Plots

function start_power(A)
    d=typeof(A)
    if (d!=Matrix{Float64})
        A=convert(Matrix{Float64},A)
    end
   return E_power=power_metod(A)
end

function power_metod(A::Matrix{Float64}, iteration=40)
    dimeansions::Int32 = size(A)[1]
    i::Int32 = size(A)[2]
    if (i < dimeansions)
        dimeansions=i
    end
    A=A[1:dimeansions,1:dimeansions]

    eigen_vector = ones(dimeansions,1)
    eigen_value = norm(eigen_vector)
     @inbounds for k in 1:iteration
        eigen_vector = A*eigen_vector./eigen_value
        eigen_value = norm(eigen_vector)
    end
    return Eigen([eigen_value], eigen_vector)
end


function start_qr(A)
    if (typeof(A)!=Matrix{Float64})
        A=convert(Matrix{Float64},A)
    end
   return E_qr=qr_q_eigen(A)
end

function qr_q_eigen(A::Matrix{T}) where {T<:Float64}
    dim::Int32 = size(A)[1]
    i::Int32 = size(A)[2]
    if (i < dim)
        dim=i
    end
    A=A[1:dim,1:dim]
    U::Matrix{Float64} = zeros(dim, dim)
     @inbounds for col in 1:dim
        U[:,col] =A[:,col]
        @inbounds for col_proj in 1:(col-1)
            U[:,col] -= (U[:,col_proj]'*A[:,col]) / (U[:,col_proj]'*U[:,col_proj]) * U[:,col_proj]
        end
    end
    @inbounds for col in 1:dim
        U[:,col] /= sqrt(sum(U[:,col].^2))
    end
    A_cpy = A

    for k in 1:70
        #QR = qr(A)
        #A = QR.Q'*A*QR.Q
        U=qr_q(A)   #niepodoba mi się ten moment, ale nie chce za bardzo tutaj ingerować, żeby się algorytm nie wykrzaczył
        A = U'*A*U

        #A = R*Q
        #?A = Q*A*Q'
    end
    eigenvalue = Vector{Float64}(undef, dim)
    eigenvector = A_cpy
    @inbounds for k in 1:dim
        eigenvalue[k] = A[k,k]
        eigenvector[:,k] = (A_cpy - Array{Float64,2}(I,dim,dim)*A[k,k])*ones(dim,1)
        eigenvector[:,k] /= max(broadcast(abs, eigenvector[:,k])...)
    end
    return eigenvalue, eigenvector
end

function qr_q(A)
    dim::Int32 = size(A)[1]
    i::Int32 = size(A)[2]
    if (i < dim)
        dim=i
    end
    A=A[1:dim,1:dim]
    U = zeros(dim, dim)
     @inbounds for col in 1:dim
        U[:,col] =(A[:,col])
        @inbounds for col_proj in 1:(col-1)
            U[:,col] -= (U[:,col_proj]'*A[:,col]) / (U[:,col_proj]'*U[:,col_proj]) * U[:,col_proj]
        end
    end
    @inbounds for col in 1:dim
        U[:,col] /= sqrt(sum(U[:,col].^2))
    end
    return U
end


function start_house(A)
    d=typeof(A)
    if (d!=Matrix{Float64})
        A=convert(Matrix{Float64},A)
    end
   return E_house=qr_house_eigen(A)
end

function qr_house_eigen(A::Matrix{Float64})
    m:: Int32,n:: Int32 = size(A)
    if (m>n)
     m=n
    end
    if (n>m)
        n=m
    end
    A=A[1:m,1:n]
    #Q = Typeof(A)(I, m,n)
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
    Q = Q'
    #R = triu(A)
    A_cpy = A
    for k in 1:70
        #QR = qr(A)
        #A = QR.Q'*A*QR.Q
        #Q = qr_q(A)
        Q = qr_house(A) #niepodoba mi się ten moment, ale nie chce za bardzo tutaj ingerować, żeby się algorytm nie wykrzaczył
        A = Q'*A*Q
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
    return eigenvalue, eigenvector
    end

    function qr_house(A)
        m:: Int64,n:: Int64 = size(A)
        if (m>n)
            m=n
        end
        if (n>m)
            n=m
        end

        A=A[1:m,1:n]

        #Q = Typeof(A)(I, m,n)
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
        Q = Q'
        #R = triu(A)
        return Q
    end




"""struct Eigen_foo{T,V,S<:AbstractMatrix,U<:AbstractVector} <: Factorization{T}
    values::U
    vectors::S
    Eigen_foo{T,V,S,U}(values::AbstractVector{V}, vectors::AbstractMatrix{T}) where {T,V,S,U} =
        new(values, vectors)
endq
Eigen_foo(values::AbstractVector{V}, vectors::AbstractMatrix{T}) where {T,V} =
    Eigen_foo{T,V,typeof(vectors),typeof(values)}(values, vectors)"""


foo = [2 -51; 6 167; 4 27]
a=start_power(foo)
b=start_qr(foo)
c=start_house(foo)





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



function eigen_qr(A, iteration=10)
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
function eigen_house(A, iteration=70)
    dim = size(A)[1]
    A_cpy = A
    for k in 1:iteration
        #QR = qr(A)
        #A = QR.Q'*A*QR.Q

        #Q = qr_q(A)
        Q = qr_house(A)
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

#qr_poiq=qr_house(foo)
#E = eigen(foo)
E_power = power_metod(foo)
E_qr = eigen_qr(foo)
#E_qra = eigen_house(foo)
#E_qriki = qr_q_eigen(foo)
#E_house1=qr_house_eigen(foo)
