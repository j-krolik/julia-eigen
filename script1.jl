using LinearAlgebra
using Plots

"""function diagonal(dim)
    tmp = zeros(dim,dim)
    for k in 1:dim
        tmp[k,k] = 1
    end
    tmp
end"""


function power_metod(A::Matrix{Int64}, iteration=40)
    dimensions = size(A)[1]
    eigen_vector = ones(dimensions,1)
    eigen_value = norm(eigen_vector)
    for k in 1:iteration
        eigen_vector = A*eigen_vector/eigen_value
        eigen_value = norm(eigen_vector)
    end
    return Eigen([eigen_value], eigen_vector)
end

function qr_q(A)
    dim = size(A)[1]
    #U::typeof(A)(undef, dim, dim)
    #proj::typeof(A)(undef, dim-1, dim)

    """    if dim>1
        A_next = A[:,2] - ()
    end"""

    U = zeros(dim, dim)
    for col in 1:dim
        U[:,col] = A[:,col]
        for col_proj in 1:(col-1)
            U[:,col] -= (U[:,col_proj]'*A[:,col]) / (U[:,col_proj]'*U[:,col_proj]) * U[:,col_proj]
        end
    end
    for col in 1:dim
        rt_mean = sqrt(sum(U[:,col].^2))
        U[:,col] /= rt_mean
    end


    #
    #A_next = rt_mean*A[:,1]
    #A_next
    U
end

function eigen_qr(A, iteration=70)
    dim = ndims(A)
    for k in 1:iteration
        #QR = qr(A)
        #A = QR.Q'*A*QR.Q
        Q = qr_q(A)
        A = Q'*A*Q
    end
    eigen_var = Vector{Float64}(undef, dim)
    for k in 1:dim
        eigen_var[k] = A[k,k]
    end
    return eigen_var
end

struct Eigen_foo{T,V,S<:AbstractMatrix,U<:AbstractVector} <: Factorization{T}
    values::U
    vectors::S
    Eigen_foo{T,V,S,U}(values::AbstractVector{V}, vectors::AbstractMatrix{T}) where {T,V,S,U} =
        new(values, vectors)
end
Eigen_foo(values::AbstractVector{V}, vectors::AbstractMatrix{T}) where {T,V} =
    Eigen_foo{T,V,typeof(vectors),typeof(values)}(values, vectors)

# A = [3 1; 1 2]
A = [0 1; 1 1]
foo = [12 -51 4; 6 167 -68; -4 24 -41]
qr_q(foo)
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
