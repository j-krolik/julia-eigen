function qr_house_old!(A,Q)
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

function qr_house_new!(A,Q)
    n::Int64,m::Int64 = size(A)
    for k in 1:n
        z = A[k:m, k]
        v = [ -sign(z[1])*norm(z) - z[1]; -z[2:end] ]
        v /= sqrt(v'*v)

        #for j in 1:n
        #    A[k:m, j] = A[k:m, j] - v*( 2*(v'*A[k:m,j]) )
        #end

            Q[k:m, :] = Q[k:m, :] - v*( 2*(v'*Q[k:m,:]) )
        
    end
    Q = Q'
    nothing
end
