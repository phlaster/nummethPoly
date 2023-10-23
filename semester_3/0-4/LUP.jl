using LinearAlgebra

function LUP(M)
    @assert ==(size(M)...) "Только квадратные матрицы!"
    N = size(M)[1]
    perm = [1:N;]
    for row in 1:N
        find_max_and_swap!(M, perm, row)
        div_under_main_diag!(M, row)
        subtract_product!(M, row)
    end
    L, U = split_LUP(M)
    return L, U, perm
end

function find_max_and_swap!(M, perm, i)
    maxval, row = findmax(abs.(M[:, i]))
    @assert abs(maxval) > 1e-16 "Невозможно найти максимальнный элемент в нулевом столбце $i!"
    M[i, :], M[row, :] = M[row, :], M[i, :]
    perm[i], perm[row] = perm[row], perm[i]
end

function div_under_main_diag!(M, i)
    M[i+1:end, i] .= M[i+1:end, i] ./ M[i,i]
end

function subtract_product!(M, i)
    N = size(M)[1]
    for j=i+1:N, k=i+1:N
        M[j,k] -= M[j,i]*M[i,k]
    end
end


function split_LUP(M)
    U = triu(M)
    L = tril(M)
    L[diagind(L)] .= 1
    return L, U
end