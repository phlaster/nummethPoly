using LinearAlgebra

rnd(x) = round.(x, digits=3)

D = diagm([7,2,1,-2])

w0 = [0.68, 0.032, 0.585, 0.61]

w = round.(w0 ./ norm(w0), digits=3)

Q = I - 2 * w * w' |> rnd

cond(Q)

A = Q' * D * Q |> rnd 

cond(A)

x = [1,2,3,4]

b = A*x |> rnd

L, U, P = lu(A)

L |> rnd
U |> rnd

(L*U)[sortperm(P), 1:4] ≈ A

sortperm(P) 