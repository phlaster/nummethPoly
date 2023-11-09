using LinearAlgebra, Statistics

L = rand(-10:10, 5,5) |> LowerTriangular |> Matrix
U = rand(-10:10, 5,5) |> UpperTriangular |> Matrix
x = rand(1:9, 5)
bL = L*x
bU = U*x
L\b

w0 = [4,5,1,6];
w = w0 ./ norm(w0);
H = I - 2w*w';
D = diagm([5,3,2,1]);
A = H'*D*H;

A[1,4] = A[4,1] = 0

C = cholesky(Hermitian(A))
C.L' * C.U'
