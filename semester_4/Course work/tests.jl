using LinearAlgebra, Statistics

w0 = [4,5,1,6];
w = w0 ./ norm(w0);
H = I - 2w*w';
D = diagm([5,3,2,1]);
A = H'*D*H;

C = cholesky(Hermitian(A))
C.L' * C.U'
