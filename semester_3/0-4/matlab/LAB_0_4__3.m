% 3б
clc;
% Размерность
n = 3;

% Матрица с нулевым определителем
A = [
    1 2 3
    4 5 6
    7 8 9];

% Вектор правой стороны
b = rand(n, 1);

% Решение
% Warning: Matrix is singular to working precision.
% Warning: Matrix is close to singular or badly scaled. 
xA = A \ b;
disp(A);
disp(['Det(A): ' num2str(det(A))]);
disp(['Cond(A): ' num2str(cond(A))]);
disp(['||rA||: ' num2str(norm(b-A*xA))]);

B = 1e8*A;
xB = B \ b;
disp(B);
disp(['Det(B): ' num2str(det(B))]);
disp(['Cond(B): ' num2str(cond(B))]);
disp(['||rB||: ' num2str(norm(b-A*xB))]);
