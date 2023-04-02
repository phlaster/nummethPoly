% Задача. Запрограммировать степенной метод,
% найти максимальное и минимальное собственные 
% числа (с.ч.)
% БАЗА (0) Создание матрицы. Нужно создать матрицу с известными с.ч
% Построение основано на свойстве подобного преобразования,
% которое не изменяет с.ч. матрицы.
% 1 способ. При помощи невырожденной матрицы В. Если есть диагональная
% матрица D (с с.ч. на диагонали), то у матрицы А=B^(-1)DB будут те же
% самые с.ч. Матрица А в общем случае не будет симметричной.
% Положительная определенность зависит от знаков элементов диагональной 
% матрицы D

% 2 способ. При помощи ортогональной матрицы Q. Если есть диагональная
% матрица D (с с.ч. на диагонали), то у матрицы А=Q'DQ будут те же самые
% с.ч. Особенность матрицы А в том, что она будет симметричной.
N = 6;
v = [1:N];
D = diag(v);% + tril(rand(N), -1);
[Q, R] = qr(rand(N));
A = Q' * D * Q;
eig(A)
csvwrite('mtr10x10_known_eig.txt',A);