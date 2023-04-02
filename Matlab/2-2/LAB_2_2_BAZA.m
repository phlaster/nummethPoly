% БАЗА (0) Создание матрицы. Нужно создать матрицу с известными с.ч.
% Построение основано на свойстве подобного преобразования, которое не
% изменяет с.ч. матрицы.


% 2 способ. Создание несимметричной матрицы при помощи ортогональной матрицы
% Q. Если есть треугольная (верхняя или нижняя) матрица B (с с.ч. на
% диагонали), то у матрицы А=Q^ТВQ будут те же самые с.ч. и матрица при
% этом получится несимметричной

N = 6;
v = [1:N];
D = diag(v) + tril(rand(N), -1);
[Q, R] = qr(rand(N));

A = Q' * D * Q