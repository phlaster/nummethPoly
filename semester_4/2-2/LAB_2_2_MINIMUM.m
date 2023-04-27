clc;
% БАЗА (0) Создание матрицы. Нужно создать матрицу с известными с.ч.
% Построение основано на свойстве подобного преобразования, которое не
% изменяет с.ч. матрицы.

% 1 способ. При помощи невырожденной матрицы В. Если есть диагональная
% матрица D (с с.ч. на диагонали), то у матрицы А=B-1
% DB будут те же самые с.ч. Матрица А в общем случае не будет симметричной.
% Положительная определенность зависит от знаков элементов диагональной
% матрицы D

% 2 способ. Создание несимметричной матрицы при помощи ортогональной матрицы
% Q. Если есть треугольная (верхняя или нижняя) матрица B (с с.ч. на
% диагонали), то у матрицы А=Q^ТВQ будут те же самые с.ч. и матрица при
% этом получится несимметричной
III = [];
PPP = [];
for N = [10, 100]
pos = 1;
v = [1:N];
D = diag(v) + tril(rand(N), -1);
[Q, R] = qr(rand(N));

% МИНИМУМ (+1) При помощи QR-алгоритма найти собственные числа матрицы
% – Создать матрицу А с действительным спектром
A = Q' * D * Q;
eig(A);
% – Одна итерация это: Получение матриц Q (ортогональную) и R (верхнюю
% треугольную) разложением матрицы А ( qr() ). Вычисление следующей матрицы
% умножением R на Q:
% [Q,R]: A=QR; A=RQ
epsilon = logspace(-1, -10, 50);
iterations = [];
prescision = [];
iters = 0;
A = hess(A);
for eps = epsilon
    while max(abs(tril(A, -1)), [], "all") > eps
        [Q, R] = qr(A);
        A = R * Q;
        iters = iters + 1;
    end
    iterations = [iterations iters];
    prescision = [prescision max(sort(diag(A)) - v')];
end
III = [III iterations'];
PPP = [PPP prescision'];
end

subplot(1, 2, 1)
semilogx(epsilon, III(:, 1))
grid on
xlabel('epsilon')
ylabel('iterations')
hold on

subplot(1, 2, 1)
semilogx(epsilon, III(:, 2))
grid on
legend('10x10', '100x100')
title("Количество шагов")

subplot(1, 2, 2)
loglog(epsilon, PPP(:, 1))
grid on
xlabel('epsilon')
ylabel('iterations')
hold on

subplot(1, 2, 2)
loglog(epsilon, PPP(:, 2))
grid on
legend('10x10', '100x100', 'Location', "northwest")
title("Фактическая точность")
xlabel('epsilon')
ylabel('prescision')